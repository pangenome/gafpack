//! Copy number estimation via negative binomial model
//!
//! Ports floco's statistical approach to Rust:
//! - MLE estimation of NB parameters via mixture model (alpha, beta)
//! - CN probability calculation for each node

use log::{debug, info, warn};
use statrs::distribution::{Discrete, NegativeBinomial};

// ─── Constants ───────────────────────────────────────────────────────────────

const OUTLIER_PCT: f64 = 0.03;      // Remove top/bottom 3% of bins
const MAX_CN_SEARCH: i32 = 100;     // Maximum CN value to consider
const MIN_EXTENSION: i32 = 1;       // Minimum steps before stopping CN search (floco: MIN_EXTENSION)
const NELDER_MEAD_ITERS: usize = 200;  // Max iterations for optimization
const ROUND_BINS: usize = 5;        // Bin rounding factor (floco: bp_step = bin_size / 5)
const EPSILON_MIX: f64 = 0.01;      // Minimum CN value in mixture model
const READ_LENGTH_SUBSAMPLE: usize = 100_000;  // Max reads for distribution fitting (floco: subset_size)

// ─── Parameter Estimation ────────────────────────────────────────────────────

/// Estimate α (mean coverage per bp at CN=1) and β (std dev) via MLE
/// Uses floco's mixture model approach: fits NB mixture with CN proportions
/// Tests each ploidy value, returns params with best likelihood
pub fn estimate_nb_params(bins: &[f64], bin_size: usize, ploidies: &[u32]) -> (f64, f64) {
    if bins.is_empty() {
        warn!("No bins for parameter estimation, using defaults");
        return (0.01, 0.01);
    }

    info!("Estimating NB parameters from {} bins", bins.len());
    debug!("Bin size: {}, testing ploidies: {:?}", bin_size, ploidies);

    // Filter outliers (top/bottom 3%)
    let filtered = filter_outliers(bins, OUTLIER_PCT);
    if filtered.is_empty() {
        warn!("All bins filtered out, using defaults");
        return (0.01, 0.01);
    }

    debug!("After outlier filtering: {} bins remain", filtered.len());

    // Convert to histogram counts (like floco's np.bincount)
    let bp_step = bin_size / ROUND_BINS;
    let max_val = filtered.iter().cloned().fold(0.0f64, f64::max);
    let max_bin_idx = (max_val / bp_step as f64).round() as usize + 1;

    let mut counts: Vec<u32> = vec![0; max_bin_idx + 1];
    for &val in &filtered {
        let idx = (val / bp_step as f64).round() as usize;
        if idx < counts.len() {
            counts[idx] += 1;
        }
    }

    debug!("Created histogram with {} bins (bp_step={})", counts.len(), bp_step);

    // Test each ploidy with mixture model
    let results: Vec<_> = ploidies
        .iter()
        .map(|&p| {
            let (ll, alpha, beta) = fit_mixture_at_ploidy(&counts, bp_step, bin_size, p);
            debug!("Ploidy {}: LL={:.2}, α={:.6}, β={:.6}", p, ll, alpha, beta);
            (-ll, alpha, beta, p)  // Return negative LL for comparison (lower is better)
        })
        .collect();

    // Find best (lowest negative log-likelihood)
    let best = results
        .iter()
        .min_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

    match best {
        Some(&(neg_ll, alpha, beta, ploidy)) => {
            info!("Best fit: ploidy={}, LL={:.2}, α={:.6}, β={:.6}", ploidy, -neg_ll, alpha, beta);
            (alpha, beta)
        }
        None => {
            warn!("No valid parameter estimates, using defaults");
            (0.01, 0.01)
        }
    }
}

/// Filter top/bottom percentile of data
fn filter_outliers(data: &[f64], pct: f64) -> Vec<f64> {
    if data.len() < 10 {
        debug!("Too few data points ({}) for outlier filtering", data.len());
        return data.to_vec();
    }

    let mut sorted = data.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

    let lo = (sorted.len() as f64 * pct) as usize;
    let hi = (sorted.len() as f64 * (1.0 - pct)) as usize;

    debug!("Filtering outliers: keeping indices {}..{} of {}", lo, hi, sorted.len());
    sorted[lo..hi].to_vec()
}

/// Fit NB mixture model at given ploidy via Nelder-Mead
/// Following floco's estimate_mean_std_at_ploidy() approach
/// Returns (negative_log_likelihood, alpha, beta)
fn fit_mixture_at_ploidy(counts: &[u32], bp_step: usize, bin_size: usize, ploidy: u32) -> (f64, f64, f64) {
    // N_CN = max(4, ploidy + 3) following floco
    let n_cn = (ploidy + 3).max(4) as usize;

    // Find mode of histogram
    let mode_idx = counts.iter()
        .enumerate()
        .max_by_key(|(_, &c)| c)
        .map(|(i, _)| i)
        .unwrap_or(1);
    let mode = (mode_idx as f64 + 0.5) * bp_step as f64;

    // Initial estimates (following floco)
    let m0 = mode / ploidy as f64;
    let v0 = m0 * m0 / 25.0;

    // Initial CN proportions: 0.1 for all except ploidy which gets 0.9
    let mut init_props: Vec<f64> = vec![0.1; n_cn];
    if (ploidy as usize) < n_cn {
        init_props[ploidy as usize] = 0.9;
    }

    debug!("Ploidy {}: mode={:.0}, init_mean={:.2}, init_var={:.2}", ploidy, mode, m0, v0);

    // Optimize using Nelder-Mead in (2 + n_cn) dimensions
    // Parameters: [mean, variance, prop_0, prop_1, ..., prop_{n_cn-1}]
    let result = nelder_mead_mixture(
        m0, v0, &init_props, ploidy as usize, counts, bp_step,
    );

    let (opt_mean, opt_var, neg_ll) = result;
    // Convert to per-bp (divide by bin_size, following floco)
    let alpha = opt_mean / bin_size as f64;
    let beta = opt_var.sqrt() / bin_size as f64;

    debug!("Ploidy {}: optimized mean={:.2}, var={:.2}", ploidy, opt_mean, opt_var);

    (neg_ll, alpha, beta)
}

/// Nelder-Mead optimization for mixture model
/// Returns (optimized_mean, optimized_variance, final_negative_log_likelihood)
fn nelder_mead_mixture(
    m0: f64,
    v0: f64,
    init_props: &[f64],
    ploidy: usize,
    counts: &[u32],
    bp_step: usize,
) -> (f64, f64, f64) {
    let n_cn = init_props.len();
    let n_params = 2 + n_cn;  // mean, variance, + proportions

    // Build initial simplex
    let mut simplex: Vec<(Vec<f64>, f64)> = Vec::with_capacity(n_params + 1);

    // Initial point
    let mut x0 = vec![m0, v0];
    x0.extend_from_slice(init_props);
    let f0 = mixture_neg_log_likelihood(&x0, counts, bp_step);
    simplex.push((x0.clone(), f0));

    // Perturbed points
    for i in 0..n_params {
        let mut xi = x0.clone();
        if i < 2 {
            xi[i] *= 1.1;  // Perturb mean/variance
        } else {
            // Perturb proportions slightly
            xi[i] = (xi[i] + 0.05).min(0.95);
        }
        let fi = mixture_neg_log_likelihood(&xi, counts, bp_step);
        simplex.push((xi, fi));
    }

    let alpha = 1.0;  // Reflection
    let gamma = 2.0;  // Expansion
    let rho = 0.5;    // Contraction
    let sigma = 0.5;  // Shrink

    for iter in 0..NELDER_MEAD_ITERS * 3 {  // More iterations for high-dim
        // Sort by function value
        simplex.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));

        let f_best = simplex[0].1;
        let f_worst = simplex[n_params].1;
        let f_second = simplex[n_params - 1].1;

        // Check convergence
        if (f_worst - f_best).abs() < 1e-6 {
            debug!("Nelder-Mead mixture converged at iteration {}", iter);
            break;
        }

        // Centroid of all but worst
        let mut centroid = vec![0.0; n_params];
        for i in 0..n_params {
            for j in 0..n_params {
                centroid[j] += simplex[i].0[j];
            }
        }
        for c in &mut centroid {
            *c /= n_params as f64;
        }

        // Reflection
        let x_r: Vec<f64> = centroid.iter()
            .zip(&simplex[n_params].0)
            .map(|(&c, &w)| c + alpha * (c - w))
            .collect();
        let x_r = clamp_params(&x_r, ploidy);
        let f_r = mixture_neg_log_likelihood(&x_r, counts, bp_step);

        if f_r < f_second && f_r >= f_best {
            simplex[n_params] = (x_r, f_r);
            continue;
        }

        if f_r < f_best {
            // Expansion
            let x_e: Vec<f64> = centroid.iter()
                .zip(&x_r)
                .map(|(&c, &r)| c + gamma * (r - c))
                .collect();
            let x_e = clamp_params(&x_e, ploidy);
            let f_e = mixture_neg_log_likelihood(&x_e, counts, bp_step);

            simplex[n_params] = if f_e < f_r { (x_e, f_e) } else { (x_r, f_r) };
            continue;
        }

        // Contraction
        let x_c: Vec<f64> = centroid.iter()
            .zip(&simplex[n_params].0)
            .map(|(&c, &w)| c + rho * (w - c))
            .collect();
        let x_c = clamp_params(&x_c, ploidy);
        let f_c = mixture_neg_log_likelihood(&x_c, counts, bp_step);

        if f_c < f_worst {
            simplex[n_params] = (x_c, f_c);
            continue;
        }

        // Shrink
        let best = simplex[0].0.clone();
        for i in 1..=n_params {
            let x_new: Vec<f64> = best.iter()
                .zip(&simplex[i].0)
                .map(|(&b, &x)| b + sigma * (x - b))
                .collect();
            let x_new = clamp_params(&x_new, ploidy);
            let f_new = mixture_neg_log_likelihood(&x_new, counts, bp_step);
            simplex[i] = (x_new, f_new);
        }
    }

    simplex.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
    let best = &simplex[0];
    (best.0[0], best.0[1], best.1)
}

/// Clamp parameters to valid ranges (following floco's bounds)
fn clamp_params(params: &[f64], ploidy: usize) -> Vec<f64> {
    let mut clamped = params.to_vec();

    // Mean: positive
    clamped[0] = clamped[0].max(1e-6);

    // Variance: must be > mean for NB
    clamped[1] = clamped[1].max(clamped[0] * 1.00001);

    // Proportions: bounded
    for i in 2..params.len() {
        let cn_idx = i - 2;
        if cn_idx == ploidy {
            // Ploidy proportion: [0.5, 1.0]
            clamped[i] = clamped[i].clamp(0.5, 1.0);
        } else {
            // Other proportions: [0.0, 0.45]
            clamped[i] = clamped[i].clamp(0.0, 0.45);
        }
    }

    clamped
}

/// Negative log-likelihood for NB mixture model
/// Following floco's MLE_NBinom function
fn mixture_neg_log_likelihood(
    params: &[f64],
    counts: &[u32],
    bp_step: usize,
) -> f64 {
    let m = params[0];
    let v = params[1].max(m * 1.00001);
    let props = &params[2..];
    let n_cn = props.len();

    // NB parameters
    let r = m.powi(2) / (v - m);
    let p = m / v;

    if r <= 0.0 || p <= 0.0 || p >= 1.0 {
        return 1e30;
    }

    // Normalize proportions to sum to 1 (log space for numerical stability)
    let prop_sum: f64 = props.iter().sum();
    if prop_sum <= 0.0 {
        return 1e30;
    }
    let log_props: Vec<f64> = props.iter().map(|&c| (c / prop_sum).ln()).collect();

    // CN-scaled r values: r_cn = max(cn, epsilon) * r
    let rs: Vec<f64> = (0..n_cn)
        .map(|cn| (cn as f64).max(EPSILON_MIX) * r)
        .collect();

    // Sum over histogram bins
    let mut total_ll = 0.0;
    for (i, &count) in counts.iter().enumerate() {
        if count == 0 {
            continue;
        }

        // Coverage value at bin center
        let cov = ((i as f64) + 0.5) * bp_step as f64;
        let k = cov.round() as u64;

        // Log-sum-exp over CN components
        let mut log_probs: Vec<f64> = Vec::with_capacity(n_cn);
        for cn in 0..n_cn {
            let r_cn = rs[cn];
            if r_cn <= 0.0 {
                log_probs.push(f64::NEG_INFINITY);
                continue;
            }

            match NegativeBinomial::new(r_cn, p) {
                Ok(nb) => {
                    let lp = log_props[cn] + nb.ln_pmf(k);
                    log_probs.push(lp);
                }
                Err(_) => {
                    log_probs.push(f64::NEG_INFINITY);
                }
            }
        }

        let lse = log_sum_exp(&log_probs);
        total_ll += count as f64 * lse;
    }

    // Return negative log-likelihood (we minimize)
    (-total_ll).min(1e30)
}

/// 3D Nelder-Mead optimization for skew-normal fitting
/// Minimizes f(x, y, z) starting from (x0, y0, z0)
fn nelder_mead_3d(x0: f64, y0: f64, z0: f64, f: impl Fn(f64, f64, f64) -> f64) -> (f64, f64, f64) {
    // Initialize 4-point simplex in 3D
    let mut simplex = [
        (x0, y0, z0, f(x0, y0, z0)),
        (x0 * 1.1, y0, z0, f(x0 * 1.1, y0, z0)),
        (x0, y0 * 1.1, z0, f(x0, y0 * 1.1, z0)),
        (x0, y0, z0 * 1.1 + 0.1, f(x0, y0, z0 * 1.1 + 0.1)),  // +0.1 in case z0=0
    ];

    let alpha = 1.0;  // Reflection
    let gamma = 2.0;  // Expansion
    let rho = 0.5;    // Contraction
    let sigma = 0.5;  // Shrink

    for iter in 0..NELDER_MEAD_ITERS * 2 {  // More iterations for 3D
        // Sort by function value
        simplex.sort_by(|a, b| a.3.partial_cmp(&b.3).unwrap_or(std::cmp::Ordering::Equal));

        let f_best = simplex[0].3;
        let f_worst = simplex[3].3;
        let f_second = simplex[2].3;

        // Check convergence
        if (f_worst - f_best).abs() < 1e-8 {
            debug!("Nelder-Mead 3D converged at iteration {}", iter);
            break;
        }

        // Centroid of best 3 points
        let x_c = (simplex[0].0 + simplex[1].0 + simplex[2].0) / 3.0;
        let y_c = (simplex[0].1 + simplex[1].1 + simplex[2].1) / 3.0;
        let z_c = (simplex[0].2 + simplex[1].2 + simplex[2].2) / 3.0;

        // Reflection
        let x_r = x_c + alpha * (x_c - simplex[3].0);
        let y_r = y_c + alpha * (y_c - simplex[3].1);
        let z_r = z_c + alpha * (z_c - simplex[3].2);
        let f_r = f(x_r.max(1e-10), y_r.max(1e-10), z_r);

        if f_r < f_second && f_r >= f_best {
            simplex[3] = (x_r, y_r, z_r, f_r);
            continue;
        }

        if f_r < f_best {
            // Expansion
            let x_e = x_c + gamma * (x_r - x_c);
            let y_e = y_c + gamma * (y_r - y_c);
            let z_e = z_c + gamma * (z_r - z_c);
            let f_e = f(x_e.max(1e-10), y_e.max(1e-10), z_e);

            simplex[3] = if f_e < f_r { (x_e, y_e, z_e, f_e) } else { (x_r, y_r, z_r, f_r) };
            continue;
        }

        // Contraction
        let x_con = x_c + rho * (simplex[3].0 - x_c);
        let y_con = y_c + rho * (simplex[3].1 - y_c);
        let z_con = z_c + rho * (simplex[3].2 - z_c);
        let f_con = f(x_con.max(1e-10), y_con.max(1e-10), z_con);

        if f_con < f_worst {
            simplex[3] = (x_con, y_con, z_con, f_con);
            continue;
        }

        // Shrink
        for i in 1..4 {
            let x_new = simplex[0].0 + sigma * (simplex[i].0 - simplex[0].0);
            let y_new = simplex[0].1 + sigma * (simplex[i].1 - simplex[0].1);
            let z_new = simplex[0].2 + sigma * (simplex[i].2 - simplex[0].2);
            simplex[i] = (x_new, y_new, z_new, f(x_new.max(1e-10), y_new.max(1e-10), z_new));
        }
    }

    simplex.sort_by(|a, b| a.3.partial_cmp(&b.3).unwrap_or(std::cmp::Ordering::Equal));
    (simplex[0].0, simplex[0].1, simplex[0].2)
}

// ─── Read Length Distribution Fitting ────────────────────────────────────────

use crate::ReadLengthParams;

/// Simple deterministic RNG (LCG) for reproducible subsampling
struct SimpleRng(u64);

impl SimpleRng {
    fn new(seed: u64) -> Self {
        Self(seed)
    }

    fn next(&mut self) -> usize {
        // Linear congruential generator
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1);
        (self.0 >> 33) as usize
    }
}

/// Fit skew-normal distribution to read lengths via MLE
/// Subsamples to READ_LENGTH_SUBSAMPLE if too many reads (like floco)
/// Returns fitted parameters (shape, loc, scale)
pub fn fit_read_length_distribution(read_lengths: &[u32]) -> ReadLengthParams {
    if read_lengths.is_empty() {
        warn!("No read lengths for distribution fitting, using defaults");
        return ReadLengthParams::new(0.0, 1000.0, 500.0);
    }

    // Subsample if too many reads (like floco's subset_size = 100000)
    let data: Vec<f64> = if read_lengths.len() > READ_LENGTH_SUBSAMPLE {
        use std::collections::HashSet;
        debug!("Subsampling {} reads to {} for distribution fitting",
               read_lengths.len(), READ_LENGTH_SUBSAMPLE);
        let mut rng = SimpleRng::new(42);  // Fixed seed for reproducibility
        let mut indices: HashSet<usize> = HashSet::with_capacity(READ_LENGTH_SUBSAMPLE);
        while indices.len() < READ_LENGTH_SUBSAMPLE {
            indices.insert(rng.next() % read_lengths.len());
        }
        indices.iter().map(|&i| read_lengths[i] as f64).collect()
    } else {
        read_lengths.iter().map(|&x| x as f64).collect()
    };

    let n = data.len() as f64;

    // Initial estimates from data moments
    let mean: f64 = data.iter().sum::<f64>() / n;
    let variance: f64 = data.iter().map(|&x| (x - mean).powi(2)).sum::<f64>() / n;
    let std_dev = variance.sqrt().max(1.0);

    // Skewness estimate for initial shape
    let skewness: f64 = data.iter()
        .map(|&x| ((x - mean) / std_dev).powi(3))
        .sum::<f64>() / n;

    // Initial shape from skewness (rough approximation)
    let init_shape = skewness.signum() * (skewness.abs() * 2.0).min(5.0);

    debug!("Read length stats: mean={:.0}, std={:.0}, skew={:.2}", mean, std_dev, skewness);

    // Fit via Nelder-Mead (minimize negative log-likelihood)
    let (opt_loc, opt_scale, opt_shape) = nelder_mead_3d(
        mean,
        std_dev,
        init_shape,
        |loc, scale, shape| -skew_normal_log_likelihood(&data, shape, loc, scale),
    );

    info!("Read length distribution fit: loc={:.0}, scale={:.0}, shape={:.2}",
          opt_loc, opt_scale, opt_shape);

    ReadLengthParams::new(opt_shape, opt_loc, opt_scale)
}

/// Log-likelihood of data under skew-normal distribution
fn skew_normal_log_likelihood(data: &[f64], shape: f64, loc: f64, scale: f64) -> f64 {
    if scale <= 0.0 {
        return f64::NEG_INFINITY;
    }

    const LN_2: f64 = std::f64::consts::LN_2;
    const LN_SQRT_2PI: f64 = 0.9189385332046727; // ln(sqrt(2*pi))

    let mut ll = data.len() as f64 * (LN_2 - LN_SQRT_2PI - scale.ln());

    for &x in data {
        let z = (x - loc) / scale;
        // -0.5 * z^2 for normal PDF part
        ll -= 0.5 * z * z;
        // ln(Phi(alpha * z)) for skew-normal adjustment
        ll += std_normal_cdf_ln(shape * z);
    }

    ll
}

/// Log of standard normal CDF (more numerically stable)
fn std_normal_cdf_ln(x: f64) -> f64 {
    // Use log1p for better precision near 0
    // Phi(x) = 0.5 * (1 + erf(x / sqrt(2)))
    let erf_val = erf(x * std::f64::consts::FRAC_1_SQRT_2);
    if erf_val > 0.0 {
        (0.5 * (1.0 + erf_val)).ln()
    } else {
        // For negative values, use: ln(Phi(x)) = ln(0.5) + ln(1 + erf(...))
        // but this can underflow. Use asymptotic expansion for very negative x
        if x < -6.0 {
            // Asymptotic: ln(Phi(x)) ≈ -x^2/2 - ln(-x) - ln(sqrt(2*pi))
            -0.5 * x * x - (-x).ln() - 0.9189385332046727
        } else {
            (0.5 * (1.0 + erf_val)).max(1e-300).ln()
        }
    }
}

/// Error function approximation (Abramowitz and Stegun)
fn erf(x: f64) -> f64 {
    let a1 =  0.254829592;
    let a2 = -0.284496736;
    let a3 =  1.421413741;
    let a4 = -1.453152027;
    let a5 =  1.061405429;
    let p  =  0.3275911;

    let sign = if x < 0.0 { -1.0 } else { 1.0 };
    let x = x.abs();

    let t = 1.0 / (1.0 + p * x);
    let y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * (-x * x).exp();

    sign * y
}

// ─── Edge Coverage Penalty ───────────────────────────────────────────────────

/// Calculate penalty for edges with insufficient supporting reads
///
/// Following floco's edge_cov_pen() logic:
/// - Expected reads = alpha * P(read_len > overlap) / 4
/// - If observed < expected.floor(), apply penalty
///
/// # Arguments
/// * `sup_reads` - Number of reads supporting (crossing) this edge
/// * `alpha` - Coverage per bp at CN=1
/// * `overlap` - Edge overlap in bp
/// * `rlen_params` - Read length distribution parameters
/// * `penalty` - Penalty value to apply (e.g., -2.0 for cheap_source)
///
/// # Returns
/// Penalty value (negative) if insufficient support, 0.0 otherwise
pub fn edge_cov_penalty(
    sup_reads: u32,
    alpha: f64,
    overlap: usize,
    rlen_params: &ReadLengthParams,
    penalty: f64,
) -> f64 {
    // Expected reads supporting edge = alpha * P(read_len > overlap) / 4
    // The /4 comes from: only ~1/4 of reads crossing a node will span the overlap
    let prob_span_overlap = rlen_params.sf(overlap as f64);
    let expected = alpha * prob_span_overlap / 4.0;

    // If we have fewer reads than expected (floor), apply penalty
    // Note: floco uses floor without max(1.0), so edges with 0 expected can have 0 support without penalty
    if (sup_reads as f64) < expected.floor() {
        penalty
    } else {
        0.0
    }
}

// ─── CN Probability Calculation ──────────────────────────────────────────────

/// Convert calibration params (alpha, beta) to NB distribution params (r, p, p0)
/// p0 is the rate for exponential distribution used for CN=0
#[inline]
fn to_nb_params(len: f64, alpha: f64, beta: f64, eps: f64) -> (f64, f64, f64) {
    let mu = alpha * len;
    let v = (beta * len).powi(2).max(mu + 1e-6);
    let r = mu.powi(2) / (v - mu);
    let p = mu / v;
    let p0 = p / (r * eps * (1.0 - p));
    (r, p, p0)
}

/// Compute log P(CN=k) for candidate CN values (aggregate coverage version)
/// Returns (lower_bound, normalized_log_probs)
///
/// # Arguments
/// * `diff_cutoff` - Stop searching when log-prob drops this much from max
///                   (floco uses 4 * |source_prob|, e.g., 4 * 20 = 80)
pub fn cn_log_probs(cov: f64, len: usize, alpha: f64, beta: f64, eps: f64, diff_cutoff: f64) -> (i32, Vec<f64>) {
    if len == 0 || alpha <= 0.0 {
        debug!("Degenerate case: len={}, alpha={}", len, alpha);
        return (0, vec![0.0]); // Degenerate case
    }

    let (r, p, p0) = to_nb_params(len as f64, alpha, beta, eps);

    // Starting estimate
    let start = ((cov / (len as f64 * alpha)).round() as i32).max(0);

    let mut probs = Vec::with_capacity(16);
    let mut lo = start;
    let mut max_prob = f64::NEG_INFINITY;

    // Search upward from start (following floco's MIN_EXTENSION logic)
    for cn in start..(start + MAX_CN_SEARCH) {
        let lp = log_prob_cn(cn, cov, r, p, p0);
        // Only stop if we've extended at least MIN_EXTENSION steps and prob is low
        if (cn - start) > MIN_EXTENSION && lp + diff_cutoff < max_prob {
            break;
        }
        max_prob = max_prob.max(lp);
        probs.push(lp);
    }

    // Search downward from start (following floco's MIN_EXTENSION logic)
    for cn in (0..start).rev() {
        let lp = log_prob_cn(cn, cov, r, p, p0);
        // Only stop if we've extended at least MIN_EXTENSION steps and prob is low
        if (start - cn) > MIN_EXTENSION && lp + diff_cutoff < max_prob {
            break;
        }
        max_prob = max_prob.max(lp);
        probs.insert(0, lp);
        lo = cn;
    }

    // Normalize via log-sum-exp
    let lse = log_sum_exp(&probs);
    for prob in &mut probs {
        *prob -= lse;
    }

    (lo, probs)
}

/// Compute log P(CN=k) using per-bin coverages (floco-compatible version)
/// Following floco's counts_to_probabs.py:19-70 exactly:
/// - NB params based on bin_size (not node length)
/// - CN=0: exponential distribution sum over bins
/// - CN>=1: sum of NB logpmf over ALL bins
///
/// # Arguments
/// * `bin_coverages` - Vector of per-bin coverage values
/// * `bin_size` - Size of each bin in bp
/// * `alpha` - Coverage per bp at CN=1
/// * `beta` - Std dev per bp
/// * `eps` - Epsilon for CN=0 sensitivity
/// * `diff_cutoff` - Stop searching when log-prob drops this much from max
///
/// # Returns
/// (lower_bound, normalized_log_probs)
pub fn cn_log_probs_bins(
    bin_coverages: &[f64],
    bin_size: usize,
    alpha: f64,
    beta: f64,
    eps: f64,
    diff_cutoff: f64,
) -> (i32, Vec<f64>) {
    let n_bins = bin_coverages.len();

    if n_bins == 0 || bin_size == 0 || alpha <= 0.0 {
        debug!("Degenerate case: n_bins={}, bin_size={}, alpha={}", n_bins, bin_size, alpha);
        return (0, vec![0.0]);
    }

    // NB params based on bin_size (NOT node length) - following floco
    let (r, p, p0) = to_nb_params(bin_size as f64, alpha, beta, eps);

    // Total coverage across all bins (for starting estimate and CN=0)
    let total_cov: f64 = bin_coverages.iter().sum();
    let mean_cov = total_cov / n_bins as f64;

    // Starting estimate based on mean bin coverage
    let start = ((mean_cov / (bin_size as f64 * alpha)).round() as i32).max(0);

    let mut probs = Vec::with_capacity(16);
    let mut lo = start;
    let mut max_prob = f64::NEG_INFINITY;

    // Search upward from start
    for cn in start..(start + MAX_CN_SEARCH) {
        let lp = log_prob_cn_bins(cn, bin_coverages, r, p, p0);
        if (cn - start) > MIN_EXTENSION && lp + diff_cutoff < max_prob {
            break;
        }
        max_prob = max_prob.max(lp);
        probs.push(lp);
    }

    // Search downward from start
    for cn in (0..start).rev() {
        let lp = log_prob_cn_bins(cn, bin_coverages, r, p, p0);
        if (start - cn) > MIN_EXTENSION && lp + diff_cutoff < max_prob {
            break;
        }
        max_prob = max_prob.max(lp);
        probs.insert(0, lp);
        lo = cn;
    }

    // Normalize via log-sum-exp
    let lse = log_sum_exp(&probs);
    for prob in &mut probs {
        *prob -= lse;
    }

    (lo, probs)
}

/// Log probability of specific CN value using per-bin coverages
/// Following floco's counts_to_probabs.py exactly
#[inline]
fn log_prob_cn_bins(cn: i32, bin_coverages: &[f64], r: f64, p: f64, p0: f64) -> f64 {
    let n_bins = bin_coverages.len();

    if cn < 1 {
        // CN=0: Exponential distribution - floco's formula:
        // -p0 * sum(bins) + n_bins * log(1 - exp(-p0))
        let total_cov: f64 = bin_coverages.iter().sum();
        -p0 * total_cov + n_bins as f64 * (1.0 - (-p0).exp()).ln()
    } else {
        // CN>=1: Sum of NB logpmf over ALL bins
        let r_cn = r * cn as f64;
        if r_cn <= 0.0 || p <= 0.0 || p >= 1.0 {
            return f64::NEG_INFINITY;
        }

        match NegativeBinomial::new(r_cn, p) {
            Ok(nb) => {
                bin_coverages.iter()
                    .map(|&cov| nb.ln_pmf(cov.round() as u64))
                    .sum()
            }
            Err(_) => f64::NEG_INFINITY,
        }
    }
}

/// Log probability of specific CN value
#[inline]
fn log_prob_cn(cn: i32, cov: f64, r: f64, p: f64, p0: f64) -> f64 {
    if cn < 1 {
        // Exponential distribution for CN=0
        -p0 * cov + (1.0 - (-p0).exp()).ln()
    } else {
        // Negative binomial for CN >= 1
        let r_cn = r * cn as f64;
        if r_cn <= 0.0 || p <= 0.0 || p >= 1.0 {
            return f64::NEG_INFINITY;
        }
        NegativeBinomial::new(r_cn, p)
            .map(|nb| nb.ln_pmf(cov.round() as u64))
            .unwrap_or(f64::NEG_INFINITY)
    }
}

/// Numerically stable log-sum-exp
#[inline]
fn log_sum_exp(xs: &[f64]) -> f64 {
    if xs.is_empty() {
        return f64::NEG_INFINITY;
    }
    let max = xs.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    if max.is_infinite() {
        return max;
    }
    max + xs.iter().map(|&x| (x - max).exp()).sum::<f64>().ln()
}

/// Simple CN call: argmax P(CN=k)
pub fn simple_cn_call(cov: f64, len: usize, alpha: f64, beta: f64, eps: f64, diff_cutoff: f64) -> u32 {
    let (lo, probs) = cn_log_probs(cov, len, alpha, beta, eps, diff_cutoff);

    probs
        .iter()
        .enumerate()
        .max_by(|a, b| a.1.partial_cmp(b.1).unwrap_or(std::cmp::Ordering::Equal))
        .map(|(i, _)| (lo + i as i32).max(0) as u32)
        .unwrap_or(0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_filter_outliers() {
        let data: Vec<f64> = (0..100).map(|x| x as f64).collect();
        let filtered = filter_outliers(&data, 0.1);
        assert!(filtered.len() < data.len());
        assert!(!filtered.contains(&0.0));
        assert!(!filtered.contains(&99.0));
    }

    #[test]
    fn test_log_sum_exp() {
        let xs = vec![1.0, 2.0, 3.0];
        let result = log_sum_exp(&xs);
        let expected = (xs.iter().map(|x| x.exp()).sum::<f64>()).ln();
        assert!((result - expected).abs() < 1e-10);
    }

    #[test]
    fn test_cn_log_probs_basic() {
        let (lo, probs) = cn_log_probs(100.0, 100, 0.5, 0.2, 0.02, 80.0);
        assert!(lo >= 0);
        assert!(!probs.is_empty());
        // Probabilities should sum to ~1 (in log space, max should be close to 0)
        let max_prob = probs.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        assert!(max_prob > -10.0); // At least one probable CN
    }

    #[test]
    fn test_cn_log_probs_bins_basic() {
        // Test bin-level CN probability calculation
        // 10 bins, each with coverage ~100 (total 1000)
        let bin_coverages = vec![100.0; 10];
        let bin_size = 100;
        let (lo, probs) = cn_log_probs_bins(&bin_coverages, bin_size, 0.5, 0.2, 0.02, 40000.0);
        assert!(lo >= 0);
        assert!(!probs.is_empty());
        // Probabilities should be normalized
        let max_prob = probs.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        assert!(max_prob > -10.0);
    }

    #[test]
    fn test_cn_log_probs_bins_empty() {
        // Empty bins should return degenerate case
        let bin_coverages: Vec<f64> = vec![];
        let (lo, probs) = cn_log_probs_bins(&bin_coverages, 100, 0.5, 0.2, 0.02, 40000.0);
        assert_eq!(lo, 0);
        assert_eq!(probs.len(), 1);
    }
}
