use flate2::read::GzDecoder;
use std::collections::HashMap;
use std::fs::File;
use std::io::{prelude::*, BufReader};
use std::path::Path;

/// Parse GFA file and extract segment information
/// Returns (segment_lengths, min_id) where segment_lengths[id - min_id] gives the length
pub fn parse_gfa(gfa_path: &str) -> std::io::Result<(Vec<usize>, usize)> {
    let path = Path::new(gfa_path);
    let mut reader = create_reader(path)?;
    let mut line = String::new();
    let mut segments_map = HashMap::new();
    let mut min_id = usize::MAX;
    let mut max_id = 0;

    loop {
        line.clear();
        let bytes_read = reader.read_line(&mut line)?;
        if bytes_read == 0 {
            break;
        }

        let line_str = line.trim();
        if !line_str.starts_with('S') {
            continue;
        }

        let mut fields = line_str.split('\t');
        let Some((id_str, seq)) = fields.next().and_then(|_type| {
            let id_str = fields.next()?;
            let seq = fields.next()?;
            Some((id_str, seq))
        }) else {
            continue;
        };

        let id = id_str.parse::<usize>().unwrap();
        min_id = min_id.min(id);
        max_id = max_id.max(id);
        segments_map.insert(id, seq.len());
    }

    let num_segments = max_id - min_id + 1;
    let mut segment_lengths = vec![0; num_segments];

    for (id, len) in segments_map {
        segment_lengths[id - min_id] = len;
    }

    Ok((segment_lengths, min_id))
}

/// Create a reader that handles compressed files
pub fn create_reader(path: &Path) -> std::io::Result<Box<dyn BufRead>> {
    let file = File::open(path)?;
    if path
        .extension()
        .is_some_and(|ext| ext == "gz" || ext == "bgz")
    {
        let decoder = GzDecoder::new(file);
        let buf_reader = BufReader::new(decoder);
        Ok(Box::new(buf_reader))
    } else {
        let buf_reader = BufReader::new(file);
        Ok(Box::new(buf_reader))
    }
}

/// Process each step in a GAF alignment line, calculating coverage for graph nodes
pub fn for_each_step(
    line: &str,
    mut callback: impl FnMut(usize, usize),
    mut get_node_len: impl FnMut(usize) -> usize,
) {
    let walk = line.split('\t').nth(5).unwrap_or("*");
    if walk == "*" {
        return;
    }

    let target_start = line
        .split('\t')
        .nth(7)
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(0);
    let target_end = line
        .split('\t')
        .nth(8)
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(0);
    let target_len = target_end - target_start;

    let fields = line
        .split('\t')
        .nth(5)
        .unwrap()
        .split(['<', '>'])
        .filter(|s| !s.is_empty())
        .filter_map(|s| s.parse::<usize>().ok())
        .enumerate()
        .collect::<Vec<(usize, usize)>>();

    let mut seen: usize = 0;
    let fields_len = fields.len();

    for (i, j) in fields {
        let mut len = get_node_len(j);
        if i == 0 {
            assert!(len >= target_start);
            len -= target_start;
        }
        if i == fields_len - 1 {
            assert!(target_len >= seen);
            len = target_len - seen;
        }
        seen += len;
        callback(j, len);
    }
}

/// Compute coverage from GAF file on a GFA graph
/// Returns coverage vector indexed by (node_id - min_id)
pub fn compute_coverage(
    gfa_path: &str,
    gaf_path: &str,
    len_scale: bool,
    weight_queries: bool,
) -> std::io::Result<(Vec<f64>, usize)> {
    let (segment_lengths, min_id) = parse_gfa(gfa_path)?;
    let num_segments = segment_lengths.len();
    let mut coverage: Vec<f64> = vec![0.0; num_segments];

    let for_each_line = |callback: &mut dyn FnMut(&str)| -> std::io::Result<()> {
        let file = File::open(gaf_path)?;
        let (reader, _compression) = niffler::get_reader(Box::new(file)).unwrap();
        let buf_reader = BufReader::new(reader);
        for line in buf_reader.lines() {
            callback(&line?);
        }
        Ok(())
    };

    if weight_queries {
        let mut query_counts: HashMap<String, usize> = HashMap::new();
        for_each_line(&mut |l: &str| {
            let fields: Vec<&str> = l.split('\t').collect();
            if fields.len() >= 4 {
                let query_key = format!("{}:{}:{}", fields[0], fields[2], fields[3]);
                *query_counts.entry(query_key).or_insert(0) += 1;
            }
        })?;

        for_each_line(&mut |l: &str| {
            let fields: Vec<&str> = l.split('\t').collect();
            let query_key = format!("{}:{}:{}", fields[0], fields[2], fields[3]);
            let count = query_counts.get(&query_key).unwrap_or(&1);

            for_each_step(
                l,
                |node_id, len| {
                    coverage[node_id - min_id] += len as f64 / *count as f64;
                },
                |node_id| segment_lengths[node_id - min_id],
            );
        })?;
    } else {
        for_each_line(&mut |l: &str| {
            for_each_step(
                l,
                |node_id, len| {
                    coverage[node_id - min_id] += len as f64;
                },
                |node_id| segment_lengths[node_id - min_id],
            );
        })?;
    }

    if len_scale {
        for (i, cov) in coverage.iter_mut().enumerate() {
            if segment_lengths[i] > 0 {
                *cov /= segment_lengths[i] as f64;
            }
        }
    }

    Ok((coverage, min_id))
}

/// Format coverage as column output (PACK format)
pub fn format_coverage_column(sample_name: &str, coverage: &[f64]) -> String {
    let mut output = String::new();
    output.push_str(&format!("##sample: {}\n", sample_name));
    output.push_str("#coverage\n");
    for &v in coverage {
        output.push_str(&format!("{}\n", v));
    }
    output
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_gfa() {
        // Test would need a sample GFA file
    }
}
