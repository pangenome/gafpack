use clap::Parser;
use flate2::read::GzDecoder;
use std::collections::HashMap;
use std::fs::File;
use std::io::{prelude::*, BufReader};
use std::path::Path;

/// Iterates through each line in a file, applying the provided callback function
///
/// # Arguments
/// * `filename` - Path to the file to read
/// * `callback` - Function to call for each line
fn for_each_line_in_file(filename: &str, mut callback: impl FnMut(&str)) {
    let file = File::open(filename).unwrap();
    let (reader, _compression) = niffler::get_reader(Box::new(file)).unwrap();
    let buf_reader = BufReader::new(reader);
    for line in buf_reader.lines() {
        callback(&line.unwrap());
    }
}

/// Process each step in a GAF alignment line, calculating coverage for graph nodes
///
/// # Arguments
/// * `line` - A GAF format alignment line
/// * `callback` - Function called for each node with (node_id, coverage_length)
/// * `get_node_len` - Function to get the length of a node by its ID
///
/// # Details
/// Parses GAF alignment lines to extract node coverage information:
/// - Handles both forward (>) and reverse (<) node traversals
/// - Adjusts coverage for partial node alignments at path ends
/// - Accumulates coverage across multi-node paths
fn for_each_step(
    line: &str,
    mut callback: impl FnMut(usize, usize),
    mut get_node_len: impl FnMut(usize) -> usize,
) {
    //eprintln!("{}", line);
    let walk = line.split('\t').nth(5).unwrap();
    if walk != "*" {
        //eprintln!("oheunotoeunthoue");
        let target_start = line.split('\t').nth(7).unwrap().parse::<usize>().unwrap();
        let target_end = line.split('\t').nth(8).unwrap().parse::<usize>().unwrap();
        let target_len = target_end - target_start;
        //eprintln!("target_len = {}", target_len);
        let fields = line
            .split('\t')
            .nth(5)
            .unwrap()
            .split(|c| c == '<' || c == '>')
            .filter(|s| !s.is_empty())
            .map(|s| s.parse::<usize>().unwrap())
            .enumerate()
            .collect::<Vec<(usize, usize)>>();
        let mut seen: usize = 0;
        let fields_len = fields.as_slice().len();
        //eprintln!("fields len = {}", fields_len);
        for (i, j) in fields {
            let mut len = get_node_len(j);
            //eprintln!("node {} len = {}", j, len);
            if i == 0 {
                //eprintln!("on first step {} {} {}", len, target_start, seen);
                assert!(len >= target_start);
                len -= target_start;
            }
            if i == fields_len - 1 {
                //eprintln!("on last step {} {} {}", len, target_end, seen);
                assert!(target_len >= seen);
                len = target_len - seen;
            }
            if i == fields_len {
                assert!(false);
            }
            //eprintln!("node {} adj len = {}", j, len);
            seen += len;
            callback(j, len);
        }
        //eprintln!("seen = {}", seen);
    }
    //eprintln!("at end");
}

/// Create a reader that handles compressed files
fn create_reader(path: &Path) -> std::io::Result<Box<dyn BufRead>> {
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

/// Simple struct to hold segment information
struct Segment {
    id: usize,
    sequence: String,
}

/// Parse GFA file and extract segments
fn parse_gfa(gfa_path: &str) -> std::io::Result<Vec<Segment>> {
    let path = Path::new(gfa_path);
    let mut reader = create_reader(path)?;
    let mut line = String::new();
    let mut segments = Vec::new();
    let mut segment_map = HashMap::new();

    loop {
        line.clear();
        let bytes_read = reader.read_line(&mut line)?;
        if bytes_read == 0 {
            break;
        }

        let line_str = line.trim();

        // Only process segment lines
        if !line_str.starts_with('S') {
            continue;
        }

        // Parse segment line format: S<tab>id<tab>sequence
        let mut fields = line_str.split('\t');
        let Some((id_str, seq)) = fields.next().and_then(|_type| {
            let id_str = fields.next()?;
            let seq = fields.next()?;
            Some((id_str, seq))
        }) else {
            continue;
        };

        // Parse segment ID
        let id = id_str.parse::<usize>().unwrap();
        segment_map.insert(id, segments.len());
        segments.push(Segment {
            id,
            sequence: seq.to_string(),
        });
    }

    // Sort segments by ID to ensure they're in order
    segments.sort_by_key(|s| s.id);

    Ok(segments)
}

/// Project a GAF alignment file into coverage over GFA graph nodes
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Input GFA pangenome graph file
    #[arg(long)]
    gfa: String,
    /// Input GAF alignment file
    #[arg(short, long)]
    gaf: String,
    /// Scale coverage values by node length
    #[arg(short, long)]
    len_scale: bool,
    /// Emit graph coverage vector in a single column
    #[arg(short, long)]
    coverage_column: bool,
    /// Weight coverage by query group occurrences
    #[arg(short = 'w', long)]
    weight_queries: bool,
}

fn main() {
    let args = Args::parse();

    // Parse GFA file
    let segments = parse_gfa(&args.gfa).unwrap();
    let num_segments = segments.len();

    // Create a map from segment ID to index for fast lookup
    let segment_id_to_index: HashMap<usize, usize> = segments
        .iter()
        .enumerate()
        .map(|(idx, seg)| (seg.id, idx))
        .collect();

    let mut coverage: Vec<f64> = vec![0.0; num_segments];

    if args.weight_queries {
        // First pass: count query occurrences
        let mut query_counts: HashMap<String, usize> = HashMap::new();
        for_each_line_in_file(&args.gaf, |l: &str| {
            let fields: Vec<&str> = l.split('\t').collect();
            if fields.len() >= 4 {
                let query_key = format!("{}:{}:{}", fields[0], fields[2], fields[3]);
                *query_counts.entry(query_key).or_insert(0) += 1;
            }
        });

        // Second pass: calculate coverage with query count adjustment
        for_each_line_in_file(&args.gaf, |l: &str| {
            let fields: Vec<&str> = l.split('\t').collect();
            let query_key = format!("{}:{}:{}", fields[0], fields[2], fields[3]);
            let count = query_counts.get(&query_key).unwrap_or(&1);

            for_each_step(
                l,
                |node_id, len| {
                    if let Some(&idx) = segment_id_to_index.get(&node_id) {
                        coverage[idx] += len as f64 / *count as f64;
                    }
                },
                |node_id| {
                    segment_id_to_index
                        .get(&node_id)
                        .map(|&idx| segments[idx].sequence.len())
                        .unwrap_or(0)
                },
            );
        });
    } else {
        // Single pass without weighting
        for_each_line_in_file(&args.gaf, |l: &str| {
            for_each_step(
                l,
                |node_id, len| {
                    if let Some(&idx) = segment_id_to_index.get(&node_id) {
                        coverage[idx] += len as f64;
                    }
                },
                |node_id| {
                    segment_id_to_index
                        .get(&node_id)
                        .map(|&idx| segments[idx].sequence.len())
                        .unwrap_or(0)
                },
            );
        });
    }

    if args.coverage_column {
        println!("##sample: {}", args.gaf);
        println!("#coverage");
        for (i, v) in coverage.into_iter().enumerate() {
            println!(
                "{}",
                if args.len_scale {
                    v / segments[i].sequence.len() as f64
                } else {
                    v
                }
            );
        }
    } else {
        print!("#sample");
        for seg in &segments {
            print!("\tnode.{}", seg.id);
        }
        println!();
        print!("{}", args.gaf);
        for (i, v) in coverage.into_iter().enumerate() {
            print!(
                "\t{}",
                if args.len_scale {
                    v / segments[i].sequence.len() as f64
                } else {
                    v
                }
            );
        }
        println!();
    }
}
