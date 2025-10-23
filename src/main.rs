use clap::Parser;
use gfa::gfa::GFA;
use std::collections::HashMap;
use std::fs::File;
use std::io::{prelude::*, BufReader};

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
    mut callback: impl FnMut(usize, usize, bool),
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
        let mut steps: Vec<(usize, bool)> = Vec::new();
        let mut current_number = String::new();
        let mut current_orientation: Option<char> = None;
        for ch in walk.chars() {
            if ch == '<' || ch == '>' {
                if let Some(orientation) = current_orientation {
                    if !current_number.is_empty() {
                        let node_id = current_number.parse::<usize>().unwrap();
                        steps.push((node_id, orientation == '<'));
                        current_number.clear();
                    }
                }
                current_orientation = Some(ch);
            } else {
                current_number.push(ch);
            }
        }
        if let Some(orientation) = current_orientation {
            if !current_number.is_empty() {
                let node_id = current_number.parse::<usize>().unwrap();
                steps.push((node_id, orientation == '<'));
            }
        }

        let mut seen: usize = 0;
        let fields_len = steps.len();
        //eprintln!("fields len = {}", fields_len);
        for (i, (node_id, is_reverse)) in steps.into_iter().enumerate() {
            let mut len = get_node_len(node_id);
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
            callback(node_id, len, is_reverse);
        }
        //eprintln!("seen = {}", seen);
    }
    //eprintln!("at end");
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
    /// Split coverage by strand when a node is visited in both orientations
    #[arg(long)]
    strand: bool,
}

fn main() {
    let args = Args::parse();
    //println!("Hello {}!", args.gfa);
    let gfa = {
        let parser = gfa::parser::GFAParser::default();
        let gfa: GFA<usize, ()> = parser.parse_file(&args.gfa).unwrap();
        gfa
    };
    //println!("{} has {} nodes", args.gfa, gfa.segments.len());
    //let mut lines = 0;
    //for_each_line_in_file(&args.gaf, |_l: &str| { lines += 1 });
    //println!("{} has {} nodes", args.gaf, lines);
    let node_count = gfa.segments.len();
    let mut coverage_forward: Vec<f64> = vec![0.0; node_count];
    let mut coverage_reverse: Vec<f64> = vec![0.0; node_count];
    let mut orientation_mask: Vec<u8> = vec![0; node_count];

    // Build orientation mask from graph paths if --strand is specified
    // This determines which orientations each node appears in within the graph structure
    if args.strand {
        for path in &gfa.paths {
            for segment_ref in &path.nodes {
                let node_id = segment_ref.segment_id;
                let idx = node_id - 1;
                if segment_ref.is_reverse {
                    orientation_mask[idx] |= 0x2; // Reverse orientation
                } else {
                    orientation_mask[idx] |= 0x1; // Forward orientation
                }
            }
        }
    }

    // Now process the GAF alignments to calculate coverage
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
                |i, j, is_reverse| {
                    let idx = i - 1;
                    let value = j as f64 / *count as f64;
                    if is_reverse {
                        coverage_reverse[idx] += value;
                    } else {
                        coverage_forward[idx] += value;
                    }
                },
                |id| gfa.segments[id - 1].sequence.len(),
            );
        });
    } else {
        // Single pass without weighting
        for_each_line_in_file(&args.gaf, |l: &str| {
            for_each_step(
                l,
                |i, j, is_reverse| {
                    let idx = i - 1;
                    let value = j as f64;
                    if is_reverse {
                        coverage_reverse[idx] += value;
                    } else {
                        coverage_forward[idx] += value;
                    }
                },
                |id| gfa.segments[id - 1].sequence.len(),
            );
        });
    }

    if args.coverage_column {
        println!("##sample: {}", args.gaf);
        println!("#coverage");
        if args.strand {
            let labels: Vec<String> = (0..node_count)
                .flat_map(|idx| {
                    let node_id = idx + 1;
                    let mask = orientation_mask[idx];
                    match mask {
                        0x3 => vec![format!("node.{}+", node_id), format!("node.{}-", node_id)],
                        0x2 => vec![format!("node.{}-", node_id)],
                        0x1 => vec![format!("node.{}+", node_id)],
                        _ => vec![format!("node.{}", node_id)],
                    }
                })
                .collect();
            println!("#nodes:\t{}", labels.join("\t"));
        }
        for idx in 0..node_count {
            let node_len = gfa.segments[idx].sequence.len() as f64;
            let forward = coverage_forward[idx];
            let reverse = coverage_reverse[idx];
            if args.strand {
                match orientation_mask[idx] {
                    0x3 => {
                        // Node appears in both orientations in the graph
                        let forward_val = if args.len_scale {
                            forward / node_len
                        } else {
                            forward
                        };
                        let reverse_val = if args.len_scale {
                            reverse / node_len
                        } else {
                            reverse
                        };
                        println!("{}", forward_val);
                        println!("{}", reverse_val);
                    }
                    0x2 => {
                        // Node appears only in reverse orientation in the graph
                        let value = if args.len_scale {
                            reverse / node_len
                        } else {
                            reverse
                        };
                        println!("{}", value);
                    }
                    0x1 => {
                        // Node appears only in forward orientation in the graph
                        let value = if args.len_scale {
                            forward / node_len
                        } else {
                            forward
                        };
                        println!("{}", value);
                    }
                    _ => {
                        // Node doesn't appear in any graph paths
                        let total = forward + reverse;
                        let value = if args.len_scale {
                            total / node_len
                        } else {
                            total
                        };
                        println!("{}", value);
                    }
                }
            } else {
                // Without --strand, combine forward and reverse coverage
                let total = forward + reverse;
                println!(
                    "{}",
                    if args.len_scale {
                        total / node_len
                    } else {
                        total
                    }
                );
            }
        }
    } else {
        // Tabular output format
        print!("#sample");
        if args.strand {
            for idx in 0..node_count {
                let node_id = idx + 1;
                match orientation_mask[idx] {
                    0x3 => {
                        print!("\tnode.{}+", node_id);
                        print!("\tnode.{}-", node_id);
                    }
                    0x2 => {
                        print!("\tnode.{}-", node_id);
                    }
                    0x1 => {
                        print!("\tnode.{}+", node_id);
                    }
                    _ => {
                        print!("\tnode.{}", node_id);
                    }
                }
            }
        } else {
            for n in 1..=node_count {
                print!("\tnode.{}", n);
            }
        }
        println!();
        print!("{}", args.gaf);
        if args.strand {
            for idx in 0..node_count {
                let node_len = gfa.segments[idx].sequence.len() as f64;
                let forward = coverage_forward[idx];
                let reverse = coverage_reverse[idx];
                match orientation_mask[idx] {
                    0x3 => {
                        // Node appears in both orientations in the graph
                        let forward_val = if args.len_scale {
                            forward / node_len
                        } else {
                            forward
                        };
                        let reverse_val = if args.len_scale {
                            reverse / node_len
                        } else {
                            reverse
                        };
                        print!("\t{}", forward_val);
                        print!("\t{}", reverse_val);
                    }
                    0x2 => {
                        // Node appears only in reverse orientation in the graph
                        let value = if args.len_scale {
                            reverse / node_len
                        } else {
                            reverse
                        };
                        print!("\t{}", value);
                    }
                    0x1 => {
                        // Node appears only in forward orientation in the graph
                        let value = if args.len_scale {
                            forward / node_len
                        } else {
                            forward
                        };
                        print!("\t{}", value);
                    }
                    _ => {
                        // Node doesn't appear in any graph paths
                        let total = forward + reverse;
                        let value = if args.len_scale {
                            total / node_len
                        } else {
                            total
                        };
                        print!("\t{}", value);
                    }
                }
            }
        } else {
            // Without --strand, combine forward and reverse coverage
            for idx in 0..node_count {
                let node_len = gfa.segments[idx].sequence.len() as f64;
                let total = coverage_forward[idx] + coverage_reverse[idx];
                print!(
                    "\t{}",
                    if args.len_scale {
                        total / node_len
                    } else {
                        total
                    }
                );
            }
        }
        println!();
    }
}
