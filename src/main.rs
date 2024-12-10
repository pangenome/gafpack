use clap::Parser;
use gfa::gfa::GFA;
use std::fs::File;
use std::io::{prelude::*, BufReader};
use std::collections::HashMap;

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
fn for_each_step(line: &str,
                 mut callback: impl FnMut(usize,usize),
                 mut get_node_len: impl FnMut(usize) -> usize) {
    //eprintln!("{}", line);
    let walk = line.split('\t').nth(5).unwrap();
    if walk != "*" {
        //eprintln!("oheunotoeunthoue");
        let target_start = line.split('\t').nth(7).unwrap().parse::<usize>().unwrap();
        let target_end = line.split('\t').nth(8).unwrap().parse::<usize>().unwrap();
        let mut target_len = target_end - target_start;
        //eprintln!("target_len = {}", target_len);
        let fields = line.split('\t').nth(5).unwrap()
            .split(|c| { c == '<' || c == '>' })
            .filter(|s| { !s.is_empty() })
            .map(|s| { s.parse::<usize>().unwrap() })
            .enumerate()
            .collect::<Vec<(usize,usize)>>();
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
            if i == fields_len-1 {
                //eprintln!("on last step {} {} {}", len, target_end, seen);
                if target_len < seen {
                    target_len = seen;
                    //eprintln!("{}", line);
                    //eprintln!("on last step {} {} {}", len, target_len, seen);
                    //assert!(false);
                }
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
    let mut coverage : Vec<f64> = vec![0.0; gfa.segments.len()];
    
    if args.weight_queries {
        // First pass: count query occurrences
        let mut query_counts: HashMap<String, usize> = HashMap::new();
        for_each_line_in_file(&args.gaf, |l: &str| {
            let query_key = l.split('\t').next().unwrap().to_string();
            *query_counts.entry(query_key).or_insert(0) += 1;
        });

        // Second pass: calculate coverage with query count adjustment
        for_each_line_in_file(&args.gaf, |l: &str| {
            let query_key = l.split('\t').next().unwrap().to_string();
            let count = query_counts.get(&query_key).unwrap_or(&1);
            
            for_each_step(
                l,
                |i, j| { coverage[i-1] += j as f64 / *count as f64; },
                |id| { gfa.segments[id-1].sequence.len()
            });
        });
    } else {
        // Single pass without weighting
        for_each_line_in_file(&args.gaf, |l: &str| {
            for_each_step(
                l,
                |i, j| { coverage[i-1] += j as f64; },
                |id| { gfa.segments[id-1].sequence.len()
            });
        });
    }

    if args.coverage_column {
        println!("##sample: {}", args.gaf);
        println!("#coverage");
        for (i, v) in coverage.into_iter().enumerate() {
            println!("{}", if args.len_scale {v as f64  / gfa.segments[i].sequence.len() as f64} else {v as f64});
        }
    } else {
        print!("#sample");
        for n in 1..gfa.segments.len()+1 {
            print!("\tnode.{}", n);
        }
        println!();
        print!("{}", args.gaf);
        for (i, v) in coverage.into_iter().enumerate() {
            print!("\t{}", if args.len_scale {v as f64  / gfa.segments[i].sequence.len() as f64} else {v as f64});
        }
        println!();
    }
}
