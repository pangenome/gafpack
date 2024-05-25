use clap::Parser;
use gfa::gfa::GFA;
use std::fs::File;
use std::io::{prelude::*, BufReader};
use std::collections::BTreeMap;

fn for_each_line_in_file(filename: &str, mut callback: impl FnMut(&str)) {
    let file = File::open(filename).unwrap();
    let reader = BufReader::new(file);
    for line in reader.lines() {
        callback(&line.unwrap());
    }
}

fn for_each_step(
    line: &str,
    mut node_callback: impl FnMut(usize, usize),
    mut edge_callback: impl FnMut(usize, usize),
    mut get_node_len: impl FnMut(usize) -> usize,
) {
    //eprintln!("{}", line);
    let walk = line.split('\t').nth(5).unwrap();
    if walk != "*" {
        //eprintln!("oheunotoeunthoue");
        let target_start = line.split('\t').nth(7).unwrap().parse::<usize>().unwrap();
        let target_end = line.split('\t').nth(8).unwrap().parse::<usize>().unwrap();
        let mut target_len = target_end - target_start;
        //eprintln!("target_len = {}", target_len);
        let fields = line.split('\t').nth(5).unwrap()
            .split(|c| c == '<' || c == '>')
            .filter(|s| !s.is_empty())
            .collect::<Vec<&str>>();
        let orientations = line.split('\t').nth(5).unwrap()
            .chars()
            .filter(|c| *c == '<' || *c == '>')
            .collect::<Vec<char>>();
        let mut seen: usize = 0;
        let fields_len = fields.as_slice().len();
        //eprintln!("fields len = {}", fields_len);
        for (i, field) in fields.iter().enumerate() {
            let j = field.parse::<usize>().unwrap();
            let mut len = get_node_len(j);
            //eprintln!("node {} len = {}", j, len);
            if i == 0 {
                //eprintln!("on first step {} {} {}", len, target_start, seen);
                assert!(len >= target_start);
                len -= target_start;
            }
            if i == fields_len - 1 {
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
            node_callback(j, len);
            if i < fields_len - 1 {
                let next_j = fields[i + 1].parse::<usize>().unwrap();
                if orientations[i] == '>' {
                    edge_callback(j, next_j);
                } else {
                    edge_callback(next_j, j);
                }
            }
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
    #[arg(short, long)]
    graph: String,
    /// Input GAF alignment file
    #[arg(short, long)]
    alignments: String,
    /// Emit node coverage vector
    #[arg(short, long)]
    node_coverage: bool,
    /// Emit edge coverage vector
    #[arg(short, long)]
    edge_coverage: bool,
    /// Scale node coverage values by node length
    #[arg(short, long)]
    len_scale: bool,
    /// Emit coverage vector in a single column
    #[arg(short, long)]
    coverage_column: bool,
}

fn main() {   
    let args = Args::parse();
    //println!("Hello {}!", args.graph);

    let gfa = {
        let parser = gfa::parser::GFAParser::default();
        let gfa: GFA<usize, ()> = parser.parse_file(&args.graph).unwrap();
        gfa
    };

    //println!("{} has {} nodes", args.graph, gfa.segments.len());
    //let mut lines = 0;
    //for_each_line_in_file(&args.alignments, |_l: &str| { lines += 1 });
    //println!("{} has {} nodes", args.alignments, lines);
    
    let mut node_coverage = vec![0; gfa.segments.len()];
    
    let mut edge_coverage: BTreeMap<(usize, usize), usize> = BTreeMap::new();
    // Initialize edge_coverage with all possible edges in the graph
    for edge in &gfa.links {
        edge_coverage.entry((edge.from_segment-1 as usize, edge.to_segment-1 as usize)).or_insert(0);
    }

    for_each_line_in_file(&args.alignments, |l: &str| {
        // get the start pos
        for_each_step(
            l,
            |i, len| { node_coverage[i-1] += len; },
            |from, to| {
                *edge_coverage.entry((from - 1, to - 1)).or_insert(0) += 1;
            },
            |id| { gfa.segments[id-1].sequence.len() });
        });

    if args.coverage_column {
        println!("##sample: {}", args.alignments);
        
        let mut header = String::from("#coverage");
        let mut header_parts = Vec::new();
        if args.node_coverage {
            header_parts.push(format!("num.nodes: {}", gfa.segments.len()));
        }
        if args.edge_coverage {
            header_parts.push(format!("num.edges: {}", edge_coverage.len()));
        }
        if !header_parts.is_empty() {
            header.push_str(&format!(" ({})", header_parts.join(", ")));
        }
        println!("{}", header);
        
        if args.node_coverage {
            for (i, v) in node_coverage.into_iter().enumerate() {
                println!("{}", if args.len_scale {v as f64  / gfa.segments[i].sequence.len() as f64} else {v as f64});
            }
        }
        
        if args.edge_coverage {
            for ((_, _), v) in edge_coverage.iter() {
                //println!("{} -> {}: {}", from + 1, to + 1, v);
                println!("{}", v);
            }
        }
    } else {
        print!("#sample");
        if args.node_coverage {
            for n in 1..=gfa.segments.len() {
                print!("\tnode.{}", n);
            }
        }
        if args.edge_coverage {
            for ((from, to), _) in edge_coverage.iter() {
                print!("\tedge.{}.{}", from + 1, to + 1);
            }
        }
        println!();
        print!("{}", args.alignments);
        if args.node_coverage {
            for (i, v) in node_coverage.into_iter().enumerate() {
                print!("\t{}", if args.len_scale {v as f64  / gfa.segments[i].sequence.len() as f64} else {v as f64});
            }
        }
        if args.edge_coverage {
            for ((_, _), v) in edge_coverage.iter() {
                print!("\t{}", v);
            }
        }
        println!();
    }
}
