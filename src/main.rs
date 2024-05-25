use clap::Parser;
use gfa::gfa::GFA;
use std::fs::File;
use std::io::{prelude::*, BufReader};
use std::collections::BTreeMap;
use log::{error};

fn for_each_line_in_file(filename: &str, mut callback: impl FnMut(&str)) {
    let file = File::open(filename).unwrap();
    let reader = BufReader::new(file);
    for line in reader.lines() {
        callback(&line.unwrap());
    }
}

fn for_each_step(
    line: &str,
    mut node_callback: Option<impl FnMut(usize, usize)>,
    mut edge_callback: Option<impl FnMut(usize, usize)>,
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
            if let Some(ref mut cb) = node_callback {
                cb(j, len);
            }
            if i < fields_len - 1 {
                let next_j = fields[i + 1].parse::<usize>().unwrap();
                if orientations[i] == '>' {
                    if let Some(ref mut cb) = edge_callback {
                        cb(j, next_j);
                    }
                } else {
                    if let Some(ref mut cb) = edge_callback {
                        cb(next_j, j);
                    }
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
    env_logger::init();
    
    let args = Args::parse();
    //println!("Hello {}!", args.graph);

    // Check if at least one of node_coverage or edge_coverage is specified
    if !args.node_coverage && !args.edge_coverage {
        error!("At least one of --node-coverage or --edge-coverage must be specified.");
        std::process::exit(1);
    }

    let gfa = {
        let parser = gfa::parser::GFAParser::default();
        let gfa: GFA<usize, ()> = parser.parse_file(&args.graph).unwrap();
        gfa
    };

    //println!("{} has {} nodes", args.graph, gfa.segments.len());
    //let mut lines = 0;
    //for_each_line_in_file(&args.alignments, |_l: &str| { lines += 1 });
    //println!("{} has {} nodes", args.alignments, lines);
    
    let mut node_coverage = if args.node_coverage {
        Some(vec![0; gfa.segments.len()])
    } else {
        None
    };

    let mut edge_coverage = if args.edge_coverage {
        let mut ec = BTreeMap::new();
        for link in &gfa.links {
            ec.entry((link.from_segment - 1, link.to_segment - 1)).or_insert(0);
        }
        Some(ec)
    } else {
        None
    };

    for_each_line_in_file(&args.alignments, |l: &str| {
        // get the start pos
        for_each_step(
            l,
            node_coverage.as_mut().map(|nc| {
                move |i, len| {
                    nc[i - 1] += len;
                }
            }),
            edge_coverage.as_mut().map(|ec| {
                move |from, to| {
                    if !ec.contains_key(&(from - 1, to - 1)) {
                        error!("Edge from {} to {} does not exist in the graph.", from, to);
                        std::process::exit(1);
                    } else {
                        *ec.entry((from - 1, to - 1)).or_insert(0) += 1;
                    }
                }
            }),
            |id| { gfa.segments[id-1].sequence.len() });
        });

    if args.coverage_column {
        println!("##sample: {}", args.alignments);
        
        let mut header = String::from("#coverage");
        let mut header_parts = Vec::new();
        if let Some(ref nc) = node_coverage {
            header_parts.push(format!("num.nodes: {}", nc.len()));
        }
        if let Some(ref ec) = edge_coverage {
            header_parts.push(format!("num.edges: {}", ec.len()));
        }
        if !header_parts.is_empty() {
            header.push_str(&format!(" ({})", header_parts.join(", ")));
        }
        println!("{}", header);
        
        if let Some(nc) = node_coverage {
            for (i, v) in nc.into_iter().enumerate() {
                println!("{}", if args.len_scale { v as f64 / gfa.segments[i].sequence.len() as f64 } else { v as f64 });
            }
        }
        
        if let Some(ec) = edge_coverage {
            for ((_, _), v) in ec.iter() {
                println!("{}", v);
            }
        }
    } else {
        print!("#sample");
        if let Some(_) = node_coverage {
            for n in 1..=gfa.segments.len() {
                print!("\tnode.{}", n);
            }
        }
        if let Some(_) = edge_coverage {
            for ((from, to), _) in edge_coverage.as_ref().unwrap().iter() {
                print!("\tedge.{}>{}", from + 1, to + 1);
            }
        }
        println!();
        print!("{}", args.alignments);
        if let Some(nc) = node_coverage {
            for (i, v) in nc.into_iter().enumerate() {
                print!("\t{}", if args.len_scale { v as f64 / gfa.segments[i].sequence.len() as f64 } else { v as f64 });
            }
        }
        if let Some(ec) = edge_coverage {
            for ((_, _), v) in ec.iter() {
                print!("\t{}", v);
            }
        }
        println!();
    }
}
