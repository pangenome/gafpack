use clap::Parser;
use gfa::gfa::GFA;
use std::fs::File;
use std::io::{prelude::*, BufReader};

fn for_each_line_in_file(filename: &str, mut callback: impl FnMut(&str)) {
    let file = File::open(filename).unwrap();
    let reader = BufReader::new(file);
    for line in reader.lines() {
        callback(&line.unwrap());
    }
}

fn for_each_step(line: &str,
                 mut callback: impl FnMut(usize,usize),
                 mut get_node_len: impl FnMut(usize) -> usize) {
    let walk = line.split('\t').nth(5).unwrap();
    let target_start = line.split('\t').nth(7).unwrap().parse::<usize>().unwrap();
    let target_end = line.split('\t').nth(8).unwrap().parse::<usize>().unwrap();
    let target_len = target_end - target_start;
    //println!("{}", line);
    //println!("target_len = {}", target_len);
    if walk != "*" {
        let fields = line.split('\t').nth(5).unwrap()
            .split(|c| { c == '<' || c == '>' })
            .filter(|s| { !s.is_empty() })
            .map(|s| { s.parse::<usize>().unwrap() })
            .enumerate()
            .collect::<Vec<(usize,usize)>>();
        let mut seen: usize = 0;
        let fields_len = fields.as_slice().len();
        //println!("fields len = {}", fields_len);
        for (i, j) in fields {
            let mut len = get_node_len(j);
            //println!("node {} len = {}", j, len);
            if i == 0 {
                //println!("on first step {} {} {}", len, target_start, seen);
                assert!(len >= target_start);
                len -= target_start;
            }
            if i == fields_len-1 {
                //println!("on last step {} {} {}", len, target_end, seen);
                if target_end < seen {
                    println!("{}", line);
                    println!("on last step {} {} {}", len, target_end, seen);
                    assert!(false);
                }
                len = target_end - seen;
            }
            if i == fields_len {
                assert!(false);
            }
            //println!("node {} adj len = {}", j, len);
            seen += len;
            callback(j, len);
        }
        //println!("seen = {}", seen);
    }
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
    /// Scale coverage values by node length
    #[arg(short, long)]
    len_scale: bool,
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
    let mut lines = 0;
    for_each_line_in_file(&args.alignments, |_l: &str| { lines += 1 });
    //println!("{} has {} nodes", args.alignments, lines);
    let mut coverage = vec![0; gfa.segments.len()];
    for_each_line_in_file(&args.alignments, |l: &str| {
        // get the start pos
        for_each_step(
            l,
            |i, j| { coverage[i-1] += j; },
            |id| { gfa.segments[id-1].sequence.len() }); });
    print!("#sample");
    for n in 1..gfa.segments.len()+1 {
        print!("\tnode.{}", n);
    }
    println!();
    print!("{}", args.alignments);
    for (i, v) in coverage.into_iter().enumerate() {
        print!("\t{}", v as f64 / gfa.segments[i].sequence.len() as f64);
    }
    println!();
}
