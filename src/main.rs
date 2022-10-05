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

fn for_each_step(line: &str, mut callback: impl FnMut(usize)) {
    let walk = line.split('\t').nth(5).unwrap();
    if walk != "*" {
        line.split('\t').nth(5).unwrap()
            .split(|c| { c == '<' || c == '>' })
            .filter(|s| { !s.is_empty() })
            .map(|s| { s.parse::<usize>().unwrap() })
            .for_each(|i| { callback(i) });
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
    for_each_line_in_file(&args.alignments, |l: &str| { for_each_step(l, |i| { coverage[i-1] += 1; }); });
    print!("#sample");
    for n in 1..gfa.segments.len()+1 {
        print!("\tnode.{}", n);
    }
    print!("\n");
    print!("{}", args.alignments);
    for v in coverage {
        print!("\t{}", v);
    }
    print!("\n");
}
