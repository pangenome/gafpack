use clap::Parser;
use gafpack::{compute_coverage, format_coverage_column};

/// Project a GAF alignment file into coverage over GFA graph nodes
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Input GFA pangenome graph file (supports .gz/.bgz compression)
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

    let (coverage, min_id) =
        compute_coverage(&args.gfa, &args.gaf, args.len_scale, args.weight_queries).unwrap();

    let num_segments = coverage.len();

    if args.coverage_column {
        print!("{}", format_coverage_column(&args.gaf, &coverage));
    } else {
        print!("#sample");
        for n in min_id..min_id + num_segments {
            print!("\tnode.{}", n);
        }
        println!();
        print!("{}", args.gaf);
        for v in coverage {
            print!("\t{}", v);
        }
        println!();
    }
}
