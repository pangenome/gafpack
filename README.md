# gafpack

A Rust tool to convert alignments from pangenome variation graphs into coverage maps, facilitating haplotype-based genotyping analysis.

## Description

gafpack processes Graph Alignment Format (GAF) files against their corresponding Graph Fragment Assembly (GFA) graphs to generate coverage statistics. This is particularly useful for:
- Analyzing read coverage across pangenome variation graphs
- Supporting haplotype-based genotyping workflows
- Quantifying alignment distribution across graph nodes

## Requirements

- Rust toolchain (cargo)
- Input files:
  - GFA format graph file
  - GAF format alignment file

## Installation

Use `cargo` to build and install:

```bash
cargo install --force --path .
```

## Usage

Basic usage requires a graph file (GFA) and alignment file (GAF):

```bash
gafpack -g graph.gfa -a alignments.gaf >coverage.tsv
```

### Options

- `-g, --graph`: Input GFA pangenome graph file (required)
- `-a, --alignments`: Input GAF alignment file (required) 
- `-l, --len-scale`: Scale coverage values by node length
- `-c, --coverage-column`: Output coverage vector in single column format

### Output Formats

1. Default (tabular):
   ```
   #sample        node.1  node.2  node.3  ...
   alignments.gaf 1.5     2.0     0.5     ...
   ```

2. Column format (-c flag):
   ```
   ##sample: alignments.gaf
   #coverage
   1.5
   2.0
   0.5
   ...
   ```

## Examples

Scale coverage by node length:
```bash
gafpack -g graph.gfa -a reads.gaf --len-scale >scaled_coverage.tsv
```

Get column format output:
```bash
gafpack -g graph.gfa -a reads.gaf --coverage-column >coverage_vector.tsv
```
