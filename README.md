# gafpack

Calculate node coverage from GAF alignments to GFA variation graphs.

This is useful for:
- Analyzing read coverage across pangenome variation graphs
- Supporting haplotype-based genotyping workflows
- Quantifying alignment distribution across graph nodes

## Install

Requires Rust 2021 edition or later. Install using:

```bash
cargo install --git https://github.com/pangenome/gafpack
```

Or build from source:

```bash
git clone https://github.com/pangenome/gafpack
cd gafpack
cargo build --release
```

## Usage

Basic usage:

```bash
gafpack --gfa graph.gfa --gaf alignments.gaf > coverage.tsv
```

## Options

- `--gfa`: Input GFA graph file (required)
- `-g, --gaf`: Input GAF alignment file (required) 
- `-l, --len-scale`: Scale coverage by node length
- `-c, --coverage-column`: Output coverage vector as single column
- `-w, --weight-queries`: Weight coverage by query occurrences

## Output Formats

1. Default (tabular):

```
#sample        node.1  node.2  node.3  ...
alignments.gaf 1.5     2.0     0.5     ...
```

2. Column format (with `-c, --coverage-column`):

```
##sample: alignments.gaf
#coverage
1.5
2.0
0.5
...
```
