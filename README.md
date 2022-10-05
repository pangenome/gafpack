# gafpack

Convert alignments to pangenome variation graphs to coverage maps useful in haplotype-based genotyping.

## building

Use `cargo` to build:

```bash
cargo install --force --path .
```

## usage

Provide a graph (in GFA) and alignments (in GAF):

```bash
gafpack -g z.gfa -a z.gaf >z.pack.tsv
```

The output gives the graph coverage vector with a header line.
