# gafpack

Calculate node coverage from GAF alignments to GFA variation graphs, with optional flow-constrained copy number estimation.

## Install

```bash
cargo install --git https://github.com/pangenome/gafpack
```

Or build from source:

```bash
git clone https://github.com/pangenome/gafpack
cd gafpack
cargo build --release
```

### Optional Features

**METIS partitioning** (for large graphs):
```bash
cargo build --release --features metis
```

**Gurobi solver** (requires Gurobi installation):
```bash
cargo build --release --features gurobi
```

Both features:
```bash
cargo build --release --features "metis,gurobi"
```

## Usage

### Coverage Mode

Basic coverage computation:

```bash
gafpack --gfa graph.gfa --gaf alignments.gaf > coverage.tsv
```

Column format output:
```bash
gafpack --gfa graph.gfa --gaf alignments.gaf -c > coverage.txt
```

### Copy Number Mode

Estimate copy number using a negative binomial model with ILP flow constraints:

```bash
gafpack --gfa graph.gfa --gaf alignments.gaf --copy-number -c > cn.txt
```

For large graphs, use partitioning:
```bash
gafpack --gfa graph.gfa --gaf alignments.gaf --copy-number \
    --partition-size 500000 -c > cn.txt
```

With METIS graph-aware partitioning (requires `metis` feature):
```bash
gafpack --gfa graph.gfa --gaf alignments.gaf --copy-number \
    --partition-size 500000 --use-metis -c > cn.txt
```

With Gurobi solver (requires `gurobi` feature and Gurobi installation):
```bash
gafpack --gfa graph.gfa --gaf alignments.gaf --copy-number \
    --solver gurobi -c > cn.txt
```

Fast ML estimation without flow constraints:
```bash
gafpack --gfa graph.gfa --gaf alignments.gaf --copy-number --no-flow -c > cn.txt
```

## Options

### General
| Option | Description |
|--------|-------------|
| `--gfa` | Input GFA graph file (supports .gz/.bgz) |
| `-g, --gaf` | Input GAF alignment file |
| `-c, --coverage-column` | Output as single column |
| `-v, --verbose` | Verbosity: 0=warn, 1=info, 2=debug |

### Coverage Mode
| Option | Description |
|--------|-------------|
| `-l, --len-scale` | Scale coverage by node length |
| `-w, --weight-queries` | Weight by query occurrences |

### Copy Number Mode
| Option | Description |
|--------|-------------|
| `--copy-number` | Enable copy number estimation |
| `--ploidy` | Background CN values for parameter estimation [default: 1,2] |
| `--bin-size` | Bin size (bp) for NB parameter estimation [default: 100] |
| `--epsilon` | CN=0 sensitivity (lower = more deletions) [default: 0.02] |
| `--no-flow` | Skip ILP, use simple ML estimation |
| `--no-dedup` | Disable read deduplication |
| `-t, --threads` | ILP solver threads (0 = auto) |

### ILP Tuning
| Option | Description |
|--------|-------------|
| `--complexity` | 1=basic, 2=+edge penalty, 3=+reverse edge penalty [default: 2] |
| `--source-prob` | Expensive super-edge penalty [default: -10000] |
| `--cheap-source` | Cheap super-edge penalty [default: -25] |
| `--prob-scale` | Coverage vs flow balance [default: 1] |
| `--diff-cutoff` | CN probability cutoff [default: 4×|source_prob|] |

### Partitioning (Large Graphs)
| Option | Description |
|--------|-------------|
| `--partition-size` | Max nodes per partition [default: 1500000] |
| `--use-metis` | Use METIS graph-aware partitioning (requires feature) |
| `--solver` | ILP solver: `highs` (default) or `gurobi` (requires feature) |

## Output Formats

### Tabular (default)
```
#sample        node.1  node.2  node.3  ...
alignments.gaf 1.5     2.0     0.5     ...
```

### Column format (`-c`)
```
##sample: alignments.gaf
#coverage
1.5
2.0
0.5
```

In copy number mode with `-c`:
```
##sample: alignments.gaf
#copy_number
2
2
1
```
