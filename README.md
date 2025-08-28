# Feature Level Report

A Rust tool for analyzing genomic alignments and generating feature-level alignment statistics.

## Overview

This tool processes alignment data with genomic features, counting aligned bases and tracking indels between query and target sequences. It's particularly useful for analyzing how genomic features (like genes or exons) align between different sequences.

## Installation

```bash
cargo build --release
```

## Usage

```bash
feature_level_report [OPTIONS]

OPTIONS:
    -i, --input <FILE>              Input file (supports .gz files)
    -m, --max-indel-size <INT>      Maximum indel size to consider as indels (default: unlimited)
    -h, --help                      Print help information
    -V, --version                   Print version information
```

## Input Format

The tool expects tab-separated input with the following fields:
- Query and target sequence information
- CIGAR string in field 12 (with "cg:Z:" prefix)
- Feature coordinates for both query and target sequences
- Feature names and strand information

## Output

Generates a tab-separated report with columns:
- Feature name and coordinates
- Query and target sequence names
- Alignment statistics (aligned bases, indels, unaligned regions)
- Separate counts for query and target sequences

## Example

```bash
# Process a gzipped alignment file
feature_level_report -i alignments.tsv.gz -m 50

# Process standard input
cat alignments.tsv | feature_level_report
```

## Author

Andrea Guarracino <aguarra1@uthsc.edu>

## License

See LICENSE file for details.