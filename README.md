# scPloidyR

Single-cell copy number calling via Hidden Markov Models.

## Overview

scPloidyR implements a Hidden Markov Model (HMM) for detecting copy number states from single-cell sequencing data using depth and B-allele frequency (BAF) observations. It supports per-chromosome fitting, reference-aware parameter learning, and parallel multi-cell processing.

## Installation

```r
# Install from GitHub
devtools::install_github("dpei/scPloidyR")
```

## Usage

scPloidyR exports two main functions:

### `call_copy_number()`

Main entry point for parallel multi-cell per-chromosome HMM copy number calling.

```r
library(scPloidyR)

results <- call_copy_number(
  depth_matrix = depth_mat,
  baf_matrix = baf_mat,
  amplicon_metadata = amp_meta,
  reference_cells = ref_cells
)
```

### `aggregate_copy_number()`

Aggregate amplicon-level copy number calls to chromosome-level via mode.

```r
chr_cn <- aggregate_copy_number(results, amplicon_metadata = amp_meta)
```

### Reference CN support

HMM supports non-diploid reference cells via the `reference_cn` parameter:

```r
# Scalar (uniform for all amplicons, default is 2)
results <- call_copy_number(..., reference_cn = 3)

# Per-chromosome (named vector)
results <- call_copy_number(..., reference_cn = c("chr10" = 3))
```

## Dependencies

- R (>= 4.0.0)
- parallel
- stats
- utils

## License

MIT
