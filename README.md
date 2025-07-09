
# HMMclone

HMMclone is an R package for inferring integer copy number profiles from pseudobulk clones derived from sparse single-cell whole genome sequencing (scWGS) data. The package uses a Hidden Markov Model (HMM) to segment the genome and call copy number states across chromosomes.

## Overview

The package provides a complete workflow for:
- Reading and counting sequencing reads from BAM files in genomic bins
- Applying GC bias correction to read counts
- Fitting Hidden Markov Models to infer copy number states
- Using the Viterbi algorithm to find the most likely copy number state sequence

## Suggested workflow

This package is designed to work with low pass single-cell whole genome sequencing with the aim of producing high resolution copy number in clones produced by aggregating information from multiple cells. The suggested workflow (and what was done in Williams et al. Nature 2025) is something along the lines of:

1. Call copy number in single cells with a larger bin size (eg 500kb). For example using [mondrian](https://github.com/mondrian-scwgs/mondrian)
2. Identify clones based on copy number across single cells. For example using a single cell phylogenetic method like MEDICC2 or a clustering algorithm such as the one available in [SIGNALS](https://github.com/shahcompbio/signals).
3. Count reads in smaller bins (eg 10kb) across cells and then aggregate the read counts for the clones
4. Normalize read counts against normal diploid cells and/or perform GC correction
5. Run HMMClone to call copy number states

## Installation

You can install the development version of HMMclone from GitHub:

```r
# Install devtools if you haven't already
if (!require("devtools")) install.packages("devtools")

# Install HMMclone
devtools::install_github("shahcompbio/HMMclone")
```

## Main Functions

### 1. Read Counting

**`count_reads()`** - Count reads per genomic bin from BAM files

```r
# Count reads in 10kb bins
read_counts <- count_reads(
  bamfile = "path/to/sample.bam",
  binsize = "10000",
  genome = "hg19",
  minMapq = 20,
  isPaired = TRUE,
  isDuplicate = FALSE
)
```

After this step, you will want to aggregate read counts across cells present in each clone.

### 2. GC Correction

**`gc_cor_hmm()`** - GC correction similar to HMMcopy approach

```r
# Apply GC correction (HMMcopy-style)
corrected_counts <- gc_cor_hmm(
  counts = read_counts$reads,
  gc = read_counts$gc,
  span1 = 0.03,
  span2 = 0.3
)
```

**`gc_cor_modal()`** - Modal regression GC correction

```r
# Apply modal GC correction
corrected_counts <- gc_cor_modal(
  counts = read_counts$reads,
  gc = read_counts$gc,
  lowess_frac = 0.2,
  q = c(0.1, 0.9)
)
```

### 3. HMM Fitting and Copy Number Calling

**`HMMclone()`** - Main function for copy number calling across multiple clones

If you already have gc corrected read counts per clone you can skip steps 1 and 2 and go directly here.

```r
# Prepare data with required columns: clone_id, chr, start, end, copy
# where 'copy' is the GC-corrected read count data

# Call copy numbers using HMM
results <- HMMclone(
  dat = your_data,
  selftransitionprob = 0.99,  # Probability of staying in same state
  maxCN = 15,                 # Maximum copy number to consider
  sd_value = 0.2,             # Standard deviation for emission probabilities
  ncores = 4,                 # Number of cores for parallel processing
  progressbar = TRUE
)
```

**`fitHMM()`** - Fit HMM to a single clone

```r
# Fit HMM to copy number data for a single clone
hmm_result <- fitHMM(
  copy = single_clone_data,
  selftransitionprob = 0.99,
  maxCN = 15,
  sd_value = 0.2,
  usecpp = TRUE  # Use C++ implementation for speed
)
```

## Advanced Usage

### Filtering Low-Quality Bins

You can exclude certain bins from inference by adding a `keep` column:

```r
# Add quality filter
your_data$keep <- your_data$mappability > 0.8 & your_data$gc > 0.3 & your_data$gc < 0.7

# Run HMM (will automatically use only bins where keep = TRUE)
results <- HMMclone(dat = your_data, ...)
```

## Output Format

The `HMMclone()` function returns a data.frame with the original data plus a `state` column containing the inferred copy number for each genomic bin:

```r
# Output columns include:
# - clone_id: Clone identifier
# - chr: Chromosome
# - start, end: Genomic coordinates
# - copy: Input copy number data (GC-corrected)
# - state: Inferred copy number state (0, 1, 2, 3, ...)
# - Additional columns from input data
```

## Performance Notes

- The package includes both R and C++ implementations of the Viterbi algorithm
- Use `usecpp = TRUE` (default) for better performance
- Multi-core processing is supported via the `ncores` parameter


## Issues and Support

For questions and issues, please file an issue on the GitHub repository.
