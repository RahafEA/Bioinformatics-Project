

# Overview

This project demonstrates a complete mini bioinformatics workflow using R programming.  
The pipeline combines:

1. RNA-Seq count normalization and preprocessing
2. DNA sequence alignment and mutation analysis
3. Proteomics quality control filtering
4. Differential protein abundance analysis and visualization


---

# Technologies & Packages Used

- R Programming Language
- RStudio
- Biostrings
- pwalign
- ggplot2
- dplyr

---

# Installation

Install required packages before running the project:

```r
install.packages("dplyr")
install.packages("ggplot2")

install.packages("BiocManager")

BiocManager::install("Biostrings")
BiocManager::install("pwalign")
```

Load libraries:

```r
library(dplyr)
library(ggplot2)
library(Biostrings)
library(pwalign)
```

---



# Part 1 — RNA-Seq Count Normalization 
Eng Ahmed(https://github.com/elmokade)


## Objective

Simulate RNA-Seq count data and apply preprocessing, filtering, normalization, and visualization techniques commonly used in transcriptomics analysis.

---

## Workflow

### 1. Simulate RNA-Seq Count Matrix

Generated a matrix of:
- 5,000 genes
- 6 samples
  - 3 Control
  - 3 Treated

```r
base_counts <- matrix(
  rpois(5000 * 6, lambda = 80),
  nrow = 5000,
  ncol = 6
)
```

---

### 2. Simulate Different Sequencing Depths

Different library sizes were simulated using scaling factors:

```r
depth_factors <- c(0.7, 1.0, 1.3, 0.8, 1.2, 1.5)
```

This mimics real RNA-Seq experiments where samples have unequal sequencing depth.

---

### 3. Gene Filtering

Genes with very low expression were removed:

```r
filtered_counts <- counts[gene_totals >= 15, ]
```

Purpose:
- Remove noise
- Improve statistical reliability
- Reduce sparsity

---

### 4. Missing Value Simulation

Five artificial missing values (`NA`) were inserted randomly:

```r
filtered_counts[random_rows[i], random_cols[i]] <- NA
```

---

### 5. Missing Value Imputation

Missing values were replaced using group-specific means:

- Control samples → Control mean
- Treated samples → Treated mean

This preserves biological group structure.

---

### 6. CPM Normalization

Formula:

```math
CPM = (Gene Count / Total Sample Counts) × 1,000,000
```

Purpose:
- Correct sequencing depth differences
- Make samples comparable

---

### 7. Log Transformation

```math
log2(counts + 1)
```

Purpose:
- Reduce skewness
- Stabilize variance
- Improve visualization

---

## Outputs

- Cleaned count matrix
- CPM-normalized matrix
- Log2-transformed matrix
- Side-by-side boxplots

---

## Key Results

| Metric | Value |
|---|---|
| Genes | 5,000 |
| Samples | 6 |
| Missing values inserted | 5 |
| Remaining NAs after imputation | 0 |

---

# Part 2 — DNA Sequence Alignment & Variant Analysis
ENG Rana (https://github.com/ranayaser)
## Objective

Analyze sequence variation between wild-type and mutant TP53 sequences.

---

## Input Files

- `TP53_WT.fasta`
- `TP53_Variant.fasta`

---

## Workflow

### 1. Read FASTA Sequences

```r
wt_seq <- readDNAStringSet("TP53_WT.fasta")[[1]]
mut_seq <- readDNAStringSet("TP53_Variant.fasta")[[1]]
```

---

### 2. Global Sequence Alignment

Performed pairwise global alignment:

```r
alignment <- pwalign::pairwiseAlignment(
  wt_seq,
  mut_seq,
  type = "global"
)
```

Purpose:
- Compare full-length sequences
- Detect mutations and mismatches

---

## Alignment Results

| Metric | Value |
|---|---|
| Alignment Type | Global |
| Alignment Score | 2318.81 |
| Sequence Identity | 99.75% |
| Sequence Length | 1,182 bp |

---

### 3. DNA → RNA Conversion

```r
wt_rna <- chartr("T", "U", as.character(wt_seq))
```

This converts DNA thymine (T) into RNA uracil (U).

---

### 4. Protein Translation

```r
variant_protein_mut <- translate(mut_rna)
```

The RNA sequence was translated into an amino acid sequence.

---

### 5. GC Content Calculation

```math
GC% = ((G + C) / (A + T + G + C)) × 100
```

Result:
- GC Content = 56.60%

---

### 6. Premature Stop Codon Detection

The translated protein sequence was scanned for internal stop codons (`*`).

Classification:
- Internal stop codon → Truncated Protein
- Stop codon only at the end → Full-Length Protein

---

## Final Biological Interpretation

The variant sequence showed:
- Very high similarity to the wild-type TP53 sequence
- No premature stop codon
- Normal full-length protein production

Conclusion:
The detected mutation likely does not severely disrupt protein integrity.

---

# Part 3 — Proteomics Data Quality Filtering
ENG Mahmoud (https://github.com/Mahmoud70-7)
## Objective

Simulate proteomics intensity data and apply biologically meaningful quality filtering.

---

## Workflow

### 1. Simulate Proteomics Matrix

Generated:
- 800 proteins
- 6 samples

```r
prot_intensities <- matrix(
  rnorm(4800, mean = 20, sd = 4),
  ncol = 6
)
```

---

### 2. Introduce Missing Values

```r
prot_intensities[sample(1:4800, 400)] <- NA
```

Missingness rate:
- 8.33%

---

### 3. Quality Filtering

Proteins were retained if they had:

- ≥2 valid values in Control
OR
- ≥2 valid values in Treated

```r
ctrl_valid <- sum(!is.na(x[1:3])) >= 2
treat_valid <- sum(!is.na(x[4:6])) >= 2
```

Purpose:
- Remove extremely sparse proteins
- Improve downstream reliability

---

### 4. Missing Value Imputation

Imputation rule:

```math
Imputation\ Value = Minimum\ Intensity / 5
```

Result:
- Minimum intensity = 6.237
- Imputed value = 1.247

Reason:
Low abundance proteins are often missing in mass spectrometry datasets.

---

## Key Results

| Metric | Value |
|---|---|
| Total Proteins | 800 |
| Missing Values | 400 |
| Missingness Rate | 8.33% |
| Retained Proteins | ~792 |
| Minimum Intensity | 6.237 |
| Imputation Value | 1.247 |

---

# Part 4 — Differential Protein Abundance Analysis
Rahaf (https://github.com/RahafEA)
## Objective

Identify proteins showing abundance differences between Control and Treated groups.

---

## Workflow

### 1. Calculate Group Means

```r
groupA_mean <- rowMeans(filtered_prot[, 1:3])
groupB_mean <- rowMeans(filtered_prot[, 4:6])
```

---

### 2. Compute Fold Change

```math
Fold\ Change = Treated\ Mean / Control\ Mean
```

Interpretation:
- FC > 1 → Upregulated
- FC < 1 → Downregulated

---

### 3. Z-score Normalization

```math
z = (x - \mu) / \sigma
```

Purpose:
- Standardize protein intensities
- Improve visualization comparability

---

### 4. Visualization

Scatter plot:
- X-axis → Control mean intensity
- Y-axis → Treated mean intensity
- Color → Fold Change

Outlier proteins appear far from the diagonal.

---

## Key Results

| Metric | Value |
|---|---|
| Example Upregulated Protein | Prot_124 (FC = 16.51) |
| Example Downregulated Protein | Prot_699 (FC = 0.097) |

---

# Biological Interpretation

Because the dataset is simulated randomly:
- No true biological signal exists
- Upregulated and downregulated proteins are approximately balanced
- Extreme fold changes mainly result from imputation effects and random variation



