# BayesMAPJoint

Experimental joint-indicator implementation of BayesMAP. This repository
retains the original samplers and adds local blocked updates for connected
SNP, gene and cell indicators.

<div align="center">

### Bayesian Multi-layer Association Mapping

*A hierarchical Bayesian framework for integrating SNPs, genes and cell annotations using individual-level or GWAS summary statistics.*

---

[![R](https://img.shields.io/badge/R-%3E%3D4.2-blue.svg)]()
[![License: GPL--3](https://img.shields.io/badge/License-GPL--3-blue.svg)]()

</div>

---

# Overview

**BayesMAP** is an R package for Bayesian multi-layer association mapping that jointly models genetic variants, genes, cell types and biological annotations within a unified hierarchical Bayesian framework.

Unlike conventional GWAS methods that analyse variants independently, BayesMAP borrows information across multiple biological layers to improve identification of causal variants and biologically relevant functional units.

BayesMAP supports three analysis modes:

- Individual-level genotype and phenotype data
- GWAS summary statistics with a dense LD matrix
- GWAS summary statistics using block-wise LD matrices for scalable genome-wide analyses

Computationally intensive MCMC updates are implemented in **Rcpp/C++**, providing substantial speed improvements while retaining the flexibility of R.

---

# Methodological framework

BayesMAP models three biological layers simultaneously
The model estimates

- SNP effects (β)
- SNP (causal) inclusion indicators (δ)
- Gene (causal) inclusion indicators (Δ)
- Cell (causal) inclusion indicators (γ)
- SNP heritability
- Genetic variance
- Residual variance
- Posterior inclusion probabilities (PIPs)

Posterior inference is performed using Gibbs sampling.

---

# Features

✔ Bayesian multi-layer association mapping

✔ Individual-level and GWAS summary-statistics analyses

✔ Dense LD implementation

✔ Block-wise LD implementation for large datasets

✔ Hierarchical Bayesian variable selection

✔ Posterior inclusion probability estimation

✔ Bayesian estimation of SNP effects

✔ Efficient Gibbs sampler

✔ Computational bottlenecks accelerated using Rcpp

✔ Example datasets included

✔ Extensive unit testing

---

# Installation

Install the latest development version from GitHub.

```r
install.packages("remotes")

remotes::install_github("mommy003/BayesMAPJoint")

library(BayesMAPJoint)
```

---

# Quick start

## 1. Individual-level analysis

```r
library(BayesMAPJoint)

data(X)
data(y)
data(L)
data(A)
data(B)

fit <- bayesmap(
    X = X,
    y = y,
    L = L,
    A = A,
    B = B,
    niter = 1000,
    verbose = TRUE
)

summary(fit)
```

---

## 2. GWAS summary statistics (dense LD)

```r
library(BayesMAPJoint)

data(X)
data(y)
data(L)
data(A)
data(B)

N <- nrow(X)

## Standardise genotype matrix
X_std <- scale(
    X,
    center = TRUE,
    scale = FALSE
)

X_std <- sweep(
    X_std,
    MARGIN = 2,
    STATS = sqrt(colMeans(X_std^2)),
    FUN = "/"
)

## Standardise phenotype
y_std <- y - mean(y)
y_std <- y_std / sqrt(mean(y_std^2))

## Construct summary statistics
LD <- crossprod(X_std) / N

bhat <- as.numeric(
    crossprod(X_std, y_std) / N
)

fit_dense <- bayesmap_summary(
    bhat = bhat,
    LD = LD,
    N = N,
    L = L,
    A = A,
    B = B,
    phenotype_variance = 1,
    niter = 5000,
    burnin = 1000,
    thin = 2,
    verbose = TRUE
)

summary(fit_dense)
```

---

## 3. GWAS summary statistics (block-wise LD)

```r
library(BayesMAPJoint)

## Divide LD matrix into blocks
block_indices <- split(
    seq_len(ncol(LD)),
    ceiling(seq_len(ncol(LD)) / 50)
)

block_indices <- unname(block_indices)

LD_blocks <- lapply(
    block_indices,
    function(index) {
        LD[index, index, drop = FALSE]
    }
)

fit_block <- bayesmap_summary_blocks(
    bhat = bhat,
    LD_blocks = LD_blocks,
    block_indices = block_indices,
    N = N,
    L = L,
    A = A,
    B = B,
    phenotype_variance = 1,
    niter = 5000,
    burnin = 1000,
    thin = 2,
    verbose = TRUE
)

summary(fit_block)
```

---

# Outputs

All fitting functions return an object containing

- Posterior samples
- SNP effects
- SNP posterior inclusion probabilities
- Gene posterior inclusion probabilities
- Cell posterior inclusion probabilities
- Estimated heritability
- Genetic variance
- Residual variance
- MCMC diagnostics

For example,

```r
summary(fit)

head(fit$pip_snp)

head(fit$pip_gene)

head(fit$pip_cell)
```

---

# Block-wise implementation

The block-wise implementation partitions the LD matrix into approximately independent genomic blocks.

Compared with the dense implementation, it

- greatly reduces memory usage
- scales to genome-wide analyses
- produces identical posterior inference when the dense LD matrix is block diagonal

For real analyses, LD blocks should correspond to biologically meaningful independent regions (e.g. chromosomes or established LD blocks) rather than arbitrary equal-sized partitions.

---

# Validation

The package has been extensively validated.

### Individual-level vs dense summary statistics

Using consistently standardised genotype and phenotype data,

- posterior means are highly concordant
- posterior inclusion probabilities are highly correlated
- parameter recovery is nearly identical

### Dense vs block-wise summary statistics

When the dense LD matrix is converted to its block-diagonal representation,

- posterior samples are identical
- posterior means are identical
- variance components are identical
- SNP, gene and cell/pathway PIPs are identical

confirming correctness of the block-wise implementation.

---

# Current package structure

```
BayesMAP
│
├── bayesmap()
│
├── bayesmap_summary()
│
├── bayesmap_summary_blocks()
│
├── summary()
│
└── Rcpp
     ├── updateBetaDelta.cpp
     ├── updateDelta.cpp
     ├── updateGamma.cpp
     ├── updateBetaDeltaSummary.cpp
     └── updateBetaDeltaSummaryBlocks.cpp
```

---

# Development status

Current version includes

- Individual-level BayesMAP
- Dense summary-statistics implementation
- Block-wise summary-statistics implementation
- Rcpp acceleration
- Unit testing
- Example datasets

Planned developments include

- Full C++ implementation
- OpenMP parallelisation
- Multi-threaded MCMC
- Prediction utilities
- Visualisation functions
- Vignettes
- CRAN submission

---

# Citation

If you use BayesMAP in your research, please cite

> **BayesMAP:.......** 

---

# License

GPL-3

---

# Authors

**Dr Md Moksedul Momin**
Email: m.momin@uq.edu.au

**Dr Jian Zeng**
Email: j.zeng@imb.uq.edu.au

---

# Acknowledgements

BayesMAP has been developed at the Institute for Molecular Bioscience, The University of Queensland, Australia.
