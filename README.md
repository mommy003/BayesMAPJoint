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
```
remotes::install_github(
  "mommy003/BayesMAPJoint",
  upgrade = "never"
)
```
Then verify:
```
library(BayesMAPJoint)

packageVersion("BayesMAPJoint")
exists("bayesmap_joint_prior")
```
