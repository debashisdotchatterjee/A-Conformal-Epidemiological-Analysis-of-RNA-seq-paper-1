# Distribution-Free Exposure–Response Inference for RNA-seq
## Conformal Epidemiological Analysis of the *airway* Dataset

This repository contains fully reproducible **R code** for the analysis presented in the manuscript:

**Distribution-Free Exposure–Response Inference for Transcriptomic Endpoints:  
A Conformal Epidemiological Analysis of RNA-seq Study**

The code implements a complete RNA-seq analysis pipeline with an explicit focus on **distribution-free statistical inference** for molecular epidemiology. Inference is obtained via **randomized conformal calibration**, layered on top of a standard RNA-seq working model, ensuring **finite-sample valid uncertainty quantification** without reliance on parametric distributional assumptions.

---

## Overview

RNA sequencing (RNA-seq) experiments are widely used in epidemiology to study how environmental or pharmacological exposures perturb molecular pathways. Standard RNA-seq workflows rely on parametric count models (e.g., Negative Binomial GLMs) whose inferential guarantees depend on correct model specification and large-sample approximations.

This repository demonstrates an alternative approach:
- A conventional RNA-seq model is used only as a **working predictor**
- **Randomized conformal inference** is applied to obtain calibrated p-values and prediction intervals
- Validity holds in finite samples under exchangeability, even if the count model is misspecified

The analysis is illustrated using the Bioconductor **`airway`** dataset, which examines the transcriptional response of human airway smooth muscle cells to dexamethasone exposure.

---

## Dataset

- **Dataset name:** `airway`
- **Source:** Bioconductor (R)
- **Design:**
  - 8 RNA-seq samples
  - 4 human airway smooth muscle cell lines
  - Paired treated vs untreated samples
  - Pharmacological exposure: dexamethasone

The dataset is publicly available and is loaded directly within the script using:

```r
library(airway)
data("airway")
