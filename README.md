# JCDR

**Joint Dimensionality Reduction for Robust Classification of High-Dimensional Data**

## Description

`JCDR` implements a **joint composite dimensionality reduction** framework tailored for high-dimensional classification tasks. Given labeled data \((X, Y)\), JCDR:

1. **Extracts two sets of ranked components** using complementary dimensionality reduction methods (e.g., PCA-style projections and PLS-style projections).
2. **Incorporates a first-moment (class-means) component**, enhancing class separation by capturing differences between per-class means.
3. **Performs a grid search over component counts**, using **k-fold cross-validation** to identify the combination of projections that **minimizes misclassification error**, yielding an optimal low-dimensional embedding and classifier.

This approach blends discriminative structural information (via PCA, PLS, and first-moment means) with robust model selection (via cross-validation). It is especially suited for high-dimensional settings where feature dimensionality greatly exceeds sample size and traditional methods risk overfitting or failing to capture class-discriminative directions.

- **PCA** ensures that directions capturing the largest variance are prioritized, providing a compact representation of high-dimensional data.  
- **PLS** projects data to maximize covariance with class labels, making it more predictive compared to variance-only methods like PCA.  
- **Cross-validation** provides principled model selection and tuning, especially for determining optimal dimensionality in supervised learning.

---

## System Requirements

- **R version**: 4.0.0 or higher (tested on R 4.2+)  
- **Operating systems**: Linux (e.g., Ubuntu 16.04+), macOS, and Windows (development and testing performed on Linux and macOS)  
- **Hardware**: A standard computer with at least 2 GB of RAM is sufficient for small datasets. For larger high-dimensional datasets, we recommend 16+ GB RAM and a multi-core CPU (≥4 cores).  
- **Dependencies**:  
  - [MASS](https://cran.r-project.org/package=MASS) (for LDA classifier)   
  - [lolR](https://cran.r-project.org/package=lolR) (for cross-validation utilities)  
  - [dplyr](https://cran.r-project.org/package=dplyr)  
  - [adabag](https://cran.r-project.org/package=adabag) (optional, AdaBoost)  
  - [class](https://cran.r-project.org/package=class) (for kNN)  
  - [randomForest](https://cran.r-project.org/package=randomForest) (optional, random forest)  

All required dependencies are installed automatically when installing **JCDR** with `devtools::install_github()`.  

---

## Installation

You can install the development version of **JCDR** from GitHub using the [`devtools`](https://cran.r-project.org/package=devtools) package:

```r
# Install devtools if not already installed
install.packages("devtools")

# Install JCDR from GitHub
devtools::install_github("riuzrj/JCDR")
```

---

## Usage

This GitHub repository contains an R package called **JCDR**, which provides the source code to perform **joint composite dimensionality reduction**. After installing the JCDR package, you can directly call the core function `spp()` to compute the projection matrix and embed high-dimensional data. For model selection, the function `spp.xval.optimal_dimselect()` can be used to perform cross-validation and choose the optimal combination and embedding dimension.

The return of `spp()` includes:
- `A`: projection matrix \([d, r]\) from original features to reduced dimensions.
- `Xr`: projected data \([n, r]\).
- `cr`: projected class centroids.
- `priors`: class prior probabilities.
- `ylabs`: unique class labels.

Here is a description of some of the important parameters:

- `X`: data matrix \([n, d]\), with `n` samples and `d` features.  
- `Y`: label vector \([n]\), with `K` unique class labels.  
- `r`: target embedding dimension. Must satisfy `r < d` and `r ≥ K`.  
- `method`: projection method. Options include `"mvi"`, `"xi"`, `"double.proj"`, and `"pca+pls"` (default).  
- `first.moment`: whether to add class-mean separation directions (`"delta"` or `FALSE`).  
- `robust.first`, `robust.second`: options for robust estimation.  

For cross-validation:

- `spp.xval.optimal_dimselect(X, Y, rs, alg, ...)`  
  - `rs`: candidate embedding dimensions to evaluate.  
  - `alg`: embedding algorithm.  
  - `k`: number of folds for CV (default `'loo'`).  
  - `classifier`: classifier for evaluation (default: `MASS::lda`).  
  - Returns: fold-wise errors, mean error per dimension, the optimal dimension, and the final model.  

---
## Examples

### Example 1: Projection with `spp()`

```r
library(JCDR)
library(mlbench)   # provides Sonar dataset

data(Sonar)
X <- as.matrix(Sonar[, 1:60])   # 60 features
Y <- Sonar[, 61]                # labels (M or R)

# Run projection into 10 dimensions
model <- spp(X, Y, r = 10, method = "pca+pls")

dim(model$A)   # projection matrix (60 x 10)
head(model$Xr) # projected data (n x 10)
```

### Example 2: Cross-validation with `spp.xval.optimal_dimselect()`

```r
library(MASS)
library(dplyr)
library(JCDR)
library(mlbench)
data(Sonar)

X <- as.matrix(Sonar[, 1:60])   # 60 features
Y <- Sonar[, 61]                # labels (M or R)

rs <- c(3, 4, 5, 10, 15, 20, 25, 30)
classifier.name <- c('lda')
cval <- spp.xval.optimal_dimselect(X, Y, rs, 
                    spp, alg.opts=list(method='pca+pls'),
                    alg.embedding="A", classifier.opts=list(),k=5)
```

## Issues

If you encounter any bugs or have feature requests, please [open an issue](https://github.com/riuzrj/JCDR/issues) on GitHub.
