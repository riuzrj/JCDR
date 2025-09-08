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



