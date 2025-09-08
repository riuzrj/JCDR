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



