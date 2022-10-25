# iRF-LOOP R Package
R Implementation of the iRF and iRF-LOOP algorithm. Please cite the following papers if you use iRF-LOOP:

```Cliff, A., Romero, J., Kainer, D., Walker, A., Furches, A., & Jacobson, D. (2019). A high-performance computing implementation of iterative random forest for the creation of predictive expression networks. Genes, 10(12), 996.```

```Walker, A. M., Cliff, A., Romero, J., Shah, M. B., Jones, P., Gazolla, J. G. F. M., ... & Kainer, D. (2022). Evaluating the performance of random forest and iterative random forest based methods when applied to gene expression data. Computational and Structural Biotechnology Journal, 20, 3372-3386.```

### What's in the package
* a Ranger-based custom implementation of iterative Random Forest (iRF)
* an implementation of iRF-LOOP built upon iRF

### The purpose of iRF-LOOP
iRF-LOOP is used to convert a matrix of features data into a network. It determines which features in the input matrix are most predictive of which other features in the input matrix, and outputs those relationships as weighted network edges. It can be applied to any sort of numeric data matrix where features are columns and samples are rows. It is very applicable to biological datasets such as gene expression, metabolomics, etc.

## The iRF-LOOP Algorithm
Given a data set of n features and m samples, iRF Leave One Out Prediction (iRF-LOOP) starts by treating one feature as the dependent variable (Y) and the remaining n − 1 features as predictors (X matrix). Using an iRF model, the importance of each feature in X, for predicting Y, is calculated. The result is a vector, of size n, of importance scores (the importance score of Y, for predicting itself is set to zero). This process is repeated for each of the n features, requiring n iRF runs which produces n vectors of importance scores. To keep importance scores on the same scale across the n runs, each run's result is normalized relative to the sum of that run's importance score vector. When all n runs are combined it represents a directed network where features are the nodes, and the ability of one feature to predict another (conditional on all other features) is represented as edges between nodes with weights defined by normalized importance scores.

![alt text](https://www.mdpi.com/genes/genes-10-00996/article_deploy/html/images/genes-10-00996-g001.png)

### Comparison to GENIE3
GENIE3 can be considered as an **RF-LOOP** algorithm. When used on a gene expression matrix, GENIE3 produces a predictive expression network (PEN) using Random Forest. With iRF-LOOP the **RF is replaced with iterative RF**. *Walker et al (2022)* showed that using iRF instead of RF when generating a PEN results in a smaller and more biologically correct network due to the feature selection and boosting process of iRF.

### A note about Speed
While this package uses Ranger (a very fast Random Forest) under the hood, it is only parallelized by threads on a single machine. This is usually fine for datasets where you have hundreds or even thousands of features as long as you have a reasonably powerful machine with lots of cores. However, if your dataset has 10s of thousands of features (e.g. a population-wide whole transcriptome dataset with 30000 genes measured in 1000 samples), then you should use the HPC implementation available here:

https://github.com/Jromero1208/RangerBasediRF
