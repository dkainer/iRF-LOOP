# iRF-LOOP
An R Implementation of the iRF and iRF-LOOP algorithms. iRF-LOOP stands for **iterative Random Forest - Leave One Out Prediction**. Please cite the following papers if you use iRF-LOOP:

>Cliff, A., Romero, J., Kainer, D., Walker, A., Furches, A., & Jacobson, D. (2019). A high-performance computing implementation of iterative random forest for the creation of predictive expression networks. Genes, 10(12), 996.

>Walker, A. M., Cliff, A., Romero, J., Shah, M. B., Jones, P., Gazolla, J. G. F. M., ... & Kainer, D. (2022). Evaluating the performance of random forest and iterative random forest based methods when applied to gene expression data. Computational and Structural Biotechnology Journal, 20, 3372-3386.

### What's in the package
* a Ranger-based custom implementation of iterative Random Forest (iRF)
* an implementation of iRF-LOOP built upon iRF

### The purpose of iRF-LOOP
iRF-LOOP is used to convert a matrix of features data into a "predictive" network. It determines which features in the input matrix are most predictive of which other features in the input matrix, and outputs those relationships as weighted network edges (in data.frame format). It can be applied to any sort of numeric data matrix where features are columns and samples are rows. It is very applicable to biological datasets such as gene expression, metabolomics, etc.

### Installation
From the R console, do this:
```
BiocManager::install("dkainer/iRF-LOOP")
```
or
```
library(devtools)
install_github("dkainer/iRF-LOOP")
```

## The iRF-LOOP Algorithm
Given a data set of n features and m samples, iRF Leave One Out Prediction (iRF-LOOP) starts by treating one feature as the dependent variable (Y) and the remaining n − 1 features as predictors (X matrix). Using an iRF model, the importance of each feature in X, for predicting Y, is calculated. The result is a vector, of size n, of importance scores (the importance score of Y, for predicting itself is set to zero). This process is repeated for each of the n features, requiring n iRF runs which produces n vectors of importance scores. To keep importance scores on the same scale across the n runs, each run's result is normalized relative to the sum of that run's importance score vector. When all n runs are combined it represents a directed network where features are the nodes, and the ability of one feature to predict another (conditional on all other features) is represented as edges between nodes with weights defined by normalized importance scores.

![alt text](https://www.mdpi.com/genes/genes-10-00996/article_deploy/html/images/genes-10-00996-g001.png)

### Comparison to GENIE3
GENIE3 can be considered as an **RF-LOOP** algorithm. When used on a gene expression matrix, GENIE3 produces a predictive expression network (PEN) using Random Forest. With iRF-LOOP the **RF is replaced with iterative RF**. *Walker et al (2022)* showed that using iRF instead of RF when generating a PEN results in a smaller and more biologically correct network due to the feature selection and boosting process of iRF. This comes at the cost of additional computational complexity.

### Comparison to Pearson's correlation (e.g Coexpression or WGCNA)
Given a feature matrix it is common to generate a network by taking pairwise correlations betwen features. This is very fast, but each model cor(fi,fj) fails to take into account the joint effects of all the other features, so it is unrealistic. It does, however, produce both positive and negative edge weights, unlike the RF-based methods.

### A note about computational speed
>**Note**:

While this package uses Ranger (a very fast Random Forest) under the hood, it is only parallelized by threads on a single machine. This is usually fine for datasets where you have hundreds or even thousands of features as long as you have a reasonably powerful machine with a good number of cores. However, if your dataset has 10s of thousands of features (e.g. a population-wide whole transcriptome dataset with 30,000 genes measured in 1000 samples) and you have access to a compute cluster, then you should use the HPC implementation available here:

https://github.com/Jromero1208/RangerBasediRF

## Examples

#### Running iRF with gene expression data
In this example we have *A.thaliana* RNAseq expression data (in TPM format) for 100 genes measured in 54 tissues.
This gives us 100 feature columns and 54 rows.
```
# load the example gene expression matrix
data(expdata)
expdata[1:10,1:10]

# use the first 99 genes to predict the expression of the 100th gene
irf <- iRF(x = expdata[, 1:99], y = expdata[, 100], iter=5, saveall=FALSE, num.trees = 500)
```
The RF from the best iteration is returned (since saveall=FALSE), while the stats from each iteration are shown in the console.
The first iteration is the equivalent of standard RF. Subsequent iterations often improve the fit of the model and 
it's Out Of Bag predictive ability while decreasing the number of active features to make a simpler model.
```
iRF iteration  1 
=================
mtry:   9.949874 
prediction error:   11.07547 
R^2:    0.6307986 
OOB cor(y,yhat):    0.825885 
Features with importance > 0: 90 
Warning: Split select weights used. Variable importance measures are only comparable for variables with equal weights.

iRF iteration  2 
=================
mtry:   9.486833 
prediction error:   7.532558 
R^2:    0.7489016 
OOB cor(y,yhat):    0.8783853 
Features with importance > 0: 79 
Warning: Split select weights used. Variable importance measures are only comparable for variables with equal weights.

iRF iteration  3 
=================
mtry:   8.888194 
prediction error:   6.890172 
R^2:    0.7703156 
OOB cor(y,yhat):    0.8855232 
Features with importance > 0: 54 
Warning: Split select weights used. Variable importance measures are only comparable for variables with equal weights.

iRF iteration  4 
=================
mtry:   7.348469 
prediction error:   6.22758 
R^2:    0.7924032 
OOB cor(y,yhat):    0.8943469 
Features with importance > 0: 45 
Warning: Split select weights used. Variable importance measures are only comparable for variables with equal weights.

iRF iteration  5 
=================
mtry:   6.708204 
prediction error:   5.969399 
R^2:    0.8010097 
OOB cor(y,yhat):    0.8956238 
Features with importance > 0: 34 
```

#### Predictive Expression Network using iRF-LOOP
```
data(expdata)
pen <- iRF_LOOP(expdata, num.trees=100)
# view the top 10 edges of highest importance
pen[order(pen$imp, decreasing = TRUE),][1:10,]

         featX     featY       imp           R2
1313 AT1G01120 AT1G01335 0.9676960 -0.005796995
7831 AT1G01710 AT1G01725 0.7650734  0.763543756
7923 AT1G01720 AT1G01650 0.7502535  0.842584195
8717 AT1G01790 AT1G01320 0.6986644  0.720483127
8517 AT1G01770 AT1G01480 0.6841163  0.590734042
7631 AT1G01700 AT1G01570 0.6773845  0.732813613
9823 AT1G01920 AT1G01510 0.6364524  0.916301947
5743 AT1G01510 AT1G01920 0.6318746  0.864956889
629  AT1G01570 AT1G01150 0.6128574  0.418503076
2120 AT1G01190 AT1G01380 0.5869559  0.531729914
```

The resulting PEN has 3349 edges (out of a possible 9900). These are directed edges - X may predict Y, but Y may not predict X with the same importance or at all.
You can see that the top edge in the PEN came from an iRF model with very poor fit (R2 ~ 0) so should be viewed with skepticism.
