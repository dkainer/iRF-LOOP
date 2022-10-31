#' @title iRF_LOOP
#' @name iRF_LOOP
#' @description iRF_LOOP (Leave One Out Prediction)
#'
#' @param xmat feature matrix (i.e. the predictors for random forest)
#' @param iter number of iterations to run in each iRF model
#' @param verbose default = FALSE
#'
#' @return a data.frame containing the edges of the resulting network. This will have 4 columns:
#' \describe{
#'   \item{featX}{predictive feature}
#'   \item{featY}{feature being predicted}
#'   \item{imp}{normalized importance score for that feature combo}
#'   \item{R2}{the fit (variance explained) of the model that generated this edge}
#' }
#' Only edges with importance > 0 are returned.
#'
#' @export
#'
#' @details Given a data set of n features and m samples, iRF Leave One Out Prediction (iRF-LOOP)
#' starts by treating one feature as the dependent variable (Y) and the
#' remaining n − 1 features as predictors (X matrix).
#' Using an iRF model, the importance of each feature in X, for predicting Y, is calculated.
#' The result is a vector, of size n, of importance scores (the importance score of Y,
#' for predicting itself is set to zero). This process is repeated for each of the n features,
#' requiring n iRF runs which produces n vectors of importance scores.
#' To keep importance scores on the same scale across the n runs, each run's result is
#' normalized relative to the sum of that run's importance score vector. When all n runs are combined
#' it represents a directed network where features are the nodes, and the ability of one feature to predict another
#' (conditional on all other features) is represented as edges between nodes with weights defined by normalized importance scores.
#'
#'     iRF-LOOP does the same overall process as **GENIE3**, but replaces the standard RF model with iRF.
#'     Walker et al (2022) demonstrated that replacing RF with iRF produces smaller and less noisy networks when
#'     applied to gene expression data. True positive edges usually rank higher in the edge list (ranked by importance score)
#'     when using iRF-LOOP rather than GENIE3.
#'
#'     iRF-LOOP can be used on any type of input data as long as it is in matrix form. However, for very large
#'     datasets where there may be tens of thousands of predictors (e.g. whole transcriptome expression data in hundreds of samples),
#'     this R package is likely to take considerable run time due to the need to run 10s of thousands of iRF models. For these
#'     situations we direct the user to the massively parallel [HPC implementation]<https://github.com/Jromero1208/RangerBasediRF>
#'
#' @references Cliff, Ashley, et al. "A high-performance computing implementation of iterative random forest for the creation of predictive expression networks." Genes 10.12 (2019): 996.
#' @references Walker AM, Cliff A, Romero J, Shah MB, Jones P, Gazolla JG, Jacobson DA, Kainer D. Evaluating the performance of random forest and iterative random forest based methods when applied to gene expression data. Computational and Structural Biotechnology Journal. 2022 Jan 1;20:3372-86.
#' @examples
#' data(expdata)
#' ## make a Predictive Expression Network
#' pen <- iRF_LOOP(expdata, num.trees=100)
#'
#' ## view top 10 edges in resulting PEN
#' ## Not run:
#' top <- pen[order(pen$imp, decreasing = TRUE),][1:10,]
#' ## note that the some top edges may come from a very poor model (low R2)
#'
#' @seealso \link[iterativeRF]{iRF}
#globalVariables("g") # we need this to make foreach work with CRAN checks
iRF_LOOP <- function(xmat, iter=3, first = NULL, last = NULL, verbose=FALSE,  ...)
{
  #for each feature in the X matrix, use it as a Y vector in iRF
  # out <- foreach::foreach(g = 1:ncol(xmat), .export = "iRF") %do%
  # {
  #     y   <- xmat[,g]
  #     x   <- xmat[,-g]
  #     irf <- iRF(x = x, y = y, iter=iter, mtry = NULL, saveall = FALSE, verbose=verbose, ...)  # must have saveall=F
  #     df  <- data.frame(featX = colnames(xmat)[-g], featY = colnames(xmat)[g], imp=irf$variable.importance, R2=irf$r.squared, row.names = NULL)
  #     df[df$imp>0,]
  # }

  out <- list()

  ## if user has declared a 'chunk' of columns to be predicted...
  if(!is.null(first) & !is.null(last) & first > 0 & last >= first)
  {
    for(g in first:last)
    {
      y   <- xmat[,g]
      x   <- xmat[,-g]
      irf <- iRF(x = x, y = y, iter=iter, mtry = NULL, saveall = FALSE, verbose=verbose, ...)  # must have saveall=F
      df  <- data.frame(featX = colnames(xmat)[-g], featY = colnames(xmat)[g], imp=irf$variable.importance, R2=irf$r.squared, row.names = NULL)
      out[[g]] <- df[df$imp>0,]
    }
  } else {
    for(g in 1:ncol(xmat))
    {
      y   <- xmat[,g]
      x   <- xmat[,-g]
      irf <- iRF(x = x, y = y, iter=iter, mtry = NULL, saveall = FALSE, verbose=verbose, ...)  # must have saveall=F
      df  <- data.frame(featX = colnames(xmat)[-g], featY = colnames(xmat)[g], imp=irf$variable.importance, R2=irf$r.squared, row.names = NULL)
      out[[g]] <- df[df$imp>0,]
    }
  }

  out <- do.call('rbind', out)
  return(out)
}
