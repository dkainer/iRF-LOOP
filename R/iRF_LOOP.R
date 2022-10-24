#' @title iRF_LOOP (Leave One Out Prediction)
#' @description iRF_LOOP (Leave One Out Prediction)
#'
#' @param xmat feature matrix (i.e. the predictors for random forest)
#' @param iter number of iterations to run in each iRF model
#' @param ntree number of trees in each random forest. Default = 20.
#' @param verbose default = FALSE
#' @param threads Parallelize it. Do not set this higher than the number of cores you have. Default = 1.
#'
#' @return a data.frame containing the edges of the resulting network. This will have 4 columns:
#' \describe{
#'   \item{featX}{predictive feature}
#'   \item{featY}{feature being predicted}
#'   \item{imp}{normalized importance score for that feature combo}
#'   \item{R2}{the fit (variance explained) of the model that generated this edge}
#' }
#'
#' @importFrom foreach %dopar%
#'
#' @export
#'
#' @references Cliff, Ashley, et al. "A high-performance computing implementation of iterative random forest for the creation of predictive expression networks." Genes 10.12 (2019): 996.
#' @references Walker AM, Cliff A, Romero J, Shah MB, Jones P, Gazolla JG, Jacobson DA, Kainer D. Evaluating the performance of random forest and iterative random forest based methods when applied to gene expression data. Computational and Structural Biotechnology Journal. 2022 Jan 1;20:3372-86.
#' @examples
#' data(expdata)
#' # make a Predictive Expression Network
#' pen <- iRF_LOOP(expdata, ntree=100)
#'
#' # view top 10 edges in resulting PEN
#' top <- pen[order(pen$imp, decreasing = TRUE),][1:10,]
#' # note that the some top edges may come from a very poor model (low R2)
#'
globalVariables("g") # we need this to make foreach work with CRAN checks
iRF_LOOP <- function(xmat, iter=3, ntree=100, verbose=FALSE, threads=1)
{
  #for each feature in the X matrix, use it as a Y vector in iRF
  cl = parallel::makeCluster(threads)
  doParallel::registerDoParallel(cl)
  out <- foreach::foreach(g = 1:ncol(xmat), .combine=rbind, .export = "iRF") %dopar%
    {
      #source("~/devwork/_SCRIPTS/iRF.R")
      y   <- xmat[,g]
      x   <- xmat[,-g]
      irf <- iRF(x = x, y = y, ntree = ntree, threads = 1, iter=iter,  mtry = NULL, saveall = FALSE)  # must have saveall=F

      df  <- data.frame(featX = colnames(xmat)[-g], featY = colnames(xmat)[g], imp=irf$variable.importance, R2=irf$r.squared, row.names = NULL)
      df[df$imp>0,]
    }
  parallel::stopCluster(cl)
  return(out)
}
