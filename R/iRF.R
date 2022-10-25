#' @title iterative Random Forest
#' @name iRF
#' @description an algorithmic advancement of Random Forest (RF),
#' which takes advantage of Random Forests ability to determine feature importance,
#' and produces a more accurate model by iteratively creating weighted forests.
#'
#' @param x Matrix. Features as columns (i.e. the predictors for random forest), samples as rows
#' @param y Numeric. The dependent variable (i.e. the variable being predicted)
#' @param iter Integer. number of iterations to run. Default = 5.
#' @param mtry Numeric. Number of variables to possibly split at in each node.
#' Default is the (rounded down) square root of the number variables with importance > 0.
#' If 0 < mtry < 1, then it will be treated as a proportion of variables with importance > 0.
#' Otherwise it is the absolute number of variables to use when splitting a node.
#' @param classification Boolean. If FALSE, then use Regression RF. Default = FALSE
#' @param saveall Boolean. if TRUE then return the final random forest from each iteration as a list of forests
#' @param usepvals Boolean. Apply permutation testing to feature importances in order to cull features based on P-values
#' @param verbose Boolean. default = TRUE
#' @param ... further parameters for Ranger::ranger function
#'
#'
#' @return One of two possible results:
#' if saveall=TRUE, returns a list of random forest objects (one per iteration).
#' if saveall=FALSE, returns the random forest object from the iteration that had the best model fit.
#' @export
#' @details
#'     Iterative Random Forest (iRF) is an algorithmic advancement of Random Forest (RF),
#'     which takes advantage of Random Forests ability to determine feature importance, and produces a more accurate
#'     model by iteratively creating weighted forests. In each iteration, the importance scores from the
#'     previous iteration are used to weight features for the current forest. The effect is akin to a lasso regression since
#'     some features have their importance boosted while others go to zero and are culled from the forest.
#'     The result is often a more accurate set of high ranked predictors from a better fit and simpler model.
#'
#'     If you set **usepvals=TRUE** then Ranger will use Altmann's method for calculating importance P-values using permutation.
#'     iRF will then cull features based on having a conservative False Discovery Rate Q > 0.2 (i.e. it sets their importance to 0).
#'     This usually results in a harsher culling of features than just allowing importances to naturally fall to zero.
#'     **WARNING**: using pvalues will be MUCH slower than not using it. So it may not be a good idea when performing iRF-LOOP
#'
#'
#' @references Basu S, Kumbier K, Brown JB et al. (2018) Iterative random forests to discover predictive and stable high-order interactions. Proc Natl Acad Sci USA 115, 1943–1948.
#' @examples
#'
#' data(expdata)
#' ## use the first 99 genes' expression to predict the expression of the 100th gene
#' irf <- iRF(x = expdata[, 1:99], y = expdata[, 100], saveall=FALSE, num.trees = 100)
#' ## view top 10 most important genes
#' ## Not run:
#' sort(irf$variable.importance, decreasing = TRUE)[1:10]
#'
#' @seealso \link[ranger]{ranger} \link[iterativeRF]{iRF_LOOP}
#'
iRF <- function(x, y, iter=5, mtry=NULL, classification=FALSE, saveall=TRUE, usepvals = FALSE, verbose=TRUE, ...)
{
  tmp <- cbind(x, Y = y)
  wt  <- rep(1/ncol(x), ncol(x)) # start with equal sample weighting per SNP
  rfs <- list()
  for(i in 1:iter)
  {

    if(is.null(mtry)){
      m = sqrt(sum(wt>0))
    } else if(mtry <= 1) {
      m = mtry*sum(wt>0)
    }

    rf <- ranger::ranger(dependent.variable.name = "Y",
                         data = tmp,
                         split.select.weights = wt,
                         mtry = m,
                         classification = classification,
                         importance = "impurity_corrected",
                         ...)

    # for the first iteration we used double as many trees and get importance p-values
    if(usepvals == TRUE)
    {
      imp.pvals <- ranger::importance_pvalues(rf, method = "altmann", formula = Y ~ ., data = tmp, num.permutations = 500)
      imp.pvals[ imp.pvals[,"importance"] < 0 , "importance" ] <- 0
      fdr <- stats::p.adjust(p = imp.pvals[,"pvalue"], method = "fdr")
      imp.pvals[ fdr > 0.2 , "importance" ] <- 0
      wt <- imp.pvals[,"importance"] / sum(imp.pvals[,"importance"])

    } else {
      imp <- rf$variable.importance
      wt  <- imp / sum(abs(imp)) # scale importances to range(0,1)
      wt[wt<0]  <- 0 # set negative weights to zero
    }

    rf$variable.importance <- wt

    if(verbose)
    {
      cat("\niRF iteration ",i,"\n")
      cat("=================\n")
      cat("mtry:  ", m, "\n")
      cat("prediction error:  ",rf$prediction.error,"\n")
      if(classification==FALSE) cat("R^2:   ",rf$r.squared,"\n")  # variance explained, computed on OOB data
      if(classification==TRUE) print(rf$confusion.matrix)

      cat("OOB cor(y,yhat):   ",stats::cor(rf$predictions,y),"\n")
      cat("Features with importance > 0:",sum(wt>0),"\n")
    }

    rfs[[i]] <- rf

    if(sum(wt>0) < max(0.01*(ncol(x)), 10))
    {
      break
    }
  }

  if(saveall==TRUE) {
    return(rfs)
  } else {
    R2 <- sapply(rfs, function(iter){iter$r.squared })  # return the iteration with max fit
    return(rfs[[which.max(R2)]])
  }

}


