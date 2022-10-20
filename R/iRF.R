

#' iterative Random Forest
#'
#' @param x feature matrix (i.e. the predictors for random forest)
#' @param y the dependent variable (i.e. the variable being predicted)
#' @param ntree number of trees in the random forest
#' @param iter number of iterations to run
#' @param classification boolean. If FALSE, then use Regression RF. Default = FALSE
#' @param threads default = 1
#' @param alwayssplits Character vector with variable names to be always selected in addition to the mtry variables tried for splitting.
#' @param mtry Number of variables to possibly split at in each node.
#' Default is the (rounded down) square root of the number variables with importance > 0.
#' If 0 < mtry < 1, then it will be treated as a proportion of variables with importance > 0.
#' Otherwise it is the absolute number of variables to use when splitting a node.
#' @param saveall if TRUE then return the final random forest from each iteration as a list of forests
#'
#' @return if saveall=TRUE, returns a list of random forest objects. Otherwise one random forest object
#' from the iteration that had the best model fit.
#' @export
#'
#' @examples
#'
#' iRF(x = xmat, y = yvec, saveall=TRUE)
#'
iRF <- function(x, y, ntree=500, iter=5, classification=F, threads=1, alwayssplits=NULL, mtry=NULL, saveall=T)
{
  tmp <- cbind(x, Y = y)
  wt  <- rep(1/ncol(x), ncol(x)) # start with equal sample weighting per SNP
  rfs <- list()
  for(i in 1:iter)
  {
    cat("\niRF iteration ",i,"\n")
    cat("=================\n")
    if(is.null(mtry)){
      m = sqrt(sum(wt>0))
    } else if(mtry <= 1) {
      m = mtry*sum(wt>0)
    }

    rf <- ranger::ranger(dependent.variable.name = "Y", data = tmp, num.trees=ntree,
                         split.select.weights = wt, classification = classification,
                         mtry = m, importance = "impurity", num.threads=threads, write.forest = T,
                         always.split.variables = alwayssplits)

    wt        <- rf$variable.importance / sum(abs(rf$variable.importance)) # scale importances to range(0,1)
    wt[wt<0]  <- 0 # set negative weights to zero

    cat("mtry:  ", m, "\n")
    cat("prediction error:  ",rf$prediction.error,"\n")
    if(classification==FALSE) cat("R^2:   ",rf$r.squared,"\n")  # variance explained, computed on OOB data
    if(classification==TRUE) print(rf$confusion.matrix)

    cat("OOB cor(y,yhat):   ",cor(rf$predictions,y),"\n")
    cat("Features with importance > 0:",sum(wt>0),"\n")

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

