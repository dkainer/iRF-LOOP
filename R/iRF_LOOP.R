
iRF_LOOP <- function(xmat, iter=3, ntree=20, verbose=F, threads=1)
{
  library(ranger)
  library(foreach)
  
  #for each feature in the X matrix, use it as a Y vector in iRF 
  cl = parallel::makeCluster(threads)
  doParallel::registerDoParallel(cl)
  out <- foreach(g = 1:ncol(xmat), .combine=rbind, .export = "iRF") %dopar%
    {
      #source("~/devwork/_SCRIPTS/iRF.R")
      y   <- xmat[,g]
      x   <- xmat[,-g]
      irf <- iRF(x = x, y = y, ntree = ntree, threads = 1, iter=iter,  mtry = NULL, saveall = F)  # must have saveall=F
      
      df  <- data.frame(featX = colnames(xmat)[-g], featY = colnames(xmat)[g], imp=irf$variable.importance, R2=irf$r.squared, row.names = NULL)
      df[df$imp>0,]
    }
  parallel::stopCluster(cl)
  return(out)
}