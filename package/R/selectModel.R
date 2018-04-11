###############################################################################
#' Model selection for time-base RT prediction using cross validation
#' @param datRTrain a data frame of retention time training set, must include
#' two columns: RT_detected.base and RT_detected.
#' @param k a numeric value specifying the fold number of the cross validation
#' @return a model object based on the best RMSE average on cross validation
#' @examples
#' libfiles <- paste(system.file("files",package="iSwathX"),
#'                              c("Lib2.txt","Lib3.txt"),sep="/")
#' datBaseLib <- readLibFile(libfiles[1], clean=TRUE, nomod=FALSE, nomc=FALSE)
#' datExtLib <- readLibFile(libfiles[2], clean=TRUE, nomod=FALSE, nomc=FALSE)
#'
#'
#' datRTrain <- getRTrain(datBaseLib, datExtLib, nomod=FALSE)
#' bestmodel <- selectModel(datRTrain)
#'
############################################################################

# source("splitData.R")

selectModel <- function(datRTrain, k=5){

	bestModel  <- NULL

  if(nrow(datRTrain) < 20) {
	warning("Size of RT training set is too small (<20) to get a good model!")
	} else {
  candidateModels <- c("lm", "svm")


  ll = splitData(datRTrain)

  rmses <- matrix(0, nrow=2, ncol=k)


  f <- as.formula("RT_detected.base ~ RT_detected")

  # best parameters for svm
  obj <- tune.svm(f, data=datRTrain, kernel="radial",
                  gamma = c(0.001,0.01, 0.1),cost=c(0.1,1,10,20))

  for(jj in 1:2){

    alg <- candidateModels[jj]

    for( ii in 1:length(ll)){
      dtrain = do.call(rbind, ll[-ii])
      dtest = ll[[ii]]

      if(alg == "lm")
        fit = get(alg)(f, data=dtrain)
      else {# svm
        fit = get(alg)(f, data=dtrain, kernel="radial", gamma=obj$best.parameters$gamma, cost=obj$best.parameters$cost)
      }
      predicted = predict(fit, dtest)

      rmses[jj,ii] <- sqrt(mean((dtest$RT_detected.base - predicted)^2))
    }
  }

  # average rmse of cross validation
  rmse.avgs <- apply(rmses, 1, mean)

  # select minum
  model.idx <- which(rmse.avgs == min(rmse.avgs))



  if(candidateModels[model.idx] == "lm")
    bestModel <- lm(f, datRTrain)

  else
    bestModel <- svm(f, datRTrain, kernel="radial", gamma=obj$best.parameters$gamma,
                                   cost=obj$best.parameters$cost)

  }
  bestModel
}




