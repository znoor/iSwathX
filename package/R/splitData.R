
##########################################################################
#' Splits data frame into arbitrary number of groups by rows randomly
#' 
#' @param dat The data frame to be split into groups
#' @param k a numeric value representing the number of groups
#' @return a list of data frames randomly split into a number of groups
##########################################################################


splitData <- function(dat, k=5){
  
  avesize <- round(nrow(dat)/k)
  lastsize <- nrow(dat) - avesize*(k-1)
  
  groupsizes <- c(rep(avesize, k-1), lastsize)    
    
    
  ids <- rep(1:k, groupsizes)
  
  
  # randomise
  set.seed(1)
  which.group <- sample(ids)
  
  
  split(dat, which.group, drop=TRUE)
}

