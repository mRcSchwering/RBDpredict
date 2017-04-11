#' addErr
#'
#' Add a prediction error estimate (produced with cvRBDfit function) to a RBDmodel (produced with RBDfit function). 
#' This enables enables easy classification based on a FDR or recall.
#' 
#' If RBPs are to be predicted from new sequences the RBDmodel object is used.
#' Without a prediction error estimate the prediction will return probabilities for each new sequence to be a RBP.
#' If a cross validation (with cvRBDfit function) was done before and added to the RBDmodel
#' a threshold can be calculated based on a FDR or recall.
#' The prediction than classifies each new sequence according to this threshold.
#'
#' @param model   RBDmodel created with RBDfit
#' @param err     RBDerr created with cvRBDfit
#' 
#' @return as argument model (ist of class RBPmodel) with following element added
#' \itemize{
#'     \item RBDerr matrix with varying rows and 4 columns: Rows represent a classification each. 
#'           Columns represent prediction errors: "TP", "FP", "TN", "FN") 
#'  }
#' 
#' @examples
#' l <- c("A", "E", "K", "S", "R", "Q", "L", "T", "N", "G", "M", "W", "D", "H", "F", "Y", "C", "I", "P", "V")
#' sequences <- sapply(sample(10:100, 500, replace = TRUE), function(x) paste(sample(l, x, replace = TRUE), collapse = ""))
#' labels <- c(rep(1, 200), rep(0, 300))
#' 
#' model <- RBDfit(sequences, labels, k = 2)
#' err <- cvRBDfit(model, sequences, labels, nfold = 3)
#' model <- addErr(model, err)
#' 
#' model
#' 
#' @export
addErr <- function(model, err) {
  
  # check
  if(class(err) != "RBDerr") stop("Argument err must be of class RBDerr")
  if(class(model)!= "RBDmodel") stop("Argument model must be of class RBDmodel")
  
  model$RBDparams$nfold <- err$RBDparams$nfold
  model$RBDerr <- err$RBDerr
  
  model
}

