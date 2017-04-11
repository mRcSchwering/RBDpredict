#' print.RBDmodel
#' 
#' summarizing the RBDmodel object
#' 
#' @export
print.RBDmodel <- function(o){
  
  cat("RBDmodel with", nrow(o$RBDfeatureExtraction), "sequences and labels\n")
  cat("Called with n =", o$RBDparams$n, ", k =", o$RBDparams$k," resulting in", o$RBDparams$feats, "features\n\n")
  
  cat("Resulting Model:\n\n")
  cat("regularization:\n")
  print(o$RBDmodel$regularization)
  cat("\nfrequencies:\n")
  print(o$RBDmodel$freqs)
  cat("\nalpha:\n")
  print(o$RBDmodel$alpha)
  cat("\nbeta:\n")
  str(o$RBDmodel$beta)
  
  if(!is.null(o$RBDerr)) cat("\nPrediction Error Estimate by", o$RBDparams$nfold, "fold Cross Validation included\n")
}


#' print.RBDerr
#' 
#' summarizing the RBDerr object
#' 
#' @export
print.RBDerr <- function(o){
  
  cat("RBDerr object\n\n")
  cat("Prediction error estimation with", o$RBDparams$n,"fold Cross Validation:\n")
  print(summary(o$RBDerr))
  cat("\nPosterior probabilities:\n")
  print(summary(o$RBDpp))
}


#' print.RBDpredict
#' 
#' summarizing the RBDpredict object
#' 
#' @export
print.RBDpredict <- function(o){
  
  cat("RBDpredict object\n")
  cat("Predicted", length(o$RBDpp), "sequences using", o$RBDparams$feats, "features\n\n")
  cat("Posterior Probabilities:\n")
  print(summary(o$RBDpp))
  
  if(!is.null(o$RBDclass)){
    cat("\nDecision boundary set to", o$RBDparams$threshold, "with", o$RBDparams$by, "set to", o$RBDparams$set, "\n")  
    cat("Classification:\n")
    cat(sum(o$RBDclass), "RBPs\t", sum(o$RBDclass == 0), "NoRBPs\n")
  }
}

