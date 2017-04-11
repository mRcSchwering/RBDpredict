#' plot.RBDerr
#' 
#' plotting posterior probability distribution of Cross Validated RBDmodel
#' 
#' @export
plot.RBDerr <- function(err){
  
  hist(err$RBDpp, breaks = 30, col = "gray50", border = "gray70", xlab = "Posterior Probability", 
       main = paste("Posterior Probability distribution after", err$RBDparams$n, "fold Cross Validation"), xlim = c(0, 1))
}


#' plot.RBDpredict
#' 
#' plotting posterior probability distribution of RBP prediction
#' 
#' @export
plot.RBDpredict <- function(pred){
  
  if(!is.null(pred$RBDclass)){
    sub <- paste("Red line indicates classification boundary with", pred$RBDparams$by, "of", pred$RBDparams$set)
  } else sub <- NULL
  
  hist(pred$RBDpp, breaks = 30, col = "gray50", border = "gray70", xlab = "Posterior Probability", 
       main = "Posterior Probability distribution of Prediction", xlim = c(0, 1), sub = sub)
  
  if(!is.null(pred$RBDclass)) abline(v = pred$RBDparams$threshold, lty = 2, col = "#db4437")
}


#' plot.RBDmodel
#' 
#' plotting 40 top ranking features
#' 
#' @export
plot.RBDmodel <- function(model){
  
  plot(model$RBDranking)
}