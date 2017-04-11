#' Err
#' 
#' compute confusion matrices from posterior probabilities and labels for a number of thresholds
#' thresholds will be distributed between 0 and 1. posterior > t will assign class as 1
#'
#' @param posterior num arr of posterior probabilities of each datapoint
#' @param labels    int arr of factor of length(posterior) of true class labels. must have {0, 1}
#' @param n         int(= 10000) for how many thresholds will be set between 0 and 1 to compute matrices.
#'                  For equal="freq" this is an estimation
#' @param equal     chr(= freq) "width" or "freq". set thresholds with equal width or equal frequency
#' 
#' @return matrix of n row and 4 cols with TP, FP, TN, FN
#' 
Err <- function(posterior, labels, n = 10000, equal = "freq"){
  
  t <- switch (equal,
               freq = Hmisc::cut2(sort(posterior), g = n, onlycuts = TRUE),
               width = seq(0, 1, length.out = n)
  )
  
  labels <- factor(labels)
  labelLevels <- levels(labels)
  
  classi <- lapply(t, function(x) as.numeric(posterior > x))
  classi <- lapply(classi, function(x) factor(x, levels = labelLevels))
  classi <- lapply(classi, function(x) table(x, labels))
  
  err <- do.call("rbind", lapply(classi, function(x) c(x[2, 2], x[2, 1], x[1, 1], x[1, 2])))
  rownames(err) <- round(t, 4)
  colnames(err) <- c("TP", "FP", "TN", "FN")
  
  err
}