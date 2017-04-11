#' cvRBDfit
#' 
#' Estimate prediction error of a RBDfit using n-fold Cross Validation.
#' 
#' 5-fold Cross Validation is a common method to estimate prediction errors.
#' This function can run for a while (roughly n-times as long as the corresponding RBDfit).
#'
#' @param model     RBDmodel created with RBDfit, this is the model for which to estimate prediction error
#' @param seq       chr arr containing amino acid sequences of proteins
#' @param labels    int arr containing labels: 0 = NoRBP, 1 = RBP
#' @param nfold     int(= 5) defining number of Cross Validations 
#'
#' @return list of class RBDerr with 3 elements
#' \itemize{
#'     \item RBDparams list of cross validation parameters
#'     \item RBDerr matrix with varying rows and 4 columns: Rows represent a classification each. 
#'           Columns represent prediction errors: "TP", "FP", "TN", "FN") 
#'     \item RBDpp num arr of posterior probabilities of each protein to be RBP 
#'  }
#' 
#' @examples
#' l <- c("A", "E", "K", "S", "R", "Q", "L", "T", "N", "G", "M", "W", "D", "H", "F", "Y", "C", "I", "P", "V")
#' sequences <- sapply(sample(10:100, 500, replace = TRUE), function(x) paste(sample(l, x, replace = TRUE), collapse = ""))
#' labels <- c(rep(1, 200), rep(0, 300))
#' 
#' model <- RBDfit(sequences, labels, k = 2)
#' err <- cvRBDfit(model, sequences, labels, nfold = 3)
#' 
#' str(err)
#' err
#' ROCplot(err, sens = 0.5)
#' plot(err)
#' 
#' @export
cvRBDfit <- function(model, seq, labels, nfold = 5) {
  
  # check
  if( class(model) != "RBDmodel" ) stop("Argument model must be of class RBDmodel")
  if( nfold < 2 ) stop("Choose reasonable value for nfold")
  if( !identical(labels, as.numeric(as.logical(labels))) ) stop("Argument labels must be of 0 and 1")

  nFeats <- model$RBDparams$feats
  k <- model$RBDparams$k_mer
  
  # feature space of 20 biological relevant amino acids
  feats <- 1:k
  feats <- lapply(feats, function(x) c("A", "E", "K", "S", "R", "Q", "L", "T", "N", "G", "M", "W", "D", "H", "F", "Y", "C", "I", "P", "V"))
  feats <- expand.grid(feats, stringsAsFactors = FALSE)
  feats <- sort(apply(feats, 1, function(x) paste(x, collapse = "")))
  
  # extract k-mers from aa sequences
  seq <- toupper(seq)
  kmers <- lapply(seq, function(x) as.vector(GetKMers(x, K = k)))
  kmers <- lapply(kmers, function(x) factor(x, levels = feats))
  kmers <- t(sapply(kmers, function(x) as.vector(table(x))))
  colnames(kmers) <- feats

  # do CV
  split <- sample(1:nfold, size = nrow(kmers), replace = TRUE)
  pred <- list()
  labs <- list()
  for( i in 1:nfold ){
    cat("\n\n\nCross Validation round", i, "\n\n")
    eval <- kmers[split == i, ]
    training <- kmers[split != i, ]
    ranks <- sda::sda.ranking(training, factor(labels[split != i], labels = c("NoRBP", "RBP"), levels = c(0, 1)))
    training <- training[, ranks[1:nFeats, 1]]
    training <- training[, order(colnames(training))]
    eval <- eval[, ranks[1:nFeats, 1]]
    eval <- eval[, order(colnames(eval))]
    trained <- sda::sda(training, factor(labels[split != i], labels = c("NoRBP", "RBP"), levels = c(0, 1)))
    pred[[i]] <- sda::predict.sda(trained, eval)$posterior
    labs[[i]] <- labels[split == i]
  }
  
  # calculate confusion matrix
  cat("\n\nCalculating Errors\n")
  pred <- do.call(rbind, pred)
  pred <- pred[, 2]
  labs <- unlist(labs)
  err <- Err(pred, labs)
  
  out <- list(RBDparams = list(nfold = nfold), RBDerr = err, RBDpp = pred)
  class(out) <- "RBDerr"
  out
}


