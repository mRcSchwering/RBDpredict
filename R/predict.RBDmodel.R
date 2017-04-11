#' predict.RBDmodel
#' 
#' @aliases predict
#'
#' Predict RBPs from a new sequence based on RBDmodel and its prediction error estimate.
#' 
#' Without a prediction error estimate (produced with cvRBDfit function) the predictor will return plain posterior probabilities
#' of each sequence being a RBP. 

#' If an estimate was added to the model (via addErr function) FDR or sens can be set (if both are set, sens is ignored).
#' In this case a classification is returned as well.
#' Note that especially for FDR the desired value cannot always be achieved.
#'
#' @param model   RBDmodel created with RBDfit
#' @param seq     chr arr of new protein sequences
#' @param sens    num(= NULL) desired sensitivtiy of classification (will be ignored if FDR is set)
#' @param FDR     num(= NULL) desired false discovery rate (1 - precision) of classification
#' 
#' @return list of class RBPprediction with 2 elements 
#' \itemize{
#'     \item RBPpp posterior probabilities of seq being RBPs
#'     \item RBPparams list of parameters used in the prediction
#' }
#' If sens or FDR was set a vector containing the classification will be added
#' \itemize{
#'     \item RBPclass classification of seq (1 = RBP, 0 = NoRBP)
#' }
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
#' newseq <- sapply(sample(10:100, 500, replace = TRUE), function(x) paste(sample(l, x, replace = TRUE), collapse = ""))
#' pred <- predict(model, newseq, sens = 0.5)
#' 
#' pred
#' plot(pred)
#' 
#' @export
predict.RBDmodel <- function(model, seq, sens = NULL, FDR = NULL) {
  
  # check
  if( !is.null(model$RBDerr) ){ cat("\nPrediction error estimate included\n\n")
  } else if( !is.null(FDR) || !is.null(sens) ) warning("No prediction error estimate included: FDR and sens are ignored.\n\n")
  if( class(model) != "RBDmodel" ) stop("Argument err must be of class RBDmodel")
  if( !is.null(FDR) && (FDR > 1 || FDR < 0) ){ stop("Please choose reasonable FDR")
  } else if( !is.null(sens) && (sens > 1 || sens < 0) ) stop("Please choose reasonable sensitivtiy")
  
  # feature space of 20 biological relevant amino acids
  feats <- 1:model$RBDparams$k_mer
  feats <- lapply(feats, function(x) c("A", "E", "K", "S", "R", "Q", "L", "T", "N", "G", "M", "W", "D", "H", "F", "Y", "C", "I", "P", "V"))
  feats <- expand.grid(feats, stringsAsFactors = FALSE)
  feats <- sort(apply(feats, 1, function(x) paste(x, collapse = "")))
  
  # extract k-mers from aa sequences
  seq <- toupper(seq)
  kmers <- lapply(seq, function(x) as.vector(GetKMers(x, K = model$RBDparams$k_mer)))
  kmers <- lapply(kmers, function(x) factor(x, levels = feats))
  kmers <- t(sapply(kmers, function(x) as.vector(table(x))))
  colnames(kmers) <- feats
  
  # reduce features
  top <- model$RBDranking[1:model$RBDparams$feats, 1]
  kmers <- kmers[, top]
  kmers <- kmers[, order(colnames(kmers))]
  
  # compute posteriors
  post <- sda::predict.sda(model$RBDmodel, kmers)$posterior

  # out
  pred <- list()
  class(pred) <- "RBDpredict"
  pred$RBDpp <- post[, 2]
  pred$RBDparams <- list(k_mer = model$RBDparams$k_mer, n = model$RBDparams$n, feats = model$RBDparams$feats)
  if( is.null(model$RBDerr) ) return(pred)
  
  # calculate precision, recall
  df <- data.frame(
    rec = model$RBDerr[, 1] / (model$RBDerr[, 1] + model$RBDerr[, 4]),
    pre = model$RBDerr[, 1] / (model$RBDerr[, 1] + model$RBDerr[, 2])
  )
  
  # threshold by sens
  if( !is.null(sens) ){
    type <- "Sensitivity"
    value <- sens
    idx <- which.min(abs(df$rec - sens))
    if( abs(df$rec[idx] - sens) > 0.01 ) warning("Sensetivity set to ", df$rec[idx])
    th <- as.numeric(rownames(model$RBDerr)[idx])
  }
  
  # threshold by FDR
  if( !is.null(FDR) ){
    type <- "FDR"
    value <- FDR
    fdr <- 1 - df$pre
    idx <- which.min(abs(fdr - FDR))
    if( abs(fdr[idx] - FDR) > 0.01 ) warning("FDR set to ", fdr[idx])
    th <- as.numeric(rownames(model$RBDerr)[idx])
  }

  # out
  pred$RBDclass <- as.integer(post[, 2] > as.numeric(th))
  pred$RBDparams$threshold <- th
  pred$RBDparams$by <- type
  pred$RBDparams$set <- value
  pred
}

