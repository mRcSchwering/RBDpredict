#' RBDfit
#' 
#' Create a model for RBP prediction based on protein sequences and labels.
#' 
#' The underlying mechanism is a shrinkage discrimminant analysis.
#' As features, the occurences of peptide combinations in each protein sequence are regarded.
#' k defines the length of these peptides. 
#' With n the number of features used in the model can be defined.
#' 
#' With small k only primitive structures can be described. With big k however the feature space becomes huge.
#' k = 3 seems to be a good trade-off with a feature space of 8000 (20^3).
#' 
#' First, the whole feature space is created.
#' Then features are ranked based on the data provided.
#' How many features the model will actually use can be defined with n.
#' With n = "all" all features are used.
#' However, most of the time it is useful to only use a subset of features.
#'
#' @param seq       chr arr containing amino acid sequences of proteins
#' @param labels    int arr containing labels: 0 = NoRBP, 1 = RBP
#' @param k         int (= 3) defining the length of peptides to be extracted
#' @param n         chr (= "Q1") int or num defining the number of features to use.
#'                  chr can be "Q1", "Q2", "Q3" for 1st 1, 2, 3 quartiles or "all" for all features.
#'                  num can be anything from 0 to 1 defining quantile of features to use.
#'                  int anything > 1 to define plain number of features to use.
#' @param model     bool (= FALSE) whether model frame should be returned as well
#' 
#' @return list of class RBPmodel 5 elements 
#' \itemize{
#'     \item RBDmodel sda object of the actual fit 
#'     \item RBDranking matrix containing information about the feature ranking
#'     \item RBDparams arguments RBPfit was called with
#'  }
#'  If model=TRUE 2 more elements are returned
#'  \itemize{
#'     \item RBDfeatureExtraction matrix: columns = features, rows = proteins, values = occurences
#'     \item RBDlabels int arr of original protein labels
#'  }
#' 
#' @examples
#' l <- c("A", "E", "K", "S", "R", "Q", "L", "T", "N", "G", "M", "W", "D", "H", "F", "Y", "C", "I", "P", "V")
#' sequences <- sapply(sample(10:100, 500, replace = TRUE), function(x) paste(sample(l, x, replace = TRUE), collapse = ""))
#' labels <- c(rep(1, 200), rep(0, 300))
#' 
#' model <- RBDfit(sequences, labels, k = 2)
#' str(model)
#' model
#' plot(model)
#' 
#' @export
RBDfit <- function(seq, labels, k = 3, n = "Q1", model = FALSE) {
  
  # check
  if( !identical(labels, as.numeric(as.logical(labels)))) stop("Argument labels must be of 0 and 1" )
  
  # feature space of 20 biological relevant amino acids
  feats <- 1:k
  feats <- lapply(feats, function(x) c("A", "E", "K", "S", "R", "Q", "L", "T", "N", "G", "M", "W", "D", "H", "F", "Y", "C", "I", "P", "V"))
  feats <- expand.grid(feats, stringsAsFactors = FALSE)
  feats <- sort(apply(feats, 1, function(x) paste(x, collapse = "")))
  
  # chose number of top features to be used
  if( class(n) == "character" ){
    nFeats <- floor(quantile(1:length(feats), switch(n, Q1 = 0.25, Q2 = 0.5, Q3 = 0.75, Q4 = , all = 1, stop("argument n invalid:", n))))
  } else if( n <= 1 && n >= 0 ){ nFeats <- floor(quantile(1:length(feats), n)) 
  } else nFeats <- if( n > length(feats) ) length(feats) else n

  # extract k-mers from aa sequences
  seq <- toupper(seq)
  kmers <- lapply(seq, function(x) as.vector(GetKMers(x, K = k)))
  kmers <- lapply(kmers, function(x) factor(x, levels = feats))
  kmers <- t(sapply(kmers, function(x) as.vector(table(x))))
  colnames(kmers) <- feats

  # sda ranking
  cat("\n\n\nRanking features\n\n")
  ranks <- sda::sda.ranking(kmers, factor(labels, labels = c("NoRBP", "RBP"), levels = c(0, 1)))
  top <- ranks[1:nFeats, 1]
  reduced <- kmers[, top]
  reduced <- reduced[, order(colnames(reduced))]
  
  # the actual fit
  cat("\n\n\nTraining a classifier\n\n")
  fit <- sda::sda(reduced, factor(labels, labels = c("NoRBP", "RBP"), levels = c(0, 1)))
  
  # out
  out = list()
  out$RBDmodel <- fit
  out$RBDranking <- ranks
  out$RBDparams <- list(k_mer = k, n = n, feats = nFeats, type = "proteins")
  if( model ){
    out$RBDfeatureExtraction <- kmers
    out$RBDlabels <- labels
  }
  class(out) <- "RBDmodel"
  out
}