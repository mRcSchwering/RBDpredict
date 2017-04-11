#' GetKMers
#'
#' Extract all K-mers of a character sequence.
#'
#' @param seq  character as input sequence
#' @param K    int defining length of K-mer
#' 
#' @return character arr of all K-mers found in seq
#' 
#' @examples
#' GetKMers(paste(letters[1:10], collapse = ""))
#' 
GetKMers <- function( seq, K = 3 ) {
  
  SP = strsplit(as.character(seq), split="")
  tmp = sapply(SP, function(x) {
    res = rep("",length(x)-K+1)
    for (i in seq_len(K)) {
      res = paste(res,x[seq_len(length(x)-K+1)+i-1],sep="")
    }
    res
  })
  tmp
}