#' ROCplot
#' 
#' Plotting a Receiver Operating Characteristic curve for RBDerr object. 
#' 
#' This is a common plot to evaluate prediction power of a predictor.
#' Correct and incorrect classifications of RBPs and NoRBPs are assessed equally.
#' Set argument sens to see Sensitivity and Specificity of at specific stringencies.
#' 
#' For plotting multiple ROC curves you can provide a list of RBDerr objects as 1st argument.
#' You must provide names for each curve with the names argument.
#' If many (>10) curves are plotted use numerical values as names. 
#' This way the legend will be a continuum.
#'
#' @param err    RBDerr object or list of RBDerr objects
#' @param sens   num(= NULL) define a sensitivity value to show as horizontal line in plot
#' @param data   bool(= FALSE) whether to just return the data.frame with precision, recall, fallout without plotting 
#' @param names  chr arr(= NULL) names for curves in case multiple curves are plotted. names is mapped to point colors.
#' @param legend chr(= "model") legend title in case multiple ROC curves are plotted, NULL for no legend, "" for no legend title
#' @param title  chr(= "ROC") title of plot, NULL for no title
#'
#' @return plot object (ggplot) or data.frame containing fallout, recall, precision
#' 
#' @examples
#' l <- c("A", "E", "K", "S", "R", "Q", "L", "T", "N", "G", "M", "W", "D", "H", "F", "Y", "C", "I", "P", "V")
#' sequences <- sapply(sample(10:100, 500, replace = TRUE), function(x) paste(sample(l, x, replace = TRUE), collapse = ""))
#' labels <- c(rep(1, 200), rep(0, 300))
#' 
#' model <- RBDfit(sequences, labels, k = 2)
#' err <- cvRBDfit(model, sequences, labels, nfold = 3)
#' ROCplot(err, sens = 0.5)
#' 
#' d <- ROCplot(err, data = TRUE)
#' summary(d)
#' plot(d$fal, d$rec)
#' 
#' err2 <- cvRBDfit(model, sequences, labels, nfold = 5)
#' ROCplot(list(err, err2), legend = "cross validation", names = c("3fold", "5fold") sens = 0.5)
#' 
#' @export
ROCplot <- function(err, sens = NULL, data = FALSE, names = NULL, legend = "model", title = "ROC"){
  
  # check
  if( class(err) != "RBDerr" && class(err) != "list" ) stop("err must be a RBDerr object or a named list of RBDerr objects")
  if( class(err) == "list" ){
    if( any(sapply(err, class) != "RBDerr") ) stop("If err is a list, it must contain RBDerr objects")
    if( is.null(names) ) stop("If err is a list of RBDerr objects, you must provide its names (names argument)")
    if( length(names) != length(err) ) stop("The number of names you provided differs from the number of elements in err")
  }
  
  # calculate recall, fallout, precision
  df <- switch(class(err),
    RBDerr = as.data.frame(err$RBDerr, stringsAsFactors = FALSE),
    list = as.data.frame(do.call(rbind, lapply(err, function(x) x$RBDerr)), stringsAsFactors = FALSE)
  )
  df <- data.frame( 
    rec = df$TP / (df$TP + df$FN), 
    fal = df$FP / (df$FP + df$TN), 
    pre = df$TP / (df$TP + df$FP)
  )
  
  if( class(err) == "list" ){
    n <- unlist(lapply(1:length(names), function(x) rep(names[x], nrow(err[[x]]$RBDerr)))) 
    df$name <- if( !any(grepl("[[:alpha:]]", n)) ) as.numeric(n) else factor(n, levels = names)
  }
  if( data ) return(df)

  # create plot
  p <- switch(class(err),
    RBDerr = ggplot2::ggplot(df, ggplot2::aes(x = fal, y = rec)),
    list = ggplot2::ggplot(df, ggplot2::aes(x = fal, y = rec, color = name)) + ggplot2::guides(color = ggplot2::guide_legend(title = legend))
  )
  p <- p + ggplot2::geom_point(size = 0.5) + ggplot2::xlab("fallout") + ggplot2::ylab("recall") +
    ggplot2::theme_classic() + ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))
  
  # sens
  if( !is.null(sens) ){
    idx <- which.min(abs(df$rec - sens))
    if( abs(df$rec[idx] - sens) > 0.01 ) warning("Sensetivity set to ", df$rec[idx])
    title <- paste(title, sep = "\n", paste("with sensitivity of", round(df$rec[idx], 2), "precision of", round(df$pre[idx], 2), "is reached"))
    p <- p + ggplot2::geom_hline(yintercept = df$rec[idx], color = "#db4437", linetype = 2)
  }
  
  # title
  if(!is.null(title)) p <- p + ggplot2::ggtitle(title)
  
  # rm legend
  if(is.null(legend)) p <- p + ggplot2::theme(legend.position = "none")
  
  p
}




