#' From \code{bn} class objects to adjacency matrix
#'
#' \code{bnclass2adj} converts a list of edges of a \code{bn} class object (e.g.
#' from the package \code{bnlearn}) into an adjacency matrix.
#'
#' @export bnclass2adj

bnclass2adj <- function(bnfit, pattern = "G", nodenames = NULL) {
  if (inherits(bnfit, what = "bn")) {
  Adj <- cbind(as.numeric(sapply(bnfit$arcs[,1], FUN = function(x) {
    strsplit(x, split = pattern)[[1]][2]})),
    as.numeric(sapply(bnfit$arcs[,2], FUN = function(x) {
      strsplit(x, split = pattern)[[1]][2]})))
  Adj <- edges2adj(Adj)
  if (is.null(nodenames)) {
    nodenames <- paste0(pattern, 1:NROW(Adj))
  }
  else {
    stopifnot(length(nodenames) == NROW(Adj))
  }
  dimnames(Adj) <- list(nodenames, nodenames)
  return(Adj)
  }
  else {
    stop("only objects of class 'bn' are currently handled.")
  }
}
