#'
#' From \code{bn} class objects to adjacency matrix
#'
#' \code{bnclass2adj} converts a list of edges of a \code{bn} class object (e.g.
#' from the package \code{bnlearn}) into an adjacency matrix.
#'
#' @param bnfit An object of class \code{bn} (e.g., from the \code{bnlearn} package)
#' @param pattern Character string used as a prefix pattern in node names. Default is "G"
#' @param nodenames Optional character vector of node names. If NULL (default), 
#'   names are generated as pattern + node number (e.g., "G1", "G2", etc.)
#'
#' @return An adjacency matrix with row and column names set to \code{nodenames}
#'
#' @export bnclass2adj

bnclass2adj <- function(bnfit, pattern = "G", nodenames = NULL) {
  if (inherits(bnfit, what = "bn")) {
  Adj <- cbind(as.numeric(sapply(bnfit$arcs[,1], FUN = function(x) {
    stringi::stri_split_fixed(x, pattern = pattern)[[1]][2]})),
    as.numeric(sapply(bnfit$arcs[,2], FUN = function(x) {
      stringi::stri_split_fixed(x, pattern = pattern)[[1]][2]})))
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
