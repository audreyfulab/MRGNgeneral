#' Plot an adjacency matrix
#'
#' This function takes an adjacency matrix representing an acyclic graph
#' and calls \link{igraph::plot.igraph} to plot the graph.
#'
#' @param adjacency A square matrix of class \code{adjacency}. This is passed to
#' \link{is.adjacency.matrix} for checking.
#'
#' @param n_vt numeric, the total number of \code{V} and \code{T} nodes
#' of the graph represented by \code{adjacency}. The default (\code{NULL}) is all. When
#' specified, \code{n_vt} must be less than or equal to \code{NCOL(adjacency)}.
#'
#' @param mode Character scalar, specifies how \code{igraph} should interpret the
#' supplied matrix. See \code{graph_from_adjacency_matrix} in the library
#' \code{igraph} for details.
#'
#' @param layout Either a function or a numeric matrix. It specifies how the
#' vertices will be placed on the plot. See \code{igraph.plotting} in the
#' library \code{igraph} for details.
#'
#' @param vertex.color,edge.arrow.size,... Additional plotting parameters.
#' See \code{igraph.plotting} in the library \code{igraph} for the complete list.
#'
#' @usage
#' # S3 method for class 'adjacency.matrix'
#' plot(adjacency, ...)
#'
#' @details
#' The function's output is random.
#' Set a random generator seed using \link{set.seed} for reproducibility.
#'
#'
#' @return Returns \code{NULL}, invisibly.
#'
#' @exportS3Method plot adjacency.matrix
#' @export plot.adjacency.matrix
# @references
#     \insertAllCited{}
#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom igraph plot.igraph
#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom igraph graph_from_adjacency_matrix
# layout_as_bipartite(), layout_as_star(), layout_as_tree(), layout_in_circle(), layout_nicely(), layout_on_grid(), layout_on_sphere(), layout_randomly(), layout_with_dh(), layout_with_fr(), layout_with_gem(), layout_with_graphopt(), layout_with_kk(), layout_with_lgl(), layout_with_mds(), layout_with_sugiyama()
#'
#' @examples
#' ### Load the network 'networka11'
#' library(MRGNgeneral)
#' data(networka11)
#'
#' ### Adjacency matrix of a subset of the network
#' adjacency <- structure(
#'   networka11$adjacency[c('V39', 'T39', 'T43', 'T52', 'W11', 'Z2', 'U28'),
#'                        c('V39', 'T39', 'T43', 'T52', 'W11', 'Z2', 'U28')],
#'   class = 'adjacency.matrix')
#' adjacency
#'
#' ### Plot the graph of the subset
#' plot (adjacency)

plot.adjacency.matrix <- function(adjacency,
                                  n_vt = NULL,
                                  mode = c("directed", "undirected", "max",
                                           "min", "upper", "lower", "plus"),
                                  layout = igraph::layout_nicely,
                                  vertex.color = 'green',
                                  edge.arrow.size = 0.2,
                                  ...) {
  stopifnot(is.adjacency.matrix(adjacency))
  number.of.nodes <- NCOL(adjacency)
  if (is.null(n_vt))
    n_vt <- number.of.nodes
  if (number.of.nodes > n_vt) {
    sub.graph.ind = 1:n_vt
    igraph.obj = igraph::graph_from_adjacency_matrix(adjacency, mode = mode[1])
    par(mfrow = c(1,2))
    igraph::plot.igraph(igraph.obj, layout = layout, edge.arrow.size = edge.arrow.size)
    igraph::plot.igraph(igraph::graph_from_adjacency_matrix(adjacency[sub.graph.ind,sub.graph.ind],
                                                            mode = mode[1]),
                        layout = layout, edge.arrow.size = edge.arrow.size,
                        vertex.color = vertex.color, ...)
  }
  else if (number.of.nodes == n_vt) {
    igraph::plot.igraph(igraph::graph_from_adjacency_matrix(adjacency, mode = mode[1]),
                        layout = layout,
                        edge.arrow.size = edge.arrow.size,
                        vertex.color = vertex.color, ...)
  }
  else {
    stop("inconsistent arguments 'adjacency' and 'n_vt'.")
  }
}
