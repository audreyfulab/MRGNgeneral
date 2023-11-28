#' Plot adjacency matrices
#'
#' This function takes any adjacency matrix representing an acyclic graph
#' and calls \link{igraph::plot.igraph} to plot the graph.
#'
#' @param Adj A square matrix of class \code{adjacency}. This is passed to
#' \link{is.adjacency.matrix} for checking.
#'
#' @param number.of.V.T numeric, the total number of \code{V} and \code{T} nodes
#' of the graph represented by \code{Adj}. The default (\code{NULL}) is all. When
#' specified, \code{number.of.V.T} must be less than or equal to \code{NCOL(Adj)}.
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
#' @return Returns \code{NULL}, invisibly.
#'
#' @exportS3Method plot adjacency.matrix
# @references
#     \insertAllCited{}
#' @import igraph
# layout_as_bipartite(), layout_as_star(), layout_as_tree(), layout_in_circle(), layout_nicely(), layout_on_grid(), layout_on_sphere(), layout_randomly(), layout_with_dh(), layout_with_fr(), layout_with_gem(), layout_with_graphopt(), layout_with_kk(), layout_with_lgl(), layout_with_mds(), layout_with_sugiyama()

plot.adjacency.matrix <- function(Adj,
                                  number.of.V.T = NULL,
                                  mode = c("directed", "undirected", "max",
                                           "min", "upper", "lower", "plus"),
                                  layout = igraph::layout_nicely,
                                  vertex.color = 'green',
                                  edge.arrow.size = 0.2,
                                  ...) {
  stopifnot(is.adjacency.matrix(Adj))
  number.of.nodes <- NCOL(Adj)
  if (is.null(number.of.V.T))
    number.of.V.T <- number.of.nodes
  if (number.of.nodes > number.of.V.T) {
    sub.graph.ind = 1:number.of.V.T
    igraph.obj = igraph::graph_from_adjacency_matrix(Adj, mode = mode[1])
    par(mfrow = c(1,2))
    igraph::plot.igraph(igraph.obj, layout = layout, edge.arrow.size = edge.arrow.size)
    igraph::plot.igraph(igraph::graph_from_adjacency_matrix(Adj[sub.graph.ind,sub.graph.ind],
                                                            mode = mode[1]),
                        layout = layout, edge.arrow.size = edge.arrow.size,
                        vertex.color = vertex.color, ...)
  }
  else if (number.of.nodes == number.of.V.T) {
    igraph::plot.igraph(igraph::graph_from_adjacency_matrix(Adj, mode = mode[1]),
                        layout = layout,
                        edge.arrow.size = edge.arrow.size,
                        vertex.color = vertex.color, ...)
  }
  else {
    stop("inconsistent arguments 'Adj' and 'number.of.V.T'.")
  }
}
