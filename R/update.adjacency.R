#' @name update.adjacency.matrix
#' @title Update an adjacency matrix
#'
#' @description This function updates the adjacency matrix of a graph based on
#' trio analysis. It is used to perform *Step* \code{3.3}, *Step* \code{3.6},
#' *Step* \code{4.3} and *Step* \code{4.6} of the \code{MRGN} algorithm.
#' This is not a user level function.
#'
#' @param adjacency numeric, an adjacency matrix
#'
#' @param n_v,n_t integer scalars, number of genetic variants (\code{n_v}), and
#' target phenotypes (\code{n_t}) in \code{adjacency}.
# Not required when \code{solve.conflicts = FALSE}.
#' @param n_q integer scalar, number of intermediate variables/common children 
#'   (\code{Q}-nodes) in \code{adjacency}. Currently not used
#'
#'
#' @param trio.set numeric, a matrix where each row indicates the column numbers
#' of variables forming a trio.
#'
#' A 3-column matrix \code{trio.set} corresponds to trios involving genetic
#' variant(s), i.e. each row of \code{trio.set} indicates the column number
#' of a genetic variant (first element), and the column numbers of two
#' expressions/phenotypes (second and third columns).
#'
#' A 4-column matrix \code{trio.set} corresponds to trios involving only
#' expressions/phenotypes, no genetic variant, i.e. each row of \code{trio.set}
#' indicates the column numbers of three phenotypes (\code{T} or \code{Q}-node),
#' and the type of triplet (last column). See \link{enumerate.triplets} for
#' details on triplet types.
#'
#' @param inferred.models a character vector of length the number of rows in
#' \code{trio.set}. Each elements of \code{inferred.models} must be one of
#' \code{M0.1}, \code{M0.2}, \code{M1.1}, \code{M1.2}, \code{M2.1}, \code{M2.2},
#' \code{M3}, \code{M4}, or \code{Other}. See Badsha and Fu (2019) for the
#' definitions of these model topologies.
#'
#' @param stringent logical, should edges absent from trio topologies \code{M1.1},
#' \code{M1.2}, \code{M2.1}, and \code{M2.2} be considered absent from the graph?
#' @param add.edges logical, should edges be added to the adjacency matrix? 
#'   Default is TRUE
#' @param solve.conflicts logical, should a resolution of conflicts be attempted?
#' If \code{FALSE}, edges are updated sequentially and conflicts are not noticed.
#' This makes the final inference depends on the order in which trios are analysed.
#'
#' If \code{TRUE} (the default), all edges are updated once. This put conflicts
#' into evidence and requires a method to deal with each type of errors.
#'
#' @param method character, only used if \code{solve.conflicts = TRUE}.
#' The method to be used to solve conflicts. The available \code{method} are:
#'
#' \code{conservative}: conflicts about the presence of an edge are solve by including
#' an edge when at least one trio analysis inferred an edge; conflicts about an
#' edge direction are solve by letting the edge undirected.
#' \code{add_methods}: add further method(s).
#'
#' @param added.edges character vector indicating all edges (and directions)
#' previously added into the network. Only used when \code{solve.conflicts = TRUE}.
#'
#' @param dropped.edges character vector indicating all edges (and directions)
#' previously dropped from the network (or inferred as non existent). Only used
#' when \code{solve.conflicts = TRUE}.
#'
#' @param cl a cluster object, created by one of the packages \code{parallel} and \code{snow}.
#' If \code{NULL}, the registered default cluster is used. Note that the latter can be
#' \code{NULL} too, in which case, no parallel computation is performed.
#' @param chunk.size integer, chunk size for parallel computation when \code{cl} is not NULL.
#'   If NULL, a default chunk size is used
#' @param ... additional arguments passed to or from other methods.
#' Currently none.
#'
#' @details \code{update.adjacency.matrix} takes an adjacency matrix, a set
#' \code{trio.set} of trios involving genetic variant(s) or not, and the underlying
#' structures under the extended Mendelian Randomization Principle, and update
#' the adjacency matrix accordingly.
#'
#' The registered default cluster is found using \code{parallel::getDefaultCluster()}.
#'
# Internal function to direct graph edges
update.adjacency.matrix <- function (adjacency,
                                     n_t,n_v,n_q, # Currently not used
                                     trio.set,
                                     inferred.models,
                                     stringent = FALSE,
                                     add.edges = TRUE,
                                     solve.conflicts = TRUE,
                                     method = "conservative",
                                     added.edges = NULL,
                                     dropped.edges = NULL,
                                     cl = NULL,
                                     chunk.size = NULL,
                                     ...) {
  # Check if 'adjacency' is a valid adjacency matrix and 'trio.set' is a matrix
  stopifnot(is.adjacency.matrix(adjacency), is.matrix(trio.set))
  # Un-comment the following lines if pushed to a user level function
  # stopifnot(is.character(inferred.models))
  # stopifnot(NROW(trio.set) != length(inferred.models))

  # Set default cluster
  if (is.null(cl))
    cl <- parallel::getDefaultCluster()

  # If 'trio.set' has four columns, then it contains triplets of only T-nodes
  if (NCOL(trio.set) == 4) {
    # In this case, the only allowed values in 'inferred.models'
    # are "M1.1", "M2.1", and "Other".
    # Un-comment the following line if pushed to a user level function
    # stopifnot(all(inferred.models %in% c("M1.1", "M2.1", "Other")))

    # Use the child function for triplets
    return(
      update_adjacency_matrix_triplet (adjacency = adjacency,
                                       n_t = n_t, n_v = n_v, n_q = n_q, # Currently not used
                                       triplet.set = trio.set,
                                       inferred.models = inferred.models,
                                       add.edges = add.edges,
                                       solve.conflicts = solve.conflicts[1],
                                       method = method[1],
                                       added.edges = added.edges,
                                       dropped.edges = dropped.edges,
                                       cl = cl, chunk.size = chunk.size,
                                       ...)
    )
  }
  # Un-comment the following lines if pushed to a user level function
  # stopifnot(NCOL(trio.set) == 3) # Ensure 'trio.set' has only three columns
  # stopifnot(all(inferred.models %in% c('M0.1', 'M0.2', 'M1.1', 'M1.2',
  #                                   'M2.1', 'M2.2', 'M3', 'M4', 'Other')))

  # When 'inferred.models' is from trio analysis
  return(
    update_adjacency_matrix_trio (adjacency = adjacency,
                                  n_t = n_t, n_v = n_v, n_q = n_q,
                                  trio.set = trio.set,
                                  inferred.models = inferred.models,
                                  stringent = stringent,
                                  add.edges = add.edges,
                                  added.edges = added.edges,
                                  dropped.edges = dropped.edges,
                                  method = method[1],
                                  cl = cl,
                                  chunk.size = chunk.size,
                                  solve.conflicts = solve.conflicts,
                                  ...)
  )
}
