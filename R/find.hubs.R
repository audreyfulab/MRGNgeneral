#'
#' Identify hubs in a network
#'
#' Find phenotypes that are hubs in a genomic network including variants,
#' phenotypes, intermediate variables, common children, and confounders.
#'
#' @param adjacency numeric (binary) matrix, adjacency matrix of the network. The nodes
#' must be ordered such that the first \code{n_v} columns of \code{adjacency} are variants,
#' the next \code{n_t} columns are phenotypes, the next \code{n_q} columns are
#' intermediate variables or common children, and the remaining columns (if any) are confounders.
#'
#' @param n_v,n_t,n_q integers, number of respectively variants (\code{V}-nodes),
#' phenotypes (\code{T}-nodes), and intermediate variables or common children
#' (\code{Q}-nodes) in the network.
#'
#' @param number.parents,number.children,number.neighbors integers, minimum number
#' of parents (respectively: children, all neighbors including parents and children)
#' required for a phenotype to be a hub. The default values identify a node as a
#' hub when the node has ten or more neighbors, irrespective of the details of the
#' numbers of parents or children.
#'
#' @param only.T logical, should only phenotypes neighbors be considered for hub
#' identification? Defaults to \code{only.T = TRUE}. Setting \code{only.T = FALSE}
#' means that all parent nodes (variants, phenotypes, intermediate variables,
#' common children, and confounders) are considered.
#'
#' @param all.conditions logical, should all the conditions for \code{number.parents},
#' \code{number.children}, and \code{number.neighbors} be simultaneously satisfied
#' for a phenotype to be a hub? The default is \code{all.conditions = FALSE}.
# Setting \code{all.conditions = FALSE} means that one condition is sufficient
# for a phenotype to be a hub.
#'
#' @param cl,chunk.size optional arguments for parallel computing, passed to
#' \link[parallel]{parLapply} (when supplied).
#'
#' @export find.hubs
#'
#' @return An object of class \code{'hubs'}, a list with two elements:
#' \describe{
#' \item{hubs}{an integer vector of the positions of the hubs in the sub-network of phenotypes.}
#' \item{number.neighbors}{a six-column matrix of where the \code{i}th row corresponds to the \code{i}th phenotype in the network;
#' and the columns indicate the numbers of parent nodes which are variants, phenotypes,
#' intermediate variables or common children, and confounders; and the numbers of child nodes
#' which are phenotypes, or intermediate variables or common children.}
#' }
#'
#' The class \code{'hubs'} has two methods: \code{update} to re-find hubs using
#' a new set of parameters (except the \code{adjacency} of the network), and
#' \code{print} to summarize the identified hubs.
#'
#' @examples
#' # Simulate a DAG
#' set.seed(21000)
#' library(MRGNgeneral)
#' data <- sample.graph.data (n_t = 100,
#'                            conf.num.vec = c(W = 50, Z = 50, U = 200, K = 0),
#'                            sample.size = 500)
#'
#' # Find hubs with at least 10 T-neighbors
#' Hubs <- find.hubs (data$adjacency,
#'                    n_v = 100, n_t = 100, n_q = 100)
#' Hubs
#'
#' # Update and print hubs with at least 20 neighbors (any type),
#' #                               including at least 10 parents.
#' Hubs <- print(Hubs,
#'               number.parents = 10,
#'               number.children = 2,
#'               number.neighbors = 20,
#'               only.T = FALSE,
#'               all.conditions = TRUE)
#'
# Check the presence of hub(s) in a network given its adjacency matrix
find.hubs <- function (adjacency,
                       n_v = 0,         # Number of variants
                       n_t = NCOL(adjacency), # Number of phenotypes
                       n_q = NCOL(adjacency) - (n_t + n_v), # Number of intermediate variables/common children
                       number.parents = if (all.conditions) 0 else Inf,   # T-node is a hub if more than 'number.parents' parents
                       number.children = if (all.conditions) 0 else Inf,  # T-node is a hub if more than 'number.children' children
                       number.neighbors = 10, # T-node is a hub if more than 10 neighbors (parents+children)
                       only.T = TRUE,
                       all.conditions = FALSE,          # Are all conditions required to declare a hub? Default is one condition is sufficient.
                       cl = NULL, chunk.size = NULL) {
  ## Save the call, not save arguments
  mcall <- match.call()
  mcall[c('n_v', 'n_t', 'n_q',
          'number.parents', 'number.children', 'number.neighbors',
          'only.T', 'all.conditions')] <- list(n_v, n_t, n_q,
                                               number.parents,
                                               number.children,
                                               number.neighbors,
                                               only.T, all.conditions)

  ## Check arguments
  stopifnot(is.numeric(adjacency))
  stopifnot(NCOL(adjacency) == NROW(adjacency))
  stopifnot(all(adjacency %in% c(0, 1)))
  stopifnot(all(diag(adjacency) == 0))
  stopifnot(is.numeric(n_v), is.numeric(n_t), is.numeric(n_q))
  stopifnot(n_t > 0)
  stopifnot(NCOL(adjacency) >= n_v + n_t + n_q)
  n_v <- n_v[1]
  n_t <- n_t[1]
  n_q <- n_q[1]
  all.conditions <- as.logical(all.conditions[1])
  number.parents <- number.parents[1]
  number.children <- number.children[1]
  number.neighbors <- number.neighbors[1]
  stopifnot(is.numeric(number.parents), is.numeric(number.children), is.numeric(number.neighbors))

  ## Labels for T-nodes
  Tlabels <- (n_v + 1):(n_v + n_t)

  ## Find for each T-node,
  ## the number of parent V, T, Q, and U-nodes, and the number of T, and Q children
  NbNeighbors <- sapply(X = Tlabels,
                        FUN = find.nb.neighbors,
                        adjacency = adjacency,
                        n_v = n_v, n_t = n_t, n_q = n_q)
  if (n_t > 1)
    NbNeighbors <- t(NbNeighbors)
  else
    NbNeighbors <- matrix(c(NbNeighbors), nrow = 1, ncol = 6)

  row.names(NbNeighbors) <- colnames(adjacency)[Tlabels]
  if (is.null(row.names(NbNeighbors)))
    row.names(NbNeighbors) <- paste0("T", 1:n_t)

  if (only.T) {
    Nbparents <- NbNeighbors[,2]
    Nbchildren <- NbNeighbors[,5]
  }
  else {
    Nbparents <- rowSums(NbNeighbors[,1:4, drop = FALSE])
    Nbchildren <- rowSums(NbNeighbors[,5:6, drop = FALSE])
  }

  ## Find T-nodes which are hubs
  if (all.conditions) {
    Hubs <- ((Nbparents >= number.parents) & (Nbchildren >= number.children)) &
      ((Nbparents + Nbchildren) >= number.neighbors)
  }
  else {
    Hubs <- ((Nbparents >= number.parents) | (Nbchildren >= number.children)) |
      ((Nbparents + Nbchildren) >= number.neighbors)
  }

  Hubs <- which(Hubs)
  #if (length(Hubs))
  #  Hubs <- Hubs + n_v

  return(structure(list(hubs = Hubs,
                        number.neighbors = NbNeighbors,
                        call = mcall),
                   class = 'hubs'))
}

# Find for a T-node with label 'label', the number of parent V, T, Q, and U-nodes, and the number of T, and Q children
# The output is a six column vector: Parent.V, Parent.T, Parent.Q, Parent.U, Child.T, Child.Q
find.nb.neighbors <- function (label, adjacency, n_v, n_t, n_q) {
  ## Find all parent nodes of 'label'
  Parents <- which(adjacency[, label] > 0)

  ## Find the number of each type of parent node
  Parents <- if (length(Parents)) {
    c(Parent.V = length(which(Parents <= n_v)),
      Parent.T = length(which((Parents > n_v) & (Parents <= (n_v + n_t)))),
      Parent.Q = length(which((Parents > (n_v + n_t)) & (Parents <= (n_v + n_t + n_q)))),
      Parent.U = length(which(Parents > (n_v + n_t + n_q))))
  }
  else {
    c(Parent.V = 0, Parent.T = 0, Parent.Q = 0, Parent.U = 0)
  }


  ## Find all child nodes
  Children <- which(adjacency[label, ] > 0)
  ## Find the number of each type of child node
  Children <- if (length(Children)) {
    c(Child.T = length(which((Children > n_v) & (Children <= (n_v + n_t)))),
      Child.Q = length(which((Children > (n_v + n_t)) & (Children <= (n_v + n_t + n_q)))))
  }
  else {
    c(Child.T = 0, Child.Q = 0)
  }

  return(c(Parents, Children))
}

# A print method for 'hubs' class objects
# ... are arguments to the generic 'print'
#' @exportS3Method print hubs
#'
print.hubs <- function(x,
                       number.parents,
                       number.children,
                       number.neighbors,
                       only.T,
                       all.conditions,
                       ...) {
  # Update x if required
  if (!all(c(missing(number.parents),
            missing(number.children),
            missing(number.neighbors),
            missing(only.T),
            missing(all.conditions)))) {
    x <- update.hubs(x,
                     number.parents = number.parents,
                     number.children = number.children,
                     number.neighbors = number.neighbors,
                     only.T = only.T,
                     all.conditions = all.conditions)
  }

  # Print 'x'
  cat("\n")
  if (x$call$only.T) {
    print(x$number.neighbors[x$hubs, c(2, 5)], ...)
  }
  else {
    print(x$number.neighbors[x$hubs,], ...)
  }
  cat("\n")
  return(invisible(x))
}

# An update method for 'hubs' class objects
#' @exportS3Method update hubs
#'
update.hubs <- function(object,
                        number.parents,
                        number.children,
                        number.neighbors,
                        only.T,
                        all.conditions,
                        ...) {
  ### Return the result in 'object' if no new argument supplied
  if (all(c(missing(number.parents),
            missing(number.children),
            missing(number.neighbors),
            missing(only.T),
            missing(all.conditions)))) {
    return(invisible(object))
  }

  ### Check supplied arguments or set default
  # number.parents
  if (missing(number.parents)) {
    ## Extract argument 'number.parents'  from 'object'
    number.parents <- object$call$number.parents
    if (is.null(number.parents)) {
      # Try to obtain 'number.parents' from the formal arguments of find.hubs
      number.parents <- formals(find.hubs)$number.parents

      # If still no value found, use 10, the default value at the time I am writing this.
      if (!is.numeric(number.parents))
        number.parents <- Inf
    }
  }
  else {
    stopifnot(is.numeric(number.parents))
    number.parents <- number.parents[1]
  }

  # number.children
  if (missing(number.children)) {
    ## Extract argument 'number.children'  from 'object'
    number.children <- object$call$number.children
    if (is.null(number.children)) {
      # Try to obtain 'number.children' from the formal arguments of find.hubs
      number.children <- formals(find.hubs)$number.children

      # If still no value found, use 10, the default value at the time I am writing this.
      if (!is.numeric(number.children))
        number.children <- Inf
    }
  }
  else {
    stopifnot(is.numeric(number.children))
    number.children <- number.children[1]
  }

  # number.neighbors
  if (missing(number.neighbors)) {
    ## Extract argument 'number.neighbors'  from 'object'
    number.neighbors <- object$call$number.neighbors
    if (is.null(number.neighbors)) {
      # Try to obtain 'number.neighbors' from the formal arguments of find.hubs
      number.neighbors <- formals(find.hubs)$number.neighbors

      # If still no value found, use 10, the default value at the time I am writing this.
      if (!is.numeric(number.neighbors))
        number.neighbors <- 10
    }
  }
  else {
    stopifnot(is.numeric(number.neighbors))
    number.neighbors <- number.neighbors[1]
  }

  # only.T
  if (missing(only.T)) {
    ## Extract argument 'only.T'  from 'object'
    only.T <- object$call$only.T
    # If 'object$call$only.T' is NULL, then the default of 'find.hubs' was used
    if (is.null(only.T)) {
      # Try to obtain only.T from the formal arguments of find.hubs
      only.T <- formals(find.hubs)$only.T

      # If still no value found, use TRUE, the default value of only.T at the time I am writing this.
      if (!is.logical(only.T))
        only.T <- TRUE
    }
  }
  else {
    stopifnot(is.logical(only.T) | is.numeric(only.T))
    only.T <- only.T[1] == 1L
  }

  # all.conditions
  if (missing(all.conditions)) {
    ## Extract argument 'all.conditions'  from 'object'
    all.conditions <- object$call$all.conditions
    if (is.null(all.conditions)) {
      # Try to obtain 'all.conditions' from the formal arguments of find.hubs
      all.conditions <- formals(find.hubs)$all.conditions

      # If still no value found, use FALSE, the default value at the time I am writing this.
      if (!is.logical(all.conditions))
        all.conditions <- FALSE
    }
  }
  else {
    stopifnot(is.logical(all.conditions) | is.numeric(all.conditions))
    all.conditions <- all.conditions[1] == 1L
  }

  ### Matrix of neighborhood information (V, T, Q, U parents and T,Q children)
  NbNeighbors <- object$number.neighbors

  ### Find number of neighbors to be considered for each T-node
  if (only.T) {
    Nbparents <- NbNeighbors[,2]
    Nbchildren <- NbNeighbors[,5]
  }
  else {
    Nbparents <- rowSums(NbNeighbors[,1:4, drop = FALSE])
    Nbchildren <- rowSums(NbNeighbors[,5:6, drop = FALSE])
  }

  ### Find T-nodes which are hubs
  if (all.conditions) {
    Hubs <- ((Nbparents >= number.parents) & (Nbchildren >= number.children)) &
      ((Nbparents + Nbchildren) >= number.neighbors)
  }
  else {
    Hubs <- ((Nbparents >= number.parents) | (Nbchildren >= number.children)) |
      ((Nbparents + Nbchildren) >= number.neighbors)
  }

  # Update object
  object$hubs <- which(Hubs)
  object$call[c('number.parents', 'number.children', 'number.neighbors',
                 'only.T', 'all.conditions')] <- list(number.parents,
                                                      number.children,
                                                      number.neighbors,
                                                      only.T, all.conditions)

  return(object)
}
