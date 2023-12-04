#'
#' Identify hubs in a network
#'
#' Find phenotypes that are hubs in a genomic network including variants,
#' phenotypes, intermediate variables, common children, and confounders.
#'
#' @param Adj numeric (binary) matrix, adjacency matrix of the network. The nodes
#' must be ordered such that the first \code{n_v} columns of \code{Adj} are variants,
#' the next \code{n_t} columns are phenotypes, the next \code{n_q} columns are
#' intermediate variables or common children, and the remaining columns (if any) are confounders.
#'
#' @param n_v,n_t,n_q integers, number of respectively variants, phenotypes, and
#' intermediate variables or common children in the network.
#'
#' @param number.parents,number.children,number.neighbors integers, minimum number
#' of parents (respectively: children, all neighbors including parents and children)
#' required for a phenotype to be a hub.
#'
#' @param all logical, should all the conditions for \code{number.parents},
#' \code{number.children}, and \code{number.neighbors} be simultaneously satisfied
#' for a phenotype to be a hub? The default is \code{all = TRUE}. Setting \code{all = FALSE}
#' means that one condition is sufficient for a phenotype to be a hub.
#'
#' @export find.hubs
#'
#' @return a list with two elements:
#' \describe{
#' \item{hubs}{an integer vector of the positions of the hubs in the sub-network of phenotypes.}
#' \item{number.neighbors}{a six-column matrix of where the \code{i}th row corresponds to the \code{i}th phenotype in the network;
#' and the columns indicate the numbers of parent nodes which are variants, phenotypes,
#' intermediate variables or common children, and confounders; and the numbers of child nodes
#' which are phenotypes, or intermediate variables or common children.}
#' }
#'
#' @examples
#' # Simulate a DAG
#' set.seed(21000)
#' library(MRGNgeneral)
#' data <- sample.graph.data (number.of.T = 100,
#'                            conf.num.vec = c(W = 50, Z = 50, U = 200, K = 0),
#'                            sample.size = 500)
#'
#' # Find hubs
#' Hubs = find.hubs (data$Adjacency,
#'                   n_v = 100, n_t = 100, n_q = 100,
#'                   all = TRUE)
#' Hubs$hubs
#' Hubs$number.neighbors[Hubs$hubs,]
#'
# Check the presence of hub(s) in a network given its adjacency matrix
find.hubs <- function (Adj,
                       n_v = 0,         # Number of variants
                       n_t = NCOL(Adj), # Number of phenotypes
                       n_q = NCOL(Adj) - n_t, # Number of intermediate variables/common children
                       number.parents = 10,   # T-node is a hub if more than 10 parents
                       number.children = 10,  # T-node is a hub if more than 10 children
                       number.neighbors = 15, # T-node is a hub if more than 15 neighbors (parents+children)
                       all = TRUE) {          # Are all conditions required to declare a hub? Default is one condition is sufficient.
  ## Check arguments
  stopifnot(is.numeric(Adj))
  stopifnot(NCOL(Adj) == NROW(Adj))
  stopifnot(all(diag(Adj) == 0))
  stopifnot(is.numeric(n_v), is.numeric(n_t), is.numeric(n_q))
  stopifnot(n_t > 0)
  stopifnot(NCOL(Adj) >= n_v + n_t + n_q)
  all <- as.logical(all[1])
  number.parents <- number.parents[1]
  number.children <- number.children[1]
  number.neighbors <- number.neighbors[1]
  stopifnot(is.numeric(number.parents), is.numeric(number.children), is.numeric(number.neighbors))

  ## Labeld for T-nodes
  Tlabels <- (n_v + 1):(n_v + n_t)

  ## Find for each T-node,
  ## the number of parent V, T, Q, and U-nodes, and the number of T, and Q children
  NbNeighbors <- sapply(X = Tlabels,
                        FUN = find.nb.neighbors,
                        Adj = Adj,
                        n_v = n_v, n_t = n_t, n_q = n_q)
  if (n_t > 1)
    NbNeighbors <- t(NbNeighbors)
  else
    NbNeighbors <- matrix(c(NbNeighbors), nrow = 1, ncol = 6)

  row.names(NbNeighbors) <- colnames(Adj)[Tlabels]
  if (is.null(row.names(NbNeighbors)))
    row.names(NbNeighbors) <- paste0("T", 1:n_t)

  Nbparents <- rowSums(NbNeighbors[,1:4, drop = FALSE])
  Nbchildren <- rowSums(NbNeighbors[,5:6, drop = FALSE])

  ## Find T-nodes which are hubs
  if (all) {
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

  return(list(hubs = Hubs,
              number.neighbors = NbNeighbors))
}

# Find for a T-node with label 'label', the number of parent V, T, Q, and U-nodes, and the number of T, and Q children
# The output is a six column vector: Parent.V, Parent.T, Parent.Q, Parent.U, Child.T, Child.Q
find.nb.neighbors <- function (label, Adj, n_v, n_t, n_q) {
  ## Find all parent nodes of 'label'
  Parents <- which(Adj[, label] > 0)

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
  Children <- which(Adj[label, ] > 0)
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
