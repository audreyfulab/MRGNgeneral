#'
#' @name enumerate.trios
#' @title Form trios involving genetic variants
#' @description Enumerate all trios involving some target genetic variants based on
#' an adjacency matrix.
#' Not a user level function.
#'
#' @export enumerate.trios
#'
#' @param i integer vector, position of the target \code{V}-nodes (genetic
#' variants) in the adjacency matrix \code{adjacency}: \code{i} must satisfy
#' \code{1 <= i <= n_v} (see the argument \code{n_v}).
# For convenience,
# \code{i} is coerced to the integer type using the function \link{as.integer}.
#' Defaults to \code{i = 1:n_v} (i.e. all genetic variants).
#'
#' @param adjacency numeric, an adjacency matrix, i.e., a binary square matrix
#' with zeros on the main diagonal.
#' The first \code{n_v} columns of \code{adjacency} represent
#' \code{V}-nodes, the next \code{n_t} columns represent \code{T}-nodes
#' (expressions/phenotypes), and the last \code{n_q} columns represent
#' \code{Q}-nodes (intermediate variables and common children).
#' Only the upper triangular part of \code{adjacency} is actually used.
#'
#' @param n_t,n_v,n_q integers, respectively: number of \code{T}-nodes (non-instrumental variables),
#' number of \code{V}-nodes (instrumental variables), and number of \code{Q}-nodes
#' (confounding variables: intermediate variables and common children).
#'
#' @param VTT logical indicating the composition of desired trios. If true (the default),
#' only trios not involving a \code{Q}-node are formed, i.e. a trio contains one
#' \code{V}-node and two \code{T}-nodes. Otherwise, only trios involving
#' \code{Q}-nodes are formed, i.e.  a trio contains one \code{V}-node,
#' one \code{T}-node, and one \code{Q}-node.
#'
#' @param cl a cluster object, created by one of the packages \code{parallel} and \code{snow}.
#' Specify \code{cl = NULL} for no parallel computation. Note that this is not the
#' genuine default. If missing, the registered default cluster (which can be
#' \code{NULL}, but may not be) is considered.
#'
#' @param chunk.size integer, number of tasks per scheduling unit during parallel computation.
#'
# enumerate.trios (i = 1:n_v,
#                  adjacency,
#                  n_t = NROW(adjacency) - n_v,
#                  n_v = length(i),
#                  n_q = NROW(adjacency) - n_t - n_v,
#                  VTT = TRUE,
#                  cl, chunk.size = NULL)
#
#' @details When \code{i} has more than one elements, trio enumeration is performed
#' for each of its unique elements (as returned by \link{unique}).
#'
#' In a trio, the variant \code{Vi} has an edge with at least one of the two
#' other variables in the trio. In addition to the variant, any trio has either
#' two different \code{T}-nodes, or one \code{T}-node and one \code{Q}-node: a
#' trio with two different \code{Q}-nodes is not allowed.
#'
#' @return A list of matrices, one for each unique element of \code{i}. If the
#' input adjacency matrix \code{adjacency} has non-\code{NULL} colnames (as returned
#' by \link{colnames}), then the elements of the list are named after the
#' corresponding columns (from the first \code{n_v} columns). Otherwise, they are
#' named as 'Vi' where 'i' are the unique elements of \code{i}.
#'
#' For each matrix, each row gives the column numbers of variables forming a
#' trio, the first giving the genetic variant (\code{i}), and the last two
#' giving phenotypes, named 'Tj', 'Tk' when \code{VTT = TRUE}, and named 'Tj',
#' 'Qk' when \code{VTT = FALSE}. If no trio can be formed for a given genetic
#' variant, the associated matrix has zero row.
#'
#' @seealso
#' \link{enumerate.new.trios} to enumerate additional trios after updating the
#' adjacency matrix \code{adjacency}.
#'
#' \link{enumerate.triplets} to enumerate trios involving no genetic
#' variant (only \code{T} and \code{Q}-nodes).
#
####################################################
# List all trios involving a target genetic variant
enumerate.trios <- function (i = 1:n_v,
                             adjacency,
                             n_t = NROW(adjacency) - n_v,
                             n_v = length(i),
                             n_q = NROW(adjacency) - n_t - n_v,
                             VTT = TRUE,
                             cl, chunk.size = NULL) {
  # Check arguments
  if (missing(i) & missing(n_v))
    stop("One of argument 'i' and 'n_v' must be specified")

  if (!missing(i)) {
    if (length(i) < 1)
      stop("'i' must be a non-NULL vector of integers.")
    i <- as.integer(unique(i))
    if (any(c(i < 1, i > n_v)))
      stop("The vector 'i' must satisfy '1 ≤ i ≤ n_v' element-wise.")
  }

  stopifnot(is.logical(VTT[1]))
  if (!VTT[1] & n_q == 0) {
    stop("'n_q > 0' is required when 'VTT = FALSE'")
  }

  # Only use the upper triangular part of adjacency
  adjacency[lower.tri(adjacency)] <- adjacency[lower.tri(t(adjacency))]

  # Labels (numbers) for all T-nodes
  n_vt <- n_v + n_t
  Tlabels <- (1 + n_v):n_vt

  # Labels (numbers) for all Q-nodes
  Qlabels <- if (!VTT[1]) (1 + n_vt):(n_vt + n_q)

  # Wrap 'enumerate.trios.i' over 'i'
  if (missing(cl))
    cl <- parallel::getDefaultCluster()
  res <- matteLapply(i,
                     FUN = enumerate.trios.i,
                     adjacency = adjacency, n_t = n_t,
                     Tlabels = Tlabels,
                     Qlabels = Qlabels,
                     cl = cl, chunk.size = chunk.size)

  # Name and return the obtained list
  Vnames <- colnames(adjacency)
  if (!is.null(Vnames)) {
    names(res) <- Vnames[i]
  }
  else {
    names(res) <- paste0('V', i)
  }
  return(res)
}

####################################################
# 'enumerate.trios.i' is the workhorse for 'enumerate.trios'
enumerate.trios.i <- function (i,
                               adjacency,
                               n_t,
                               Tlabels,
                               Qlabels = NULL) {
  # Binary vector indicating T-nodes associated with 'Vi'
  Tlabeli <- adjacency[i, Tlabels]

  # number of T nodes that have an edge with 'Vi'
  pi <- sum(Tlabeli)

  # Terminate if 'Vi' has no edge with T-nodes
  if (pi == 0)
    return(NULL)

  if (is.null(Qlabels)) { # If trios NOT involving Q-nodes are desired
    # Doublets of T-nodes both associated with 'Vi'
    Apairs <- if (pi > 1) {
      t(combn(Tlabels[as.logical(Tlabeli)], 2))
    }

    # Pairs involving a T-node non-associated to Vi (with an associated one)
    Npairs <- if (pi < n_t) {
      form.doublets(Tlabeli = Tlabeli, Tlabels = Tlabels,
                    adjacency = adjacency)
    }

    # Bind the two groups of doublets
    res <- rbind(Apairs, Npairs)

  }
  else { # If trios involving Q-nodes are desired
    # Pairs involving a T-node and a Q-node
    res <- form.doublets_TQ (Tlabeli = Tlabeli,
                             Tlabels = Tlabels,
                             Qlabels = Qlabels,
                             adjacency = adjacency)

  }

  # If there is in fact no doublet, return a matrix with zero row
  if (NROW(res) == 0)
    res <- matrix(NA, ncol = 3, nrow = 0)
  else { # Otherwise, format the doublets and get trios
    # Order the doublets if many
    if (NROW(res) > 1) {
      # Using order (a, b); thanks to Bandita
      res <- res[order(res[,1], res[,2]), , drop = FALSE]
    }

    # Get trios
    res <- cbind(i, res)
  }

  # Names columns and return
  colnames(res) <- c("Vi", "Tj", if (is.null(Qlabels)) "Tk" else "Qk")
  return(res)
}

####################################################
# A routine for 'enumerate.trios.i': form all doublets involving each a T-node
# associated with Vi and a T-node non associated with 'Vi'
form.doublets <- function(Tlabeli, Tlabels, adjacency) {
  # Get a list of doublets involving each T-node
  res <- lapply(Tlabels[as.logical(Tlabeli)], FUN = form.doublets.j,
                Tlabeli = Tlabeli, Tlabels=Tlabels, adjacency=adjacency)

  # Unlist the result
  res <- unlist(res, recursive = TRUE)

  # 'unvec' the result if any
  if (!is.null(res))
    return(unique(matrix(res, ncol = 2, byrow = TRUE)))

  # Return a matrix with zero row if no result
  return(matrix(NA, ncol = 2, nrow = 0))
}

# A child function for the T-node numbered 'j'
form.doublets.j <- function(j, Tlabeli, Tlabels, adjacency) {
  # Binary vector indicating T-nodes associated with 'Tj'
  assoc <- adjacency[j, Tlabels]

  # Eliminating T-nodes already in Tlabeli (T-nodes directly associated with 'Vi')
  # to avoid duplicating T-nodes already accounted for in pairwise combinations
  assoc <- assoc * (1 - Tlabeli)

  # Terminate if no new association exists
  if (sum(assoc) == 0)
    return(NULL)

  # Form doublets (columns)
  resj <- rbind(j, Tlabels[as.logical(assoc)])

  # Sort each column
  # resj <- apply(resj, MARGIN = 2, FUN = sort) # (useful???)
  # Guess not useful, maybe counter productive

  # Return doublets involving 'Tj' (applying 'vec')
  return(c(resj))
}

form.doublets_TQ <- function(Tlabeli, Tlabels, Qlabels, adjacency) {
  # Get a list of doublets involving each T-node
  res <- lapply(Tlabels[as.logical(Tlabeli)], FUN = form.doublets_TQ.j,
                Qlabels = Qlabels,
                adjacency=adjacency)

  # Unlist the result
  res <- unlist(res, recursive = TRUE)

  # 'unvec' the result if any
  if (!is.null(res))
    return(unique(matrix(res, ncol = 2, byrow = TRUE)))

  # Return a matrix with zero row if no result
  return(matrix(NA, ncol = 2, nrow = 0))

}

# A child function for the T-node numbered 'j'
form.doublets_TQ.j <- function(j, Qlabels, adjacency) {
  # Binary vector indicating Q-nodes associated with 'Tj'
  assoc <- adjacency[j, Qlabels]

  # Terminate if no new association exists
  if (sum(assoc) == 0)
    return(NULL)

  # Form doublets (columns)
  resj <- rbind(j, Qlabels[as.logical(assoc)])

  # Sort each column
  # resj <- apply(resj, MARGIN = 2, FUN = sort) # (useful???)
  # Guess not useful, maybe counter productive

  # Return doublets involving 'Tj' (applying 'vec')
  return(c(resj))
}
