#' @name enumerate.triplets
#' @title List all (update-able) trios involving only phenotypes
#' @description Enumerate update-able triplets involving only expressions/phenotypes
#'  (no genetic variant) based on an adjacency matrix. Not a user level function.
#'
#' @export enumerate.triplets
#'
#' @param Adj numeric, a binary matrix indicating adjacency and possibly edge
#' directions (i.e. a non-directed, or a partially directed graph). All the
#' \code{p} columns of \code{Adj} represent \code{T}-nodes (phenotypes). As such,
#' \code{Adj} will in general be an extract of a larger adjacency matrix (last
#' \code{p} rows and columns in this package). See Section \code{Details} for more
#' on the interpretation of the values \code{0/1} in \code{Adj} for this function.
#'
#' @param p,q,r integers, respectively: number of \code{T}-nodes (non-instrumental variables),
#' number of \code{V}-nodes (instrumental variables), and number of \code{Q}-nodes
#' (confounding variables: intermediate variables and common children).
#'
#' @param TTT logical indicating the composition of desired trios. If true (the default),
#' only trios not involving a \code{Q}-node are formed, i.e. a trio contains one
#' three \code{T}-nodes. Otherwise, only trios involving \code{Q}-nodes are
#' formed, i.e.  a trio contains two \code{T}-nodes, and one \code{Q}-node.
#'
#' @param return.3 logical, should triplets of type (3) be returned as an
#' attribute named \code{triplet.3} (of class \code{matrix}) for the output?
#' Defaults to \code{FALSE}. See Section \code{Details} for more on type (1)
#' triplets.
#'
#' @param cl a cluster object, created by one of the packages \code{parallel} and \code{snow}.
#' If \code{NULL}, the registered default cluster is used. Note that the latter can be
#' \code{NULL} too, in which case, no parallel computation is performed.
#'
#' @param chunk.size integer, number of tasks per scheduling unit during parallel computation.
#'
#' @details The binary matrix \code{Adj} represent a non-directed, or a partially
#' directed graph. For two nodes \code{Ti} and \code{Tj}, the matrix elements
#' \code{a_ij} and \code{a_ji} define both the presence of an edge and the
#' direction of the edge if any: an edge is present between nodes \code{Ti}
#' and \code{Tj} if \code{a_ij + a_ji > 0}. If \code{a_ij + a_ji = 2}, then the
#' edge is bi-directed (or equivalently, undirected). Otherwise, if \code{a_ij = 1},
#' then the edge goes from \code{Ti} to \code{Tj}, and if \code{a_ji = 1}, the
#' edge goes from \code{Tj} to \code{Ti}.
#'
#' The function only returns trios for which graph models can be further
#' identified under the extended interpretation of the Principle of Mendelian
#' Randomization (PMR). For trios of strictly \code{T}-nodes with at
#' least one undirected edge, two general situations are possible: (A) the trio
#' has two edges; and (B) the trio has three edges. In situation (B), no further
#' inference is possible, and such triplets are *not* returned.
#'
#' In situation (A), we have three possibilities: (1) two undirected edges (\code{T1-T2-T3});
#' (2) one directed edge, and the parent \code{T}-node has only one edge (\code{T1->T2-T3});
#' and (3) one directed edge, and the parent \code{T}-node has two edges (\code{T1<-T2-T3}).
#'
#' In case (1), one of the edges of the triplet can potentially be directed using
#' a conditional independence test between the two \code{T}-nodes not directly
#' related by an edge. Here, the ability to direct the second edge depends on
#' the test result. Such triplets are returned in the outputed matrix.
#'
#' In case (2), the parent \code{T}-node can be treated as a genetic variant.
#' The second edge of the triplet can be directed using a conditional independence
#' test between the parent \code{T}-node and the non directly-related \code{T}-node.
#' Such triplets are returned in the outputed matrix.
#'
#' In case (3), no further inference is possible under PMR. In the absence of an
#' additional genetic variant linked to one or many of the \code{T}-nodes, we
#' have three different but Markov equivalent graph candidates, i.e., these
#' graphs share the same set of conditional and marginal independence relations
#' (see Kvamme and Fu, 2024). Such triplets are NOT returned in the outputed
#' matrix. They are nevertheless returned as an attribute named \code{triplet.3}
#' (of class \code{matrix}) for the matrix output if \code{return.3 = TRUE}
#' (note that the default is \code{return.3 = FALSE}).
#'
#' In the outputed matrix, the columns are named 'Ti', 'Tj', and 'Tk', and each
#' row indicates three \code{T}-nodes whose edges can potentially be updated
#' (i.e. directed). The \code{T}-nodes are arranged such that 'Tj' has two
#' edges (\code{Ti-Tj-Tk}). This way, any further inference only requires regressing 'Ti' (or 'Tk')
#' on 'Tj' and 'Tk' (or 'Ti'), and confounding variables, if any. If a triplet
#' is of type \code{2}, then 'Ti' is the parent node (i.e. 'Ti' points to 'Tj').
#'
#' @return A list with two elements:
#'
#' \item{triplets}{an \code{n} by \code{3} matrix. Each row of the output gives
#' the indices of non-instrumental variables forming a triplet with two edges.
#' The columns of the matrix are named 'Ti', 'Tj', and 'Tk'. The \code{T}-nodes
#' are arranged such that "Tj" has two edges. If no update-able triplet can be
#' formed, the returned matrix has zero row.}
#'
#' \item{types}{a \code{n}-vector of integers. The \code{i}th element indicates the
#' type of triplet corresponding to the i^th row the matrix \code{triplets}.}
#'
#' If \code{return.3 = TRUE}, then the output matrix \code{triplets} has an
#' attribute named \code{triplet.3} giving triplets of type (3). See Section
#' \code{Details} for more on type (3) triplets. The attribute \code{triplet.3}
#' is a \code{matrix}, and it has zero row if no triplet of type (3) is
#' found.
#'
#' @seealso \link{enumerate.trios}.
#'
#
####################################################
# List update-able triplets involving strictly T-nodes
enumerate.triplets <- function (Adj,
                                p,
                                q,
                                r = NROW(Adj) - p - q,
                                TTT = TRUE,
                                return.3 = FALSE,
                                cl, chunk.size = NULL) {
  # Check arguments
  if (missing(p))
    stop("Argument 'p' must be specified")

  stopifnot(is.logical(TTT[1]))
  if (!TTT[1] & r == 0) {
    stop("'r > 0' is required when 'TTT = FALSE'")
  }

  # Set default cluster
  if (is.null(cl))
    cl <- parallel::getDefaultCluster()

  # Keep only T-nodes (discard V-nodes) and Q-nodes if any
  Adj <- Adj[(q + 1):(q + p + if (TTT[1]) 0 else r),
             (q + 1):(q + p + if (TTT[1]) 0 else r)]

  # Adjacency matrix of the un-directed graph
  UndirAdj <- ((Adj + t(Adj)) > 0) + 0

  # Labels (numbers) for all T-nodes
  Tlabels <- 1:p

  # Labels (numbers) for all Q-nodes if required
  Qlabels <- if (!TTT[1]) (p+1):(p+r)

  if (TTT[1]) {
    # An indicator for T,Q-nodes with at least two edges
    Nedges <- rowSums(UndirAdj)
    margins <- Nedges >= 2

    # List all triplets involving a T-node with at least two edges
    triplets <- matteLapply (X = Tlabels[margins],
                             FUN = find.triplets.i,
                             Adj = UndirAdj,
                             Tlabels = Tlabels,
                             cl = cl, chunk.size = chunk.size)
  }
  else {
    # An indicator for Q-nodes with at least one edge with a T-node
    Nedges <- rowSums(UndirAdj[Qlabels, Tlabels])
    margins <- Nedges >= 1

    # List all triplets involving a Q-node, and with at least two edges
    triplets <- matteLapply (X = Qlabels[margins],
                             FUN = find.Qtriplets.i,
                             Adj = UndirAdj,
                             Tlabels = Tlabels,
                             cl = cl, chunk.size = chunk.size)
  }
  triplets <- do.call('rbind', triplets)

  # Terminate if no strictly T,Q-nodes triplet can be formed
  if (length(dim(triplets)) < 2) {
    return(triplet.null (return.3))
  }

  # Get update-able triplets (and types of triplet) if any
  triplets <- t(matteApply(triplets,
                           MARGIN = 1,
                           FUN = get.triplet.type.i,
                           Adj = Adj,
                           Tlabels = Tlabels,
                           TTT = TTT[1],
                           Qlabels = Qlabels,
                           cl = cl, chunk.size = chunk.size))

  # An indicator to remove null rows
  keep <- rowSums(triplets) > 0

  # Terminate if no strictly T-nodes triplet can be formed
  if (!any(keep)) {
    return(triplet.null (return.3))
  }

  # Remove null rows and extract 'types'
  triplets <- triplets[keep, , drop = FALSE]
  colnames(triplets) <- c("Ti", "Tj", "Tk", "types")
  types <- triplets[,4]

  # Indicator of type (1) and (2) triplets
  keep <- types < 3

  # Return
  if (return.3[1]) {
    triplet.3 <- if (!all(keep))
      triplets[!keep, 1:3, drop = FALSE]
    else
      triplet.null (FALSE)$triplets
    triplets <- triplets[keep, 1:3, drop = FALSE]
    attr(triplets, 'triplet.3') <- triplet.3
  }
  else {
    triplets <- if (any(keep))
      triplets[keep, 1:3, drop = FALSE]
    else
      triplet.null (FALSE)$triplets
  }
  return(list(triplets = triplets + q, types = if (any(keep)) types[keep]))
}

####################################################
# A routine check if a triplet is update-able
get.triplet.type.i <- function(x, Adj, Tlabels,
                               TTT = TRUE,
                               Qlabels = NULL) {
  res <- find.triplet.type (Aijk = Adj[x, x],
                            labelsijk = c(Tlabels, if (!TTT) Qlabels)[x],
                            TTT = TTT)

  if (is.null(res$types))
    return(rep(0, 4))
  return(c(res$triplets, res$types))
}
# Sub routine for 'get.triplet.type.i'
# 'Aijk' is a 3 by 3 binary matrix
# 'labelsijk' is 3-vector of numeric labels (numbers of the 3 T-nodes)
find.triplet.type <- function (Aijk, labelsijk, TTT = TRUE) {
  # Terminate if no edge is present
  if(!sum(Aijk)) {
    return(triplet.null (FALSE))
  }

  # Presence of edges per node
  n_edges <- pmax(rowSums(Aijk), colSums(Aijk))

  # Terminate if any node is isolated from the two others
  if(any(n_edges == 0)) {
    return(triplet.null (FALSE))
  }

  # Count edge indicators
  Edge_c <- Aijk[lower.tri(Aijk)] + Aijk[upper.tri(Aijk)]

  # Terminate if all present edges are already directed
  if(all(Edge_c <= 1)) {
    return(triplet.null (FALSE))
  }

  # Determine the presence of edge(s)
  #p_e <- (Edge_c > 0) + 0
  # Determine the number of edge(s)
  #n_e <- sum(p_e)
  # Terminate if the number of edges in the trio is not exactly 2
  if (sum(Edge_c > 0) != 2) {
    return(triplet.null (FALSE))
  }

  if (TTT) {
    # Which one of T1, T2 and T3 is the mediator?
    if (all(Edge_c[1:2] > 0)) { # T1 is the mediator
      # Do we have two undirected edges from/to T1?
      if (all(Edge_c[1:2] == 2)) {
        # Form a row matrix with node numbers
        triplets <- matrix(labelsijk[c(2, 1, 3)], nrow = 1, ncol = 3)
        types <- 1 # T2--T1--T3
      }
      else {
        if (Edge_c[1] == 1) { # If T1-T2 is directed
          # Form a row matrix with node numbers
          triplets <- matrix(labelsijk[c(2, 1, 3)], nrow = 1, ncol = 3)

          if (Aijk[2,1] == 1)
            types <- 2 # T2->T1--T3
          else
            types <- 3 # T2<-T1--T3
        }
        else {  # If T1-T3 is directed
          # Form a row matrix with node numbers
          triplets <- matrix(labelsijk[c(3, 1, 2)], nrow = 1, ncol = 3)

          if (Aijk[3,1] == 1)
            types <- 2 # T3->T1--T2
          else
            types <- 3 # T3<-T1--T2
        }
      }

    }
    else if(all(Edge_c[2:3] > 0)) { # T3 is the mediator
      # Do we have two undirected edges from/to T3?
      if (all(Edge_c[1:2] == 2)) {
        # Form a row matrix with node numbers
        triplets <- matrix(labelsijk[c(1, 3, 2)], nrow = 1, ncol = 3)
        types <- 1 # T1--T3--T2
      }
      else {
        if (Edge_c[2] == 1) { # If T1-T3 is directed
          # Form a row matrix with node numbers
          triplets <- matrix(labelsijk[c(1, 3, 2)], nrow = 1, ncol = 3)

          if (Aijk[1,3] == 1)
            types <- 2 # T1->T3--T2
          else
            types <- 3 # T1<-T3--T2
        }
        else {  # If T2-T3 is directed
          # Form a row matrix with node numbers
          triplets <- matrix(labelsijk[c(2, 3, 1)], nrow = 1, ncol = 3)

          if (Aijk[2,3] == 1)
            types <- 2 # T2->T3--T1
          else
            types <- 3 # T2<-T3--T1
        }
      }
    }
    else { # T2 is the mediator
      # Do we have two undirected edges from/to T2?
      if (all(Edge_c[c(1,3)] == 2)) {
        # Form a row matrix with node numbers
        triplets <- matrix(labelsijk, nrow = 1, ncol = 3)
        types <- 1 # T1--T2--T3
      }
      else {
        if (Edge_c[1] == 1) { # If T1-T2 is directed
          # Form a row matrix with node numbers
          triplets <- matrix(labelsijk, nrow = 1, ncol = 3)

          if (Aijk[1,2] == 1)
            types <- 2 # T1->T2--T3
          else
            types <- 3 # T1<-T2--T3
        }
        else {  # If T2-T3 is directed
          # Form a row matrix with node numbers
          triplets <- matrix(labelsijk[c(3, 2, 1)], nrow = 1, ncol = 3)

          if (Aijk[3,2] == 1)
            types <- 2 # T3->T2--T1
          else
            types <- 3 # T3<-T2--T1
        }
      }
    }

    colnames(triplets) <- c("Ti", "Tj", "Tk")
    return(list(triplets = triplets, types = types))

    ############################################################################
    ################### Fossil code; NOT RUN; Remove  ###################
    # Indicator for: number of edges per node = 2 or not
    {
      n_2e <- n_edges == 2

    # Put the node with two edges in the first position if it is not yet
    if(!n_2e[1]) {
      labelsijk <- c(labelsijk[n_2e], labelsijk[!n_2e])
    }

    # Form a row matrix with node numbers
    triplets <- matrix(labelsijk, nrow = 1, ncol = 3)
    colnames(triplets) <- c("Ti", "Tj", "Tk")

    # Number of directed edges
    n_d <- sum(Edge_c == 1) # = 0 or 1 (case n_d = 2 already ruled out)

    # If n_d = 0, return triplet of type (1)
    if (n_d == 0) {
      return(list(triplets = triplets, types = 1))
    }

    # If n_d == 1, determine if we have a triplet of type (2) or (3)
    # Test if the node with two edges is the parent node
    if (sum(Aijk[n_2e,]) == 2) { # type = '3'
      return(list(triplets = triplets, types = 3))
    }

    # Otherwise, type = '2'
    if(!n_2e[1]) { # Reorganize the matrix to match 'labelsijk'
      i <- sum(n_2e * (1:3))
      Aijk <- cbind(Aijk[,n_2e], Aijk[, !n_2e])
      Aijk <- rbind(Aijk[n_2e,], Aijk[!n_2e,])
    }


    # Ensure that 'Tj' is the parent node (put it in position 2)
    if (sum(Aijk[,2]) > 0)
      triplets[,2:3] <- triplets[,3:2]

    return(list(triplets = triplets, types = 2))
  }
    ############################################################################

  }
  else {
    # These trios are here ordered as [Qi, Tj, Tk]

    # Is Q a mediator?
    if (all(Edge_c[1:2] > 0)) {# Q is a mediator
      # Do we have two undirected edges from/to Q?
      if (all(Edge_c[1:2] == 2)) {
        # Form a row matrix with node numbers
        triplets <- matrix(labelsijk[c(2, 1, 3)], nrow = 1, ncol = 3)
        types <- 1 # T1--Q--T2
      }
      else {
        if (Edge_c[1] == 1) { # If Q-T1 is directed
          # Form a row matrix with node numbers
          triplets <- matrix(labelsijk[c(2, 1, 3)], nrow = 1, ncol = 3)

          if (Aijk[2,1] == 1)
            types <- 2 # T1->Q--T2
          else
            types <- 3 # T1<-Q--T2
        }
        else {  # If Q-T2 is directed
          # Form a row matrix with node numbers
          triplets <- matrix(labelsijk[c(3, 1, 2)], nrow = 1, ncol = 3)

          if (Aijk[3,1] == 1)
            types <- 2 # T2->Q--T1
          else
            types <- 3 # T2<-Q--T1
        }
      }
      colnames(triplets) <- c("Ti", "Qj", "Tk")
    }
    else if (all(Edge_c[1] == 0)) { # Q is not a mediator and we have T2-Q edge (not have T1-Q)
      # Do we have two undirected edges from/to T2?
      if (all(Edge_c[2:3] == 2)) {
        # Form a row matrix with node numbers: T1--T2--Q
        triplets <- matrix(labelsijk[c(2, 3, 1)], nrow = 1, ncol = 3)
        colnames(triplets) <- c("Ti", "Tj", "Qk")
        types <- 1 # T1--T2--Q
      }
      else {
        if (Edge_c[3] == 1) { # If T2-T1 is directed
          # Form a row matrix with node numbers: T1--T2--Q
          triplets <- matrix(labelsijk[c(2, 3, 1)], nrow = 1, ncol = 3)
          colnames(triplets) <- c("Ti", "Tj", "Qk")

          if (Aijk[2,3] == 1)
            types <- 2 # T1->T2--Q
          else
            types <- 3 # T1<-T2--Q
        }
        else {  # If T2-Q is directed
          # Form a row matrix with node numbers: Q--T2--T1
          triplets <- matrix(labelsijk[c(1, 3, 2)], nrow = 1, ncol = 3)
          colnames(triplets) <- c("Qi", "Tj", "Tk")

          if (Aijk[1,3] == 1)
            types <- 2 # Q->T2--T1
          else
            types <- 3 # Q<-T2--T1
        }
      }
    }
    else { # Q is not a mediator and we have T1-Q edge
      # Do we have two undirected edges from/to T1?
      if (all(Edge_c[c(1,3)] == 2)) {
        # Form a row matrix with node numbers: T2--T1--Q
        triplets <- matrix(labelsijk[c(3, 2, 1)], nrow = 1, ncol = 3)
        colnames(triplets) <- c("Ti", "Tj", "Qk")
        types <- 1 # T2--T1--Q
      }
      else {
        if (Edge_c[3] == 1) { # If T2-T1 is directed
          # Form a row matrix with node numbers: T2--T1--Q
          triplets <- matrix(labelsijk[c(3, 2, 1)], nrow = 1, ncol = 3)
          colnames(triplets) <- c("Ti", "Tj", "Qk")

          if (Aijk[3,2] == 1)
            types <- 2 # T2->T1--Q
          else
            types <- 3 # T2<-T1--Q
        }
        else {  # If T1-Q is directed
          # Form a row matrix with node numbers: Q--T1--T2
          triplets <- matrix(labelsijk[c(1, 2, 3)], nrow = 1, ncol = 3)
          colnames(triplets) <- c("Qi", "Tj", "Tk")

          if (Aijk[1,2] == 1)
            types <- 2 # Q->T1--T2
          else
            types <- 3 # Q<-T1--T2
        }
      }
    }

    return(list(triplets = triplets, types = types))
  }
}
####################################################

# Parallelized version of 'enumerate.triplets'
enumerate.triplets.parallel <- function (Adj, return.3 = FALSE, cl,
                                         chunk.size = NULL) {

  # Labels (numbers) for all T-nodes
  p <- NCOL(Adj)
  Tlabels <- 1:p

  # Adjacency matrix of the un-directed graph
  UndirAdj <- ((Adj + t(Adj)) > 0) + 0

  # An indicator for T-nodes with at least two edges
  Nedges <- rowSums(UndirAdj)
  margins <- Nedges >= 2

  # List all triplets involving a T-node with at least two edges
  triplets <- parallel::parLapply (cl = cl,
                                   X = Tlabels[margins], fun = find.triplets.i,
                                   Adj = UndirAdj, Tlabels = Tlabels)
  triplets <- do.call('rbind', triplets)

  # Terminate if no strictly T-nodes triplet can be formed
  if (length(dim(triplets)) < 2) {
    return(triplet.null (return.3))
  }

  # Get update-able triplets (and types of triplet) if any
  triplets <- t(parallel::parApply (cl = cl,
                                    X = triplets,
                                    MARGIN = 1, FUN = function(x) {
    res <- find.triplet.type (Aijk = Adj[x, x],
                              labelsijk = Tlabels[x])
    if (is.null(res$types))
      return(rep(0, 4))
    return(c(res$triplets, res$types))
  }))

  # An indicator to remove null rows
  keep <- rowSums(triplets) > 0

  # Terminate if no strictly T-nodes triplet can be formed
  if (!any(keep)) {
    return(triplet.null (return.3))
  }

  # Remove null rows and extract 'types'
  triplets <- triplets[keep, , drop = FALSE]
  colnames(triplets) <- c("Ti", "Tj", "Tk", "types")
  types <- triplets[,4]

  # Indicator of type (2) and (3) triplets
  keep <- types > 1

  # Return
  if (return.3[1]) {
    triplet.3 <- if (!all(keep))
      triplets[!keep, 1:3, drop = FALSE]
    else
      triplet.null (FALSE)$triplets
    triplets <- triplets[keep, 1:3, drop = FALSE]
    attr(triplets, 'triplet.3') <- triplet.3
  }
  else {
    triplets <- if (any(keep)) {
      triplets[keep, 1:3, drop = FALSE]
    }
    else
      triplet.null (FALSE)$triplets
  }
  return(list(triplets = triplets, types = if (any(keep)) types[keep]))
}

# List update-able triplets involving strictly T-nodes
# Initial version of enumerate.triplets
# Uses a greedy search
# Kept to test the effectiveness of 'enumerate.triplets'
enumerate.triplets.greedy <- function (Adj, return.3 = FALSE) {
  # Labels (numbers) for all T-nodes
  p <- NCOL(Adj)
  Tlabels <- 1:p

  # An indicator to exclude isolated T-nodes
  margins <- (rowSums(Adj) + colSums(Adj)) > 0

  # Number of non completely isolated T-nodes
  nnodes <- sum(margins)

  # Terminate if no strictly T-nodes triplet can be formed
  if (nnodes < 3) {
    return(triplet.null (return.3))
  }

  # Special case of only one possible triplet
  if (nnodes == 3) {
    # Get triplet (and type of triplet) if any
    triplets <- find.triplet.type (Aijk = Adj[margins, margins],# Extract sub-adjacency matrix
                                   labelsijk = Tlabels[margins])# Pick the 3 labels

    # Return a 0 row matrix if triplet is of type (1)
    if (triplets$types == 1) {
      # Save 'triplet.3' as an attribute if asked for
      if (return.3[1]) {
        triplet.3 <- triplets$triplets
        triplets <- triplet.null (FALSE)
        attr(triplets$triplets, 'triplet.3') <- triplet.3
        return(triplets)
      }
      return(triplet.null (FALSE))
    }

    # Give a 0 row matrix as the attribute 'triplet.3' if required
    if (return.3[1]) {
      attr(triplets$triplets, 'triplet.3') <- triplet.null (FALSE)$triplets
    }

    # Return the triplet
    return(triplets)
  }

  # List all non-isolated triplets
  triplets <- combn(Tlabels[margins], 3)

  # Get update-able triplets (and types of triplet) if any
  triplets <- t(apply(triplets, MARGIN = 2, FUN = function(x) {
    res <- find.triplet.type (Aijk = Adj[x, x],
                              labelsijk = Tlabels[x])
    if (is.null(res$types))
      return(rep(0, 4))
    return(c(res$triplets, res$types))
  }))

  # An indicator to remove null rows
  keep <- rowSums(triplets) > 0

  # Terminate if no strictly T-nodes triplet can be formed
  if (!any(keep)) {
    return(triplet.null (return.3))
  }

  # Remove null rows and extract 'types'
  triplets <- triplets[keep, , drop = FALSE]
  colnames(triplets) <- c("Ti", "Tj", "Tk", "types")
  types <- triplets[,4]

  # Indicator of type (2) and (3) triplets
  keep <- types > 1

  # Return
  if (return.3[1]) {
    triplet.3 <- if (!all(keep))
      triplets[!keep, 1:3, drop = FALSE]
    else
      triplet.null (FALSE)$triplets
    triplets <- triplets[keep, 1:3, drop = FALSE]
    attr(triplets, 'triplet.3') <- triplet.3
  }
  else {
    triplets <- if (any(keep)) {
      triplets[keep, 1:3, drop = FALSE]
    }
    else
      triplet.null (FALSE)$triplets
  }
  return(list(triplets = triplets, types = if (any(keep)) types[keep]))
}

####################################################
# A routine to return an empty matrix of triplets
# 'return.3' is a logical scalar
triplet.null <- function (return.3) {
  triplets <- matrix(NA, ncol = 3, nrow = 0)
  colnames(triplets) <- c("Ti", "Tj", "Tk")
  if (return.3[1]) {
    attr(triplets, 'triplet.3') <- triplets
  }
  return(list(triplets = triplets, types = NULL))
}

####################################################
# Routine for 'enumerate.triplets' (which enumerates triplets sequentially)
# Find triplets involving a particular T-node
# This function is a slightly modified version of 'enumerate.trios.i'
find.triplets.i <- function(i, Adj, Tlabels) {
  # Binary vector indicating T-nodes associated with 'Ti'
  Tlabelj <- Adj[i, Tlabels]

  # number of T-nodes that have an edge with 'Ti'
  pi <- sum(Tlabelj)

  # return a matrix with zero row if 'Ti' does not have edges with at least 2 other T-nodes
  if (pi < 2) {
    res <- matrix(NA, ncol = 3, nrow = 0)
  }
  else {
    # Doublets of T-nodes both associated with 'Ti'
    res <- t(combn(Tlabels[as.logical(Tlabelj)], 2))
    # Order the doublets if many
    if (NROW(res) > 1) {
      # Using order (a, b); thanks to Bandita
      res <- res[order(res[,1], res[,2]), , drop = FALSE]
    }
    # Get triplets
    res <- cbind(i, res)
  }

  # Names columns and return
  colnames(res) <- c("Ti", "Tj", "Tk")
  return(res)
}

# This function is very close to 'enumerate.trios.i' (with Qi instead of Vi)
# Should find a way to avoid this duplicate
find.Qtriplets.i <- function(i, Adj, Tlabels) {
  # Binary vector indicating T-nodes associated with 'Qi'
  Tlabeli <- Adj[i, Tlabels]

  # number of T-nodes that have an edge with 'Qi'
  pi <- sum(Tlabeli)

  # return a matrix with zero row if 'Qi' does not have edges with at least one other T-node
  if (pi == 0) {
    res <- matrix(NA, ncol = 3, nrow = 0)
  }
  else {
    # Doublets of T-nodes both associated with 'Qi'
    Apairs <- if (pi > 1) {
      t(combn(Tlabels[as.logical(Tlabeli)], 2))
    }

    # Pairs involving a T-node non-associated to Qi (with an associated one)
    Npairs <- if (pi < length(Tlabels)) {
      form.doublets(Tlabeli = Tlabeli, Tlabels = Tlabels, Adj = Adj)
    }

    # Bind the two groups of doublets
    res <- rbind(Apairs, Npairs)

    # Order the doublets if many
    if (NROW(res) > 1) {
      # Using order (a, b); thanks to Bandita
      res <- res[order(res[,1], res[,2]), , drop = FALSE]
    }

    # Get triplets
    res <- cbind(i, res)
  }

  # Names columns and return
  colnames(res) <- c("Qi", "Tj", "Tk")
  return(res)
}
