# Determine how reliable an inferred edge is
# Reliability of edges inferred present by MRGNgeneral
# object = of class MRGN
get.reliability.MRGN <- function (object, only.edges = FALSE) {
  p <- eval(object$call$p)
  q <- eval(object$call$q)
  r <- eval(object$call$r)
  nb.nodes <- p+q+r
  if (nb.nodes < 3)
    stop("Graph of three or more nodes required.")
  Tindex <- (q+1):(p+q)

  Adj <- as.matrix(object$Adj[1:nb.nodes, 1:nb.nodes])
  nb.mediator <- nb.children <- Adj[Tindex, Tindex]

  # matrix (lower tri) of number of direction between two nodes
  AugAdj <- nb.mediator + t(nb.mediator)
  AugAdj[upper.tri(AugAdj)] <- 0

  if (only.edges) {
    # Directed edges
    DirEdges <- which(AugAdj == 1, arr.ind = TRUE) + q
    nb.dir <- NROW(DirEdges)
    RelDir <- apply(DirEdges,
                    MARGIN = 1,
                    FUN = edge.reliability,
                    Adj = Adj)
    if (nb.dir > 1)
      RelDir <- t(RelDir)
    else
      RelDir <- matrix(c(RelDir), nrow = 1)

    # Undirected edges
    UndirEdges <- which(AugAdj == 2, arr.ind = TRUE) + q
    nb.undir <- NROW(UndirEdges)
    RelUnDir <- apply(UndirEdges,
                      MARGIN = 1,
                      FUN = edge.reliability,
                      Adj = Adj)
    if (nb.undir > 1)
      RelUnDir <- t(RelUnDir)
    else
      RelUnDir <- matrix(c(RelUnDir), nrow = 1)

    # Fill the matrices nb.mediator and nb.children
    nb.mediator[AugAdj == 1] <- RelDir[,1]
    nb.mediator[AugAdj == 2] <- RelUnDir[,1]

    nb.children[AugAdj == 1] <- RelDir[,2]
    nb.children[AugAdj == 2] <- RelUnDir[,2]
  }
  else {
    # Directed edges
    Edges <- which(lower.tri(AugAdj), arr.ind = TRUE) + q
    nb.edges <- NROW(Edges)
    RelEdges <- apply(Edges,
                    MARGIN = 1,
                    FUN = edge.reliability,
                    Adj = Adj)
    if (nb.edges > 1)
      RelEdges <- t(RelEdges)
    else
      RelEdges <- matrix(c(RelEdges), nrow = 1)

    # Fill the matrices nb.mediator and nb.children
    nb.mediator[lower.tri(AugAdj)] <- RelEdges[,1]
    nb.children[lower.tri(AugAdj)] <- RelEdges[,2]
  }

  # Ensure symmetry
  nb.mediator <- nb.mediator + t(nb.mediator)
  nb.children <- nb.children + t(nb.children)

  return(list(nb.mediator = nb.mediator,
              nb.children = nb.children))

}

# Returns the number of mediators and the number of common children for two T-nodes
edge.reliability <- function (x, Adj) {
  ## Find mediators
  # T-nodes pointing to T1
  Tj <- which(Adj[,x[1]] == 1)
  mediators <- if (length(Tj)) {
    # Those of Tj that T2 is pointing to
    which(Adj[x[2], Tj] == 1)
  }

  # T-nodes pointing to T2
  Tk <- which(Adj[,x[2]] == 1)
  mediators <- c(mediators,
                 if (length(Tk)) {
                   # Those of Tk that T1 is pointing to
                   which(Adj[x[1], Tk] == 1)
                 })
  mediators <- unique(mediators)

  ## Find common children
  T1 <- which(Adj[x[1],] == 1)
  T2 <- which(Adj[x[2],] == 1)

  children <- intersect(T1, T2)

  return(c(length(mediators), length(children)))
}
