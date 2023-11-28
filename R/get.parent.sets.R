### Build the list of true confounding variables for T-nodes
### From an adjacency matrix
### Refer to MRGN for notation (p,q,r,u)
### offset is the number of variables in Adj before true confounding variables
get.true.conf.set <- function(Adj, p, u, q = 0, r = 0, offset = p + q + r) { # offset = nb.VTnodes + r
  if (u == 0) {
    return(list())
  }
  if (all(dim(Adj) == c(u, p))) {
    KU_T.Adj <- Adj
  }
  else {
    KU_T.Adj <- Adj[(p + q + 1):(p + q + u), (q + 1):(p + q)]
  }
  ConfSet <- lapply(1:p, FUN = function(j) {
    TK <- KU_T.Adj[,j] * (1:u)
    if (any((Kj <- TK > 0)))
      return(TK[Kj] + offset) # True confounders for T node j
    NULL
  })
  names(ConfSet) <- colnames(KU_T.Adj)
  return(ConfSet)
}

### Build the list of parent V-nodes for T-nodes
### From an adjacency matrix
get.true.variant.set <- function(Adj, p, q) {
  if (q == 0) {
    return(list())
  }
  if (all(dim(Adj) == c(q, p))) {
    V_T.Adj <- Adj
  }
  else {
    V_T.Adj <- Adj[1:q, (q + 1):(p + q)]
  }
  VSet <- lapply(1:p, FUN = function(j) {
    VT <- V_T.Adj[,j] * (1:q)
    if (any((Vj <- VT > 0)))
      return(VT[Vj])
    NULL
  })
  names(VSet) <- colnames(V_T.Adj)
  return(VSet)
}

### Build the list of parent T-nodes for each T-node
### From an adjacency matrix
get.true.parent.genes <- function(Adj, q = 0, p = NCOL(Adj) - q) {
  if (p == 0) {
    return(list())
  }
  Adj <- Adj[(q+1):(p+q), (q+1):(p+q)]
  TSet <- lapply(1:p, FUN = function(j) {
    TT <- Adj[,j] * (1:p)
    if (any((Tj <- TT > 0)))
      return(TT[Tj] + q)
    NULL
  })
  names(TSet) <- colnames(Adj)
  return(TSet)
}
