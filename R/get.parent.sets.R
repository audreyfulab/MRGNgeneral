### Build the list of true confounding variables for T-nodes
### From an adjacency matrix
### Refer to MRGN for notation (n_t,n_v,n_wz,n_u); n_wz = n_q (the true n_q)
### offset is the number of variables in Adj before true confounding variables
# Get confounder set
get.true.conf.set <- function(Adj, n_v = 0, n_t, n_wz = 0, n_u,
                              Upool = TRUE,
                              offset = n_t + n_v + if (Upool) n_wz else 0) {
  if (Upool) {
    if (n_u == 0) {
      return(list())
    }

    if (all(dim(Adj) == c(n_u, n_t))) {
      KU_T.Adj <- Adj
    }
    else {
      KU_T.Adj <- Adj[(n_t + n_v + n_wz + 1):(n_t + n_v + n_wz + n_u), (n_v + 1):(n_t + n_v)]
    }

    ConfSet <- lapply(1:n_t, FUN = function(j) {
      TK <- KU_T.Adj[,j] * (1:n_u)
      if (any((Kj <- TK > 0)))
        return(TK[Kj] + offset) # True confounders for T node j
      NULL
    })

    names(ConfSet) <- colnames(KU_T.Adj)
  }
  else {
    if (n_wz == 0) {
      return(list())
    }

    if (all(dim(Adj) == c(n_wz, n_t))) {
      Q_T.Adj <- Adj
    }
    else {
      Adj <- ((Adj + t(Adj)) > 0) + 0
      Q_T.Adj <- Adj[(n_t + n_v + 1):(n_t + n_v + n_wz), (n_v + 1):(n_t + n_v)]
    }

    ConfSet <- lapply(1:n_t, FUN = function(j) {
      TQ <- Q_T.Adj[,j] * (1:n_wz)
      if (any((Qj <- TQ > 0)))
        return(TQ[Qj] + offset) # True W,Z-nodes for T node j
      NULL
    })
    names(ConfSet) <- colnames(Q_T.Adj)

  }

  return(ConfSet)
}

### Build the list of parent V-nodes for T-nodes
### From an adjacency matrix
get.true.variant.set <- function(Adj, n_t, n_v) {
  if (n_v == 0) {
    return(list())
  }
  if (all(dim(Adj) == c(n_v, n_t))) {
    V_T.Adj <- Adj
  }
  else {
    V_T.Adj <- Adj[1:n_v, (n_v + 1):(n_t + n_v)]
  }
  VSet <- lapply(1:n_t, FUN = function(j) {
    VT <- V_T.Adj[,j] * (1:n_v)
    if (any((Vj <- VT > 0)))
      return(VT[Vj])
    NULL
  })
  names(VSet) <- colnames(V_T.Adj)
  return(VSet)
}

### Build the list of parent T-nodes for each T-node
### From an adjacency matrix
get.true.parent.genes <- function(Adj, n_v = 0, n_t = NCOL(Adj) - n_v) {
  if (n_t == 0) {
    return(list())
  }
  Adj <- Adj[(n_v+1):(n_t+n_v), (n_v+1):(n_t+n_v)]
  TSet <- lapply(1:n_t, FUN = function(j) {
    TT <- Adj[,j] * (1:n_t)
    if (any((Tj <- TT > 0)))
      return(TT[Tj] + n_v)
    NULL
  })
  names(TSet) <- colnames(Adj)
  return(TSet)
}
