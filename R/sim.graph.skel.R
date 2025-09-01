# Generate a graph skeleton and the effects

#' @importFrom MRGN gen.conf.coefs

sim.graph.effects <- function(Adj, b.snp, b.med,
                               n_t,
                               n_v.t = 1, # numeric vector of the numbers of V-nodes per T-node, or its average
                               family.n_v = NULL, # distribution of the numbers of V-nodes per T-node, n_v.t is then an average
                               conf.coef.ranges,
                              nb.T.per.KU = 2,
                              neg.freq, degree, method = "er") {

  #random topology for T
  Adj.sub = methods::as(pcalg::randDAG(n = n_t, d = degree, method = method),"matrix")
  #replace with 1's
  Adj.sub = abs(Adj.sub) > 0

  # Add V topology
  if (is.numeric(n_v.t)) {
    n_v.t = rep(n_v.t, length.out = n_t)
  }
  else {
    stop("'n_v.t' must be numeric")
  }
  if (!is.null(family.n_v)) {
    stop("not yet implemented")
  }
  n_v = sum(n_v.t)
  if (n_v > 0) {
    if (all(n_v.t == 1)) {
      # Case one variant per T-node
      Adj.sub = rbind(diag(n_v),
                      Adj.sub)
    }
    else if (all(n_v.t == n_v.t[1])) {
      # Case of a constant number (>1) of variants per T-node
      Adj.sub = rbind(dirprod(diag(n_t),
                              rep(1, n_v.t[1])),
                      Adj.sub)
    }
    else {
      # General case where the number of variants varies from one T-node to another one
      Adj.sub = rbind(sapply(1:n_t,
                             FUN = Adj.sub.V,
                             n_t = n_t,
                             n_v.t = n_v.t),
                      Adj.sub)
    }
    Adj.sub = cbind(matrix(0, nrow = n_v + n_t, ncol = n_v),
                    Adj.sub)
  }

  #storing info
  V.idx = if (n_v) 1:n_v
  T.idx = (1:n_t)+n_v

  #create a second matrix for storing the simulation effects
  coefs.sub = Adj.sub
  #get the total number of V and T edges
  number.of.edges.V = if (n_v) sum(Adj.sub[V.idx,]) else 0
  number.of.edges.T = sum(Adj.sub[T.idx,])

  #replace edges with their simulated effects
  if (n_v)
    coefs.sub[V.idx,][coefs.sub[V.idx,]==1] = sample(b.snp, number.of.edges.V, replace = TRUE)
  coefs.sub[T.idx,][coefs.sub[T.idx,]==1] = sample(b.med, number.of.edges.T, replace = TRUE)
  #convert to igraph and get to the topological ordering
  topo.order = igraph::topo_sort(igraph::graph_from_adjacency_matrix(Adj.sub))
  Adj[1:(n_v+n_t), 1:(n_v+n_t)] = coefs.sub

  #handle confounders
  letter.id = c("W", "Z", "U", "K")
  for(i in 3:4){
    loc = which(grepl(letter.id[i], colnames(Adj)))
    for(j in loc){
      weight = gen.conf.coefs(n.effects = 1, coef.range.list = conf.coef.ranges[[i]], neg.freq = neg.freq)
      Adj[j, sample(T.idx, size = nb.T.per.KU, replace = F)] = weight
    }
  }

  #handle intermediate variables
  loc.int = which(grepl(letter.id[1], colnames(Adj)))
  for(i in loc.int){
    weights = gen.conf.coefs(n.effects = 2, coef.range.list = conf.coef.ranges[[1]], neg.freq = neg.freq)
    child.T = sample(topo.order[-c(1:(n_v+1))], 1)
    loc.in.topo.order = which(topo.order==child.T)
    poss.parents = c(topo.order[(n_v+1):(loc.in.topo.order-1)])
    if(length(poss.parents)>1){
      parent.T = sample(poss.parents, 1)
    }else{
      parent.T = poss.parents
    }
    Adj[i, child.T] = weights[1]
    Adj[parent.T, i] = weights[2]
  }

  #handle common children
  loc.cc = which(grepl(letter.id[2], colnames(Adj)))
  for(i in loc.cc){
    weight = gen.conf.coefs(n.effects = 1, coef.range.list = conf.coef.ranges[[2]], neg.freq = neg.freq)
    Adj[sample(T.idx, 2, replace = F), i] = weight
  }

  return(Adj)

}

Adj.sub.V <- function(k, n_t, n_v.t) {
  n_v.t = rep(n_v.t, length.out = n_t)
  n_v <- cumsum(n_v.t[1:n_t])
  Vcol <- numeric(n_v[n_t])
  if (n_v.t[k] == 0) {
    return(Vcol)
  }
  if (k == 1) {
    Vcol[1:n_v[1]] <- 1
    return(Vcol)
  }
  Vcol[(n_v[k - 1] + 1):n_v[k]] <- 1
  return(Vcol)
}

dirprod <- function (x, y)
{
  if (is.vector(x)) {
    A <- matrix(x)
  }
  else {
    if (is.matrix(x)) {
      A <- x
    }
    else {
      stop("Argument x is not a matrix or vector")
    }
  }
  if (is.vector(y)) {
    B <- matrix(y)
  }
  else {
    if (is.matrix(y)) {
      B <- y
    }
    else {
      stop("Argument y is not a matrix or vector")
    }
  }
  matrices <- list()
  for (i in 1:nrow(A)) {
    matrices[[i]] <- list()
    for (j in 1:ncol(A)) {
      matrices[[i]][[j]] <- A[i, j] * B
    }
  }
  row.matrices <- list()
  for (i in 1:nrow(A)) {
    row.matrices[[i]] <- matrices[[i]][[1]]
    if (ncol(A) > 1) {
      for (j in 2:ncol(A)) {
        row.matrices[[i]] <- cbind(row.matrices[[i]],
                                   matrices[[i]][[j]])
      }
    }
  }
  C <- row.matrices[[1]]
  if (length(row.matrices) > 1) {
    for (i in 2:length(row.matrices)) {
      C <- rbind(C, row.matrices[[i]])
    }
  }
  return(C)
}

