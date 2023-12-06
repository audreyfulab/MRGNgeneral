#' sample.graph.data
#'
#' @description Different from \code{sim.graph.data} in the way \code{V}-nodes are
#' included in the network skeleton: here there are just added (not simulated).
#' The number of \code{V}-nodes is specified for each \code{T}-node, or sampled
#' from a count distribution. Each \code{V}-node is only related to only one \code{T}-node.
#'
#' @export sample.graph.data
#'
# Import 'find.parents' from MRGN when MRGN export it.
# @importFrom MRGN find.parents
#'
#' @importFrom igraph graph_from_adjacency_matrix
#'
#' @importFrom igraph topo_sort
#'
#' @importFrom stats rnorm
#' @importFrom stats na.omit
#'
#' @importFrom pcalg randDAG
# 'graph_type' specifies a limited number of graph structure
# one of "small-world", "scale-free" or "random graph"
# 'method' is passed to 'pcalg::randDAG',
# one of "regular", "watts", "er", "power", "bipartite", "barabasi", "geometric", or "interEr"
#'
# Important change on Oct.-12-2023: reordered variables as V, T, W, Z, U, K (old order: V, T, K, U, W, Z)
# The order of elements of 'conf.num.vec' has been accordingly changed.

sample.graph.data <- function (number.of.T,
                            number.of.V.T = 1, # numeric vector of the numbers of V-nodes per T-node, or its average
                            family.number.of.V.T = NULL, # distribution of the numbers of V-nodes per T-node, number.of.V.T is then an average
                            conf.num.vec = rep(0, 4), # c("W","Z", "U","K")
                            graph_type = "scale-free", # not used when 'method' is supplied
                            method,
                            degree = 3,
                            connection_prob = 0.05, # not used now
                            mixed = "None", # not used now
                            theta = .5,
                            b0.1 = 0,
                            b.snp = 1,
                            b.med = 0.8 * b.snp,
                            sd.1 = 1,
                            neg.freq = 0.5,
                            conf.coef.ranges = list(W = c(0.15,0.5),
                                                    Z = c(1, 1.5),
                                                    U = c(0.15,0.5),
                                                    K = c(0.01, 0.1)),
                            scale = FALSE,
                            sample.size) {

  # Total number of nodes
  if (is.numeric(number.of.V.T)) {
    number.of.V.T = rep(number.of.V.T, length.out = number.of.T)
  }
  else {
    stop("'number.of.V.T' must be numeric")
  }
  if (!is.null(family.number.of.V.T)) {
    stop("not yet implemented")
  }
  number.of.V <- sum(number.of.V.T)
  number.nodes <-  number.of.V + number.of.T + sum(conf.num.vec)

  # Initialize the adjacency matrix
  A <- matrix(0, nrow = number.nodes, ncol = number.nodes)

  # Create all confounding variable names
  if(sum(conf.num.vec)>0){
    conf.node.names = NULL
    letter.id = c("W", "Z", "U", "K")
    nz.letter.id = letter.id[which(conf.num.vec>0)]
    conf.num.vec2 = conf.num.vec[which(conf.num.vec>0)]
    for(i in 1:length(conf.num.vec2)){
      conf.node.names = append(conf.node.names,
                               paste0(nz.letter.id[i], c(1:conf.num.vec2[i])))
    }
  }
  else
    conf.node.names = NULL
  row.names(A) = colnames(A) = c(paste0("V", c(1:number.of.V)),
                                 paste0("T", c(1:number.of.T)),
                                 conf.node.names)

  # Create the 'method' argument for 'get.custom.graph' starting from 'graph_type'
  if (missing(method)) {
    method <- switch(graph_type,
                     `scale-free` = "power",
                     `small-world` = "watts",
                     `random graph` = "er")
  }
  stopifnot(method %in% c("regular",
                          "watts",
                          "er",
                          "power",
                          "bipartite",
                          "barabasi",
                          "geometric",
                          "interEr"))

  # Generate the graph skeleton and the effects
  A = sim.graph.effects(Adj = A,
                         b.snp = b.snp,
                         b.med = b.med,
                         number.of.V.T = number.of.V.T,
                         number.of.T = number.of.T,
                         conf.coef.ranges = conf.coef.ranges,
                         neg.freq = neg.freq,
                         degree = degree,
                         method = method)

  # Save the effects
  effects.adj = A

  # Convert effects matrix to adjacency matrix
  A[A!=0] = 1

  # Create a graph object
  igraph.obj = igraph::graph_from_adjacency_matrix(A)
  graph.attr <- list(adjacency = A, effects.adj = effects.adj, igraph.obj = igraph.obj)

  # Initialize a data matrix
  X = as.data.frame(matrix(0, nrow = sample.size, ncol = dim(graph.attr$adjacency)[2]))
  colnames(X) = colnames(graph.attr$adjacency)

  # Get the topological ordering
  topo.order = colnames(graph.attr$adjacency)[as.vector(igraph::topo_sort(graph.attr$igraph.obj))]

  # Simulate data for each node
  for(i in 1:length(topo.order)){
    #generate V nodes in topo order
    location = match(topo.order[i], colnames(graph.attr$adjacency))
    parent.list = find.parents(Adjacency = graph.attr$adjacency, location = location)

    if(grepl("V", topo.order[i])){
      X[,location] = c(sample(c(0, 1, 2), size = sample.size, replace = TRUE,
                              prob = c((1 -theta)^2, 2 * theta * (1 - theta), theta^2)))
    }else{
      #catch nodes with no parents of any kind
      if(sum(unlist(lapply(parent.list, is.na)))==6){
        X[, location] = stats::rnorm(n = sample.size, mean = b0.1, sd = sd.1)
      }else{
        #simulate all other types of nodes according to parental list, topo.order, and the effects adj
        coefs = graph.attr$effects.adj[stats::na.omit(unlist(parent.list)), topo.order[i]]
        Cond.vals = as.matrix(X[, stats::na.omit(unlist(parent.list))])

   #     print(topo.order[i])

#        if (topo.order[i] == "T24") {
#          browser()
#       }

        X[, location] = stats::rnorm(n = sample.size, mean = (if (scale) 0 else b0.1) + Cond.vals %*% coefs, sd = sd.1)
        if (scale) {
          X[, location] = sd.1 * X[, location] / sd(X[, location]) + b0.1
        }
      }
    }
  }

  return(list(data = X,
              dims = list(p = number.of.T,
                          q = number.of.V,
                          r1 = conf.num.vec[1],
                          r2 = conf.num.vec[2],
                          u1 = conf.num.vec[3],
                          u2 = conf.num.vec[4]),
              sd.1 = sd.1,
              Adjacency = graph.attr$adjacency,
              Effects = graph.attr$effects.adj,
              igraph = graph.attr$igraph.obj))
}

# Internal function of MRGN (not exported from MRGN) copied on 04 December, 2023
find.parents <- function (Adjacency, location){
  #define letter identifier for node types
  letter.id = c("V","T","K","U","W","Z")
  all.parents = row.names(Adjacency)[which(Adjacency[, location] == 1)]
  #allocate each type of parent indexes to list by type
  parent.list = lapply(letter.id,function(x,y,z){if(any(grepl(x,y))) match(y[which(grepl(x,y))], z) else NA},
                       y = all.parents, z = colnames(Adjacency))
  names(parent.list)=c("V","T","K","U","W","Z")

  return(parent.list)
}
