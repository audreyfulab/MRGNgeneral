
#' @name sample.graph.data
#' @title Sample graph data with genetic variants
#' @description Different from \code{sim.graph.data} in the way \code{V}-nodes are
#' included in the network skeleton: here there are just added (not simulated).
#' The number of \code{V}-nodes is specified for each \code{T}-node, or sampled
#' from a count distribution. Each \code{V}-node is only related to only one \code{T}-node.
#'
#' @param n_t Integer, number of target phenotypes (T-nodes) in the network. Default is 100
#'
#' @param n_v.t Numeric vector of the numbers of V-nodes per T-node, or a scalar for 
#'   the average/fixed number. Default is 1
#'
#' @param family.n_v Distribution family for sampling the number of V-nodes per T-node 
#'   when \code{n_v.t} should be treated as an average. Default is NULL (not yet implemented)
#'
#' @param conf.num.vec Named numeric vector specifying the numbers of different confounder 
#'   types: W (intermediate variables), Z (common children), U (confounders), K (confounders 
#'   with specific properties), and I (independent confounders). 
#'   Default is c(W = 0, Z = 0, U = 200, K = 0, I = 100)
#'
#' @param graph_type Character, type of graph structure for T-nodes. One of "scale-free", 
#'   "small-world", or "random graph". Ignored when \code{method} is supplied. 
#'   Default is "scale-free"
#'
#' @param method Character, method for generating the random DAG, passed to 
#'   \code{pcalg::randDAG}. One of "regular", "watts", "er", "power", "bipartite", 
#'   "barabasi", "geometric", or "interEr". If missing, determined from \code{graph_type}
#'
#' @param degree Integer, average degree (number of connections) per node in the graph. 
#'   Default is 3
#'
#' @param theta Numeric in (0, 1), minor allele frequency for simulating V-nodes (variants). 
#'   Default is 0.4
#'
#' @param b0 Numeric, intercept/baseline value for simulating node values. Default is 0
#'
#' @param b.snp Numeric vector of length 2, range \code{c(min, max)} for effect sizes of 
#'   V-nodes (SNPs/variants). Default is c(-0.5, 0.5)
#'
#' @param b.med Numeric vector of length 2, range \code{c(min, max)} for effect sizes of 
#'   T-nodes (mediators/phenotypes). Default is c(-0.8, 0.8)
#'
#' @param sigma Numeric, standard deviation of random noise added to node values. Default is 0.1
#'
#' @param neg.freq Numeric in (0, 1), frequency of negative effects. Default is 0.5
#'
#' @param conf.coef.ranges Named list specifying ranges of coefficient values for different 
#'   confounder types (W, Z, U, K). Each element should be a vector c(min, max). 
#'   Default is list(W = c(0.15, 0.5), Z = c(1, 1.5), U = c(0.15, 0.5), K = c(0.15, 0.5))
#'
#' @param scale Logical, should simulated values be scaled to have marginal standard 
#'   deviation equal to \code{sigma}? Default is FALSE
#'
#' @param sample.size Integer, number of samples (observations/individuals) to generate
#'
#' @param seed scalar integer, seed for reproducibility of random number generation.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{data}: A data frame with simulated values for all nodes
#'   \item \code{dims}: A list with named counts of each node type (n_v, n_t, n_w, n_z, n_u, n_k, n_i)
#'   \item \code{b0}: The baseline/intercept value used
#'   \item \code{sigma}: The standard deviation used
#'   \item \code{adjacency}: The adjacency matrix of the generated graph
#'   \item \code{effects}: The effects adjacency matrix with coefficient values
#'   \item \code{igraph}: An igraph object of the generated graph
#' }
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
#
# Important change on Oct.-12-2023: reordered variables as V, T, W, Z, U, K (old order: V, T, K, U, W, Z)
# The order of elements of 'conf.num.vec' has been accordingly changed.
# Important change on Dec-06-2023: Adding I-nodes: independent confounding variables

sample.graph.data <- function (n_t = 100,
                               n_v.t = 1, # numeric vector of the numbers of V-nodes per T-node, or its average
                               family.n_v = NULL, # distribution of the numbers of V-nodes per T-node, n_v.t is then an average
                               conf.num.vec = c(W = 0, Z = 0, U = 200, K = 0, I = 100), # c("W","Z", "U","K", "I")
                               graph_type = "scale-free", # ignored when 'method' is supplied
                               method,
                               degree = 3, # connection_prob = 0.05, mixed = "None", # not used now
                               theta = .4,
                               b0 = 0,
                               b.snp = c(-0.5, 0.5), # b.v
                               b.med = c(-0.8, 0.8), # b.t
                               sigma = 0.1,
                               neg.freq = 0.5,
                               conf.coef.ranges = list(W = c(0.15, 0.5),
                                                       Z = c(1, 1.5),
                                                       U = c(0.15, 0.5),
                                                       K = c(0.15, 0.5)),
                               scale = FALSE,
                               sample.size,
                               seed = NULL) {

  # Save random generator state to restitute it later
  if(!is.null(seed)) {
    saved.seed <- .Random.seed
    set.seed(seed)
  }
  # Total number of nodes, except independent confounding variables
  if (is.numeric(n_v.t)) {
    n_v.t <- rep(n_v.t, length.out = n_t)
  }
  else {
    stop("'n_v.t' must be numeric")
  }
  if (!is.null(family.n_v)) {
    stop("not yet implemented")
  }
  n_v <- sum(n_v.t)
  stopifnot(is.numeric(conf.num.vec))
  if (length(conf.num.vec) == 5) {
    number.of.I <- conf.num.vec[5]
    conf.num.vec <- conf.num.vec[1:4]
  }
  else if (length(conf.num.vec) == 4) {
    number.of.I <- 0
  }
  else {
    stop("'conf.num.vec' have four or five element")
  }

  number.nodes <- n_v + n_t + sum(conf.num.vec)

  # Initialize the adjacency matrix
  A <- matrix(0, nrow = number.nodes, ncol = number.nodes)

  # Create all confounding variable names
  conf.node.names <- NULL
  if(sum(conf.num.vec)>0) {
    letter.id <- c("W", "Z", "U", "K")
    nz.letter.id <- letter.id[which(conf.num.vec>0)]
    conf.num.vec2 <- conf.num.vec[which(conf.num.vec>0)]
    for(i in 1:length(conf.num.vec2)){
      conf.node.names <- append(conf.node.names,
                                paste0(nz.letter.id[i], c(1:conf.num.vec2[i])))
    }
  }

  rownames(A) <- colnames(A) <- c(paste0("V", 1:n_v),
                                  paste0("T", 1:n_t),
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
  A <- sim.graph.effects(Adj = A,
                         b.snp = b.snp,
                         b.med = b.med,
                         n_v.t = n_v.t,
                         n_t = n_t,
                         conf.coef.ranges = conf.coef.ranges,
                         neg.freq = neg.freq,
                         degree = degree,
                         method = method)

  # Save the effects
  effects.adj <- A

  # Convert effects matrix to adjacency matrix
  A[A!=0] <- 1

  # Create a graph object
  igraph.obj <- igraph::graph_from_adjacency_matrix(A)
  graph.attr <- list(adjacency = A, effects.adj = effects.adj, igraph.obj = igraph.obj)

  # Initialize a data matrix
  X <- as.data.frame(matrix(0, nrow = sample.size, ncol = dim(graph.attr$adjacency)[2]))
  colnames(X) <- colnames(graph.attr$adjacency)

  # Get the topological ordering
  topo.order <- colnames(graph.attr$adjacency)[as.vector(igraph::topo_sort(graph.attr$igraph.obj))]

  # Simulate data for each node
  for(i in 1:length(topo.order)){
    #generate V nodes in topo order
    location <- match(topo.order[i], colnames(graph.attr$adjacency))
    parent.list <- find.parents(Adjacency = graph.attr$adjacency, location = location)

    if(grepl("V", topo.order[i])){
      X[,location] <- c(sample(c(0, 1, 2), size = sample.size, replace = TRUE,
                               prob = c((1 -theta)^2, 2 * theta * (1 - theta), theta^2)))
    }else{
      #catch nodes with no parents of any kind
      if(sum(unlist(lapply(parent.list, is.na)))==6){
        X[, location] <- stats::rnorm(n = sample.size, mean = b0, sd = sigma)
      }else{
        #simulate all other types of nodes according to parental list, topo.order, and the effects adj
        coefs <- graph.attr$effects.adj[stats::na.omit(unlist(parent.list)), topo.order[i]]
        Cond.vals <- as.matrix(X[, stats::na.omit(unlist(parent.list))])

        X[, location] <- stats::rnorm(n = sample.size, mean = (if (scale) 0 else b0) + Cond.vals %*% coefs, sd = sigma)

        # Scale each column to have marginal standard deviation sigma
        if (scale) {
          X[, location] <- sigma * X[, location] / sd(X[, location]) + b0
        }
      }
    }
  }

  # Add independent confounding variables if required
  if (number.of.I) {
    Imat <- matrix(rnorm(sample.size * number.of.I, mean = 0, sd = sigma),
                   nrow = sample.size, ncol = number.of.I)

    colnames(Imat) <- paste0("I", 1:number.of.I)

    X <- cbind(X, Imat)
  }

  # Restitute random generator state
  if(!is.null(seed)) {
    .Random.seed <- saved.seed
  }

  return(list(data = X,
              dims = list(n_v = c(V=as.numeric(n_v)),
                          n_t = c(T=as.numeric(n_t)),
                          n_w = c(W=as.numeric(conf.num.vec[1])),
                          n_z = c(Z=as.numeric(conf.num.vec[2])),
                          n_u = c(U=as.numeric(conf.num.vec[3])),
                          n_k = c(K=as.numeric(conf.num.vec[4])),
                          n_i = c(I=as.numeric(number.of.I))),
              b0 = b0,
              sigma = sigma,
              adjacency = graph.attr$adjacency,
              effects = graph.attr$effects.adj,
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
