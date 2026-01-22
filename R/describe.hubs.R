#'
#' Empirical assessment of hubs in inferred graphs
#'
#' Empirical description of hubs in inferred graph with known adjacency.
#' The function find hubs in a Direct Acyclic Graph (DAG), describe their correlation
#' structure based on a random sample for each node in the network, and summarize
#' confounder selection and network inference results for the hubs.
#'
#' @param adjacency numeric (binary) matrix, adjacency matrix of a DAG. The nodes
#' must be ordered such that the first \code{n_v} columns of \code{adjacency} are variants,
#' the next \code{n_t} columns are phenotypes, the next \code{n_q} columns are
#' intermediate variables or common children, and the next \code{n_u} columns (if any) are confounders.
#'
#' Note that \code{adjacency} must be square, binary with only zeros on the main diagonal.
#' Only the first \code{n_v+n_t+n_q+n_u} rows and columns are used.
#'
#' @param n_v,n_t,n_q,n_u integers, describe the types of nodes in the network
#' numbers of variants (\code{V}-nodes), phenotypes (\code{T}-nodes),
#' intermediate variables/common children (\code{Q}-nodes), and confounders (\code{U}-nodes).
#'
#' @param number.parents,number.children,number.neighbors integer arguments used to define a hub. See \link{find.hubs}.
#'
#' @param only.T,all.conditions logical arguments used to define a hub. See \link{find.hubs}.
#'
#' @param data \code{data.frame} object with a random sample for each node in
#' the network \code{adjacency} matrix to be used to describe the network.
#' @param conf.sets An object of class \code{'conf.sets'} containing the selected 
#'   confounders for each phenotype. See \link{get.conf.sets}
#'
#' @param inferred.adjacency Numeric (binary) matrix, the adjacency matrix of the 
#'   inferred network to be evaluated
#'
#' @param ... Additional arguments passed to \code{method} when \code{method} is a function
#'
#' @param method either a \code{character} or a \code{function} indicating the
#' association measure to use. When \code{method} is a character, One of "pearson"
#'
#' @param method either a \code{character} or a \code{function} indicating the
#' association measure to use. When \code{method} is a character, One of "pearson"
#' (default), "kendall", or "spearman" (can be abbreviated). A function \code{method}
#' must take in an input matrix \code{x} and return a corresponding partial
#' correlation matrix. Note that only the absolute value of the results are used.
#'
#' @param combine.nodes a function or a character giving the name of a function
#' (e.g. \code{combine.hubs} = 'mean') to be used to combine signals (absolute
#' association measures) over parent or child \code{T}-nodes (of a hub), if many.
#'
#' @param combine.hubs character or function, either \code{'largest'} (for computing the association
#' measures for the largest hub), a function or the name of a function (e.g. \code{combine.hubs} = 'mean')
#' to be used to combine signals over hubs, if many.
#'
#' @param neighbors character, when \code{combine.hubs = 'largest'} (otherwise
#' ignored), which neighbors are to be considered to select the largest hub?
#' One of \code{'parents'}, \code{'children'} and \code{all} (i.e. parents and
#' children).
#'
#' @param verbose logical, should warning messages from running \code{method} be
#' allowed? Defaults to \code{TRUE}.
#'
#' @param cl,chunk.size optional arguments for parallel computing, passed to
#' \link[parallel]{parLapply} (when supplied).
#'
#' @details
#' The function \code{describe.hubs} finds hub \code{T}-nodes in the true network
#' (by passing \code{adjacency} to \link{find.hubs}), and then computes for each hub:
#' \itemize{
#' \item the strength of the signal from parent \code{T}-nodes,
#' typically absolute Pearson's partial correlation to each parent \code{T}-node, given
#' all parents of the hub (including other \code{T}-nodes, parent \code{V} and
#' \code{U}-nodes), and then average them over the parents, if many ( *i*- \code{Parent.cor});
#' \item the strength of the signal from parent \code{T}-nodes
#' to each parent given the set of selected confounders for the hub (including
#' selected \code{V} and \code{T}-nodes, and \code{C}-nodes selected as
#' \code{U}-nodes), and then average them over the parents, if many ( *ii*- \code{Parent.cor.hat});
#' \item the strength of the signal to each child \code{T}-node
#' given all other parents of the child (including other \code{T}-nodes, parent \code{V} and
#' \code{U}-nodes), and then average them over the child \code{T}-nodes, if many ( *iii*- \code{Child.cor});
#' \item the strength of the signal to each child \code{T}-node
#' given the union of the sets of selected confounders for the hub and its child,
#' and then average them over the child \code{T}-nodes, if many ( *iv*- \code{Child.cor.hat});
#' \item recall and precision of confounder selection
#' (considering the set of true parent \code{V}, \code{T} and \code{U}-nodes as
#' the truth) ( *v*- \code{conf.rec} and *vi*- \code{conf.pre});
#' \item the recall of parent \code{T}-nodes, and the
#' related precision ( *vii*- \code{Parent.rec} and *viii*- \code{Parent.pre}); and
#' \item the recall of child \code{T}-nodes by network inference (in \code{inferred.adjacency}), and the
#' related precision ( *ix*- \code{Child.rec} and *x*- \code{Child.pre}).
#' }
#'
#' These 10 measures (*i* to *x*) are returned in a 10-column matrix named \code{raw.measures},
#' with one row for each found hub. If no parent or child node is found for a hub,
#' the corresponding measures are \code{NA}s.
#'
#' The association measure used when \code{method} is a character is \link[ppcor]{pcor}.
#' A named function \code{combine.hubs} must be available in the calling environment,
#' and return a numeric scalar output for a numeric vector input. For two elements
#' \code{a} and \code{b}, \code{combine.hubs(c(a, b))} should be in the closed \code{[a, b]}.
#' This is checked for \code{a = 0} and \code{b = 1}, and the function stops if
#' this fails.
#'
#' @return An object of class \code{'summary.hubs'} which is a list with the following elements:
#' \describe{
#' \item{measure}{ a numeric 10-element vector: column wise combinations of the
#' \code{raw.measures} matrix. The combination of the elements of a column of
#' \code{raw.measures} uses the function specified by the argument \code{combine.hubs}.}
#' \item{raw.measure}{ Either a numeric 10-column matrix with each row giving the
#' equivalent of \code{measure} per hub, or \code{NULL} when no hub is found or
#' when \code{combine.hubs = 'largest'} (see Section 'Details' for more on
#' \code{raw.measures}).}
#' \item{number.parents}{ a numeric vector with length the number of rows of \code{raw.measure},
#' giving the number of parent \code{T}-nodes for each found hub, or \code{NA} if no hub was found.}
#' \item{number.parents}{ a numeric vector with length the number of rows of \code{raw.measure},
#' giving the number of child \code{T}-nodes for each found hub, or \code{NA} if no hub was found.}
#' \item{hubs}{ an object of class \code{'hubs'} as returned by \link{find.hubs}.}
#' \item{call}{ the matched call to \code{describe.hubs}.}
#' }
#'
#' The class\code{'summary.hubs'} has a \code{print} method which shows a summary
#' of hubs' association with their neighbors, and recall and precision of
#' confounder selection and network inference.
#'
#' @export describe.hubs
#' @exportS3Method print summary.hubs
#'
#' @importFrom ppcor pcor
#'
#' @examples
#' ## Describe hubs in the 'networkA11' data
#' library (MRGNgeneral)
#' data ('networkA11')
#' data ('confsetsA11')
#' data ('mrgninferA11')
#'
#' # Using the average signal over parent/child nodes
#' Hubs <- describe.hubs (adjacency = networkA11$adjacency,
#'                        data = networkA11$data,
#'                        inferred.adjacency = mrgninferA11$adjacency,
#'                        conf.sets = confsetsA11,
#'                        n_t = 100, n_v = 100, n_q = 100, n_u = 200,
#'                        combine.nodes = mean,
#'                        combine.hubs = mean,
#'                        verbose = FALSE)
#' Hubs
#'
#' # Using the maximum signal over parent/child nodes instead of the mean
#' Hubs <- describe.hubs (adjacency = networkA11$adjacency,
#'                        data = networkA11$data,
#'                        inferred.adjacency = mrgninferA11$adjacency,
#'                        conf.sets = confsetsA11,
#'                        n_t = 100, n_v = 100, n_q = 100, n_u = 200,
#'                        combine.nodes = max,
#'                        combine.hubs = mean,
#'                        verbose = FALSE)
#' Hubs

# A function to:
# Return partial correlation given true confounders and V-nodes (done)
# Return partial correlation given selected confounders and V-nodes
# Return recall for conf selection for each hub
# Return network inference recall and precision for each hub
# ... are passed to 'method' if a function.
describe.hubs <- function (adjacency,                 # adjacency matrix of the true network
                         n_v = 0,                   # Number of variants
                         n_t = NCOL(adjacency),     # Number of phenotypes
                         n_q = NCOL(adjacency) - (n_t + n_v),       # Number of intermediate variables/common children
                         n_u = NCOL(adjacency) - (n_t + n_v + n_q), # Number of confounders
                         data,                      # Data frame to be used to describe the network
                         conf.sets,                 # 'conf.sets' object
                         inferred.adjacency,        # adjacency matrix of the inferred network
                         number.parents = if (all.conditions) 0 else Inf,   # T-node is a hub if more than 'number.parents' parents
                         number.children = if (all.conditions) 0 else Inf,  # T-node is a hub if more than 'number.children' children
                         number.neighbors = 10,     # T-node is a hub if more than 10 neighbors (parents+children)
                         only.T = TRUE, all.conditions = FALSE,          # Are all conditions required to declare a hub? Default is one condition is sufficient.
                         method = c('pearson', 'kendall', 'spearman'),
                         combine.nodes = c('mean', 'max'),
                         combine.hubs = c('mean', 'max', 'largest'),
                         neighbors = c('all', 'parents', 'children'),
                         verbose = TRUE,
                         cl = NULL, chunk.size = NULL, # Used to parallelize computation (if cl is not NULL) over hubs, if many
                         ...) { # Argument passed to method
  ## Save the call, not save arguments
  mcall <- match.call()

  ## Check arguments
  stopifnot(is.numeric(adjacency))
  stopifnot(NCOL(adjacency) == NROW(adjacency))
  stopifnot(all(adjacency %in% c(0, 1)))
  stopifnot(all(diag(adjacency) == 0))

  stopifnot(is.numeric(inferred.adjacency))
  stopifnot(NCOL(inferred.adjacency) == NROW(inferred.adjacency))
  stopifnot(all(inferred.adjacency %in% c(0, 1)))
  stopifnot(all(diag(inferred.adjacency) == 0))

  stopifnot(is.numeric(n_v), is.numeric(n_t), is.numeric(n_q), is.numeric(n_u))
  stopifnot(n_t > 0)
  n_v <- n_v[1]
  n_t <- n_t[1]
  n_q <- n_q[1]
  n_u <- n_u[1]
  stopifnot(NCOL(adjacency) >= n_v + n_t + n_q + n_u)

  data <- as.data.frame(data)
  stopifnot(is.data.frame(data))

  stopifnot(is.conf.sets(conf.sets))
  if (is.null(conf.sets$confounders))
    conf.sets$confounders <- mapply (function(x,y,z) c(x, y, z), conf.sets$Vconfounders, conf.sets$Tconfounders, conf.sets$Uconfounders, SIMPLIFY = FALSE)

  stopifnot(is.character(neighbors))
  neighbors <- neighbors[1]
  stopifnot(neighbors %in% c('all', 'parents', 'children'))

  stopifnot(is.character(method) | is.function(method))
  if (is.character(method)) {
    method <- method[1]
    stopifnot(method %in% c('pearson', 'kendall', 'spearman'))
    method.fun <- function (x) {
      ppcor::pcor (x = x, method = method)$estimate
    }
  }
  else {
    method.fun <- function (x) {
      method (x, ...)
    }
  }

  stopifnot(is.character(combine.nodes) | is.function(combine.nodes))
  if (is.character(combine.nodes)) {
    combine.nodes <- combine.nodes[1]
    combine.nodes <- get(combine.nodes, mode = "function", envir = parent.frame())
  }

  stopifnot(is.character(combine.hubs) | is.function(combine.hubs))
  Largest <- FALSE
  if (is.character(combine.hubs)) {
    combine.hubs <- combine.hubs[1]
    if (!identical(tolower(combine.hubs), 'largest')) {
      combine.hubs <- get(combine.hubs, mode = "function", envir = parent.frame())
    }
    else {
      Largest <- TRUE
    }
  }

  ### Find hubs in the network
  hubs <- find.hubs(adjacency,
                    n_v = n_v, n_t = n_t, n_q = n_q,
                    number.parents = number.parents,
                    number.children = number.children,
                    number.neighbors = number.neighbors,
                    only.T = only.T, all.conditions = all.conditions,
                    cl = cl, chunk.size = chunk.size)
  nb.hubs <- length(hubs$hubs)
  only.T <- hubs$call$only.T

  ### Stop if no hub
  if (nb.hubs == 0) {
    # write code for returved value when there is no hub

    # Build output to return
    out <- structure(list(measure = rep(NA, 10),
                          number.parents = NA,
                          number.children = NA,
                          raw.measures = NULL,
                          hubs = hubs,
                          call = mcall),
                     class = 'summary.hubs')

    return(out)
  }

  ### Count parents and children per hub
  if (only.T) {
    Nbparents <- hubs$number.neighbors[hubs$hubs,2]
    Nbchildren <- hubs$number.neighbors[hubs$hubs,5]
  }
  else {
    if (nb.hubs == 1) {
      Nbparents <- sum(hubs$number.neighbors[hubs$hubs,1:4])
      Nbchildren <- sum(hubs$number.neighbors[hubs$hubs,5:6])
    }
    else {
      Nbparents <- rowSums(hubs$number.neighbors[hubs$hubs,1:4, drop = FALSE])
      Nbchildren <- rowSums(hubs$number.neighbors[hubs$hubs,5:6, drop = FALSE])
    }
  }

  ### Sets of nodes
  Vindex <- if (n_v > 0) 1:n_v
  Tindex <- (n_v + 1):(n_v + n_t)
  Uindex <- if (n_u > 0) (n_v + n_t + n_q + 1):(n_v + n_t + n_q + n_u)

  ### Compute correlations for hubs
  if (Largest) {
    ### Count all considered neighbors per hub
    switch(neighbors,
           all = {
             NbNeighbors <- Nbparents + Nbchildren
           },
           parents = {
             NbNeighbors <- Nbparents
           },
           children = {
             NbNeighbors <- Nbchildren
           })

    # Find the T-node hub to keep
    keephubs <- which.max(NbNeighbors)

    # Compute parent-hub and hub-children associations (correlations) for the hub
    measure <- get.hub.measures (n_v + hubs$hubs[keephubs],
                                 adjacency = adjacency,
                                 n_v = n_v, n_t = n_t, n_q = n_q, n_u = n_u,
                                 data = data,
                                 conf.sets = conf.sets,
                                 inferred.adjacency = inferred.adjacency,
                                 Vindex = Vindex, Tindex = Tindex, Uindex = Uindex,
                                 method = method.fun,
                                 combine.nodes = combine.nodes,
                                 verbose = verbose)

    # Build output to return
    out <- structure(list(measure = measure,            # a vector of two elements
                          number.parents = Nbparents,   # number of parents used to compute the first element of \code{measure}
                          number.children = Nbchildren, #
                          raw.measures = NULL,
                          hubs = hubs,
                          call = mcall),
                     class = 'summary.hubs')

  }
  else {
    ##
    raw.measures <- matteSapply(n_v + hubs$hubs,
                                FUN = get.hub.measures,
                                adjacency = adjacency,
                                n_v = n_v, n_t = n_t, n_q = n_q, n_u = n_u,
                                data = data,
                                conf.sets = conf.sets,
                                inferred.adjacency = inferred.adjacency,
                                Vindex = Vindex, Tindex = Tindex, Uindex = Uindex,
                                method = method.fun,
                                combine.nodes = combine.nodes,
                                verbose = verbose,
                                cl = cl, chunk.size = chunk.size)

    if (nb.hubs > 1) {
      # Test the function 'combine.hubs'
      rr <- combine.hubs (0:1)
      stopifnot((length(rr) == 1) & is.numeric(rr) & !is.na(rr) & (0 <= rr) & (rr <= 1))

      measure <- apply(raw.measures,
                       MARGIN = 1,
                       FUN = function(x) combine.hubs (na.omit(x)))

      # Build output to return
      out <- structure(list(measure = measure,            # a vector of 10 elements
                            number.parents = Nbparents,
                            number.children = Nbchildren,
                            raw.measures = t(raw.measures),
                            hubs = hubs,
                            call = mcall),
                       class = 'summary.hubs')

    }
    else {
      # Build output to return
      out <- structure(list(measure = raw.measures,
                            number.parents = Nbparents,
                            number.children = Nbchildren,
                            raw.measures = raw.measures,
                            hubs = hubs,
                            call = mcall),
                       class = 'summary.hubs')
    }
  }

  # Save some call elements used for printing 'summary.hubs' objects
  out$args <- list(method = method, combine.nodes = combine.nodes,
                   combine.hubs = combine.hubs, neighbors = neighbors)

  return(out)
}

## Compute descriptive and performance measures for hub Tj
get.hub.measures <- function (j, # column number of a T-node hub in adjacency
                              adjacency,
                              n_v, n_t, n_q, n_u,
                              data,
                              conf.sets,
                              inferred.adjacency,
                              Vindex = if (n_v > 0) 1:n_v,
                              Tindex = (n_v + 1):(n_v + n_t),
                              Uindex = if (n_u > 0) (n_v + n_t + n_q + 1):(n_v + n_t + n_q + n_u),
                              method, combine.nodes,
                              verbose = TRUE) {
  stopifnot(j >= n_v, j <= n_v + n_t)
  #### Parents - Hub measures
  ### Find hub T-parent nodes
  Tparents <- which(adjacency[Tindex, j] == 1)
  nb.parents <- length(Tparents)

  ### Find hub parents that are not T-nodes, nor Q-nodes
  if (nb.parents == 0) {
    ParentCorr <- HatParentCorr <- NA
    conf.perf <- ParentInf.perf <- c(NA, NA)
  }
  else {
    ## make Tparents column index f
    Tparents <- Tparents + n_v

    ## Compute correlation given true V and U parent nodes
    # Find all V,U parents of Tj
    VUparents <- if (n_u) which(adjacency[Uindex, j] == 1)
    if (length(VUparents))
      VUparents <- VUparents + n_v + n_t + n_q

    VUparents <- c(if (n_v) which(adjacency[Vindex, j] == 1),
                   VUparents)

    # Get the measure
    if (verbose) {
      ParentCorr <- method (data[, c(j, Tparents, VUparents), drop = FALSE])[1, 2:(nb.parents + 1)]
    }
    else {
      ParentCorr <- catch.conditions({
        method (data[, c(j, Tparents, VUparents), drop = FALSE])[1, 2:(nb.parents + 1)]
      })$value
    }
    ParentCorr <- abs(ParentCorr)

    ## Compute correlation given selected confounders
    if (verbose) {
      HatParentCorr <- sapply(Tparents, FUN = hatcorfun,
                              j = j, conf.sets = conf.sets,
                              data = data, n_v = n_v, method = method)

    }
    else {
      HatParentCorr <- catch.conditions({
        sapply(Tparents, FUN = hatcorfun,
               j = j, conf.sets = conf.sets,
               data = data, n_v = n_v, method = method)
      })$value
    }
    HatParentCorr <- abs(HatParentCorr)

    ## Compute recall and precision of confounder selection for Tj
    conf.perf <- c(
      mean(c(Tparents, VUparents) %in% conf.sets$confounders[[j - n_v]]), # recall
      mean(conf.sets$confounders[[j - n_v]] %in% c(Tparents, VUparents))  # precision
    )

    ## Compute recall and precision of network inference for Tj hub (Tj + parents)
    hub_parent_adj <- rbind(c(0, adjacency[j, Tparents]), #Tj -> parents
                            cbind(adjacency[Tparents, j], matrix(0, nrow = nb.parents, ncol = nb.parents)))

    inferred_hub_parent_adj <- rbind(c(0, inferred.adjacency[j, Tparents]), #Tj -> parents
                            cbind(inferred.adjacency[Tparents, j], matrix(0, nrow = nb.parents, ncol = nb.parents)))
    dimnames(hub_parent_adj) <- dimnames(inferred_hub_parent_adj) <-
      list(paste0("G", 0:nb.parents), paste0("G", 0:nb.parents))
    ParentInf.perf <- RecallPrecision(g1 = as (hub_parent_adj, 'graphNEL'),
                                      g2 = as (inferred_hub_parent_adj, 'graphNEL'),
                                      GV = 0, includeGV = FALSE,
                                      edge.presence = 1, edge.direction = 0.5)
    ParentInf.perf <- c(
      ParentInf.perf$Recall, # recall
      ParentInf.perf$Precision # precision
    )

    ## Combine over all parent T-nodes of Tj
    if (nb.parents > 1) {
      ParentCorr <- combine.nodes (ParentCorr)
      HatParentCorr <- combine.nodes (HatParentCorr)
    }
  }

  #### Hub - Children measures
  ### Find hub child-parent nodes
  Tchildren <- which(adjacency[j, Tindex] == 1)
  nb.children <- length(Tchildren)

  ### Find hub parents that are not T-nodes, nor Q-nodes
  if (nb.children == 0) {
    ChildCorr <- HatChildCorr <- NA
    ChildInf.perf <- c(NA, NA)
  }
  else {
    ## make Tchildren column index
    Tchildren <- Tchildren + n_v

    ## For each child node of Tj, find other parents (V,T,U) and condition the correlation with Tj on them
    if (verbose) {
      ChildCorr <- sapply(Tchildren, FUN = childcorfun,
                          j = j, adjacency = adjacency,
                          n_v = n_v, n_t = n_t, n_q = n_q,
                          Vindex = Vindex, Tindex = Tindex, Uindex = Uindex,
                          data = data, method = method)
    }
    else {
      ChildCorr <- catch.conditions({
        sapply(Tchildren, FUN = childcorfun,
               j = j, adjacency = adjacency,
               n_v = n_v, n_t = n_t, n_q = n_q,
               Vindex = Vindex, Tindex = Tindex, Uindex = Uindex,
               data = data, method = method)
      })$value
    }
    ChildCorr <- abs(ChildCorr)

    ## Compute correlation given selected confounders
    if (verbose) {
      HatChildCorr <- sapply(Tchildren, FUN = hatcorfun,
                             j = j, conf.sets = conf.sets,
                             data = data, n_v = n_v, method = method)
    }
    else {
      HatChildCorr <- catch.conditions({
        sapply(Tchildren, FUN = hatcorfun,
               j = j, conf.sets = conf.sets,
               data = data, n_v = n_v, method = method)
      })$value
    }
    HatChildCorr <- abs(HatChildCorr)

    ## Compute recall and precision of network inference for Tj hub (Tj + children)
    hub_child_adj <- rbind(c(0, adjacency[j, Tchildren]), #Tj -> children
                            cbind(adjacency[Tchildren, j], matrix(0, nrow = nb.children, ncol = nb.children)))

    inferred_hub_child_adj <- rbind(c(0, inferred.adjacency[j, Tchildren]), #Tj -> parents
                                     cbind(inferred.adjacency[Tchildren, j], matrix(0, nrow = nb.children, ncol = nb.children)))
    dimnames(hub_child_adj) <- dimnames(inferred_hub_child_adj) <-
      list(paste0("G", 0:nb.children), paste0("G", 0:nb.children))
    ChildInf.perf <- RecallPrecision(g1 = as (hub_child_adj, 'graphNEL'),
                                     g2 = as (inferred_hub_child_adj, 'graphNEL'),
                                     GV = 0, includeGV = FALSE,
                                     edge.presence = 1, edge.direction = 0.5)
    ChildInf.perf <- c(
      ChildInf.perf$Recall,
      ChildInf.perf$Precision
    )

    ## Combine over all child T-nodes of Tj
    if (nb.children > 1) {
      ChildCorr <- combine.nodes (ChildCorr)
      HatChildCorr <- combine.nodes (HatChildCorr)
    }
  }

  return(c(Parent.cor = ParentCorr,
           Parent.cor.hat = HatParentCorr,
           Child.cor = ChildCorr,
           Child.cor.hat = HatChildCorr,
           conf.rec = conf.perf[1],
           conf.pre = conf.perf[2],
           Parent.rec = ParentInf.perf[1],
           Parent.pre = ParentInf.perf[2],
           Child.rec = ChildInf.perf[1],
           Child.pre = ChildInf.perf[2]))
}

#### Internal to compute correlation given true parent nodes
## For each child node of Tj, find other parents (V,T,U) and condition the correlation with Tj on them
childcorfun <- function (k, j, adjacency,
                         n_v, n_t, n_q,
                         Vindex, Tindex, Uindex,
                         data, method) {

  # Find all V,U,T parents of the child Tk
  Uk <- which(adjacency[Uindex, k] == 1)
  if (length(Uk))
    Uk <- Uk + n_v + n_t + n_q
  Tk <- which(adjacency[Tindex, k] == 1)
  if (length(Tk))
    Tk <- Tk + n_v
  Vk <- which(adjacency[Vindex, k] == 1)
  VUTk <- c(Vk, Tk, Uk)

  # Have the parent Tj in first position
  VUTk <- unique(c(j, VUTk))

  # Get the measure
  return(method (data[, c(k, VUTk), drop = FALSE])[1, 2])
}

#### Internal to compute correlation given selected confounders
hatcorfun <- function (k, j, conf.sets, data, n_v, method) {
  # Extract confounders selected for Tj
  Tjconfs <- conf.sets$confounders[[j - n_v]]

  # Extract confounders selected for Tk
  Tkconfs <- conf.sets$confounders[[k - n_v]]

  # Get the measure
  method (data[, unique(c(j, k, Tjconfs, Tkconfs)), drop = FALSE])[1, 2]
}

# Print an object of class 'summary.hubs'
#' @export
print.summary.hubs <- function (x, digits = max(3, getOption("digits") - 3),
                                ...) {
  # Print the call tha generated 'x'
  if (!is.null(x$call)) {
    cat("\n")
    cat("Call: ")
    print(x$call)
  }

  # Print 'x'
  cat("\n")
  nb.hubs <- length(x$hubs$hubs)
  if (nb.hubs == 0) {
    cat(" No hub found.\n")
    cat("\n")
    return(invisible(x))
  }

  cat(  paste0('Number of hubs:                    ', nb.hubs))
  cat(paste0('\nSignal from parents:               ', format(x$measure[2], digits = digits), '        (vs truth: ', format(x$measure[1], digits = digits), ')'))
  cat(paste0('\nSignal to children:                ', format(x$measure[4], digits = digits), '        (vs truth: ', format(x$measure[3], digits = digits), ')'))
  cat(paste0('\nConfounder (recall, precision):   (', format(x$measure[5], digits = digits), ', ', format(x$measure[6], digits = digits), ')'))
  cat(paste0('\nParent edges (recall, precision): (', format(x$measure[7], digits = digits), ', ', format(x$measure[8], digits = digits), ')'))
  cat(paste0('\nChild edges (recall, precision):  (', format(x$measure[9], digits = digits), ', ', format(x$measure[10], digits = digits), ')\n'))
  method <- x$call$method
  if (is.null(method))
    method <- 'pearson'
  combine.hubs <- x$call$combine.hubs
  if (is.null(combine.hubs))
    combine.hubs <- 'largest'

  cat(paste0(' (combine.hubs: ', as.character(combine.hubs), ') \n'))
  cat("\n")

  if (is.null(x$raw.measures)) {
    return(invisible(x))
  }

  cat('Parents: \n')
  df <- data.frame(Nb.parents = as.numeric(x$hubs$number.neighbors[x$hubs$hubs, 2]),
                   `True.Signal` = as.numeric(x$raw.measures[,1]),
                   `Used.Signal` = as.numeric(x$raw.measures[,2]),
                   Recall = as.numeric(x$raw.measures[,7]),
                   Precision = as.numeric(x$raw.measures[,8]),
                   Conf.recall = as.numeric(x$raw.measures[,5]),
                   Conf.precision = as.numeric(x$raw.measures[,6])
                   )
  rownames(df) <- names(x$hubs$hubs)
  print(df[df[,1] > 0,], digits = digits, ...)

  cat('\nChildren: \n')
  df <- data.frame(Nb.children = as.numeric(x$hubs$number.neighbors[x$hubs$hubs, 5]),
                   `True.Signal` = as.numeric(x$raw.measures[,3]),
                   `Used.Signal` = as.numeric(x$raw.measures[,4]),
                   Recall = as.numeric(x$raw.measures[,9]),
                   Precision = as.numeric(x$raw.measures[,10])
  )
  rownames(df) <- names(x$hubs$hubs)
  print(df[df[,1] > 0,], digits = digits, ...)

  combine.nodes <- x$call$combine.nodes
  if (is.null(combine.nodes))
    combine.nodes <- 'mean'
  cat(paste0(' (measure: ', as.character(method),
             ', combine.nodes: ', as.character(combine.nodes), ')'))

  cat("\n")
  return(invisible(x))
}

