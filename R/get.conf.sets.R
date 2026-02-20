
# Procedure to be tested
# Select U-nodes (partial correlations given all V-nodes [or just one or two from eQTL analysis per T-node])
# Select T-nodes (partial correlations given selected U-nodes and one or two V-nodes from eQTL analysis per T-node)
# Select V-nodes (partial correlations given selected U-nodes, and selected T-nodes)

#' Select confounding variables for a genomic network
#'
#' Build a set of confounding variables for genes in a genomic network including
#' genes (e.g. expression values), genetic variants, and confounding variables (e.g. sex, age,
#' PC scores from whole-genome expression) by testing marginal or partial associations
#' of each potential confounding variable and each gene. For each gene, the function
#' finds variants and/or other genes (in the network) that can confound the
#' regulatory mechanisms between a particular gene and another target gene.
#'
#' @param data \code{matrix} or \code{dataframe} of size \code{N} samples by \code{m} variables
#' where the first \code{n_v} columns represents genetic variants, the following
#' \code{n_t} columns represent target genes (e.g. expression values), and the next \code{n_c} columns
#' represent potential confounding variables (e.g. sex, age, other genes).
# If \code{conditional.vars} is specified, then all the columns indexed in \code{conditional.vars} must also
# be present in \code{data}.
#' Hence \code{m} must equal \code{n_t+n_v+n_c} or more.
#'
#' @param scale.data logical scalar, should the data be scaled before processing?
#' Scaling refers here to standardization, \code{i.e.} subtracting the sample mean
#' and dividing by the sample standard deviation as achieved by \link{scale}. Only
#' genes and confounders are scaled, variants are not scaled even when
#' \code{scale.data = TRUE} (the default).
#'
#' @param n_v,n_t,n_c numeric scalars, numbers of respectively \code{T}-nodes (genes/expression),
#' \code{V}-nodes (variants) and \code{C}-nodes (candidate confounding variables) in the
#' genomic network. The number of columns of \code{data} must be \code{n_t+n_v+n_c} or more.
#'
#' @param skip.V,skip.T,skip.C logical, should the selection of respectively
#' \code{V}-nodes, \code{T}-nodes, or \code{C}-nodes be skipped? If \code{TRUE},
#' the corresponding selection is skipped and an empty list is returned. This can be
#' useful when a dataset is expanded with only one or two types of nodes, and there
#' is no need to select again confounding variables from the set(s) that remained unchanged.
#'
#' @param T.measure,C.measure characters, indicating the association measures to
#' be used for different selections. One of \code{"partial"} for conditional
#' Pearson's correlation given \code{V}-nodes, and \code{"marginal"} for correlation
#' coefficient. The measure used for \code{V}-nodes selection is always the
#' marginal correlation coefficient.
#'
#' @param FDRcontrol,V.FDRcontrol,T.FDRcontrol,C.FDRcontrol,Q.FDRcontrol characters indicating
#' the FDR control methods to be used for different selections.
#' One of \code{"none"}, \code{"qvalue"} (see \link[qvalue]{qvalue}), or
#' \code{"bonferroni"} (see \link[stats]{p.adjust}).
#' If any of \code{T.FDRcontrol}, \code{C.FDRcontrol}, \code{V.FDRcontrol}, \code{Q.FDRcontrol} is missing,
#' \code{FDRcontrol} is used, otherwise \code{FDRcontrol} is ignored.
#'
#' @param adjust_by,V.adjust_by,T.adjust_by,C.adjust_by,Q.adjust_by character indicating the
#' adjustment scheme for tests. One of \code{"none"} (no adjustment is desired),
#' \code{"individual"} (adjust p-values for each gene separately),
#' and \code{"all"} (adjust all p-values for all genes at once).
#' If any of \code{T.adjust_by}, \code{C.adjust_by}, \code{V.adjust_by}, \code{Q.adjust_by} is missing,
#' \code{adjust_by} is used, otherwise \code{adjust_by} is ignored.
#'
#' @param T.filter logical, should \code{T}-node selection attempt to filter out
#' potential child \code{T}-nodes from the selected pool for each \code{T}-node?
#' Defaults to \code{FALSE}.
#'
#' @param fdr,lambda,pi0.method,alpha See \link[MRGN]{get.conf.matrix}
#' (\code{selection_fdr} is used there for \code{fdr}).
#'
#' @param parallel logical, should computations be parallelized?
#'
#' @param cl,chunk.size a cluster object (\code{cl}) for parallel computations,
#' and a numeric scalar (\code{chunk.size}) to schedule parallel tasks.
#' Only used if \code{parallel = TRUE}.
#'
#' @param blocksize integer, block size for computing correlation matrices in chunks
#'   using \link[propagate]{bigcor}. Default is 100. Smaller values use less memory
#'   but may be slower
#'
#' @param verbose integer, verbosity level. 0 (default) for silent operation,
#'   higher values for more detailed progress messages
#'
#' @param seed integer, random seed for reproducible results in parallel computing.
#'   If NULL (default), no seed is set
#' @param save.list (logical) if TRUE the output is saved as a \code{.RData}
#' object (default = FALSE).
#'
#' @param save.path string specifying the path name of the output list to be save
#' as a \code{.RData} structure. Only used if \code{save.list = TRUE}.
#'
#' @details
# This function is a wrapper for \link[MRGN]{get.conf.matrix}.
#' For a graph \code{G(V, T, C)} where \code{V} is the set of genetic variants,
#' \code{T} is the set of gene expressions, and \code{C} is the set of candidate
#' confounding variables in the network; \code{get.conf.sets} performs
#' three independent tasks.
#' \describe{
#'   \item{V-node selection:}{ the *marginal* association between each gene (\code{T}-node) and
#'   each variant (\code{V}-node) is tested using a *t*-test. If none of \code{V.adjust_by}
#'   and \code{adjust_by} is \code{'none'}, then either all the tests involving a
#'   particular \code{T}-node (when \code{adjust_by = 'individual'}), or all the tests
#'   involving all \code{T}-nodes (when \code{adjust_by = 'all'}) are adjusted for
#'   multiple comparison using the method \code{V.FDRcontrol}. All \code{V}-nodes
#'   significantly associated with a particular \code{T}-node are then retained as
#'   selected confounders for the \code{T}-node.}
#'   \item{T-node selection:}{ the process is similar to *V-node selection* except
#'   that here, the associations between each \code{T}-node and all other \code{T}-nodes
#'   in the network are tested, and the association measure between two \code{T}-nodes can
#'   be partial (the default) or marginal correlation (specified via \code{T.measure}).
#'   When \code{T.measure = "partial"}, partial correlations are computed given all
#'   \code{V}-nodes in the network, and *z*-tests are used instead of *t*-tests.
#'   When sample size is smaller than V-nodes + 3, selected V-nodes are used instead.}
#'   \item{C-node selection:}{ this step has two sub-steps. The first sub-step is a pre-selection
#'   which is similar to *T-node selection*, except that the association of each \code{T}-node
#'   with each \code{C}-node is tested. The results is a pool of confounding variables
#'   that can include true confounders (\code{U}-nodes) on one hand (\code{U}-nodes
#'   are marginally uncorrelated with the \code{V}-nodes in the network), and
#'   intermediate variables (\code{W}-node) and common children (\code{Z}-node)
#'   on the other hand (\code{W} and \code{Z}-nodes are marginally correlated
#'   with some \code{V}-nodes). The second sub-step thus tests for marginal
#'   association between all \code{V}-nodes and each pre-selected \code{C}-node,
#'   to assign the later as a confounder (\code{U}-node) or not (\code{W} or
#'   \code{Z}-node).}
#' }
#'
#' Correlation matrices are computed using \link[propagate]{bigcor} which is
#' more efficient in large datasets \insertCite{Spiess2018propagate}{MRGNgeneral}.
#' Partial correlations are obtained using \link[ppcor]{pcor} \insertCite{kim2015ppcor}{MRGNgeneral}.
#' The \code{C}-node selection step follows
#' \insertCite{yang2017identifying}{MRGNgeneral}.
# \insertCite{yang2017identifying;textual}{MRGNgeneral}
#'
#' @return an object of class \code{'conf.sets'}, i.e. a named list with elements:
#' \describe{
#' \item{\code{Vconfounders}}{a list of length \code{n_t} giving the vector of
#' \code{V}-nodes selected as confounders for each of the \code{n_t}
#' \code{T}-nodes in the network.}
#' \item{\code{Tconfounders}}{a list of length \code{n_t} giving the vector of
#' \code{T}-nodes selected as confounders for each of the \code{n_t}
#' \code{T}-nodes in the network.}
#' \item{\code{Uconfounders}}{a list of length \code{n_t} giving the vector of
#' \code{U}-nodes selected as confounders for each of the \code{n_t}
#' \code{T}-nodes in the network.}
#' \item{\code{confounders}}{a list of length \code{n_t} obtained as
#' \code{lbind(Vconfounders, Tconfounders, Uconfounders)}.}
#' \item{\code{UWZconfounders}}{a list of length \code{n_t} giving the vector of
#' \code{U,W,Z}-nodes selected as confounders for each of the \code{n_t}
#' \code{T}-nodes in the network.}
#' \item{\code{WZconfounders}}{a list of length \code{n_v} giving the vector of
#' \code{W,Z}-nodes identified as associated to each of the \code{n_v}
#' \code{V}-nodes in the network.}
#' \item{\code{UWZindices}}{a vector giving the pool of all selected
#' \code{U,W,Z}-nodes.}
#' \item{\code{WZindices}}{a vector giving the pool of all selected
#' \code{W,Z}-nodes.}
#' \item{\code{raw}}{a list of the raw results (as returned by
#' \link[MRGN]{get.conf.matrix}) for each selection step.}
#' \item{\code{time}}{a list of the CPU time (as returned by
#' \link{system.time}) for each selection step.}
#' \item{new.order}{ \code{NULl} (not returned by this function),
#' but can be a numeric vector of new column-orders used to alter other slots
#' (except \code{time}), see \link{reorder_conf_sets}.}
#' }
#'
#' @export get.conf.sets
#'
# @importFrom MRGN get.conf.matrix
#' @importFrom MRGN p.from.parcor
#' @importFrom MRGN p.from.cor
#' @importFrom MRGN adjust.q
#' @importFrom stats p.adjust
#' @importFrom ppcor pcor
#'
#' @seealso \link{assess.conf.selection} to evaluate the performance of the
#' selection procedure given the adjacency matrix of the true network.
#
# \link[MRGN]{get.conf.matrix}.
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#' ## Simulate some data with 20 phenotypes
#' set.seed(167)
#' net20data <- sample.graph.data (n_t = 20,
#'                                 n_v.t = 1,
#'                                 family.n_v = NULL,
#'                                 conf.num.vec = c(W = 10, Z = 10,
#'                                                  U = 40, K = 0, I = 20),
#'                                 graph_type = "scale-free",
#'                                 degree = 3,
#'                                 theta = .4,
#'                                 b0 = 0,
#'                                 b.snp = c(-0.5, 0.5),
#'                                 b.med = c(-0.8, 0.8),
#'                                 sigma = 0.1,
#'                                 neg.freq = 0.5,
#'                                 conf.coef.ranges = list(W = c(0.4, 0.5),
#'                                                         Z = c(1, 1.5),
#'                                                         U = c(0.4, 0.5),
#'                                                         K = c(0.01, 0.1)),
#'                                 scale = FALSE,
#'                                 sample.size = 100)
#'
#' ## Performing confounding variable selection
#' confnet20 <- get.conf.sets(data = net20data$data,
#'                            n_v = net20data$dims$n_v,
#'                            n_t = net20data$dims$n_t,
#'                            n_c = NCOL(net20data$data) - net20data$dims$n_v - net20data$dims$n_t,
#'                            blocksize = 10,
#'                            T.measure = 'partial',
#'                            C.measure = 'partial',
#'                            FDRcontrol = 'qvalue',
#'                            adjust_by = 'individual',
#'                            alpha = 0.01,
#'                            fdr = 0.05,
#'                            lambda = 0.05,
#'                            pi0.method = 'smoother')
#'
#' ## Recall and precision of the selection procedure
#' Perf <- assess.conf.selection (confnet20,
#'                                adjacency = net20data$adjacency,
#'                                n_v = net20data$dims$n_v,
#'                                n_t = net20data$dims$n_t,
#'                                n_w = net20data$dims$n_w,
#'                                n_z = net20data$dims$n_z,
#'                                n_u = net20data$dims$n_u)
#' Perf$recall
#' Perf$precision
#'

get.conf.sets <- function (data,
                           scale.data = TRUE,
                           n_v,                      # Number of Variants
                           n_t,                      # Number of Genes
                           n_c,                      # Number of variables in the pool of candidate confounding variables (to select from)
                           skip.V = FALSE,           # Skip V-node selection?
                           skip.T = FALSE,
                           skip.C = FALSE,
                           T.measure = c("partial", "marginal"), # Selection measure
                           C.measure = T.measure,
                           blocksize = 100,
                           FDRcontrol = c("qvalue", "bonferroni", "none"),
                           V.FDRcontrol = FDRcontrol,
                           T.FDRcontrol = FDRcontrol,
                           C.FDRcontrol = FDRcontrol,
                           Q.FDRcontrol = FDRcontrol,
                           adjust_by = c("individual", "all", "none"),
                           V.adjust_by = adjust_by,
                           T.adjust_by = adjust_by,
                           C.adjust_by = adjust_by,
                           Q.adjust_by = adjust_by,
                           T.filter = FALSE,
                           fdr = 0.05, # Only used if FDRcontrol = 'qvalue'
                           lambda = 0.05,
                           pi0.method = c("smoother", "boostrap"),
                           alpha = 0.01, # Not used if FDRcontrol = 'qvalue'
                           parallel = FALSE,
                           cl = parallel::getDefaultCluster(),
                           chunk.size = NULL, # scalar number; number of invocations of fun or FUN in one chunk; a chunk is a unit for scheduling.
                           verbose = 0,
                           save.list = FALSE, # Only used if 'save.path' is not missing
                           save.path = "/path/to/save/location/",
                           seed = NULL) { # seed for reproducible results in parallel computing
  ### Save the call
  mcall <- match.call()

  ### Check arguments
  eval(check.get.conf.sets.args())

  ### Get pools of variables (T.pool, C.pool, V.pool)
  eval(get.candidate.pools())

  ### Select V-nodes
  VTime <- system.time({
    if (skip.V | n_v == 0) {
      ## Return empty list if no V-node
      Vraw <- NULL
      Vconfounders <- vector(mode = "list", length = n_t)
    }
    else {
      if (verbose) {
        cat("\n        * selecting 'V-nodes' using marginal correlations ... \n")
      }

      ### Call 'get.conf.matrix' on 'V.pool'
      Vraw <- catch.conditions({
        get.conf.matrix (data = data[, T.pool, drop = FALSE],
                         cov.pool = data[, V.pool, drop = FALSE],
                         measure = "correlation",
                         conditional.vars = NULL,
                         blocksize = blocksize,
                         selection_fdr = fdr,
                         adjust_by = V.adjust_by,
                         apply.qval = V.apply.qval,
                         lambda = lambda,
                         pi0.method = pi0.method,
                         alpha = alpha,
                         verbose = max(0, verbose - 1),
                         cl = cl, chunk.size = chunk.size)
      })$value

      ## Check the result: error or not?
      if (any(class(Vraw) %in% c("simpleError", "error", "condition"))) {
        if (identical(Vraw$message, "wrong sign in 'by' argument")) {
          if (verbose) {
            print(Vraw$message)
            cat("\n 'get.conf.matrix' failled ...  \n")
          }
        }

        # Change pi0.method to 'bootstrap' if it was 'smoother' and 'get.conf.matrix' failed
        if (identical(pi0.method, "smoother")) {
          if (verbose) {
            print(Vraw)
            cat("\n pi0.meth = 'smoother' failed, trying pi0.meth = 'bootstrap' \n")
          }

          Vraw <- catch.conditions({
            get.conf.matrix (data = data[, T.pool, drop = FALSE],
                             cov.pool = data[, V.pool, drop = FALSE],
                             measure = "correlation",
                             conditional.vars = NULL,
                             blocksize = blocksize,
                             selection_fdr = fdr,
                             adjust_by = V.adjust_by,
                             apply.qval = V.apply.qval,
                             lambda = lambda,
                             pi0.method = "bootstrap",
                             alpha = alpha,
                             verbose = max(0, verbose - 1),
                             cl = cl, chunk.size = chunk.size)
          })$value

          ## Return an empty n_t list if 'bootstrap' method also failed
          if (any(class(Vraw) %in% c("simpleError", "error", "condition"))) {

            Vraw$sig.asso.covs <- vector(mode = "list", length = n_t)

          }

        }
        else {

          # Return an empty n_t list
          Vraw$sig.asso.covs <- vector(mode = "list", length = n_t)

        }
      }

      Vconfounders <- Vraw$sig.asso.covs
    }
  })

  ### Select T-nodes
  TTime <- system.time({
    if (skip.T) {
      ## Return empty list if no V-node
      Traw <- NULL
      Tconfounders <- vector(mode = "list", length = n_t)
    }
    else {
    if (verbose) {
      cat(paste0("\n        * selecting 'T-nodes' using ",
                 if (n_v > 0 & identical(T.measure, "partial")) 'partial' else 'marginal',
                 " correlations ... \n"))
    }
    ## get.conf.Tset
    if (n_v > 0 & identical(T.measure, "partial")) {

      ## Call get.conf.Tset
      Traw <- catch.conditions({
        get.conf.Tset (data = data, # data matrix
                       T.pool = T.pool,
                       V.pool = V.pool,
                       n_t = n_t,
                       n_v = n_v,
                       FDRcontrol = T.FDRcontrol,
                       adjust_by = T.adjust_by,
                       T.filter = T.filter,
                       Vconfounders = Vconfounders,
                       fdr = fdr,
                       lambda = lambda,
                       pi0.method = pi0.method,
                       alpha = alpha,
                       parallel = parallel,
                       cl = cl,
                       chunk.size = chunk.size,
                       verbose = max(0, verbose - 1))
      })$value

      ## Check the result: error or not?
      if (any(class(Traw) %in% c("simpleError", "error", "condition"))) {
        if (identical(Traw$message, "wrong sign in 'by' argument")) {
          if (verbose) {
            print(Traw$message)
            cat("\n 'get.conf.Tset' failled ...  \n")
          }
        }

        # Change pi0.method to 'bootstrap' if it was 'smoother' and 'get.conf.matrix' failed
        if (identical(pi0.method, "smoother")) {
          if (verbose) {
            print(Traw)
            cat("\n pi0.meth = 'smoother' failed, trying pi0.meth = 'bootstrap' \n")
          }

          Traw <- catch.conditions({
            get.conf.Tset (data = data, # data matrix
                           T.pool = T.pool,
                           V.pool = V.pool,
                           n_t = n_t,
                           n_v = n_v,
                           FDRcontrol = T.FDRcontrol,
                           adjust_by = T.adjust_by,
                           T.filter = T.filter,
                           Vconfounders = Vconfounders,
                           fdr = fdr,
                           lambda = lambda,
                           pi0.method = "bootstrap",
                           alpha = alpha,
                           parallel = parallel,
                           cl = cl,
                           chunk.size = chunk.size,
                           verbose = max(0, verbose - 1))
          })$value

          ## Return an empty n_t list if 'bootstrap' method also failed
          if (any(class(Traw) %in% c("simpleError", "error", "condition"))) {

            Traw$sig.asso.covs <- vector(mode = "list", length = n_t)

          }

        }
        else {
          # Return an empty n_t list
          Traw$sig.asso.covs <- vector(mode = "list", length = n_t)
        }
      }
    }
    else {
      if (verbose) {
        cat("\n        * selecting 'T-nodes' using marginal correlations ... \n")
      }

      ### Call 'get.conf.matrix' on 'T.pool'
      Traw <- catch.conditions({
        get.conf.matrix (data = data[, T.pool, drop = FALSE],
                         cov.pool = data[, T.pool, drop = FALSE],
                         measure = "correlation",
                         conditional.vars = NULL,
                         blocksize = blocksize,
                         selection_fdr = fdr,
                         adjust_by = T.adjust_by,
                         apply.qval = T.apply.qval,
                         lambda = lambda,
                         pi0.method = pi0.method,
                         alpha = alpha,
                         verbose = max(0, verbose - 1),
                         cl = cl, chunk.size = chunk.size)
      })$value

      ## Check the result: error or not?
      if (any(class(Traw) %in% c("simpleError", "error", "condition"))) {
        if (identical(Traw$message, "wrong sign in 'by' argument")) {
          if (verbose) {
            print(Traw$message)
            cat("\n 'get.conf.matrix' failled ...  \n")
          }
        }

        # Change pi0.method to 'bootstrap' if it was 'smoother' and 'get.conf.matrix' failed
        if (identical(pi0.method, "smoother")) {
          if (verbose) {
            print(Traw)
            cat("\n pi0.meth = 'smoother' failed, trying pi0.meth = 'bootstrap' \n")
          }

          Traw <- catch.conditions({
            get.conf.matrix (data = data[, T.pool, drop = FALSE],
                             cov.pool = data[, T.pool, drop = FALSE],
                             measure = "correlation",
                             conditional.vars = NULL,
                             blocksize = blocksize,
                             selection_fdr = fdr,
                             adjust_by = T.adjust_by,
                             apply.qval = T.apply.qval,
                             lambda = lambda,
                             pi0.method = "bootstrap",
                             alpha = alpha,
                             verbose = max(0, verbose - 1),
                             cl = cl, chunk.size = chunk.size) # $sig.asso.covs
          })$value

          ## Return an empty n_t list if 'bootstrap' method also failed
          if (any(class(Traw) %in% c("simpleError", "error", "condition"))) {

            Traw$sig.asso.covs <- vector(mode = "list", length = n_t)

          }

        }
        else {
          # Return an empty n_t list
          Traw$sig.asso.covs <- vector(mode = "list", length = n_t)
        }
      }
    }

    Tconfounders <- Traw$sig.asso.cov

    # Add n_t + n_v to the column indices returned by 'get.conf.matrix' to indicate 'U-nodes' in data
    Tconfounders <- lapply(1:n_t,
                           FUN = function(j) {
                             x <- Tconfounders[[j]]
                             if (is.null(x))
                               return(NULL)

                             # index of T-j selected for T_j
                             keep <- x != j

                             if (all(!keep))
                               return(NULL)

                             return(x[keep] + n_v)

                           })
    }
  })

  ### Select U, W and Z-nodes
  UWZTime <- system.time({
    if (n_c == 0 | skip.C) {
      ## Return NULL vector if no candidate W, Z-node or no V-node to distinguish U-nodes from Z,W-nodes
      UWZraw <- UWZconfounders <- NULL

      # Pool of candidate U,W,Z-nodes
      UWZindices <- NULL

      ### Number of U,W,Z-nodes
      n_uwz <- 0
    }
    else {
      ## Sample size
      n_samples <- NROW(data)

      ## Check if we need fallback approach
      use_C_fallback <- (n_v > 0) && (n_samples <= n_v + 3) && identical(C.measure, "partial")

      if (verbose) {
        cat(paste0("\n        * selecting 'U,W,Z-nodes' using ",
                   if (n_v > 0 & identical(C.measure, "partial")) 'partial' else 'marginal',
                   " correlations ... \n"))
      }

      if (n_v > 0 & identical(C.measure, "partial")) {

        if (use_C_fallback) {
          ## Use fallback approach
          if (verbose) {
            cat(paste0("            NOTE: Sample size (", n_samples,
                       ") <= V-nodes + 3 (", n_v + 3,
                       "). Using selected V-nodes instead of all V-nodes for T-C partial correlations.\n"))
          }

          ## Create grid of T-C pairs
          iterable <- expand.grid(T_idx = 1:n_t, C_idx = 1:n_c)

          ## Compute correlations with fallback
          results <- matteApply(X = iterable,
                                MARGIN = 1,
                                FUN = pcorTCwithFallback,
                                T.pool = T.pool,
                                C.pool = C.pool,
                                V.pool = V.pool,
                                Vconfounders = Vconfounders,
                                data = data,
                                n_samples = n_samples,
                                simplify = FALSE,
                                cl = cl, chunk.size = chunk.size)

          ## Extract correlations, methods, and conditioning set sizes
          r_vec <- sapply(results, function(x) x$r)
          methods_used <- sapply(results, function(x) x$method)
          n_cond_vec <- sapply(results, function(x) x$n_cond)
          n_used_vec <- sapply(results, function(x) x$n_used)

          ## Report method usage
          if (verbose) {
            method_counts <- table(methods_used)
            cat("            Method usage for T-C partial correlations:\n")
            for (m in names(method_counts)) {
              pct <- round(100 * method_counts[m] / length(methods_used), 1)
              cat(paste0("              - ", m, ": ", method_counts[m], " pairs (", pct, "%)\n"))
            }
          }

          ## Compute p-values based on method used
          p_vec <- mapply(function(r, nc, method, nu) {
            if (method == "marginal") {
              ## t-test for marginal correlation
              if (abs(r) >= 1) return(if (abs(r) == 1) 0 else NA)
              t_stat <- r * sqrt((nu - 2) / (1 - r^2))
              2 * (1 - pt(abs(t_stat), df = nu - 2))
            } else {
              ## z-test for partial correlation
              if (abs(r) >= 1) return(if (abs(r) == 1) 0 else NA)
              z <- (sqrt(nu - nc - 3) / 2) * log((1 + r) / (1 - r))
              2 * (1 - pnorm(abs(z)))
            }
          }, r_vec, n_cond_vec, methods_used, n_used_vec)

          ## Reshape to matrices: expand.grid gives T_idx varying fastest
          ## Result is n_t rows x n_c cols, need to transpose to n_c x n_t
          r.mat <- matrix(r_vec, nrow = n_t, ncol = n_c, byrow = FALSE)
          p.mat <- matrix(p_vec, nrow = n_t, ncol = n_c, byrow = FALSE)
          r.mat <- t(r.mat)  # n_c x n_t (covariates as rows, T-nodes as columns)
          p.mat <- t(p.mat)  # n_c x n_t

          ## Apply FDR control
          if (verbose) {
            if (C.apply.qval) {
              cat(paste0("            applying qvalue correction to control the FDR at ", fdr, "\n"))
            } else if (C.adjust_by != "none") {
              cat(paste0("            applying ", C.FDRcontrol, " correction with threshold ", alpha, "\n"))
            }
          }

          if (C.apply.qval) {
            ## qvalue correction with smoother -> bootstrap retry
            switch(C.adjust_by,
                   individual = {
                     qsig.list <- catch.conditions({
                       apply(p.mat, MARGIN = 2,
                             FUN = adjust.q,
                             fdr = fdr,
                             lambda = lambda,
                             pi0.meth = pi0.method)
                     })$value

                     ## Retry with bootstrap if smoother failed
                     if (any(class(qsig.list) %in% c("simpleError", "error", "condition")) &&
                         identical(pi0.method, "smoother")) {
                       if (verbose)
                         cat("            pi0.meth = 'smoother' failed for T-C qvalue, trying 'bootstrap'\n")
                       qsig.list <- catch.conditions({
                         apply(p.mat, MARGIN = 2,
                               FUN = adjust.q,
                               fdr = fdr,
                               lambda = lambda,
                               pi0.meth = "bootstrap")
                       })$value
                     }

                     ## Extract results or fall back to alpha threshold
                     if (any(class(qsig.list) %in% c("simpleError", "error", "condition"))) {
                       if (verbose)
                         cat("            qvalue correction failed, falling back to raw p-values with alpha threshold\n")
                       sig.mat <- p.mat <= alpha
                       q.mat <- matrix(NA, nrow = n_c, ncol = n_t)
                     } else {
                       sig.mat <- sapply(qsig.list, FUN = function(x) x$significant)
                       q.mat <- sapply(qsig.list, function(x) x$qvalue)
                     }
                   },
                   all = {
                     qsig.all <- catch.conditions({
                       adjust.q(as.vector(p.mat),
                                fdr = fdr,
                                lambda = lambda,
                                pi0.meth = pi0.method)
                     })$value

                     ## Retry with bootstrap if smoother failed
                     if (any(class(qsig.all) %in% c("simpleError", "error", "condition")) &&
                         identical(pi0.method, "smoother")) {
                       if (verbose)
                         cat("            pi0.meth = 'smoother' failed for T-C qvalue, trying 'bootstrap'\n")
                       qsig.all <- catch.conditions({
                         adjust.q(as.vector(p.mat),
                                  fdr = fdr,
                                  lambda = lambda,
                                  pi0.meth = "bootstrap")
                       })$value
                     }

                     ## Extract results or fall back to alpha threshold
                     if (any(class(qsig.all) %in% c("simpleError", "error", "condition"))) {
                       if (verbose)
                         cat("            qvalue correction failed, falling back to raw p-values with alpha threshold\n")
                       sig.mat <- p.mat <= alpha
                       q.mat <- matrix(NA, nrow = n_c, ncol = n_t)
                     } else {
                       sig.mat <- matrix(qsig.all$significant, nrow = n_c, ncol = n_t, byrow = FALSE)
                       q.mat <- matrix(qsig.all$qvalue, nrow = n_c, ncol = n_t, byrow = FALSE)
                     }
                   },
                   none = {
                     q.mat <- matrix(NA, nrow = n_c, ncol = n_t)
                     sig.mat <- p.mat <= alpha
                   })
            p.adj.mat <- matrix(NA, nrow = n_c, ncol = n_t)
          }
          else {
            ## Bonferroni or other methods
            switch(C.adjust_by,
                   individual = {
                     p.adj.mat <- apply(p.mat, MARGIN = 2,
                                        FUN = stats::p.adjust,
                                        method = C.FDRcontrol)
                     sig.mat <- p.adj.mat <= alpha
                     q.mat <- matrix(NA, nrow = n_c, ncol = n_t)
                   },
                   all = {
                     p.adj.mat <- matrix(stats::p.adjust(as.vector(p.mat), method = C.FDRcontrol),
                                         nrow = n_c, ncol = n_t, byrow = FALSE)
                     sig.mat <- p.adj.mat <= alpha
                     q.mat <- matrix(NA, nrow = n_c, ncol = n_t)
                   },
                   none = {
                     p.adj.mat <- q.mat <- matrix(NA, nrow = n_c, ncol = n_t)
                     sig.mat <- p.mat <= alpha
                   })
          }

          ## Extract significant covariates for each T-node
          sig.asso.covs <- lapply(1:n_t, FUN = function(j) {
            which(sig.mat[, j])
          })

          ## Create UWZraw structure for compatibility
          UWZraw <- list(sig.asso.covs = sig.asso.covs,
                         pvalues = p.mat,
                         qvalues = q.mat,
                         cors = r.mat,
                         sig = sig.mat,
                         adj.p = p.adj.mat)
        }
        else {
          ## Original approach: use get.conf.matrix with all V-nodes
          UWZraw <- catch.conditions({
            get.conf.matrix(data = data[, T.pool, drop = FALSE],
                            cov.pool = data[, C.pool, drop = FALSE],
                            measure = "partial_corr",
                            conditional.vars = data[, V.pool, drop = FALSE],
                            blocksize = blocksize,
                            selection_fdr = fdr,
                            adjust_by = C.adjust_by,
                            apply.qval = C.apply.qval,
                            lambda = lambda,
                            pi0.method = pi0.method,
                            alpha = alpha,
                            verbose = max(0, verbose - 1),
                            cl = cl, chunk.size = chunk.size)
          })$value

          ## Check the result: error or not?
          if (any(class(UWZraw) %in% c("simpleError", "error", "condition"))) {
            if (identical(UWZraw$message, "wrong sign in 'by' argument")) {
              if (verbose) {
                print(UWZraw$message)
                cat("\n 'get.conf.matrix' failled ...  \n")
              }
            }

            # Change pi0.method to 'bootstrap' if it was 'smoother' and 'get.conf.matrix' failed
            if (identical(pi0.method, "smoother")) {
              if (verbose) {
                print(UWZraw)
                cat("\n pi0.meth = 'smoother' failed, trying pi0.meth = 'bootstrap' \n")
              }

              UWZraw <- catch.conditions({
                get.conf.matrix(data = data[, T.pool, drop = FALSE],
                                cov.pool = data[, C.pool, drop = FALSE],
                                measure = "partial_corr",
                                conditional.vars = data[, V.pool, drop = FALSE],
                                blocksize = blocksize,
                                selection_fdr = fdr,
                                adjust_by = C.adjust_by,
                                apply.qval = C.apply.qval,
                                lambda = lambda,
                                pi0.method = "bootstrap",
                                alpha = alpha,
                                verbose = max(0, verbose - 1),
                                cl = cl, chunk.size = chunk.size)
              })$value

              ## Return an empty n_t list if 'bootstrap' method also failed
              if (any(class(UWZraw) %in% c("simpleError", "error", "condition"))) {

                UWZraw$sig.asso.covs <- vector(mode = "list", length = n_t)

              }

            }
            else {

              # Return an empty n_t list
              UWZraw$sig.asso.covs <- vector(mode = "list", length = n_t)

            }
          }
        }
      }
      else {
        ## Marginal correlation path
        if (verbose) {
          cat("\n        * selecting 'U,W,Z-nodes' using marginal correlations ... \n")
        }

        UWZraw <- catch.conditions({
          get.conf.matrix(data = data[, T.pool, drop = FALSE],
                          cov.pool = data[, C.pool, drop = FALSE],
                          measure = "correlation",
                          conditional.vars = NULL,
                          blocksize = blocksize,
                          selection_fdr = fdr,
                          adjust_by = C.adjust_by,
                          apply.qval = C.apply.qval,
                          lambda = lambda,
                          pi0.method = pi0.method,
                          alpha = alpha,
                          verbose = max(0, verbose - 1),
                          cl = cl, chunk.size = chunk.size)
        })$value

        ## Error handling
        if (any(class(UWZraw) %in% c("simpleError", "error", "condition"))) {
          if (identical(pi0.method, "smoother")) {
            if (verbose) {
              print(UWZraw)
              cat("\n pi0.meth = 'smoother' failed, trying pi0.meth = 'bootstrap' \n")
            }

            UWZraw <- catch.conditions({
              get.conf.matrix(data = data[, T.pool, drop = FALSE],
                              cov.pool = data[, C.pool, drop = FALSE],
                              measure = "correlation",
                              conditional.vars = NULL,
                              blocksize = blocksize,
                              selection_fdr = fdr,
                              adjust_by = C.adjust_by,
                              apply.qval = C.apply.qval,
                              lambda = lambda,
                              pi0.method = "bootstrap",
                              alpha = alpha,
                              verbose = max(0, verbose - 1),
                              cl = cl, chunk.size = chunk.size)
            })$value

            if (any(class(UWZraw) %in% c("simpleError", "error", "condition"))) {
              UWZraw$sig.asso.covs <- vector(mode = "list", length = n_t)
            }
          }
          else {
            UWZraw$sig.asso.covs <- vector(mode = "list", length = n_t)
          }
        }
      }

      # Extract selected confounding variables
      UWZconfounders <- UWZraw$sig.asso.covs

      # Add n_t + n_v to the column indices returned by 'get.conf.matrix' to indicate 'U-nodes' in data
      UWZconfounders <- lapply(UWZconfounders,
                             FUN = function(x) {
                               if (!is.null(x)) x + n_t + n_v
                             })

      # Extract the unique indices of all selected confounding variables: candidate U,W,Z-nodes
      UWZindices <- unique(unlist(UWZconfounders, recursive = TRUE))

      ### Number of U,W,Z-nodes (selected)
      n_uwz <- length(UWZindices)

      # Add n_t + n_v to the column indices returned by 'get.conf.matrix' to indicate 'Z,W-nodes' in data
      if (n_uwz > 0)
        UWZindices <- sort(UWZindices)
    }
  })

  ### Filtering: separate W,Z- nodes from U-nodes
  WZTime <- system.time({
    if (n_uwz == 0 | skip.C) {
      ## Return empty list if no candidate Z,W-node
      WZraw <- WZindices <- NULL
      WZconfounders <- vector(mode = "list", length = n_v)
      Uconfounders <- vector(mode = "list", length = n_t)

      # Number of W,Z-nodes, number of U-nodes
      n_wz <- n_u <- 0
    }
    else {
      ## Determine which V-nodes to use for W,Z identification
      ## When sample size <= n_v + 3, use union of selected V-nodes for consistency with fallback approach
      n_samples <- NROW(data)
      use_selected_V_for_WZ <- (n_v > 0) && (n_samples <= n_v + 3) && identical(C.measure, "partial")

      if (use_selected_V_for_WZ) {
        ## Use union of selected V-nodes
        V_for_WZ <- sort(unique(unlist(Vconfounders)))
        if (length(V_for_WZ) == 0) {
          ## No V-nodes were selected for any T-node, skip W,Z identification
          WZraw <- WZindices <- NULL
          WZconfounders <- vector(mode = "list", length = n_v)
          Uconfounders <- UWZconfounders  # All UWZ become U-nodes
          n_wz <- 0
          n_u <- n_uwz
        }
      } else {
        V_for_WZ <- V.pool
      }

      if (verbose) {
        cat("        * selecting 'W,Z-nodes' using marginal correlations ... \n")
        if (use_selected_V_for_WZ && length(V_for_WZ) > 0) {
          cat(paste0("            NOTE: Sample size (", n_samples,
                     ") <= V-nodes + 3 (", n_v + 3,
                     "). Using ", length(V_for_WZ),
                     " unique selected V-nodes instead of all ", n_v, " V-nodes.\n"))
        }
      }

      ## Only proceed if we have V-nodes for W,Z identification
      if (length(V_for_WZ) > 0) {
        ### Call 'get.conf.matrix' on 'UWZindices'
        WZraw <- catch.conditions({
          get.conf.matrix (data = data[, V_for_WZ, drop = FALSE],
                           cov.pool = data[, UWZindices, drop = FALSE],
                           measure = "correlation",
                           conditional.vars = NULL,
                           blocksize = blocksize,
                           selection_fdr = fdr,
                           adjust_by = Q.adjust_by,
                           apply.qval = Q.apply.qval,
                           lambda = lambda,
                           pi0.method = pi0.method,
                           alpha = alpha,
                           verbose = max(0, verbose - 1),
                           cl = cl, chunk.size = chunk.size)
        })$value

        ## Check the result: error or not?
        if (any(class(WZraw) %in% c("simpleError", "error", "condition"))) {
          if (identical(WZraw$message, "wrong sign in 'by' argument")) {
            if (verbose) {
              print(WZraw$message)
              cat("\n 'get.conf.matrix' failled ...  \n")
            }
          }

          # Change pi0.method to 'bootstrap' if it was 'smoother' and 'get.conf.matrix' failed
          if (identical(pi0.method, "smoother")) {
            if (verbose) {
              print(WZraw)
              cat("\n pi0.meth = 'smoother' failed, trying pi0.meth = 'bootstrap' \n")
            }

            WZraw <- catch.conditions({
              get.conf.matrix (data = data[, V_for_WZ, drop = FALSE],
                               cov.pool = data[, UWZindices, drop = FALSE],
                               measure = "correlation",
                               conditional.vars = NULL,
                               blocksize = blocksize,
                               selection_fdr = fdr,
                               adjust_by = Q.adjust_by,
                               apply.qval = Q.apply.qval,
                               lambda = lambda,
                               pi0.method = "bootstrap",
                               alpha = alpha,
                               verbose = max(0, verbose - 1),
                               cl = cl, chunk.size = chunk.size)
            })$value

            ## Return an empty list if 'bootstrap' method also failed
            if (any(class(WZraw) %in% c("simpleError", "error", "condition"))) {

              WZraw$sig.asso.covs <- vector(mode = "list", length = length(V_for_WZ))

            }

          }
          else {

            # Return an empty list
            WZraw$sig.asso.covs <- vector(mode = "list", length = length(V_for_WZ))

          }
        }
      }

      ## Only process WZ results if W,Z identification was performed
      if (length(V_for_WZ) > 0) {
        WZconfounders <- WZraw$sig.asso.covs

        # Adjust WZconfounders to indicate column numbers in data
        WZconfounders <- lapply(WZconfounders,
                                FUN = function(x) {
                                  if (!is.null(x)) UWZindices[x]
                                })

        # When using selected V-nodes, remap WZconfounders to full V.pool length (n_v)
        # so that names(WZconfounders) <- nm.variants works correctly downstream
        if (use_selected_V_for_WZ) {
          WZconf_full <- vector(mode = "list", length = n_v)
          for (k in seq_along(V_for_WZ)) {
            WZconf_full[[ V_for_WZ[k] ]] <- WZconfounders[[k]]
          }
          WZconfounders <- WZconf_full
        }

        # Extract the unique indices of all selected confounding variables: candidate U,W,Z-nodes
        WZindices <- unique(unlist(WZconfounders, recursive = TRUE))

        ### Number of W,Z-nodes (selected)
        n_wz <- length(WZindices)

        # Add n_t + n_v to the column indices returned by 'get.conf.matrix' to indicate 'Z,W-nodes' in data
        if (n_wz > 0)
          WZindices <- sort(WZindices)

        # Number of U-nodes
        n_u <- n_uwz - n_wz

        # Exclude W,Z-nodes from UWZconfounders to get Uconfounders
        Uconfounders <- lapply (UWZconfounders,
                                FUN = function(x) {
                                  Uconf <- setdiff(x, WZindices)
                                  return(if (length(Uconf)) Uconf )
                                })
      }

    }
  })

  ### Names/labels for V and T-nodes
  nm.variants <- colnames(data[, V.pool, drop = FALSE])
  if (is.null(nm.variants))
    nm.variants <- paste0("V", 1:length(V.pool))

  nm.genes <- colnames(data[, T.pool, drop = FALSE])
  if (is.null(nm.genes))
    nm.genes <- paste0("T", 1:length(T.pool))

  ### Name elements of each list
  names(Vconfounders) <- nm.genes
  names(Tconfounders) <- nm.genes
  names(Uconfounders) <- nm.genes
  if (!is.null(UWZconfounders))
    names(UWZconfounders) <- nm.genes
  names(WZconfounders) <- nm.variants

  ### Bind all selected nodes in one list
  out <- mapply (function(x,y,z) c(x, y, z), Vconfounders, Tconfounders, Uconfounders, SIMPLIFY = FALSE)

  ### Returning each type of nodes if required
  out <- list(confounders = out,
              Vconfounders = Vconfounders,
              Tconfounders = Tconfounders,
              Uconfounders = Uconfounders,
              UWZconfounders = UWZconfounders,
              WZconfounders = WZconfounders,
              UWZindices = UWZindices,
              WZindices = WZindices,
              raw = list(Vraw = Vraw,
                         Traw = Traw,
                         UWZraw = UWZraw,
                         WZraw = WZraw),
              time = list(VTime = VTime,
                          TTime = TTime,
                          UWZTime = UWZTime,
                          WZTime = WZTime),
              call = mcall)

  out <- structure(out, class = "conf.sets")

  ### Save output if required
  if(save.list == TRUE & !missing(save.path)) {
    if (verbose)
      cat(paste0("\n          saving output at ", save.path, ".RData"))
    save(out, file = paste0(save.path,".RData"))
  }

  return(out)

}

check.get.conf.sets.args <- function () {
  expression({
    ### Check arguments
    FDRcontrol <- FDRcontrol[1]
    stopifnot(FDRcontrol %in% c("none", "qvalue", p.adjust.methods))
    adjust_by <- adjust_by[1]
    stopifnot(adjust_by %in% c("all", "individual", "none"))

    C.adjust_by <- C.adjust_by[1]
    stopifnot(C.adjust_by %in% c("all", "individual", "none"))
    C.FDRcontrol <- C.FDRcontrol[1]
    stopifnot(C.FDRcontrol %in% c("none", "qvalue", p.adjust.methods))
    C.apply.qval <- (C.FDRcontrol == "qvalue") & (C.adjust_by != "none")

    T.adjust_by <- T.adjust_by[1]
    stopifnot(T.adjust_by %in% c("all", "individual", "none"))
    T.FDRcontrol <- T.FDRcontrol[1]
    stopifnot(T.FDRcontrol %in% c("none", "qvalue", p.adjust.methods))
    T.apply.qval <- (T.FDRcontrol == "qvalue") & (T.adjust_by != "none")

    V.adjust_by <- V.adjust_by[1]
    stopifnot(V.adjust_by %in% c("all", "individual", "none"))
    V.FDRcontrol <- V.FDRcontrol[1]
    stopifnot(V.FDRcontrol %in% c("none", "qvalue", p.adjust.methods))
    V.apply.qval <- (V.FDRcontrol == "qvalue") & (V.adjust_by != "none")

    Q.adjust_by <- Q.adjust_by[1]
    stopifnot(Q.adjust_by %in% c("all", "individual", "none"))
    Q.FDRcontrol <- Q.FDRcontrol[1]
    stopifnot(Q.FDRcontrol %in% c("none", "qvalue", p.adjust.methods))
    Q.apply.qval <- (Q.FDRcontrol == "qvalue") & (Q.adjust_by != "none")

    T.filter <- as.logical(T.filter)
    T.filter <- T.filter[1]
    stopifnot(is.logical(T.filter))

    ### If no data at all
    if (is.null(data)) {
      out <- list(confounders = list(),
                  Vconfounders = list(),
                  Tconfounders = list(),
                  Uconfounders = list(),
                  UWZconfounders = list(),
                  WZconfounders = list(),
                  UWZindices = list(),
                  WZindices = list(),
                  raw = list(Vraw = list(),
                             Traw = list(),
                             UWZraw = list(),
                             WZraw = list()),
                  time = list(VTime = system.time(""),
                              TTime = system.time(""),
                              UWZTime = system.time(""),
                              WZTime = system.time(""))
                  )

      ### Save output if required
      if(save.list == TRUE & !missing(save.path)) {
        if (verbose)
          cat(paste0("\n          saving output at ", save.path, ".RData"))
        save(out, file = paste0(save.path,".RData"))
      }

      return(out)
    }

    nvars <- NCOL(data)
    if (missing(n_c))
      n_c <- nvars - n_t - n_v
    else {
      stopifnot(n_t + n_v + n_c <= nvars)
    }
    m <- n_t + n_v + n_c
    stopifnot(m > 0, c(n_t, n_v, n_c) >= 0)
    stopifnot(all(c(n_t, n_v, n_c, m) <= nvars))

    # If no T-node to do selection for
    if (nvars == n_v) {
      out <- list(confounders = list(),
                  Vconfounders = list(),
                  Tconfounders = list(),
                  Uconfounders = list(),
                  UWZconfounders = list(),
                  WZconfounders = list(),
                  UWZindices = list(),
                  WZindices = list(),
                  raw = list(Vraw = list(),
                             Traw = list(),
                             UWZraw = list(),
                             WZraw = list()),
                  time = list(VTime = system.time(""),
                              TTime = system.time(""),
                              UWZTime = system.time(""),
                              WZTime = system.time(""))
                  )

      ### Save output if required
      if(save.list == TRUE & !missing(save.path)) {
        if (verbose)
          cat(paste0("\n          saving output at ", save.path, ".RData"))
        save(out, file = paste0(save.path,".RData"))
      }

      return(out)
    }

    # Only keep the useful data in 'data'
    data <- data[,1:m, drop = FALSE]

    if (scale.data) {
      if (n_v > 0)
        data[,-c(1:n_v)] <- scale(data[,-c(1:n_v)], center = TRUE, scale = TRUE)
      else
        data <- scale(data, center = TRUE, scale = TRUE)
    }

    ## Check correlation measure
    stopifnot(C.measure %in% c("marginal", "partial"))
    C.measure <- C.measure[1]
    stopifnot(T.measure %in% c("marginal", "partial"))
    T.measure <- T.measure[1]

    stopifnot(pi0.method %in% c("smoother", "boostrap"))
    pi0.method <- pi0.method[1]

    # Adjust blocksize to avoid bigcor errors (only when user explicitly provided it)
    if (!missing(blocksize))
      blocksize <- min(blocksize, n_t, n_v)

    # Parallelize computations or not?
    stopifnot(is.logical(parallel))
    if (parallel[1]) {
      if (is.null(cl)) {
        # Use option cl.cores to choose an appropriate cluster size
        cl <- parallel::makeCluster(getOption("cl.cores", 2L))
      }
      # Set up cluster with library paths, seed, and required packages
      setup_cluster(cl, packages = c("MRGN", "ppcor", "propagate"), seed = seed)
      if (!is.null(chunk.size)) {
        stopifnot(is.numeric(chunk.size))
        chunk.size <- chunk.size[1]
        stopifnot(chunk.size > 0)
        chunk.size <- ceiling(chunk.size)
      }
    }
    else {
      cl <- NULL
    }
  })
}
# p.adjust.methods bonferroni

# Generate expression to extract candidate counfounder pools
get.candidate.pools <- function () {
  expression({

    # V-nodes
    if (n_v > 0) {
      V.pool <- 1:n_v
    }
    else
      V.pool <- NULL

    # T-nodes
    T.pool <- (n_v+1):(n_v+n_t)

    # U-nodes
    if (n_c > 0) {
      C.pool <- (n_v+n_t+1):(n_v+n_t+n_c)
    }
    else
      C.pool <- NULL

  })
}


#' Compute T-T partial correlation with fallback logic
#'
#' When sample size is too small for partial correlation with all V-nodes,
#' this function falls back to using selected V-nodes (union for each pair),
#' or marginal correlation if selected V-nodes are still too many.
#'
#' @param jk vector of two T-node column indices in data
#' @param V.pool vector of V-node column indices
#' @param Vconfounders list of selected V-nodes for each T
#' @param data the data matrix
#' @param n_v number of V-nodes
#' @param n_samples sample size
#' @return list with r (correlation), method used, and n_cond (number of conditioning variables)
pcorTTwithFallback <- function(jk, V.pool, Vconfounders, data, n_v, n_samples) {

  # Get indices in Vconfounders list (T-nodes are indexed after V-nodes)
  i_idx <- jk[1] - n_v
  j_idx <- jk[2] - n_v

  # Step 1: Try all V-nodes if sample size is sufficient
  if (n_samples > length(V.pool) + 3) {
    mat <- na.omit(cbind(data[, c(jk, V.pool)]))
    if (nrow(mat) > length(V.pool) + 3) {
      r <- ppcor::pcor(mat)$estimate[1, 2]
      return(list(r = r, method = "all_V", n_cond = length(V.pool), n_used = nrow(mat)))
    }
  }

  # Step 2: Try selected V-nodes (union of V-nodes for both T-nodes)
  V_i <- Vconfounders[[i_idx]]
  V_j <- Vconfounders[[j_idx]]
  V_selected <- union(V_i, V_j)

  if (length(V_selected) > 0 && n_samples > length(V_selected) + 3) {
    mat <- na.omit(cbind(data[, c(jk, V_selected)]))
    if (nrow(mat) > length(V_selected) + 3) {
      r <- ppcor::pcor(mat)$estimate[1, 2]
      return(list(r = r, method = "selected_V", n_cond = length(V_selected), n_used = nrow(mat)))
    }
  }

  # Step 3: Fall back to marginal correlation
  r <- cor(data[, jk[1]], data[, jk[2]], use = "complete.obs")
  n_used <- sum(complete.cases(data[, c(jk[1], jk[2])]))
  return(list(r = r, method = "marginal", n_cond = 0, n_used = n_used))
}


#' Compute T-C partial correlation with fallback logic
#'
#' When sample size is too small for partial correlation with all V-nodes,
#' this function falls back to using selected V-nodes for the T-node,
#' or marginal correlation if selected V-nodes are still too many.
#'
#' @param idx vector of two indices: T index (1 to n_t) and C index (1 to n_c)
#' @param T.pool vector of T-node column indices in data
#' @param C.pool vector of C-node column indices in data
#' @param V.pool vector of V-node column indices in data
#' @param Vconfounders list of selected V-nodes for each T
#' @param data the data matrix
#' @param n_samples sample size
#' @return list with r (correlation), method used, and n_cond (number of conditioning variables)
pcorTCwithFallback <- function(idx, T.pool, C.pool, V.pool, Vconfounders, data, n_samples) {

  t_idx <- idx[1]  # Index in T.pool (1 to n_t)
  c_idx <- idx[2]  # Index in C.pool (1 to n_c)

  t_col <- T.pool[t_idx]
  c_col <- C.pool[c_idx]

  # Step 1: Try all V-nodes if sample size is sufficient
  if (n_samples > length(V.pool) + 3) {
    mat <- na.omit(cbind(data[, t_col], data[, c_col], data[, V.pool]))
    if (nrow(mat) > length(V.pool) + 3) {
      r <- ppcor::pcor(mat)$estimate[1, 2]
      return(list(r = r, method = "all_V", n_cond = length(V.pool), n_used = nrow(mat)))
    }
  }

  # Step 2: Try selected V-nodes for this T-node
  V_selected <- Vconfounders[[t_idx]]

  if (length(V_selected) > 0 && n_samples > length(V_selected) + 3) {
    mat <- na.omit(cbind(data[, t_col], data[, c_col], data[, V_selected]))
    if (nrow(mat) > length(V_selected) + 3) {
      r <- ppcor::pcor(mat)$estimate[1, 2]
      return(list(r = r, method = "selected_V", n_cond = length(V_selected), n_used = nrow(mat)))
    }
  }

  # Step 3: Fall back to marginal correlation
  r <- cor(data[, t_col], data[, c_col], use = "complete.obs")
  n_used <- sum(complete.cases(data[, c(t_col, c_col)]))
  return(list(r = r, method = "marginal", n_cond = 0, n_used = n_used))
}


# Select T-nodes with fallback for high-dimensional settings
get.conf.Tset <- function (data,
                           T.pool,
                           V.pool,
                           n_t = length(T.pool),
                           n_v = length(V.pool),
                           FDRcontrol = "qvalue",
                           adjust_by = 'individual',
                           T.filter = FALSE,
                           Vconfounders = NULL,
                           fdr = 0.05,
                           lambda = 0.05,
                           pi0.method = "smoother",
                           alpha = 0.01,
                           parallel = FALSE,
                           cl = NULL,
                           chunk.size = NULL,
                           verbose = 0,
                           save.list = FALSE,
                           save.path = "/path/to/location") {

  ## Center data columns
  data[, c(T.pool, V.pool)] <- scale(data[, c(T.pool, V.pool), drop = FALSE],
                                     center = TRUE, scale = FALSE)

  ## Sample size
  n_samples <- NROW(data)

  ## Number of V-nodes selected for each T-node
  N_Vi <- if (is.null(Vconfounders)) 0 else sapply(Vconfounders, length)

  ## Check if we need fallback approach (sample size too small for all V-nodes)
  use_fallback <- (n_v > 0) && (n_samples <= n_v + 3)

  if (use_fallback && verbose) {
    cat(paste0("            NOTE: Sample size (", n_samples,
               ") <= V-nodes + 3 (", n_v + 3,
               "). Using selected V-nodes instead of all V-nodes for T-T partial correlations.\n"))
  }

  if (T.filter & any(N_Vi)) {
    ## Asymmetric matrix case (T.filter logic)
    ## Create a grid of (x,y) pairs (i.e. pairs of indices of columns in data)
    xygrid <- expand.grid(x = T.pool, y = T.pool)
    nonzero <- xygrid[,1] != xygrid[,2]

    if (use_fallback) {
      ## Fallback approach: for each pair (T_j, T_k), condition on
      ## setdiff(V_k, V_j), i.e. T_k's instruments not shared with T_j.
      ## This mirrors the non-fallback setdiff(V.pool, V_j) at pair-specific
      ## scale, producing an asymmetric partial correlation matrix.
      ## Note: any(N_Vi) is guaranteed TRUE by the outer condition at entry.
      results <- matteApply(X = xygrid[nonzero, , drop = FALSE],
                            MARGIN = 1,
                            FUN = pcorTTgivenVjWithFallback,
                            Vconfounders = Vconfounders,
                            data = data,
                            n_v = n_v,
                            n_samples = n_samples,
                            simplify = FALSE,
                            chunk.size = chunk.size, cl = cl)

      ## Extract correlations and methods
      r.mat <- numeric(length = NROW(xygrid))
      r.mat[nonzero] <- sapply(results, function(x) x$r)
      methods_used <- sapply(results, function(x) x$method)
      n_cond_vec <- sapply(results, function(x) x$n_cond)
      n_used_vec <- sapply(results, function(x) x$n_used)

      ## Report method usage
      if (verbose) {
        method_counts <- table(methods_used)
        cat("            Method usage for T-T partial correlations (T.filter mode):\n")
        for (m in names(method_counts)) {
          pct <- round(100 * method_counts[m] / length(methods_used), 1)
          cat(paste0("              - ", m, ": ", method_counts[m], " pairs (", pct, "%)\n"))
        }
        if (all(methods_used == "marginal")) {
          cat("            NOTE: all pairs cascaded to marginal correlation.\n")
        }
      }

      ## Compute p-values based on method used for each pair
      pvalues <- numeric(length = NROW(xygrid)) + 1
      pvalues[nonzero] <- mapply(function(r, nc, method, nu) {
        if (method == "marginal") {
          ## t-test for marginal correlation
          if (abs(r) >= 1) return(if (abs(r) == 1) 0 else NA)
          t_stat <- r * sqrt((nu - 2) / (1 - r^2))
          2 * (1 - pt(abs(t_stat), df = nu - 2))
        } else {
          ## z-test for partial correlation
          if (abs(r) >= 1) return(if (abs(r) == 1) 0 else NA)
          z <- (sqrt(nu - nc - 3) / 2) * log((1 + r) / (1 - r))
          2 * (1 - pnorm(abs(z)))
        }
      }, r.mat[nonzero], n_cond_vec, methods_used, n_used_vec)

      ## Reorganize into matrices
      r.mat <- matrix(r.mat, nrow = n_t, ncol = n_t, byrow = FALSE)
      p.mat <- matrix(pvalues, nrow = n_t, ncol = n_t, byrow = FALSE)
    }
    else {
      ## Original approach: use all V-nodes minus V_j for each pair
      r.mat <- numeric(length = NROW(xygrid))
      r.mat[nonzero] <- matteApply(X = xygrid[nonzero, , drop = FALSE],
                                   MARGIN = 1,
                                   FUN = pcorTTgivenVj,
                                   V.pool = V.pool,
                                   Vset = Vconfounders,
                                   data = data, n_v = n_v,
                                   chunk.size = chunk.size, cl = cl)

      ## Perform Pearson tests and extract p-values
      pvalues <- numeric(length = NROW(xygrid)) + 1
      pvalues[nonzero] <- p.from.parcor(r.mat[nonzero], n = n_samples, S = N_Vi)$pvalue

      ## Reorganize into matrices
      r.mat <- matrix(r.mat, nrow = n_t, ncol = n_t, byrow = FALSE)
      p.mat <- matrix(pvalues, nrow = n_t, ncol = n_t, byrow = FALSE)
    }
  }
  else {
    ## Symmetric matrix case
    ## Create a grid of (x,y) pairs - only lower triangle
    xygrid <- t(combn(T.pool, m = 2))

    if (use_fallback) {
      ## Use fallback approach with selected V-nodes or marginal correlation
      results <- matteApply(X = xygrid,
                            MARGIN = 1,
                            FUN = pcorTTwithFallback,
                            V.pool = V.pool,
                            Vconfounders = Vconfounders,
                            data = data,
                            n_v = n_v,
                            n_samples = n_samples,
                            simplify = FALSE,
                            chunk.size = chunk.size, cl = cl)

      ## Extract correlations and methods
      pcorrs <- sapply(results, function(x) x$r)
      methods_used <- sapply(results, function(x) x$method)
      n_cond_vec <- sapply(results, function(x) x$n_cond)
      n_used_vec <- sapply(results, function(x) x$n_used)

      ## Report method usage
      if (verbose) {
        method_counts <- table(methods_used)
        cat("            Method usage for T-T partial correlations:\n")
        for (m in names(method_counts)) {
          pct <- round(100 * method_counts[m] / length(methods_used), 1)
          cat(paste0("              - ", m, ": ", method_counts[m], " pairs (", pct, "%)\n"))
        }
      }

      ## Compute p-values based on method used for each pair
      pvalues <- mapply(function(r, nc, method, nu) {
        if (method == "marginal") {
          ## t-test for marginal correlation
          if (abs(r) >= 1) return(if (abs(r) == 1) 0 else NA)
          t_stat <- r * sqrt((nu - 2) / (1 - r^2))
          2 * (1 - pt(abs(t_stat), df = nu - 2))
        } else {
          ## z-test for partial correlation
          if (abs(r) >= 1) return(if (abs(r) == 1) 0 else NA)
          z <- (sqrt(nu - nc - 3) / 2) * log((1 + r) / (1 - r))
          2 * (1 - pnorm(abs(z)))
        }
      }, pcorrs, n_cond_vec, methods_used, n_used_vec)

    }
    else {
      ## Original approach: partial correlations given all V-nodes
      pcorrs <- matteApply(X = xygrid,
                           MARGIN = 1,
                           FUN = pcorTTgivenV,
                           V.pool = V.pool,
                           data = data,
                           chunk.size = chunk.size, cl = cl)

      ## Compute p-values
      pvalues <- p.from.parcor(pcorrs, n = n_samples, S = n_v)$pvalue
    }

    ## Reorganize into matrices
    r.mat <- p.mat <- diag(n_t)
    r.mat[lower.tri(r.mat, diag = FALSE)] <- pcorrs
    p.mat[lower.tri(p.mat, diag = FALSE)] <- pvalues

    ## Make the matrices symmetric
    r.mat <- r.mat + t(r.mat)
    diag(r.mat) <- 0
    p.mat <- p.mat + t(p.mat)
    diag(p.mat) <- 1
  }

  ## Perform FDR control
  switch(FDRcontrol,
         none = {
           p.adj.mat <- qsig.mat <- q.mat <- matrix(NA, nrow = n_t, ncol = n_t)
           sig.mat <- p.mat <= alpha
         },
         qvalue = {
           if (verbose) {
             cat(paste0("            applying qvalue correction to control the FDR at ", fdr, "\n"))
           }
           switch(adjust_by,
                  individual = {
                    qsig.mat <- apply(p.mat, MARGIN = 2,
                                      FUN = adjust.q,
                                      fdr = fdr,
                                      lambda = lambda,
                                      pi0.meth = pi0.method)
                    sig.mat <- sapply(qsig.mat, FUN = function(x) x$significant)
                    q.mat <- sapply(qsig.mat, function(x) x$qvalue)
                  },
                  all = {
                    if (T.filter & any(N_Vi)) {
                      ## Asymmetric case: use all non-diagonal p-values
                      nondiag <- row(p.mat) != col(p.mat)
                      pvals_nondiag <- p.mat[nondiag]

                      qsig.mat <- adjust.q(pvals_nondiag,
                                           fdr = fdr,
                                           lambda = lambda,
                                           pi0.meth = pi0.method)

                      sig.mat <- matrix(FALSE, nrow = n_t, ncol = n_t)
                      sig.mat[nondiag] <- qsig.mat$significant

                      q.mat <- matrix(1, nrow = n_t, ncol = n_t)
                      q.mat[nondiag] <- qsig.mat$qvalue
                    }
                    else {
                      ## Symmetric case: use lower triangle only
                      pvals_lower <- p.mat[lower.tri(p.mat, diag = FALSE)]
                      qsig.mat <- adjust.q(pvals_lower,
                                           fdr = fdr,
                                           lambda = lambda,
                                           pi0.meth = pi0.method)

                      ## Restructure into matrices
                      sig.mat <- diag(n_t)
                      sig.mat[lower.tri(sig.mat, diag = FALSE)] <- qsig.mat$significant
                      sig.mat <- (sig.mat + t(sig.mat)) > 0
                      diag(sig.mat) <- FALSE

                      q.mat <- diag(n_t)
                      q.mat[lower.tri(q.mat, diag = FALSE)] <- qsig.mat$qvalue
                      q.mat <- q.mat + t(q.mat)
                      diag(q.mat) <- 1
                    }
                  },
                  none = {
                    qsig.mat <- q.mat <- matrix(NA, nrow = n_t, ncol = n_t)
                    sig.mat <- p.mat <= alpha
                  },
                  stop(paste0("Input to argument 'adjust_by' = ", adjust_by,
                              " is not recognized. Use one of 'individual', 'all', or 'none'"))
           )
           p.adj.mat <- matrix(NA, nrow = n_t, ncol = n_t)
         },
         {
           ## Bonferroni or other p.adjust methods
           if (verbose) {
             cat(paste0("            applying ", FDRcontrol, " correction with threshold ", alpha, "\n"))
           }
           switch(adjust_by,
                  individual = {
                    p.adj.mat <- apply(p.mat, MARGIN = 2,
                                       FUN = stats::p.adjust,
                                       method = FDRcontrol)
                    sig.mat <- apply(p.adj.mat, 2, function(x) x <= alpha)
                    q.mat <- matrix(NA, nrow = n_t, ncol = n_t)
                  },
                  all = {
                    if (T.filter & any(N_Vi)) {
                      ## Asymmetric case: use all non-diagonal p-values
                      nondiag <- row(p.mat) != col(p.mat)
                      pvals_nondiag <- p.mat[nondiag]

                      p.adj.nondiag <- stats::p.adjust(pvals_nondiag, method = FDRcontrol)

                      p.adj.mat <- matrix(1, nrow = n_t, ncol = n_t)
                      p.adj.mat[nondiag] <- p.adj.nondiag

                      sig.mat <- p.adj.mat <= alpha
                      q.mat <- matrix(NA, nrow = n_t, ncol = n_t)
                    }
                    else {
                      ## Symmetric case: use lower triangle only
                      pvals_lower <- p.mat[lower.tri(p.mat, diag = FALSE)]
                      p.adj.lower <- stats::p.adjust(pvals_lower, method = FDRcontrol)

                      p.adj.mat <- diag(n_t)
                      p.adj.mat[lower.tri(p.adj.mat, diag = FALSE)] <- p.adj.lower
                      p.adj.mat <- p.adj.mat + t(p.adj.mat)
                      diag(p.adj.mat) <- 1

                      sig.mat <- p.adj.mat <= alpha
                      q.mat <- matrix(NA, nrow = n_t, ncol = n_t)
                    }
                  },
                  none = {
                    p.adj.mat <- q.mat <- matrix(NA, nrow = n_t, ncol = n_t)
                    sig.mat <- p.mat <= alpha
                  },
                  stop(paste0("Input to argument 'adjust_by' = ", adjust_by,
                              " is not recognized. Use one of 'individual', 'all', or 'none'"))
           )
         })

  ## Extract significant covariates
  if (verbose) {
    cat("            selecting significant covariates... \n")
  }
  sig.asso.covs <- lapply(1:NCOL(sig.mat), FUN = function(j) {
    which(sig.mat[, j])
  })

  out.list <- list(sig.asso.covs = sig.asso.covs,
                   pvalues = p.mat,
                   qvalues = q.mat,
                   cors = r.mat,
                   sig = sig.mat,
                   adj.p = p.adj.mat)

  if (save.list == TRUE) {
    message(paste0("saving output at ", save.path, ".RData"))
    save(out.list, file = paste0(save.path, ".RData"))
  }

  return(out.list)
}


# Get T-T partial correlation given all V-nodes (calls ppcor::pcor)
pcorTTgivenV <- function (jk, V.pool, data) {
  tryCatch({
    mat <- na.omit(cbind(data[, c(jk, V.pool)]))
    ppcor::pcor(mat)$estimate[1, 2]
  }, error = function(e) {
    warning("pcorTTgivenV failed for pair (", jk[1], ", ", jk[2], "): ", e$message)
    NA
  })
}

# Get T-T partial correlation given selected V-nodes
pcorTTgivenVj <- function (jk, V.pool, Vset, data, n_v) {

  Vjk <- setdiff(V.pool, Vset[[jk[1] - n_v]]) # Remove V-nodes selected for Tk from the pool of V-nodes

  tryCatch({
    mat <- na.omit(cbind(data[, c(jk, Vjk)]))
    ppcor::pcor(mat)$estimate[1, 2]
  }, error = function(e) {
    warning("pcorTTgivenVj failed for pair (", jk[1], ", ", jk[2], "): ", e$message)
    NA
  })
}

#' T-T partial correlation with fallback for T.filter mode
#'
#' Computes partial correlation between two T-nodes using pair-specific
#' selected V-nodes while excluding T_j's instruments (T.filter philosophy).
#' For pair (T_j, T_k), the conditioning set is setdiff(V_k, V_j), i.e.
#' T_k's selected instruments that are not shared with T_j. Falls back to
#' marginal correlation when the conditioning set is empty or sample size
#' is insufficient.
#'
#' @param jk vector of two column indices (T_j, T_k) in data
#' @param Vconfounders list of selected V-node indices for each T-node
#' @param data data matrix
#' @param n_v number of V-nodes
#' @param n_samples sample size
#' @return list with r (correlation), method (string), n_cond (conditioning set size)
pcorTTgivenVjWithFallback <- function(jk, Vconfounders, data, n_v, n_samples) {

  ## T-node indices (1 to n_t)
  t_j_idx <- jk[1] - n_v  # row node
  t_k_idx <- jk[2] - n_v  # column node

  ## Get V-nodes selected for T_j and T_k
  V_j <- Vconfounders[[t_j_idx]]
  V_k <- Vconfounders[[t_k_idx]]
  if (is.null(V_j)) V_j <- integer(0)
  if (is.null(V_k)) V_k <- integer(0)

  ## Conditioning set: T_k's instruments not shared with T_j
  ## This is the pair-specific analog of setdiff(V.pool, V_j)
  V_cond <- setdiff(V_k, V_j)

  if (length(V_cond) > 0 && n_samples > length(V_cond) + 3) {
    mat <- na.omit(cbind(data[, jk], data[, V_cond]))
    if (nrow(mat) > length(V_cond) + 3) {
      r <- ppcor::pcor(mat)$estimate[1, 2]
      return(list(r = r, method = "selected_V_k_minus_V_j", n_cond = length(V_cond), n_used = nrow(mat)))
    }
  }
  ## Fall back to marginal correlation
  r <- cor(data[, jk[1]], data[, jk[2]], use = "complete.obs")
  n_used <- sum(complete.cases(data[, c(jk[1], jk[2])]))
  return(list(r = r, method = "marginal", n_cond = 0, n_used = n_used))
}
