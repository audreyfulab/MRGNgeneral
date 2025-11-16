
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
#' @param FDRcontrol,V.FDRcontrol,T.FDRcontrol,C.FDRcontrol characters indicating
#' the FDR control methods to be used for different selections.
#' One of \code{"none"}, \code{"qvalue"} (see \link[qvalue]{qvalue}), or
#' \code{"bonferroni"} (see \link[stats]{p.adjust}).
#' If any of \code{T.FDRcontrol}, \code{C.FDRcontrol}, \code{V.FDRcontrol} is missing,
#' \code{FDRcontrol} is used, otherwise \code{FDRcontrol} is ignored.
#'
#' @param adjust_by,V.adjust_by,T.adjust_by,C.adjust_by,Q.FDRcontrol character indicating the
#' adjustment scheme for tests. One of \code{"none"} (no adjustment is desired),
#' \code{"individual"} (adjust p-values for each gene separately),
#' and \code{"all"} (adjust all p-values for all genes at once).
#' If any of \code{T.adjust_by}, \code{C.adjust_by}, \code{V.adjust_by} is missing,
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
#'   \code{V}-nodes in the network, and *z*-tests are used instead of *t*-tests.}
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
#' (except \code{time}), see \link{reorder.conf.sets}.}
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
      if (verbose) {
        cat(paste0("\n        * selecting 'U,W,Z-nodes' using ", C.measure[1], " correlations ... \n"))
      }

      ### Call 'get.conf.matrix' on 'C.pool'
      UWZraw <- catch.conditions({
        get.conf.matrix (data = data[, T.pool, drop = FALSE],
                         cov.pool = data[, C.pool, drop = FALSE],
                         measure = if (n_v > 0 & identical(C.measure, "partial")) "partial_corr" else "correlation",
                         conditional.vars = if (n_v > 0 & identical(C.measure, "partial")) data[, V.pool, drop = FALSE],
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
            get.conf.matrix (data = data[, T.pool, drop = FALSE],
                             cov.pool = data[, C.pool, drop = FALSE],
                             measure = if (n_v > 0 & identical(C.measure, "partial")) "partial_corr" else "correlation",
                             conditional.vars = if (n_v > 0 & identical(C.measure, "partial")) data[, V.pool, drop = FALSE],
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
      if (verbose) {
        cat("        * selecting 'W,Z-nodes' using marginal correlations ... \n")
      }

      ### Call 'get.conf.matrix' on 'UWZindices'
      WZraw <- catch.conditions({
        get.conf.matrix (data = data[, V.pool, drop = FALSE],
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
            get.conf.matrix (data = data[, V.pool, drop = FALSE],
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

          ## Return an empty n_t list if 'bootstrap' method also failed
          if (any(class(WZraw) %in% c("simpleError", "error", "condition"))) {

            WZraw$sig.asso.covs <- vector(mode = "list", length = n_t)

          }

        }
        else {

          # Return an empty n_t list
          WZraw$sig.asso.covs <- vector(mode = "list", length = n_t)

        }
      }

      WZconfounders <- WZraw$sig.asso.covs

      # Adjust WZconfounders to indicate column numbers in data
      WZconfounders <- lapply(WZconfounders,
                              FUN = function(x) {
                                if (!is.null(x)) UWZindices[x]
                              })

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

# Select T-nodes
get.conf.Tset <- function (data, # data matrix
                           T.pool,
                           V.pool,
                           n_t = length(T.pool),    # Number of Genes
                           n_v = length(V.pool),    # Number of Variants
                           FDRcontrol = "qvalue",
                           adjust_by = 'individual',
                           T.filter = FALSE,
                           Vconfounders = NULL,
                           fdr = 0.05,
                           lambda = 0.05,
                           pi0.method = "smoother",
                           alpha = 0.01, # Used if FDRcontrol = 'none'
                           parallel = FALSE,   # Use this in the CODE
                           cl = NULL,
                           chunk.size = NULL, # scalar number; number of invocations of fun or FUN in one chunk; a chunk is a unit for scheduling.
                           verbose = 0,
                           save.list = FALSE,
                           save.path="/path/to/location") {

  ## Center data columns to be used (only useful when using inverse updating to compute artial correlations)
  data[, c(T.pool, V.pool)] <- scale(data[, c(T.pool, V.pool), drop = FALSE],
                                     center = TRUE, scale = FALSE)

  ## Number of V-nodes selected for each T-node
  N_Vi <- if (is.null(Vconfounders)) 0 else sapply(Vconfounders, length)

  if (T.filter & any(N_Vi)) { ## the matrix of correlations may be asymmetric in this case
    ## Create a grid of (x,y) pairs (i.e. pairs of indices of columns in data)
    xygrid <- expand.grid(x = T.pool, y = T.pool)

    ## Get the partial correlations
    r.mat <- numeric(length = NROW(xygrid)) # initialize to zeros
    nonzero <- xygrid[,1] != xygrid[,2]      # Keep zeros for corr(Tj, Tj)
    r.mat[nonzero] <- matteApply (X = xygrid[nonzero, , drop = FALSE],
                          MARGIN = 1,
                          FUN = pcorTTgivenVj,
                          V.pool = V.pool,
                          Vset = Vconfounders,
                          data = data, n_v = n_v,
                          chunk.size = chunk.size, cl = cl)

    ## Perform Pearson tests and extract n_t values
    pvalues <- numeric(length = NROW(xygrid)) + 1
    pvalues[nonzero] <- p.from.parcor(r.mat[nonzero], n = NROW(data), S = N_Vi)$pvalue

    ## Reorganize into matrices
    r.mat <- matrix(r.mat, nrow = n_t, ncol = n_t, byrow = FALSE)
    p.mat <- matrix(pvalues, nrow = n_t, ncol = n_t, byrow = FALSE)
  }
  else { ## Here we have a symmetric matrix of correlations
    ## Create a grid of (x,y) pairs (i.e. pairs of indices of columns in data)
    ## Here, we only consider indices (i, j) such that j < i
    xygrid <- t(combn(T.pool, m = 2))

    ## Get the partial correlations (lower triangle elements)
    pcorrs <- matteApply (X = xygrid,
                          MARGIN = 1,
                          FUN = pcorTTgivenV,
                          V.pool = V.pool,
                          data = data,
                          chunk.size = chunk.size, cl = cl)

    ## Perform Pearson tests and extract n_t values
    pvalues <- p.from.parcor(pcorrs, n = NROW(data), S = n_v)$pvalue

    ## Reorganize into matrices
    r.mat <- p.mat <- diag(n_t)
    r.mat[lower.tri(r.mat, diag = FALSE)] <- pcorrs # (lower triangle)
    p.mat[lower.tri(p.mat, diag = FALSE)] <- pvalues # (lower triangle)

    ## Make the matrices symmetric
    r.mat <- r.mat + t(r.mat) # (upper triangle)
    diag(r.mat) <- 0  # Setting cor(T_j, T_j) to zero to avoid selecting T_j as a confounder for T_j
    p.mat <- p.mat + t(p.mat)
    diag(p.mat) <- 1  # Setting pvalue = 1 for cor(T_j, T_j) to avoid selecting T_j as a confounder for T_j
  }

  ## Perform q value correction, if required (code adapted from Jarred's 'get.conf.matrix' function)
  switch(FDRcontrol,
         none =  {
           p.adj.mat <- qsig.mat <- q.mat <- matrix(NA, nrow = n_t, ncol = n_t)
           sig.mat <- p.mat <= alpha
         },
         qvalue = {
           #qvalue correction
           if (verbose) {
             cat(paste0("            applying qvalue correction to control the FDR at ", fdr, "\n"))
           }
           #adjustment by columns or by all pvalues
           switch(adjust_by,
                  individual = {
                    qsig.mat <- apply(p.mat, MARGIN = 2,
                                     FUN = adjust.q,
                                     fdr = fdr,
                                     lambda = lambda,
                                     pi0.meth = pi0.method)

                    #significance matrix (binary matrix)
                    sig.mat <- sapply(qsig.mat, FUN = function(x) x$significant)

                    #qvalue matrix
                    q.mat <- sapply(qsig.mat, function(x) x$qvalue)
                  },
                  all = {
                    if (T.filter & any(N_Vi)) {
                      qsig.mat <- adjust.q (pvalues[nonzero],
                                            fdr = fdr,
                                            lambda = lambda,
                                            pi0.meth = pi0.method)

                      sig.mat <- nonzero
                      sig.mat[nonzero] <- qsig.mat$significant
                      q.mat <- numeric(NROW(xygrid)) + 1
                      q.mat[nonzero] <- qsig.mat$qvalue

                      qsig.mat$significant
                      qsig.mat$qvalue[nonzero] <- qsig.mat0$qvalue; rm(qsig.mat0)

                      #restructure results into significance matrix of snps X covariates
                      sig.mat <- matrix(sig.mat, nrow = n_t, ncol = n_t, byrow = FALSE)
                      q.mat <- matrix(q.mat, nrow = n_t, ncol = n_t, byrow = FALSE)

                    }
                    else {
                      qsig.mat <- adjust.q (pvalues,
                                            fdr = fdr,
                                            lambda = lambda,
                                            pi0.meth = pi0.method)

                      #restructure results into significance matrix of snps X covariates
                      sig.mat <- diag(n_t)
                      sig.mat[lower.tri(sig.mat, diag = FALSE)] <- qsig.mat$significant # (lower triangle)
                      sig.mat <- (sig.mat + t(sig.mat)) > 0
                      diag(sig.mat) <- FALSE
                      #restructure results into qvalue matrix of snps X covariates
                      q.mat <- diag(n_t)
                      q.mat[lower.tri(q.mat, diag = FALSE)] <- qsig.mat$qvalue # (lower triangle)
                      q.mat <- q.mat + t(q.mat)
                      diag(q.mat) <- 1
                    }
                  },
                  none = {
                    qsig.mat <- q.mat <- matrix(NA, nrow = n_t, ncol = n_t)
                    #significance matrix (binary matrix)
                    sig.mat <- p.mat <=  alpha
                  },
                  stop(paste0('Input to argument \'adjust_by \' = ', adjust_by, ' is not recognized. use one of \'individual\', \' all \', or \'none\' '))
           )
           #empty matrix
           p.adj.mat <- matrix(NA, nrow = n_t, ncol = n_t)
         },
         {
           #bonferroni correction
           if (verbose) {
             cat(paste0("            applying Bonferroni correction with threshold ", alpha, "\n"))
           }
           #adjustment by columns or by all pvalues
           switch(adjust_by,
                  individual = {
                    #adjust by columns
                    #matrix of adjusted pvalues
                    p.adj.mat <- apply(p.mat,
                                      MARGIN = 2,
                                      FUN = stats::p.adjust,
                                      method = FDRcontrol)

                    #significance matrix (binary matrix)
                    sig.mat <- apply(p.adj.mat, 2, function(x) x <= alpha)
                    #empty matrix
                    q.mat <- matrix(NA, nrow = n_t, ncol = n_t)

                  },
                  all = {
                    if (T.filter & any(N_Vi)) {
                      # Adjust pvalues
                      p.adj.mat <- numeric(NROW(xygrid)) + 1
                      p.adj.mat[nonzero] <- stats::p.adjust(pvalues[nonzero],
                                                            method = FDRcontrol)

                      #restructure results into significance matrix of snps X covariates
                      p.adj.mat <- matrix(p.adj.mat, nrow = n_t, ncol = n_t, byrow = FALSE)

                      #significance matrix (binary matrix)
                      sig.mat <- p.adj.mat <= alpha

                      #empty matrix
                      q.mat <- matrix(NA, nrow = nrow(p.mat), ncol = ncol(p.mat))
                    }
                    else {
                      #adjust by all pvalues
                      #matrix of adjusted pvalues
                      p.adj.mat <- stats::p.adjust(pvalues,
                                                   method = FDRcontrol)

                      #significance matrix (binary matrix)
                      sig.mat <- diag(n_t)
                      sig.mat[lower.tri(sig.mat, diag = FALSE)] <- p.adj.mat # (lower triangle)
                      sig.mat <- sig.mat + t(sig.mat)
                      diag(sig.mat) <- 1
                      p.adj.mat <- sig.mat
                      sig.mat <- p.adj.mat <= alpha

                      #empty matrix
                      q.mat <- matrix(NA, nrow = nrow(p.mat), ncol = ncol(p.mat))
                    }
                  },
                  none = {
                    p.adj.mat <- q.mat <- matrix(NA, nrow = n_t, ncol = n_t)
                    sig.mat <- p.mat <= alpha
                  },
                  stop(paste0('Input to argument \'adjust_by \' = ', adjust_by, ' is not recognized. use one of \'individual\', \' all \', or \'none\' '))
           )
         })


  #-------------------obtain-final-list-of-significant-covs----------------
  #find the covs that correlated with every column of the matrix
  if (verbose) {
    cat("            selecting significant covariates... \n")
  }
  #extract significant covariates:
  sig.asso.covs <- lapply(1:NCOL(sig.mat), FUN = function(j) {
    x <- sig.mat[,j]
    which(x)
  })


  out.list <- list(sig.asso.covs = sig.asso.covs,
                   pvalues = p.mat,
                   qvalues = q.mat,
                   cors = r.mat,
                   sig = sig.mat,
                   adj.p = p.adj.mat)


  if(save.list == TRUE) {
    message(paste0("saving output at ", save.list, ".RData"))
    save(out.list, file = paste0(save.path,".RData"))
  }

  return(out.list)

}

# Get T-T partial correlation given all V-nodes (calls ppcor::pcor)
pcorTTgivenV <- function (jk, V.pool, data) {
  ppcor::pcor(cbind(data[, c(jk, V.pool) ]))$estimate[1,2]
}

# Get T-T partial correlation given selected V-nodes
pcorTTgivenVj <- function (jk, V.pool, Vset, data, n_v) {

  Vjk <- setdiff(V.pool, Vset[[jk[1] - n_v]]) # Remove V-nodes selected for Tk from the pool of V-nodes

  ppcor::pcor(cbind(data[, c(jk, Vjk) ]))$estimate[1,2]
}
