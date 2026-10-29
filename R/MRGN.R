
# Keep the function simple for now
# Assume we know (as input) the set of confounding variables for each T node
# Know the type of each confounding variable (actual confounder, or one of intermediate or common child variable)
# Dealing with missing values ???

#' @name MRGN
#'
#' @title Infer a causal network using the MRGN algorithm
#'
#' @description Infer a causal network/graph with directed and undirected edges from observational data.
#'
#' @param data \code{data.frame}
#'
#' @param Qlabels labels indicating columns of \code{data} which are
#' \code{Q}-nodes (common children or intermediate variables).
#'
#' @param n_v,n_t,n_q,n_u integers, numbers of respectively \code{T}-nodes and \code{V}-nodes
#' in \code{data}.
#'
#' @param confounders a list of length \code{n_t} giving the set of known/selected
#' confounders for each phenotype, as returned by \link{get.conf.sets}.
#'
#' @param adjacency input adjacency matrix (binary)
#'
#' @param scale.data Logical, should the data be scaled (standardized) before processing?
#'   Default is TRUE. Only genes and confounders are scaled; variants are not scaled
#'
#' @param is.CNA Logical, indicating if a genetic variant is a copy number alteration (CNA).
#'   Default is FALSE
#'
#' @param alpha Numeric in (0, 0.5), type I error rate for individual tests. Default is 0.01
#'
#' @param use.perm Logical, should permutation tests be used? See \link[MRGN]{infer.trio}.
#'   Default is TRUE
#'
#' @param gamma Numeric, significance level for permutation tests when \code{use.perm = TRUE}.
#'   Default is 0.05
#'
#' @param nperms Integer, number of permutations to use when \code{use.perm = TRUE}.
#'   Default is 10000
#'
#' @param FDRcontrol Character, multiple-testing correction applied to the trio
#'   (and triplet) analyses at each iteration. One of "BH" (Benjamini-Hochberg),
#'   "bonferroni" (or any other method accepted by \link[stats]{p.adjust}),
#'   "qvalue", or "none". Default is "BH" (since 2026-10-22; previously "bonferroni")
#'
#' @param fdr Numeric in (0, 0.5), false discovery rate for all trio analyses at one iteration:
#'   the rejection threshold for q-values (\code{FDRcontrol = "qvalue"}) or for
#'   adjusted p-values (\code{FDRcontrol = "BH"}, \code{"bonferroni"}, etc.).
#'   Not used if \code{FDRcontrol = "none"}. Default is 0.05
#'
#' @param lambda Numeric, tuning parameter for pi0 estimation in q-value calculation.
#'   If NULL (default), initialized using \code{lambda.step}
#'
#' @param lambda.step Numeric in (0, 0.1), step size used to initialize \code{lambda}
#'   when \code{lambda = NULL}. Default is 0.05
#'
#' @param pi0.meth Character, method for pi0 estimation. One of "bootstrap" or "smoother".
#'   Default is "bootstrap"
#'
#' @param analyse.trios Logical, should Step 3 (VTT trios, followed by the
#'   \code{VTTthresh} majority vote; see \link{MRGN.VTT}) be run? Default is TRUE.
#'   With \code{FALSE}, Step 4 (if requested) runs directly on \code{adjacency}.
#'   Added 2026-10-22.
#'
#' @param analyse.triplets Logical, should triplets (trios with no genetic variants) be analyzed
#'   (Step 4; see \link{MRGN.TTT})? Default is TRUE
#'
#' @param solve.conflicts Logical, should conflicts in edge directions be resolved?
#'   Default is TRUE
#'
#' @param method Character, method for resolving conflicts when \code{solve.conflicts = TRUE}.
#'   Default is "conservative". Passed to \link{update.adjacency.matrix}
#'
#' @param maxiter Integer, maximum number of iterations for the iterative algorithm.
#'   If NULL (default), no limit is imposed
#'
#' @param parallel Logical, should computations be parallelized? Default is FALSE
#'
#' @param verbose Integer, verbosity level. 0 (default) for minimal output, higher values
#'   for more detailed progress messages
#'
#' @param seed Integer, random seed for reproducible results in parallel computing.
#'   If NULL (default), no seed is set
#'

#' @param stringent logical, should edges not appearing in \code{'M1'} and
#' \code{'M2'} trio structures be dropped during adjacency matrix update?
#' Passed to \link{update.adjacency.matrix}.
#'
#' @param VTTthresh Numeric in \eqn{[0, 1]}. Defaults to \code{0.3}: the T-T
#'   block of the final adjacency matrix is replaced by a majority-vote
#'   procedure over all VTT trios. A T-T edge is retained only if the fraction
#'   of trios supporting it (via models \code{M1.1}, \code{M1.2}, \code{M2.1},
#'   \code{M2.2}, \code{M4}, or \code{Other.1}) strictly exceeds
#'   \code{VTTthresh}, and its direction is determined by majority vote among
#'   directional models. Set to \code{0} to run the original MRGN algorithm
#'   unchanged instead. Values outside \eqn{[0, 1]} trigger a warning and
#'   revert to the default \code{0.3}.
#'
#' @param TTTthresh Numeric in \eqn{[0, 1]}. Default \code{0}: no Step 4 majority
#'   vote. When in \eqn{(0, 1]}, the direction of each \code{Tj-Tk} edge tested by
#'   "candidate trio (2)" triplets is re-decided by majority vote after Step 4; see
#'   \link{MRGN.TTT}. Never removes an edge. Re-added on 2026-10-29.
#'
#' @param cl,chunk.size optional arguments for parallel computing, passed to
#' \link[parallel]{parLapply} (when supplied).
#'
#' @param V.trioconfounders.filter logical, should each trio's confounder set be
#' dynamically filtered, right before that trio's analysis, to exclude any
#' candidate confounder correlated with that trio's OWN V-node (instrument)?
#' Unlike the (separate, \link{get.conf.sets}-level) Tcandidates filter, which
#' statically restricts confounder selection per T-node using that T-node's
#' overall set of selected V-confounders, this filter is recomputed fresh for
#' every trio, using the specific V anchoring that one trio. Defaults to
#' \code{TRUE} (since 2026-10-22; previously \code{FALSE}). Set to \code{FALSE}
#' to leave confounder sets completely unchanged. Added 2026-09-17.
#'
#' @param V.filter.cutoff numeric, the \code{|cor(Vi, confounder)|} cutoff
#' used when \code{V.trioconfounders.filter = TRUE}. Ignored otherwise. Defaults
#' to \code{0.2} (since 2026-10-22; previously \code{0.35}). Added 2026-09-17.
#'
#' @details
#' Since 2026-10-22, \code{MRGN} checks its arguments once and then calls the
#' workers behind \link{MRGN.VTT} (Step 3) and \link{MRGN.TTT} (Step 4) in turn.
#' Steps 3 and 4 keep separate iteration counters, each limited by \code{maxiter}.
#' The returned \code{niter} is \code{1 +} (extra Step 3 iterations) \code{+}
#' (extra Step 4 iterations). The output also contains \code{VTT.adjacency.prevote}
#' and \code{VTT.adjacency}, the skeleton after Step 3 before and after the
#' \code{VTTthresh} majority vote (\code{NULL} if \code{analyse.trios = FALSE}).
#'
#' @seealso \link{MRGN.VTT}, \link{MRGN.TTT}; \link[MRGN]{infer.trio} for inferring edges in small networks of one
#'  genetic variant and only two genes.

#' @export
#' @import parallel
#' @importFrom MRGN class.vec
#' @importFrom MRGN infer.trio
#' @importFrom MRGN adjust.q
#' @importFrom methods as is
#' @importFrom stats cor cov cov2cor p.adjust.methods pnorm sd
#' @importFrom utils combn installed.packages
#'
#'
#' @examples
#' ## Run MRGN on the 'networkA11' data
#' library (MRGNgeneral)
#' data ('networkA11')  # .rda structure with simulated data
#' data ('confsetsA11') #  obtained from get.conf.sets() on 'networkA11'
#'
#' # Build a graph skeleton
#' Adjacency0 <- get.initial.skeleton (data = networkA11$data,
#'                                     n_v = networkA11$dims$n_v,
#'                                     n_t = networkA11$dims$n_t,
#'                                     threshold_v = 0.2,
#'                                     threshold_m = 0.05,
#'                                     conf.sets = confsetsA11)
#'
#'
#' # Run MRGN: running time ~ 12 minutes on MAC OS Ventura 13.6.3 (16 GB memory, 8 cores).
#' \dontrun{
#' MRGNfit <- MRGN(data = networkA11$data,
#'                 n_v = networkA11$dims$n_v,
#'                 n_t = networkA11$dims$n_t,
#'                 Qlabels = confsetsA11$WZindices,
#'                 n_q = length(confsetsA11$WZindices),
#'                 n_u = networkA11$dims$n_w + networkA11$dims$n_z +
#'                   networkA11$dims$n_u + networkA11$dims$n_k +
#'                   networkA11$dims$n_i - length(confsetsA11$WZindices),
#'                 adjacency = Adjacency0,
#'                 confounders = confsetsA11$confounders,
#'                 alpha = 0.01,
#'                 FDRcontrol = 'bonferroni',
#'                 fdr = 0.05,
#'                 verbose = TRUE)
#' }

# We have four blocks of columns in data: V, T, (I&C) and confounders

# Post filtering based on marginal tests ???


# Keep track of all added edges (directions) so that
# at each iteration, I can check for conflicts with the current adds
# but also previous adds

# Use one 'control' argument, a list of (all?) control arguments (which ones?)

# To do: make dropped.edges and added.edges lists, so as to save the edges dropped/added at each iteration

MRGN <- function (data, # input n-by-m data matrix: 'n_v' Variants, 'n_t' Phenotypes, n_q Intermediate&Common_Child and 'n_u' confounders
                  scale.data = TRUE,
                  Qlabels = NULL, # column number or names of W,Z-nodes in data
                  n_t, n_v, n_q = length(Qlabels),
                  n_u = m - n_t - n_v - n_q,
                  confounders = list(), # a list of length n_t giving the set of confounders for each phenotype
                  adjacency, # a (n_t+n_v+n_q)-by-(n_t+n_v+n_q) adjacency matrix of an undirected acyclic graph (zeros on diagonal).
                  is.CNA = FALSE, # logical indicating if a genetic variant is a copy number alteration.
                  alpha = 0.01, # in (0, .5) type I error rate for individual tests, in the open (0, .5)
                  use.perm = FALSE, # See infer.trio
                  gamma = 0.05, #
                  nperms = 10000, #
                  ## CHANGED (2026-10-22): default is now "BH" (Benjamini-Hochberg via
                  ## stats::p.adjust, rejection at adjusted p <= fdr), previously "bonferroni".
                  FDRcontrol = c("BH", "bonferroni", "qvalue", "none"),
                  fdr = 0.05, # in (0, .5) False discovery rate for all trio analyses at one iteration, in the open (0, .5)
                  lambda = NULL, # The value of the tuning parameter for n_v-value estimation
                  lambda.step = 0.05, # in (0, .1); Used to initialize "lambda" when lambda = NULL
                  pi0.meth = c("bootstrap", "smoother"),
                  stringent = TRUE,
                  analyse.trios = TRUE,      # NEW (2026-10-22): run Step 3 (VTT)?
                  analyse.triplets = TRUE,
                  solve.conflicts = TRUE, # Just used for testing, to be removed?
                  method = "conservative",
                  maxiter = NULL,
                  parallel = FALSE,
                  cl = parallel::getDefaultCluster(),
                  chunk.size = NULL, # scalar number; number of invocations of fun or FUN in one chunk; a chunk is a unit for scheduling.
                  verbose = 0L,
                  VTTthresh = 0.3,
                  TTTthresh = 0,   # NEW (2026-10-29): optional Step 4 majority vote, default off
                  ## ---- NEW (2026-09-17): per-trio, dynamic confounder filter.
                  ## For each [Vi,Tj,Tk] (or [Vi,Tj,Qk]) trio, right before it goes
                  ## into trio analysis, excludes from that trio's confounder set
                  ## any candidate confounder correlated with THIS trio's own Vi
                  ## (|cor| >= V.filter.cutoff). Unlike the Tcandidates filter in
                  ## get.conf.sets() (static, per-T-node, based on that T-node's
                  ## OVERALL selected V-confounders), this is computed fresh per
                  ## trio, using the specific Vi anchoring that one trio -- V/T/UWZ
                  ## confounder selection in get.conf.sets() itself is completely
                  ## unaffected. CHANGED (2026-10-22): default TRUE (was FALSE).
                  ## See analyse.trio.i() in R/analyse.trio.R for the actual filter.
                  V.trioconfounders.filter = TRUE,
                  V.filter.cutoff = 0.2,   # CHANGED (2026-10-22): was 0.35
                  seed = NULL) { # seed for reproducible results in parallel computing
  # ============================================================================
  # Save the call to MRGN
  # ----------------------------------------------------------------------------
  mrgncl <- match.call()

  par.setup <- .mrgn.setup.parallel(parallel = parallel, cl = cl,
                                    chunk.size = chunk.size, seed = seed)
  cl <- par.setup$cl
  chunk.size <- par.setup$chunk.size

  # ============================================================================
  # Check arguments & Pre-process some inputs (done once for both steps)
  # ----------------------------------------------------------------------------
  eval(check.mrgn.args())
  VTTthresh <- .mrgn.check.VTTthresh(VTTthresh)
  TTTthresh <- .mrgn.check.TTTthresh(TTTthresh)

  # ============================================================================
  # Since 2026-10-22, Step 3 (VTT trios + VTTthresh majority vote) and Step 4
  # (TTT triplets) live in R/MRGN.VTT.R and R/MRGN.TTT.R; MRGN() calls their
  # internal workers in turn. Each step can also be run on its own with the
  # exported MRGN.VTT() / MRGN.TTT().
  # ----------------------------------------------------------------------------
  vtt <- ttt <- NULL
  added.edges <- dropped.edges <- NULL

  ### Step 3: VTT trio analysis (+ VTTthresh majority vote)
  if (analyse.trios) {
    vtt <- .MRGN.VTT.core(data = data, adjacency = adjacency,
                          n_v = n_v, n_t = n_t, n_q = n_q, Qlabels = Qlabels,
                          confounders = confounders,
                          alpha = alpha, FDRcontrol = FDRcontrol, fdr = fdr,
                          lambda = lambda, lambda.step = lambda.step, pi0.meth = pi0.meth,
                          use.perm = use.perm, gamma = gamma, is.CNA = is.CNA, nperms = nperms,
                          stringent = stringent, solve.conflicts = solve.conflicts,
                          method = method, maxiter = maxiter, VTTthresh = VTTthresh,
                          V.trioconfounders.filter = V.trioconfounders.filter,
                          V.filter.cutoff = V.filter.cutoff,
                          verbose = verbose, cl = cl, chunk.size = chunk.size)
    adjacency <- vtt$adjacency
    added.edges <- vtt$added.edges
    dropped.edges <- vtt$dropped.edges
  }

  ### Step 4: TTT triplet analysis
  if (analyse.triplets) {
    ttt <- .MRGN.TTT.core(data = data, adjacency = adjacency,
                          n_v = n_v, n_t = n_t, n_q = n_q,
                          confounders = confounders,
                          alpha = alpha, FDRcontrol = FDRcontrol, fdr = fdr,
                          lambda = lambda, lambda.step = lambda.step, pi0.meth = pi0.meth,
                          stringent = stringent, solve.conflicts = solve.conflicts,
                          method = method, maxiter = maxiter,
                          added.edges = added.edges, dropped.edges = dropped.edges,
                          TTTthresh = TTTthresh,
                          verbose = verbose, cl = cl, chunk.size = chunk.size)
    adjacency <- ttt$adjacency
    added.edges <- ttt$added.edges
    dropped.edges <- ttt$dropped.edges
  }

  # Indicator of the presence of undirected edges in the final network
  Not.fully.directed <- any(rowSums((adjacency + t(adjacency)) == 2) > 0)

  # Total number of passes, counted as before the split: 1 + extra Step 3
  # iterations + extra Step 4 iterations
  niter <- 1 + (if (!is.null(vtt)) vtt$niter - 1 else 0) +
    (if (!is.null(ttt)) ttt$niter - 1 else 0)

  # Report on the presence of undirected edges in the final network
  if (verbose) {
    if (Not.fully.directed)
      cat("\n      # the final network is not fully directed. \n")
    else
      cat("\n      # the final network is fully directed. \n")
  }

  structure(
    list(adjacency = structure(adjacency, class = 'adjacency.matrix'),
         VTTtrio.set = vtt$trio.set,
         iter.nb.VTTtrios = vtt$iter.nb.trios,
         VTTtrio.analysis = vtt$trio.analysis,
         VTQtrio.set = vtt$Qtrio.set,
         iter.nb.VTQtrios = vtt$iter.nb.Qtrios,
         VTQtrio.analysis = vtt$Qtrio.analysis,
         TTTtrio.set = ttt$triplet.set,
         iter.nb.TTTtrios = ttt$iter.nb.triplets,
         TTTtrio.analysis = ttt$triplet.analysis,
         TTQtrio.set = ttt$Qtriplet.set,
         iter.nb.TTQtrios = ttt$iter.nb.Qtriplets,
         TTQtrio.analysis = ttt$Qtriplet.analysis,
         added.edges = added.edges,
         dropped.edges = dropped.edges,
         input.adjacency = structure(Adj0, class = 'adjacency.matrix'),
         Qlabels = Qlabels, niter = niter,
         VTT.adjacency.prevote = vtt$adjacency.prevote,   # NEW (2026-10-22)
         TTT.adjacency.prevote = ttt$adjacency.prevote,   # NEW (2026-10-29): Step 4 result before the optional TTTthresh vote
         VTT.adjacency = vtt$adjacency,                   # NEW (2026-10-22): skeleton entering Step 4
         cl = cl, chunk.size = chunk.size,
         call = mrgncl),
    class = "MRGN")
}

# A brief print method for MRGN output (class 'MRGN')
# Issue with the count of edges: not consistent with graphNEL result
# when there are some undirected edges
#' @export
print.MRGN <- function (x, TTonly = FALSE, format = "adjacency",
                        digits = max(3, getOption("digits") - 3),
                        ...) {
  n_t <- eval(x$call$n_t)
  n_v <- eval(x$call$n_v)
  n_q <- if (TTonly) 0 else eval(x$call$n_q)
  adjacency <- as.matrix(x$adjacency[1:(n_t+n_v+n_q), 1:(n_t+n_v+n_q)])

  switch(format,
         graphNEL = {
           print(as(adjacency, "graphNEL"), digits = digits, ...)
         },
         igraph = { # igraph object
           print(igraph::graph_from_adjacency_matrix(adjacency), digits = digits, ...)
         },
         {
           cat("\nCall:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"),
               "\n\n", sep = "")

           #AugAdj <- adjacency + t(adjacency)
           #AugAdj[upper.tri(AugAdj)] <- 0

           #DirEdges <- which(AugAdj == 1, arr.ind = TRUE)
           #UndirEdges <- which(AugAdj == 2, arr.ind = TRUE)

           #NROW(DirEdges)
           #NROW(UndirEdges)

           EdgeDirCount <- adjacency[lower.tri(adjacency)] + adjacency[upper.tri(adjacency)]
           unidir <- EdgeDirCount == 1
           nb.uni <- sum(unidir)
           bidir <- EdgeDirCount == 2
           nb.bi <- sum(bidir)
           nb.edges <- nb.uni + nb.bi
           cat(paste0(" MRGN Inferred a Direct Acyclic Graph with: \n"))
           cat(paste0("    ", n_t+n_v+n_q, " nodes \n"))
           cat(paste0("        ", n_v, " variants (V-nodes) \n"))
           cat(paste0("        ", n_t, " genes (T-nodes) \n"))
           if (n_q > 0)
             cat(paste0("        ", n_q, " intermadiate variables or common children (Q-nodes) \n"))
           if (nb.uni > 0 & nb.bi > 0)
             cat(paste0("    ", nb.edges, " edges \n"))
           if (nb.uni > 0) {
             if (nb.bi > 0)
               cat(paste0("        ", nb.uni, " directed edges \n"))
             else
               cat(paste0("    ", nb.uni, " directed edges \n"))
           }
           if (nb.bi > 0) {
             if (nb.uni > 0)
               cat(paste0("        ", nb.bi, " undirected edges \n"))
             else
               cat(paste0("    ", nb.bi, " undirected edges \n"))
           }

         })
}
