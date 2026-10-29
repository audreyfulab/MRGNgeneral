#' Step 3 of MRGN: VTT trio analysis
#'
#' Run Step 3 of the MRGN algorithm on its own: form \code{[Vi,Tj,Tk]} (and
#' \code{[Vi,Tj,Qk]}) trios from the current skeleton, infer each trio's causal
#' structure, update the adjacency matrix, and iterate as long as new edges are
#' added; then (if \code{VTTthresh > 0}) replace the \code{T}-\code{T} block by
#' the majority vote over all VTT trios. This is exactly the VTT part of
#' \link{MRGN}, split out on 2026-10-22 so that it can be called separately,
#' e.g. to inspect the VTT result before running \link{MRGN.TTT} on it.
#'
#' @inheritParams MRGN
#'
#' @return an object of class \code{'MRGN.VTT'}, a list with elements
#' \describe{
#' \item{\code{adjacency}}{the adjacency matrix after Step 3 and the majority vote.}
#' \item{\code{adjacency.prevote}}{the adjacency matrix after Step 3, before the
#' \code{VTTthresh} majority vote (identical to \code{adjacency} when \code{VTTthresh = 0}).}
#' \item{\code{trio.set}, \code{trio.analysis}, \code{iter.nb.trios}}{all
#' \code{[Vi,Tj,Tk]} trios analysed, their analysis results, and the number of
#' trios formed at each iteration.}
#' \item{\code{Qtrio.set}, \code{Qtrio.analysis}, \code{iter.nb.Qtrios}}{same for
#' \code{[Vi,Tj,Qk]} trios (\code{NULL} if \code{n_q = 0}).}
#' \item{\code{added.edges}, \code{dropped.edges}}{edges added/dropped so far;
#' passed on to \link{MRGN.TTT} when its \code{adjacency} is an \code{'MRGN.VTT'} object.}
#' \item{\code{niter}}{number of Step 3 passes (1 + number of extra iterations).}
#' \item{\code{Not.fully.directed}}{logical, are there undirected edges left?}
#' \item{\code{input.adjacency}, \code{Qlabels}, \code{call}}{the (formatted)
#' input skeleton, \code{Q}-node labels, and the matched call.}
#' }
#'
#' @seealso \link{MRGN}, \link{MRGN.TTT}
#'
#' @export MRGN.VTT
MRGN.VTT <- function (data,
                      scale.data = TRUE,
                      Qlabels = NULL,
                      n_t, n_v, n_q = length(Qlabels),
                      n_u = m - n_t - n_v - n_q,
                      confounders = list(),
                      adjacency,
                      is.CNA = FALSE,
                      alpha = 0.01,
                      use.perm = FALSE,
                      gamma = 0.05,
                      nperms = 10000,
                      FDRcontrol = c("BH", "bonferroni", "qvalue", "none"),
                      fdr = 0.05,
                      lambda = NULL,
                      lambda.step = 0.05,
                      pi0.meth = c("bootstrap", "smoother"),
                      stringent = TRUE,
                      solve.conflicts = TRUE,
                      method = "conservative",
                      maxiter = NULL,
                      parallel = FALSE,
                      cl = parallel::getDefaultCluster(),
                      chunk.size = NULL,
                      verbose = 0L,
                      VTTthresh = 0.3,
                      V.trioconfounders.filter = TRUE,
                      V.filter.cutoff = 0.2,
                      seed = NULL) {
  mcall <- match.call()

  par.setup <- .mrgn.setup.parallel(parallel = parallel, cl = cl,
                                    chunk.size = chunk.size, seed = seed)
  cl <- par.setup$cl
  chunk.size <- par.setup$chunk.size

  eval(check.mrgn.args())
  VTTthresh <- .mrgn.check.VTTthresh(VTTthresh)

  out <- .MRGN.VTT.core(data = data, adjacency = adjacency,
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

  out$adjacency <- structure(out$adjacency, class = 'adjacency.matrix')
  out$adjacency.prevote <- structure(out$adjacency.prevote, class = 'adjacency.matrix')
  out$input.adjacency <- structure(Adj0, class = 'adjacency.matrix')
  out$Qlabels <- Qlabels
  out$call <- mcall
  structure(out, class = "MRGN.VTT")
}

# Internal worker for Step 3 (VTT trios + VTTthresh majority vote). Assumes the
# arguments have already been checked/formatted by check.mrgn.args() -- called
# by MRGN() (which checks once for both steps) and by MRGN.VTT().
# Body moved verbatim from MRGN() on 2026-10-22 (only the bookkeeping of
# Step 4 variables was removed and 'adjacency.prevote' added).
.MRGN.VTT.core <- function (data, adjacency, n_v, n_t, n_q, Qlabels, confounders,
                            alpha, FDRcontrol, fdr, lambda, lambda.step, pi0.meth,
                            use.perm, gamma, is.CNA, nperms,
                            stringent, solve.conflicts, method, maxiter, VTTthresh,
                            V.trioconfounders.filter, V.filter.cutoff,
                            verbose = 0L, cl = NULL, chunk.size = NULL) {
  ### Save the current Adjacency matrix
  old.adjacency <- adjacency

  # Initialize some outputs as empty vector/lists
  trio.analysis <- NULL
  dropped.edges <- added.edges <- NULL # list()
  Qtrio.set <- Qtrio.analysis <- NULL

  # Initialize iteration-wise counters/indicators
  new.edges <- NULL # indicator of trios with new added edges at an iteration (length = NROW(trio.set))
  new_edge_added <- FALSE # scalar indicator of at least one new edge added # new_edge_added <- any(new.edges)
  nb.new_edges <- 0 # total number of new edges added into the network at an iteration
  iter <- 1

  ### Step 3: List trios involving genetic variants and direct related edges
  # ----------------------------------------------------------------------------
  ## Step 3.1: Identify target trios, i.e. generate trio list & Bind trios identifiers as rows
  ## a [Vi,Tj,Tk] trio is formed exactly when Vi–Tj is an existing edge, and Tj–Tk is an existing edge
  ## (either because Vi–Tk also exists directly, or because Tk connects to Vi only indirectly through Tj)
  if (verbose) {
    cat("\n      # Generating an initial list of trios with genetic variants. \n")
  }
  trio.set <- enumerate.trios(i = 1:n_v, adjacency = adjacency,
                              n_t = n_t, n_v = n_v, n_q = n_q,
                              VTT = TRUE,
                              cl = cl,
                              chunk.size = chunk.size)
  trio.set <- do.call('rbind', trio.set)

  # Number of initial trios
  nb.trios <- NROW(trio.set)

  if (verbose) {
    cat(paste0("          * ",
               if (nb.trios == 0) "zero trio"
               else if (nb.trios == 1) "one [Vi,Tj,Tk] trio"
               else paste0(nb.trios, " [Vi,Tj,Tk] trios"),
               " generated. \n"))
  }

  ## Step 3.2-3.3: Determine the directed structure for each trio and update adjacency
  if (nb.trios > 0) {
    ## Step 3.2: Determine the directed structure for each trio
    if (verbose)
      cat("      # Inferring the causal network for each [Vi,Tj,Tk] trio. \n")
    trio.analysis <- analyse.trio.set (trio.set, nb.trios = nb.trios, Qlabels = Qlabels,
                                       alpha = alpha, FDRcontrol = FDRcontrol,
                                       fdr = fdr, lambda = lambda, lambda.step = lambda.step, pi0.meth = pi0.meth,
                                       data = data, confounders = confounders, n_t = n_t, n_v = n_v,
                                       use.perm = use.perm, gamma = gamma, is.CNA = is.CNA,
                                       nperms = nperms, verbose = verbose,
                                       V.trioconfounders.filter = V.trioconfounders.filter, V.filter.cutoff = V.filter.cutoff,
                                       cl = cl, chunk.size = chunk.size)

    ## Step 3.3: Update the adjacency matrix
    if (verbose) {
      cat(paste0("      # Updating the adjacency matrix: adding edge ",
                 if (nb.trios == 1) "direction for one [Vi,Tj,Tk] trio. \n"
                 else "directions for ", nb.trios, " [Vi,Tj,Tk] trios. \n"))
    }
    update.object <- update.adjacency.matrix (adjacency = adjacency, n_t = n_t, n_v = n_v, n_q = n_q,
                                              trio.set = trio.set,
                                              inferred.models = trio.analysis$Inferred.Model,
                                              stringent = stringent, add.edges = TRUE,
                                              solve.conflicts = solve.conflicts,
                                              method = method, cl = cl, chunk.size = chunk.size)
    # Update dropped.edges and added.edges
    added.edges <- update.object$added.edges
    dropped.edges <- update.object$dropped.edges

    # Report conflicts, if any
    if (verbose & solve.conflicts) {
      if (update.object$nb.conflicts)
        cat(paste0("          * ",
                   if (update.object$nb.conflicts == 1) "one conflict detected"
                   else paste0(update.object$nb.conflicts, " conflicts detected"),
                   " and handled with 'method = ", method, "'. \n"))
    }

    # Logical vector indicating trios in which new edges were added
    nb.new_edges <- update.object$nb.new_edges
    new.edges <- (update.object$new.edges == 1) & (nb.new_edges > 0)

    # Check if new edges were added from [Vi,Tj,Tk] trio analysis
    new_edge_added <- any(new.edges)
    if (verbose & new_edge_added) {
      cat(paste0("          * ",
                 if (nb.new_edges == 1) "one new edge"
                 else paste0(nb.new_edges, " new edges"),
                 " added into the network. \n"))
    }

    # Save the new adjacency
    adjacency <- update.object$adjacency
  }

  ## Step 3.4: Generating an initial list of trios involving a genetic variant and a Q-node
  if (n_q > 0) {
    if (verbose) {
      cat("\n      # Generating an initial list of [Vi,Tj,Qk] trios. \n")
    }
    Qtrio.set <- enumerate.trios(i = 1:n_v,
                                 adjacency = adjacency,
                                 n_t = n_t, n_v = n_v,
                                 n_q = n_q,
                                 VTT = FALSE,
                                 cl = cl, chunk.size = chunk.size)
    Qtrio.set <- do.call('rbind', Qtrio.set)

    # Number of initial trios
    nb.Qtrios <- NROW(Qtrio.set)

    if (verbose) {
      cat(paste0("          * ",
                 if (nb.Qtrios == 0) "zero trio"
                 else if (nb.Qtrios == 1) "one [Vi,Tj,Qk] trio"
                 else paste0(nb.Qtrios, " [Vi,Tj,Qk] trios"),
                 " generated. \n"))
    }

    if (nb.Qtrios > 0) {
      # Step 3.5: Determine the directed structure for each trio
      if (verbose)
        cat("      # Inferring the causal network for each trio. \n")
      Qtrio.analysis <- analyse.trio.set (Qtrio.set, nb.trios = nb.Qtrios, Qlabels = Qlabels,
                                          alpha = alpha, FDRcontrol = FDRcontrol,
                                          fdr = fdr, lambda = lambda, lambda.step = lambda.step, pi0.meth = pi0.meth,
                                          data = data, confounders = confounders, n_t = n_t, n_v = n_v,
                                          use.perm = use.perm, gamma = gamma, is.CNA = is.CNA,
                                          nperms = nperms, verbose = verbose,
                                          V.trioconfounders.filter = V.trioconfounders.filter, V.filter.cutoff = V.filter.cutoff,
                                          cl = cl, chunk.size = chunk.size)

      ## Step 3.6: Update the adjacency matrix
      if (verbose) {
        cat(paste0("      # Updating the adjacency matrix: adding edge ",
                   if (nb.Qtrios == 1) "direction for one [Vi,Tj,Qk] trio. \n"
                   else "directions for ", nb.Qtrios, " [Vi,Tj,Qk] trios. \n"))
      }
      Qupdate.object <- update.adjacency.matrix (adjacency = adjacency,
                                                 n_t = n_t, n_v = n_v, n_q = n_q,
                                                 trio.set = Qtrio.set,
                                                 inferred.models = Qtrio.analysis$Inferred.Model,
                                                 stringent = stringent, add.edges = FALSE,
                                                 solve.conflicts = solve.conflicts, method = method,
                                                 added.edges = added.edges, dropped.edges = dropped.edges,
                                                 cl = cl, chunk.size = chunk.size)
      dropped.edges <- Qupdate.object$dropped.edges

      # Update adjacency
      adjacency <- Qupdate.object$adjacency
    }
  }
  # ----------------------------------------------------------------------------

  ## Iterate Step 3 each time new edges are added (and the network is not fully directed)
  # Indicator of the presence of undirected edges in the updated network
  Not.fully.directed <- any(rowSums((adjacency + t(adjacency)) == 2) > 0)

  # Number of [Vi,Tj,Tk] trios generated at each iteration
  iter.nb.trios <- nb.trios

  # Number of [Vi,Tj,Qk] trios generated at each iteration
  iter.nb.Qtrios <- if (n_q > 0) nb.Qtrios

  # Initialize the number of new trios to nb.trios
  nb.new.trios <- nb.trios

  while ((new_edge_added & Not.fully.directed) & iter < maxiter) {
    ## Step 3.1: Form additional trios
    if (verbose) {
      cat("\n      # Generating additional [Vi,Tj,Tk] trios. \n")
    }
    new.trio.set <- enumerate.new.trios (trio.set = trio.set,
                                         new.edges = c(rep(FALSE, nb.trios - nb.new.trios),
                                                       new.edges),
                                         old.adjacency = old.adjacency, adjacency = adjacency,
                                         n_t = n_t, n_v = n_v, n_q = n_q,
                                         VTT = TRUE,
                                         cl = cl, chunk.size = chunk.size)

    # Stop if no new trio generated
    nb.new.trios <- NROW(new.trio.set)
    if (nb.new.trios == 0) {
      new_edge_added <- FALSE
      if (verbose) {
        cat("          * no additional [Vi,Tj,Tk] trio found.")
      }

      break
    }

    if (verbose) {
      cat(paste0("          * ",
                 if (nb.new.trios == 1) "one new [Vi,Tj,Tk] trio"
                 else paste0(nb.new.trios, " new [Vi,Tj,Tk] trios"),
                 " generated. \n"))
    }

    # Save the current Adjacency matrix (at iter)
    old.adjacency <- adjacency

    # iteration counter
    iter <- iter + 1

    ## Step 3.2: Determine the directed structure for each trio
    if (verbose) {
      cat(paste0("      # Inferring the causal network for each new [Vi,Tj,Tk] trio (iteration ",
                 iter, "). \n"))
    }
    new.trio.analysis <- analyse.trio.set (new.trio.set, nb.trios = nb.new.trios, Qlabels = Qlabels,
                                           alpha = alpha, FDRcontrol = FDRcontrol,
                                           fdr = fdr, lambda = lambda,
                                           lambda.step = lambda.step, pi0.meth = pi0.meth,
                                           data = data, confounders = confounders, n_t = n_t, n_v = n_v,
                                           use.perm = use.perm, gamma = gamma, is.CNA = is.CNA,
                                           nperms = nperms, verbose = verbose,
                                           V.trioconfounders.filter = V.trioconfounders.filter, V.filter.cutoff = V.filter.cutoff,
                                           cl = cl, chunk.size = chunk.size)

    ## Step 3.3: Update the adjacency matrix
    if (verbose) {
      cat(paste0("      # Updating the adjacency matrix: adding edge ",
                 if (nb.new.trios == 1) "direction for one [Vi,Tj,Tk] trio. \n"
                 else "directions for ", nb.new.trios, " [Vi,Tj,Tk] trios. \n"))
    }
    update.object <- update.adjacency.matrix (adjacency = adjacency,
                                              n_t = n_t, n_v = n_v, n_q = n_q,
                                              trio.set = new.trio.set,
                                              inferred.models = new.trio.analysis$Inferred.Model,
                                              stringent = stringent, add.edges = TRUE,
                                              solve.conflicts = solve.conflicts, method = method,
                                              added.edges = added.edges,
                                              dropped.edges = dropped.edges,
                                              cl = cl, chunk.size = chunk.size)

    # Update dropped.edges and added.edges
    added.edges <- update.object$added.edges
    dropped.edges <- update.object$dropped.edges

    # Report conflicts, if any
    if (verbose & solve.conflicts) {
      if (update.object$nb.conflicts)
        cat(paste0("          * ",
                   if (update.object$nb.conflicts == 1) "one conflict detected"
                   else paste0(update.object$nb.conflicts, " conflicts detected"),
                   " and handled with 'method = ", method, "'. \n"))
    }

    # Logical vector indicating trios in which new edges were added
    nb.new_edges <- update.object$nb.new_edges
    new.edges <- (update.object$new.edges == 1) & (nb.new_edges > 0)

    # Check if new edges were added from [V,T1,T2] trio analysis
    new_edge_added <- any(new.edges)
    if (verbose & new_edge_added) {
      cat(paste0("          * ",
                 if (nb.new_edges == 1) "one new edge"
                 else paste0(nb.new_edges, " new edges"),
                 " added into the network. \n"))
    }

    # Save the new adjacency
    adjacency <- update.object$adjacency

    ## Step 3.4: Generating a list of trios involving a genetic variant and a Q-node
    if (n_q > 0) {
      if (verbose) {
        cat("\n      # Generating additional [Vi,Tj,Qk] trios. \n")
      }

      new.Qtrio.set <- enumerate.new.trios (trio.set = new.trio.set,
                                            new.edges = new.edges,
                                            old.adjacency = old.adjacency, adjacency = adjacency,
                                            n_t = n_t, n_v = n_v, n_q = n_q,
                                            VTT = FALSE,
                                            Qtrio.set = Qtrio.set,
                                            cl = cl, chunk.size = chunk.size)

      # Number of initial trios
      nb.new.Qtrios <- NROW(new.Qtrio.set)

      if (verbose) {
        cat(paste0("          * ",
                   if (nb.new.Qtrios == 0) "zero new trio"
                   else if (nb.new.Qtrios == 1) "one new [Vi,Tj,Qk] trio"
                   else paste0(nb.new.Qtrios, " new [Vi,Tj,Qk] trios"),
                   " generated. \n"))
      }

      if (nb.new.Qtrios > 0) {
        # Step 3.5: Determine the directed structure for each trio
        if (verbose)
          cat("      # Inferring the causal network for each trio. \n")
        new.Qtrio.analysis <- analyse.trio.set (new.Qtrio.set, nb.trios = nb.new.Qtrios, Qlabels = Qlabels,
                                                alpha = alpha, FDRcontrol = FDRcontrol,
                                                fdr = fdr, lambda = lambda, lambda.step = lambda.step, pi0.meth = pi0.meth,
                                                data = data, confounders = confounders, n_t = n_t, n_v = n_v,
                                                use.perm = use.perm, gamma = gamma, is.CNA = is.CNA,
                                                nperms = nperms, verbose = verbose,
                                                V.trioconfounders.filter = V.trioconfounders.filter, V.filter.cutoff = V.filter.cutoff,
                                                cl = cl, chunk.size = chunk.size)

        ## Step 3.6: Update the adjacency matrix
        if (verbose) {
          cat(paste0("      # Updating the adjacency matrix: adding edge ",
                     if (nb.new.Qtrios == 1) "direction for one [Vi,Tj,Qk] trio. \n"
                     else "directions for ", nb.new.Qtrios, " [Vi,Tj,Qk] trios. \n"))
        }
        Qupdate.object <- update.adjacency.matrix (adjacency = adjacency,
                                                   n_t = n_t, n_v = n_v, n_q = n_q,
                                                   trio.set = new.Qtrio.set,
                                                   inferred.models = new.Qtrio.analysis$Inferred.Model,
                                                   stringent = stringent, add.edges = FALSE,
                                                   solve.conflicts = solve.conflicts,
                                                   method = method, dropped.edges = dropped.edges,
                                                   cl = cl, chunk.size = chunk.size)
        dropped.edges <- Qupdate.object$dropped.edges

        # Update adjacency
        adjacency <- Qupdate.object$adjacency

        # Save the additional trios in 'Qtrio.set'
        Qtrio.set <- rbind(Qtrio.set, new.Qtrio.set)

        # Save the additional trios analysis result in 'Qtrio.analysis'
        Qtrio.analysis <- rbind(Qtrio.analysis, new.Qtrio.analysis) # (not useful at this point)
      }

      # Save the number of additional trios generated
      iter.nb.Qtrios <- c(iter.nb.Qtrios, nb.new.Qtrios)
    }

    # Save the number of additional trios generated
    iter.nb.trios <- c(iter.nb.trios, nb.new.trios)

    # Save the additional trios in 'trio.set'
    trio.set <- rbind(trio.set, new.trio.set)

    # Save the additional trios analysis result in 'trio.analysis'
    trio.analysis <- rbind(trio.analysis, new.trio.analysis) # (not useful at this point)

    # Update the total number of trios analysed so far
    nb.trios <- nb.trios + nb.new.trios

    # Indicator of the presence of undirected edges in the updated network
    Not.fully.directed <- any(rowSums((adjacency + t(adjacency)) == 2) > 0)
  }

  # Keep the Step 3 skeleton before the vote (returned for comparison)
  adjacency.prevote <- adjacency

  # Apply VTTthresh majority vote on T-T block before Step 4 so that Step 4
  # propagates from the majority-vote skeleton rather than the raw iterative skeleton.
  if (VTTthresh > 0 && !is.null(trio.set) && !is.null(trio.analysis)) {
    Tidx <- (n_v + 1L):(n_v + n_t)
    adjacency <- .majority_adj(adjacency = adjacency,
                               ts        = trio.set,
                               ana       = trio.analysis,
                               Tidx      = Tidx,
                               thresh    = VTTthresh)
    Not.fully.directed <- any(rowSums((adjacency + t(adjacency)) == 2) > 0)
  }

  list(adjacency = adjacency,
       adjacency.prevote = adjacency.prevote,
       trio.set = trio.set,
       iter.nb.trios = iter.nb.trios,
       trio.analysis = trio.analysis,
       Qtrio.set = Qtrio.set,
       iter.nb.Qtrios = iter.nb.Qtrios,
       Qtrio.analysis = Qtrio.analysis,
       added.edges = added.edges,
       dropped.edges = dropped.edges,
       niter = iter,
       Not.fully.directed = Not.fully.directed)
}

# VTTthresh majority vote over VTT trios (moved from R/MRGN.R on 2026-10-22)
.majority_adj <- function(adjacency, ts, ana, Tidx,        # Tidx = which rows/cols of adjacency are T-nodes
                          keep = c("M1.1","M1.2","M2.1","M2.2","M4","Other.1"),
                          thresh = 0.5, clear.tt = TRUE) {
  # md is the classified model for every trio, one entry per trio (e.g.
  # "M1.1"). E.g. idx <- which(key == "103.165") finds all 13 trios that
  # tested the T3-T65 pair, and md[idx] gives their 13 classifications
  # (2 "M1.1", 11 "Other.1"). These per-trio labels are what the existence
  # check (md %in% keep, below) and the direction vote (further below) both
  # work from.
  md <- ana$Inferred.Model
  ok <- !is.na(md); ts <- ts[ok, , drop = FALSE]; md <- md[ok]  # drop trios with no classification (NA) before using md/ts below

  A <- adjacency                          # a copy
  if (clear.tt) A[Tidx, Tidx] <- 0        # VTT: wipe the whole T-T block first; only pairs the vote approves below get restored. clear.tt = FALSE skips this and leaves existing edges alone (used by the optional TTTthresh vote in MRGN.TTT(), which only re-directs edges)

  lo <- pmin(ts[,2], ts[,3]); hi <- pmax(ts[,2], ts[,3])  # canonical (lo,hi) so the same pair gets the same key regardless of column order
  key <- paste0(lo, ".", hi)              # one key per Tj-Tk pair; pools all trios that tested it, whatever the third node (Vi or Ti) was

  # For each Tj-Tk pair, compute what fraction of its trios have a model in
  # 'keep' (md %in% keep is TRUE/FALSE per trio; averaging TRUE/FALSE within
  # a pair gives the fraction that were TRUE). Keep the pair only if that
  # fraction is > thresh -- not all trios need to agree, just enough of them.
  ex <- tapply(md %in% keep, key, mean)
  keep.pairs <- names(ex)[ex > thresh]

  for (k in keep.pairs) {           # for each pair that passed the existence check above
    idx <- which(key == k)          # all trios that tested this specific pair
    a <- lo[idx[1]]; b <- hi[idx[1]] # the two node ids for this pair (same for every trio in idx)
    dir <- md[idx] %in% c("M1.1","M1.2","M2.1","M2.2")  # which of this pair's trios have an opinion on direction (M4/Other.1 don't)
    if (any(dir)) {
      di <- idx[dir]; m <- md[di]
      # translate each directional trio's model into an actual edge direction:
      # M1.1/M2.2 mean col2 -> col3; M1.2/M2.1 mean col3 -> col2 (see update_adj_mat_trio.R)
      from <- ifelse(m %in% c("M1.1","M2.2"), ts[di,2], ts[di,3])
      to   <- ifelse(m %in% c("M1.1","M2.2"), ts[di,3], ts[di,2])
      ab <- sum(from == a & to == b); ba <- sum(from == b & to == a)  # count votes for a->b vs b->a
    } else ab <- ba <- 0   # no directional trios at all -> treated as a tie below (undirected)
    # majority direction wins, and the losing direction is cleared explicitly;
    # a tie (including ab == ba == 0) leaves the edge undirected
    if      (ab > ba) { A[a,b] <- 1; A[b,a] <- 0 }
    else if (ba > ab) { A[b,a] <- 1; A[a,b] <- 0 }
    else            { A[a,b] <- 1; A[b,a] <- 1 }
  }
  A
}
