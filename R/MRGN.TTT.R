#' Step 4 of MRGN: TTT triplet analysis
#'
#' Run Step 4 of the MRGN algorithm on its own: form \code{[Ti,Tj,Tk]} (and
#' \code{[Ti,Tj,Qk]}) triplets with no genetic variant from the current
#' skeleton, infer their causal structure, update the adjacency matrix, and
#' iterate on triplets classified as \code{"Other"} while new directions can be
#' propagated. This is exactly the TTT part of \link{MRGN}, split out on
#' 2026-10-22 so that it can be called separately, e.g. on the result of
#' \link{MRGN.VTT}, or directly on an initial skeleton.
#'
#' @inheritParams MRGN
#'
#' @param adjacency either an adjacency matrix (as for \link{MRGN}), or an
#' object of class \code{'MRGN.VTT'} (or \code{'MRGN'}), whose \code{adjacency}
#' is then used, together with its \code{added.edges} and \code{dropped.edges}
#' (unless these are supplied explicitly).
#'
#' @param added.edges,dropped.edges edges already added/dropped by earlier
#' steps, used by \link{update.adjacency.matrix} to check conflicts. Default
#' \code{NULL} (none), or taken from \code{adjacency} when it is an
#' \code{'MRGN.VTT'} or \code{'MRGN'} object.
#'
#' @param TTTthresh Numeric in \eqn{[0, 1]}. Default \code{0}: no majority vote,
#'   Step 4 runs as its sequential pass alone. When in \eqn{(0, 1]}, after Step 4
#'   the direction of every \code{Tj-Tk} edge tested by "candidate trio (2)"
#'   triplets (\code{Ti -> Tj - Tk}) is re-decided by majority vote over all such
#'   triplets: the edge is directed only if the fraction of triplets voting
#'   \code{M1.1} or \code{M2.1} (as opposed to \code{"Other"}) strictly exceeds
#'   \code{TTTthresh}, by the majority among those votes (\code{M1.1}:
#'   \code{Tj -> Tk}; \code{M2.1}: \code{Tk -> Tj}; a tie leaves it undirected).
#'   Unlike \code{VTTthresh}, this never removes an edge. Values outside
#'   \eqn{[0, 1]} trigger a warning and revert to \code{0}. Re-added on
#'   2026-10-29 (first added 2026-07-26, removed 2026-08-19).
#'
#' @return an object of class \code{'MRGN.TTT'}, a list with elements
#' \code{adjacency}, \code{adjacency.prevote} (before the optional
#' \code{TTTthresh} vote; identical to \code{adjacency} when \code{TTTthresh = 0}),
#' \code{triplet.set}, \code{iter.nb.triplets},
#' \code{triplet.analysis}, \code{Qtriplet.set}, \code{iter.nb.Qtriplets},
#' \code{Qtriplet.analysis}, \code{added.edges}, \code{dropped.edges},
#' \code{niter} (number of Step 4 passes, counted separately from Step 3),
#' \code{Not.fully.directed}, \code{input.adjacency}, \code{Qlabels} and \code{call}.
#'
#' @seealso \link{MRGN}, \link{MRGN.VTT}
#'
#' @export MRGN.TTT
MRGN.TTT <- function (data,
                      scale.data = TRUE,
                      Qlabels = NULL,
                      n_t, n_v, n_q = length(Qlabels),
                      n_u = m - n_t - n_v - n_q,
                      confounders = list(),
                      adjacency,
                      alpha = 0.01,
                      FDRcontrol = c("BH", "bonferroni", "qvalue", "none"),
                      fdr = 0.05,
                      lambda = NULL,
                      lambda.step = 0.05,
                      pi0.meth = c("bootstrap", "smoother"),
                      stringent = TRUE,
                      solve.conflicts = TRUE,
                      method = "conservative",
                      maxiter = NULL,
                      added.edges = NULL,
                      dropped.edges = NULL,
                      TTTthresh = 0,
                      parallel = FALSE,
                      cl = parallel::getDefaultCluster(),
                      chunk.size = NULL,
                      verbose = 0L,
                      seed = NULL) {
  mcall <- match.call()

  ## Accept the output of MRGN.VTT() (or MRGN()) as 'adjacency'
  if (inherits(adjacency, c("MRGN.VTT", "MRGN"))) {
    prev.step <- adjacency
    if (is.null(added.edges))
      added.edges <- prev.step$added.edges
    if (is.null(dropped.edges))
      dropped.edges <- prev.step$dropped.edges
    adjacency <- prev.step$adjacency
    rm(prev.step)
  }

  par.setup <- .mrgn.setup.parallel(parallel = parallel, cl = cl,
                                    chunk.size = chunk.size, seed = seed)
  cl <- par.setup$cl
  chunk.size <- par.setup$chunk.size

  eval(check.mrgn.args())
  TTTthresh <- .mrgn.check.TTTthresh(TTTthresh)

  out <- .MRGN.TTT.core(data = data, adjacency = adjacency,
                        n_v = n_v, n_t = n_t, n_q = n_q,
                        confounders = confounders,
                        alpha = alpha, FDRcontrol = FDRcontrol, fdr = fdr,
                        lambda = lambda, lambda.step = lambda.step, pi0.meth = pi0.meth,
                        stringent = stringent, solve.conflicts = solve.conflicts,
                        method = method, maxiter = maxiter,
                        added.edges = added.edges, dropped.edges = dropped.edges,
                        TTTthresh = TTTthresh,
                        verbose = verbose, cl = cl, chunk.size = chunk.size)

  out$adjacency <- structure(out$adjacency, class = 'adjacency.matrix')
  out$adjacency.prevote <- structure(out$adjacency.prevote, class = 'adjacency.matrix')
  out$input.adjacency <- structure(Adj0, class = 'adjacency.matrix')
  out$Qlabels <- Qlabels
  out$call <- mcall
  structure(out, class = "MRGN.TTT")
}

# Internal worker for Step 4 (TTT triplets). Assumes the arguments have already
# been checked/formatted by check.mrgn.args() -- called by MRGN() and MRGN.TTT().
# Body moved verbatim from MRGN() on 2026-10-22, with three changes:
#  - Step 4 has its own iteration counter 'iter' (before the split it continued
#    Step 3's counter, so Step 3 iterations used up Step 4's 'maxiter' budget);
#  - the 'analyse.triplets' switch is handled by the caller;
#  - iter.nb.triplets now records nb.new.triplets (was nb.new.trios, a Step 3
#    variable -- bookkeeping only, no effect on the inferred network).
.MRGN.TTT.core <- function (data, adjacency, n_v, n_t, n_q, confounders,
                            alpha, FDRcontrol, fdr, lambda, lambda.step, pi0.meth,
                            stringent, solve.conflicts, method, maxiter,
                            added.edges = NULL, dropped.edges = NULL,
                            TTTthresh = 0,
                            verbose = 0L, cl = NULL, chunk.size = NULL) {
  # Initialize outputs
  triplet.set <- triplet.analysis <- NULL
  Qtriplet.set <- Qtriplet.analysis <- NULL
  iter.nb.triplets <- iter.nb.Qtriplets <- NULL

  # Step 4 iteration counter (separate from Step 3 since 2026-10-22)
  iter <- 1

  # Indicator of the presence of undirected edges in the input network
  Not.fully.directed <- any(rowSums((adjacency + t(adjacency)) == 2) > 0)

  ### Step 4: List trios involving no genetic variant and direct related edges
  # ----------------------------------------------------------------------------

  if (all(c(Not.fully.directed, n_t + n_q >= 3))) { # n_t >= 3 means at least one triplet can be formed
    ## Step 4.1: Form [Ti,Tj,Tk]  trios
    if (verbose) {
      cat("\n      # Generating an initial list of triplets involving only T-nodes. \n")
    }
    triplet.set <- enumerate.triplets (adjacency,
                                       n_t = n_t,
                                       n_v = n_v,
                                       n_q = n_q,
                                       TTT = TRUE,
                                       cl = cl,
                                       chunk.size = chunk.size)

    # Number of initial triplets
    nb.triplets <- length(triplet.set$types)
    if (!nb.triplets) {
      if (verbose) {
        cat("          * no updateable [Ti,Tj,Tk] trio found. \n")
      }
      triplet.set <- triplet.analysis <- NULL
    }
    else {
      triplet.set <- cbind(triplet.set$triplets, type = triplet.set$types)

      if (verbose) {
        cat(paste0("          * ",
                   if (nb.triplets == 1) "one updateable [Ti,Tj,Tk] trio"
                   else paste0(nb.triplets, " updateable  [Ti,Tj,Tk] trios"),
                   " generated. \n"))
      }

      # Step 4.2: Determine the directed structure for each triplet
      if (verbose)
        cat("      # Inferring the causal network for each [Ti,Tj,Tk] trio. \n")
      triplet.analysis <- analyse.triplet.set (triplet.set, nb.triplets = nb.triplets,
                                               alpha = alpha, FDRcontrol = FDRcontrol, fdr = fdr,
                                               lambda = lambda, lambda.step = lambda.step,
                                               pi0.meth = pi0.meth,
                                               data = data, confounders = confounders,
                                               n_t = n_t, n_v = n_v, verbose = verbose,
                                               cl = cl, chunk.size = chunk.size)

      # Step 4.3: Update the adjacency matrix
      if (verbose) {
        cat(if (nb.triplets == 1) paste0("      # Updating the adjacency matrix: adding edge ",
                                         "direction for one [Ti,Tj,Tk] trio. \n")
            else paste0("      # Updating the adjacency matrix: adding edge directions for ",
                        nb.triplets, " [Ti,Tj,Tk] trios. \n"))
      }
      update.object <- update.adjacency.matrix (adjacency = adjacency,
                                                n_v = n_v,
                                                trio.set = triplet.set,
                                                inferred.models = triplet.analysis$Inferred.Model,
                                                stringent = stringent,
                                                add.edges = FALSE,
                                                solve.conflicts = solve.conflicts,
                                                method = method,
                                                added.edges = added.edges,
                                                dropped.edges = dropped.edges,
                                                cl = cl, chunk.size = chunk.size)
      dropped.edges <- update.object$dropped.edges

      # Update adjacency
      adjacency <- update.object$adjacency
    }

    # Step 4.4: Form trio with T-nodes and Q-nodes
    if (n_q > 0) {
      if (verbose) {
        cat("\n      # Generating an initial list of [Ti,Tj,Qk] trios. \n")
      }
      Qtriplet.set <- enumerate.triplets(adjacency = adjacency,
                                         n_t = n_t, n_v = n_v,
                                         n_q = n_q,
                                         TTT = FALSE,
                                         cl = cl, chunk.size = chunk.size)

      # Number of initial trios
      nb.Qtriplets <- length(Qtriplet.set$types)
      if (!nb.Qtriplets) {
        if (verbose) {
          cat("          * no updateable [Ti,Tj,Qk] trio found. \n")
        }
        Qtriplet.set <- Qtriplet.analysis <- NULL
      }
      else {
        Qtriplet.set <- cbind(Qtriplet.set$triplets, type = Qtriplet.set$types)

        if (verbose) {
          cat(paste0("          * ",
                     if (nb.Qtriplets == 1) "one updateable [Ti,Tj,Tk] trio"
                     else paste0(nb.Qtriplets, " updateable  [Ti,Tj,Tk] trios"),
                     " generated. \n"))
        }

        # Step 4.5: Determine the directed structure for each trio
        if (verbose)
          cat("      # Inferring the causal network for each [Ti,Tj,Qk] trio. \n")
        Qtriplet.analysis <- analyse.triplet.set (Qtriplet.set, nb.triplets = nb.Qtriplets,
                                                  alpha = alpha, FDRcontrol = FDRcontrol, fdr = fdr,
                                                  lambda = lambda, lambda.step = lambda.step,
                                                  pi0.meth = pi0.meth,
                                                  data = data, confounders = confounders,
                                                  n_t = n_t, n_v = n_v, verbose = verbose,
                                                  cl = cl, chunk.size = chunk.size)

        # Step 4.3: Update the adjacency matrix
        if (verbose) {
          cat(if (nb.Qtriplets == 1) paste0("      # Updating the adjacency matrix: adding edge ",
                                            "direction for one [Ti,Tj,Qk] trio. \n")
              else paste0("      # Updating the adjacency matrix: adding edge directions for ",
                          nb.Qtriplets, " [Ti,Tj,Qk] trios. \n"))
        }
        Qupdate.object <- update.adjacency.matrix (adjacency = adjacency,
                                                   n_v = n_v,
                                                   trio.set = Qtriplet.set,
                                                   inferred.models = Qtriplet.analysis$Inferred.Model,
                                                   stringent = stringent,
                                                   add.edges = FALSE,
                                                   solve.conflicts = solve.conflicts,
                                                   method = method,
                                                   added.edges = added.edges,
                                                   dropped.edges = dropped.edges,
                                                   cl = cl, chunk.size = chunk.size)
        dropped.edges <- Qupdate.object$dropped.edges

        # Update adjacency
        adjacency <- Qupdate.object$adjacency
      }
    }

    ## Iterate Step 4
    # Indicator of the presence of undirected edges in the updated network
    Not.fully.directed <- any(rowSums((adjacency + t(adjacency)) == 2) > 0)

    # Repeat step 4: re-evaluate the type of each triplet with two edges, and one undirected.
    # Re-interpret the result of triplet analysis, and update adjacency.
    # Number of [Vi,Tj,Tk] trios generated at each iteration
    iter.nb.triplets <- nb.triplets
    new.triplet.set <- triplet.set
    new.triplet.analysis <- triplet.analysis
    other.triplets <- if (nb.triplets) new.triplet.analysis$Inferred.Model == "Other" else FALSE
    new_directions <- TRUE
    if (n_q > 0) {
      # Number of [Vi,Tj,Qk] trios generated at each iteration
      iter.nb.Qtriplets <- if (n_q > 0) nb.Qtriplets
      new.Qtriplet.set <- Qtriplet.set
      new.Qtriplet.analysis <- Qtriplet.analysis
      other.Qtriplets <- if (nb.Qtriplets) new.Qtriplet.analysis$Inferred.Model == "Other" else FALSE
    }
    else {
      other.Qtriplets <- FALSE
    }

    Not.stop <- Not.fully.directed & (any(other.Qtriplets) | any(other.triplets)) # & FALSE

    while(Not.stop & iter < maxiter) {

      if (any(other.triplets)) {
        ## Step 4.1: Extract [Ti,Tj,Tk] trios to be updated
        if (verbose) {
          cat("\n      # Identifying updateable [Ti,Tj,Tk] trios. \n")
        }

        # Save old triplet types
        old.types <- new.triplet.set[other.triplets, 4]

        # Find new triplet types based on the updated adjacency
        new.triplet.set <- get.new.triplet.type (triplet.set = new.triplet.set[other.triplets,, drop = FALSE],
                                                 adjacency = adjacency,
                                                 cl = cl,
                                                 chunk.size = chunk.size)
        keep.new.triplets <- (new.triplet.set[,4] <= 2) & (new.triplet.set[,4] != old.types)
        new.triplet.set <- new.triplet.set[keep.new.triplets,,drop = FALSE]

        # Number of initial triplets
        nb.new.triplets <- sum(keep.new.triplets)

        # Stop if no new updateble trio found
        if (nb.new.triplets == 0) {
          new_directions <- FALSE
          if (verbose) {
            cat("          * no additional [Ti,Tj,Tk] trio found.")
          }
          other.triplets <- FALSE
        }

        if (nb.new.triplets > 0) {
          if (verbose) {
            cat(paste0("          * ",
                       if (nb.new.triplets == 1) "one updateable [Ti,Tj,Tk] trio"
                       else paste0(nb.new.triplets, " updateable  [Ti,Tj,Tk] trios"),
                       " found \n"))
          }

          # Step 4.2: Extract previous triplet analysis results for kept triplets
          new.triplet.analysis <- new.triplet.analysis[other.triplets,, drop = FALSE]
          new.triplet.analysis <- new.triplet.analysis[keep.new.triplets,, drop = FALSE]

          # Step 4.3: Update the adjacency matrix
          if (verbose) {
            cat(if (nb.new.triplets == 1) paste0("      # Updating the adjacency matrix: adding edge ",
                                                 "direction for one [Ti,Tj,Tk] trio. \n")
                else paste0("      # Updating the adjacency matrix: adding edge directions for ",
                            nb.new.triplets, " [Ti,Tj,Tk] trios. \n"))
          }
          update.object <- update.adjacency.matrix (adjacency = adjacency,
                                                    n_v = n_v,
                                                    trio.set = new.triplet.set,
                                                    inferred.models = new.triplet.analysis$Inferred.Model,
                                                    stringent = stringent,
                                                    add.edges = FALSE,
                                                    solve.conflicts = solve.conflicts,
                                                    method = method,
                                                    added.edges = added.edges,
                                                    dropped.edges = dropped.edges,
                                                    cl = cl, chunk.size = chunk.size)
          dropped.edges <- update.object$dropped.edges

          # Update adjacency
          adjacency <- update.object$adjacency

          # Indicator of the presence of triplets of type "Other"
          other.triplets <- new.triplet.analysis$Inferred.Model == "Other"

          # Save the number of additional trios generated
          iter.nb.triplets <- c(iter.nb.triplets, nb.new.triplets) # FIX (2026-10-22): was nb.new.trios

          # Save the additional trios in 'trio.set'
          triplet.set <- rbind(triplet.set, new.triplet.set)

          # Save the additional trios analysis result in 'trio.analysis'
          triplet.analysis <- rbind(triplet.analysis, new.triplet.analysis) # (not useful at this point)

        }

      }

      if (n_q > 0 & any(other.Qtriplets)) {
        ## Step 3.4: Extract [Ti,Tj,Qk] trios to be updated
        if (verbose) {
          cat("\n      # Identifying updateable [Ti,Tj,Qk] trios. \n")
        }

        # Save old triplet types
        old.Qtypes <- new.Qtriplet.set[other.Qtriplets, 4]

        # Find new triplet types based on the updated adjacency
        new.Qtriplet.set <- get.new.triplet.type (triplet.set = new.Qtriplet.set[other.Qtriplets,, drop = FALSE],
                                                  adjacency = adjacency,
                                                  cl = cl,
                                                  chunk.size = chunk.size)
        keep.new.Qtriplets <- (new.Qtriplet.set[,4] <= 2) & (new.Qtriplet.set[,4] != old.Qtypes)
        new.Qtriplet.set <- new.Qtriplet.set[keep.new.Qtriplets,,drop = FALSE]

        # Number of initial triplets
        nb.new.Qtriplets <- sum(keep.new.Qtriplets)

        # Stop if no new updateble trio found
        if (nb.new.Qtriplets == 0) {
          new_directions <- FALSE
          if (verbose) {
            cat("          * no additional [Ti,Tj,Qk] trio found.")
          }

          other.Qtriplets <- FALSE
        }

        if (nb.new.Qtriplets > 0) {
          if (verbose) {
            cat(paste0("          * ",
                       if (nb.new.Qtriplets == 1) "one updateable [Ti,Tj,Qk] trio"
                       else paste0(nb.new.Qtriplets, " updateable  [Ti,Tj,Qk] trios"),
                       " found. \n"))
          }

          # Step 4.2: Extract previous triplet analysis results for kept triplets
          new.Qtriplet.analysis <- new.Qtriplet.analysis[other.Qtriplets,, drop = FALSE]
          new.Qtriplet.analysis <- new.Qtriplet.analysis[keep.new.Qtriplets,, drop = FALSE]

          # Step 4.3: Update the adjacency matrix
          if (verbose) {
            cat(if (nb.new.Qtriplets == 1) paste0("      # Updating the adjacency matrix: adding edge ",
                                                  "direction for one [Ti,Tj,Qk] trio. \n")
                else paste0("      # Updating the adjacency matrix: adding edge directions for ",
                            nb.new.Qtriplets, " [Ti,Tj,Qk] trios. \n"))
          }
          update.object <- update.adjacency.matrix (adjacency = adjacency,
                                                    n_v = n_v,
                                                    trio.set = new.Qtriplet.set,
                                                    inferred.models = new.Qtriplet.analysis$Inferred.Model,
                                                    stringent = stringent,
                                                    add.edges = FALSE,
                                                    solve.conflicts = solve.conflicts,
                                                    method = method,
                                                    added.edges = added.edges,
                                                    dropped.edges = dropped.edges,
                                                    cl = cl, chunk.size = chunk.size)
          dropped.edges <- update.object$dropped.edges

          # Update adjacency
          adjacency <- update.object$adjacency

          # Save the additional trios in 'Qtrio.set'
          Qtriplet.set <- rbind(Qtriplet.set, new.Qtriplet.set)

          # Save the additional trios analysis result in 'Qtrio.analysis'
          Qtriplet.analysis <- rbind(Qtriplet.analysis, new.Qtriplet.analysis) # (not useful at this point)

          # Indicator of the presence of triplets of type "Other"
          other.Qtriplets <- new.Qtriplet.analysis$Inferred.Model == "Other"
        }

        # Save the number of additional trios generated
        iter.nb.Qtriplets <- c(iter.nb.Qtriplets, nb.new.Qtriplets)
      }

      # iteration counter
      iter <- iter + 1

      # Indicator of the presence of undirected edges in the updated network
      Not.fully.directed <- any(rowSums((adjacency + t(adjacency)) == 2) > 0)

      Not.stop <- Not.fully.directed & any(c(other.triplets, other.Qtriplets))
    }

  }

  # Keep the Step 4 skeleton before the optional TTT vote (returned for comparison)
  adjacency.prevote <- adjacency

  ## ---- Optional TTTthresh majority vote (re-added 2026-10-29, default off) ----
  ## Same rule as the TTTthresh option added on 2026-07-26 (commit b2e7ddc) and
  ## removed on 2026-08-19 (commit 82b2fe6): for "candidate trio (2)" triplets
  ## [Ti,Tj,Tk] (Ti -> Tj established, Tj - Tk undirected), all triplets testing
  ## the same Tj-Tk pair are pooled; if the fraction voting M1.1 or M2.1 (i.e. not
  ## "Other") strictly exceeds TTTthresh, Tj-Tk is directed by the majority of
  ## those votes (M1.1: Tj -> Tk, M2.1: Tk -> Tj; a tie leaves it undirected).
  ## Unlike VTTthresh, this never removes an edge (clear.tt = FALSE): Tj-Tk is
  ## known to exist, so pairs at or below the threshold, or never tested by a
  ## type-2 triplet, are left as Step 4's sequential pass produced them.
  if (TTTthresh > 0 && !is.null(triplet.set) && !is.null(triplet.analysis)) {
    type2 <- triplet.set[, 4] == 2
    if (any(type2)) {
      Tidx <- (n_v + 1L):(n_v + n_t)
      adjacency <- .majority_adj(adjacency = adjacency,
                                 ts        = triplet.set[type2, , drop = FALSE],
                                 ana       = triplet.analysis[type2, , drop = FALSE],
                                 Tidx      = Tidx,
                                 keep      = c("M1.1", "M2.1"),
                                 thresh    = TTTthresh,
                                 clear.tt  = FALSE)
      Not.fully.directed <- any(rowSums((adjacency + t(adjacency)) == 2) > 0)
    }
  }

  list(adjacency = adjacency,
       adjacency.prevote = adjacency.prevote,
       triplet.set = triplet.set,
       iter.nb.triplets = iter.nb.triplets,
       triplet.analysis = triplet.analysis,
       Qtriplet.set = Qtriplet.set,
       iter.nb.Qtriplets = iter.nb.Qtriplets,
       Qtriplet.analysis = Qtriplet.analysis,
       added.edges = added.edges,
       dropped.edges = dropped.edges,
       niter = iter,
       Not.fully.directed = Not.fully.directed)
}
