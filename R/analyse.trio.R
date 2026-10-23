#'
#' Infer the causal networks for a set trios
#'
#' This function takes a set of trios including each a genetic variant,
#' phenotypes, and optionally confounding variables, and performs the MRGN inference.
#' This a wrapper for the function \link[MRGN]{infer.trio} of the package \code{MRGN}.
#' \code{analyse.trio.set} returns a \code{data.frame} object with one row for each input trio.
#'
#' @param trio.set numeric matrix with three column, each row giving the numbers
#' of the columns of \code{data} defining a target trio.
#'
#' @param nb.trios number of trios to analyse. If \code{nb.trios} is less than
#' the number of rows of \code{trio.set}, only the first \code{nb.trios} trios are
#' analysed.
#'
#' @param Qlabels numeric vector of column numbers indicating the columns of \code{data}
#' that correspond to common children or intermediate variables in the general network
#' from which trios are formed.
#'
#' @param alpha significance level for individual tests, i.e. the rejection threshold
#' for each Wald tests. Passed to \link[MRGN]{infer.trio}
#'
#' @param FDRcontrol character, method for p-value adjustment for multiple
#' comparisons. Either \code{'qvalue'}, or any method accepted by \link[stats]{p.adjust}.
#'
#' @param fdr,lambda,pi0.meth arguments passed to \link[MRGN]{adjust.q}.
#'
#' @param lambda.step used to define the default value of \code{lambda} when the latter is \code{NULL}.
#'
#' @param use.perm,gamma,is.CNA,nperms,verbose arguments passed to \link[MRGN]{infer.trio}.
#'
#' @param data \code{data.frame} object where columns indexed in \code{trio.set},
#' \code{Qlabels}, \code{confounders} are taken from for trio analysis.
#' The first \code{n_t} columns are \code{V}-nodes (variants), the next \code{n_t}
#' columns are \code{T}-nodes (genes), and the remaining columns are
#' \code{Q}-nodes (common children and intermediate variables), and
#' \code{U}-nodes (confounders).
#'
#' @param confounders \code{n_t}-list of column numbers indexing confounders for each \code{T}-node in \code{data}.
#'
#' @param n_v,n_t integers, numbers of respectively \code{T}-nodes and \code{V}-nodes
#' in \code{data}.
#'
#' @param cl,chunk.size optional arguments for parallel computing, passed to
#' \link[parallel]{parLapply} (when supplied).
#'
#' @importFrom MRGN infer.trio
#' @importFrom MRGN adjust.q
#'
#' @export analyse.trio.set
#'
#'
# Wrap 'analyse.trio.i' over a set of trios
analyse.trio.set <- function(trio.set, nb.trios = NROW(trio.set), Qlabels = NULL,
                             alpha, FDRcontrol, fdr, lambda, lambda.step,
                             pi0.meth = "bootstrap",
                             data, confounders, n_v,n_t,
                             use.perm, gamma, is.CNA,
                             nperms, verbose, cl = NULL, chunk.size = NULL,
                             ## ---- NEW (2026-09-17): threaded through to
                             ## analyse.trio.i() -- see there for details.
                             ## Default TRUE since 2026-10-22 (matches MRGN()).
                             V.trioconfounders.filter = TRUE,
                             V.filter.cutoff = 0.2) {   # was 0.35 before 2026-10-22
  switch(FDRcontrol,
         none = { # No correction
           all.stats <- matteApply(trio.set,
                                   MARGIN = 1,
                                   FUN = analyse.trio.i,
                                   Qlabels = Qlabels,
                                   data = data, confounders = confounders, n_t = n_t, n_v = n_v,
                                   use.perm = use.perm, gamma = gamma, is.CNA = is.CNA,
                                   alpha = alpha, nperms = nperms, half1 = FALSE, verbose = verbose > 1,
                                   V.trioconfounders.filter = V.trioconfounders.filter, V.filter.cutoff = V.filter.cutoff,
                                   cl = cl, chunk.size = chunk.size)
           all.stats <- do.call('rbind', all.stats)
         },
         qvalue = { # qvalue method
           # Half trio analysis
           all.stats <- matteApply(trio.set,
                                   MARGIN = 1,
                                   FUN = analyse.trio.i,
                                   Qlabels = Qlabels,
                                   data = data, confounders = confounders, n_t = n_t, n_v = n_v,
                                   use.perm = use.perm, gamma = gamma, is.CNA = is.CNA,
                                   nperms = nperms, half1 = TRUE, verbose = verbose > 1,
                                   V.trioconfounders.filter = V.trioconfounders.filter, V.filter.cutoff = V.filter.cutoff,
                                   cl = cl, chunk.size = chunk.size)
           all.stats <- do.call('rbind', all.stats)
           ## FIX (2026-10-07): select by name, not position 8 -- MRGN >= 1.3.9's
           ## infer.trio() inserted coef11/coef12/coef21/coef22 columns before
           ## Inferred.Model, which silently shifted what a fixed numeric index
           ## pointed to. See R/analyse.trio.i()'s half1 extraction below for the
           ## matching fix on the upstream side.
           Inferred.Model0 <- all.stats[,"Inferred.Model"]

           # FDR control: apply qvalue correction
           if (is.null(lambda)) {
             lambda <- seq(min(c(as.vector(as.matrix(all.stats[, 1:6])), lambda.step), na.rm = T),
                           max(c(as.vector(as.matrix(all.stats[, 1:6])), lambda.step), na.rm = T), lambda.step)
           }

           Qvalues <- catch.conditions (adjust.q (as.vector(as.matrix(all.stats[, 1:6])), fdr = fdr,
                                                  lambda = lambda, pi0.meth = pi0.meth))$value
           if (any(class(Qvalues) %in% c("simpleError", "error", "condition"))) {
             if (identical(pi0.meth, "smoother")) {
               if (verbose) {
                 cat("      #         * q value correction with 'pi0.meth = smoother' failled, trying 'pi0.meth = bootstrap'  \n")
               }
               Qvalues <- catch.conditions (adjust.q (as.vector(as.matrix(all.stats[, 1:6])), fdr = fdr,
                                                      lambda = lambda, pi0.meth = "bootstrap"))$value
             }
             if (any(class(Qvalues) %in% c("simpleError", "error", "condition"))) {
               if ((ll <- length(lambda)) > 1) {
                 Qvalues <- catch.conditions (adjust.q (as.vector(as.matrix(all.stats[, 1:6])),
                                                        fdr = fdr, lambda = lambda[ll],
                                                        pi0.meth = pi0.meth))$value
                 if (any(class(Qvalues) %in% c("simpleError", "error", "condition"))) {
                   Qvalues <- catch.conditions (adjust.q (as.vector(as.matrix(all.stats[, 1:6])),
                                                          fdr = fdr, lambda = lambda[1],
                                                          pi0.meth = pi0.meth))$value
                 }
               }

               #if (any(class(Qvalues) %in% c("simpleError", "error", "condition"))) {
               #  Qvalues <- catch.conditions (adjust.q (as.vector(as.matrix(all.stats[, 1:6])),
               #                                         fdr = fdr, lambda = 0.05,
               #                                         pi0.meth = pi0.meth))$value
               #}

               if (any(class(Qvalues) %in% c("simpleError", "error", "condition"))) {
                 if (verbose) {
                   cat("            -- 'qvalue::qvalue' failed. q value correction skiped. \n")
                 }
                 Qvalues <- list(qvalue = as.vector(as.matrix(all.stats[, 1:6])),
                                 significant = as.vector(as.matrix(all.stats[, 1:6])) <= alpha)
               }
             }
           }

           all.stats <- cbind(matrix(Qvalues$significant, nrow = nb.trios, byrow = FALSE),
                              matrix(Qvalues$qvalue, nrow = nb.trios, byrow = FALSE),
                              Minor.freq = all.stats[,7])
           all.stats = as.data.frame(all.stats)
           colnames(all.stats) <- c("b11","b12", "b21","b22", "V1:T1", "V1:T2", "pb11",
                                    "pb12", "pb21","pb22","pV1:T1","pV1:T2", "Minor.freq")

           # Infer model structure (Second half trio analysis)
           # Replace NAs in significance columns with FALSE to guard against class.vec NA crash
           all.stats[, 1:6][is.na(all.stats[, 1:6])] <- FALSE
           all.stats$Inferred.Model <- c(matteApply(all.stats, MARGIN = 1, FUN = class.vec,
                                                    cl = cl, chunk.size = chunk.size))

           # Also return the inferred model with no correction
           all.stats$Inferred.Model0 <- Inferred.Model0

         },
         { # BH (MRGN() default since 2026-10-22), Bonferroni and all other stats::p.adjust methods
           # Half trio analysis
           all.stats <- matteApply(trio.set,
                                   MARGIN = 1,
                                   FUN = analyse.trio.i,
                                   Qlabels = Qlabels,
                                   data = data, confounders = confounders, n_t = n_t, n_v = n_v,
                                   use.perm = use.perm, gamma = gamma, is.CNA = is.CNA,
                                   nperms = nperms, half1 = TRUE, verbose = verbose > 1,
                                   V.trioconfounders.filter = V.trioconfounders.filter, V.filter.cutoff = V.filter.cutoff,
                                   cl = cl, chunk.size = chunk.size)
           all.stats <- do.call('rbind', all.stats)
           ## FIX (2026-10-07): select by name, not position 8 -- see matching
           ## fix in the qvalue branch above for the full explanation.
           Inferred.Model0 <- all.stats[,"Inferred.Model"]

           # p.adjust correction (BH, Bonferroni, ...), pooled over all trios' p-values
           p.adj <- stats::p.adjust(as.vector(as.matrix(all.stats[, 1:6])), method=FDRcontrol)
           p.adj <- matrix(p.adj, nrow = nb.trios, byrow = FALSE)
           all.stats <- cbind(p.adj <= fdr, p.adj,
                              Minor.freq = all.stats[,7])
           all.stats = as.data.frame(all.stats)
           colnames(all.stats) <- c("b11","b12", "b21","b22", "V1:T1", "V1:T2", "pb11",
                                    "pb12", "pb21","pb22","pV1:T1","pV1:T2", "Minor.freq")

           # Infer model structure (Second half trio analysis)
           # Replace NAs in significance columns with FALSE to guard against class.vec NA crash
           all.stats[, 1:6][is.na(all.stats[, 1:6])] <- FALSE
           all.stats$Inferred.Model <- c(matteApply(all.stats, MARGIN = 1, FUN = class.vec,
                                                    cl = cl, chunk.size = chunk.size))

           # Also return the inferred model with no correction
           all.stats$Inferred.Model0 <- Inferred.Model0

         })

  ### Distinguish 'Other.1' and 'Other.2' trio structures
  ## FIX (2026-10-07): select/assign by name, not position 14. Under the
  ## "none" branch, all.stats keeps the FULL width returned by MRGN::infer.trio()
  ## (no half1 truncation), so with MRGN >= 1.3.9 (which inserts coef11/coef12/
  ## coef21/coef22 before Inferred.Model), column 14 there is "coef11" (numeric),
  ## not "Inferred.Model" (character) -- the old position-based check/assignment
  ## silently compared/wrote the wrong column and this reclassification never
  ## fired. Columns 1:13 (b11..Minor.freq) are unaffected by the new columns and
  ## still safe to pass to class.other() positionally.
  Others <- all.stats[, "Inferred.Model"] == "Other"
  if (any(Others)) {
    all.stats[Others, "Inferred.Model"] <- c(matteApply(all.stats[Others, 1:13, drop = FALSE],
                                          MARGIN = 1, FUN = class.other,
                                          cl = cl, chunk.size = chunk.size))
  }

  return(all.stats)

}

# Extract the correct confounder set for a trio and Call infer.trio for trio analysis
analyse.trio.i <- function (col.indices, # indicate column number in the adjacency matrix
                            data, n_v, n_t,
                            confounders, # indicate column number in the data fame
                            Qlabels = NULL, # NULL means that 'col.indices' also indicate positions in data (not just in the Adj matrix, i.e. columns of data are ordered so that they match those of Adj)
                            use.perm = TRUE, gamma = 0.05, is.CNA = FALSE,
                            alpha = 0.01, nperms = 10000, half1 = FALSE, verbose = FALSE,
                            ## ---- NEW (2026-09-17): per-trio, dynamic confounder
                            ## filter by correlation with THIS trio's own Vi. Default
                            ## TRUE since 2026-10-22 (was FALSE; matches MRGN()).
                            V.trioconfounders.filter = TRUE,
                            V.filter.cutoff = 0.2) {   # was 0.35 before 2026-10-22

  # Replace Adj column index by data column index if any of the two T-nodes is a Q-node
  n_vt <- n_t + n_v
  if (!is.null(Qlabels)) {
    if (col.indices[2] > n_vt) {
      col.indices[2] <- Qlabels[col.indices[2] - n_vt]
    }
    if (col.indices[3] > n_vt) {
      col.indices[3] <- Qlabels[col.indices[3] - n_vt]
    }
  }

  # For each trio, take the union of the confounders for the two T nodes
  conf.set.trio <- c(if (col.indices[2] <= n_vt) confounders[[col.indices[2] - n_v]], # Take confounder index for T-node
                     if (col.indices[3] <= n_vt) confounders[[col.indices[3] - n_v]]) # Result is NULL if Q-node

  #conf.set.trio <- unique(conf.set.trio)

  # Find and remove duplicates (V-nodes, T-nodes), if any
  #VTduplicates <- conf.set.trio %in% col.indices
  #if (any(VTduplicates))
  #  conf.set.trio <- conf.set.trio[!VTduplicates]

  ### ---- BEGIN NEW STEP (2026-09-17): dynamic, per-trio confounder filter ----
  ## Unlike the Tcandidates filter (get.conf.sets.R, 2026-09-08, which restricts
  ## the T-T confounder POOL once, statically, per T-node, based on that T-node's
  ## own OVERALL set of selected V-confounders), this filters the confounder set
  ## ALREADY ASSEMBLED for THIS SPECIFIC trio (conf.set.trio, the union of Tj's
  ## and Tk's selected confounders as usual -- confounder selection in
  ## get.conf.sets() itself is completely untouched), based on correlation with
  ## col.indices[1] = Vi, the ACTUAL V-node anchoring this one trio -- not Tj's
  ## or Tk's V-confounders in general, which may differ from Vi. A candidate
  ## confounder correlated with this trio's own instrument (Vi) is excluded from
  ## THIS trio's regression only; other trios sharing the same Tj/Tk but a
  ## different Vi are filtered independently, using their own Vi.
  if (V.trioconfounders.filter && length(conf.set.trio) > 0) {
    Vi_col <- col.indices[1]
    cors_to_Vi <- sapply(conf.set.trio, function(cidx) {
      cor(data[, Vi_col], data[, cidx], use = "pairwise.complete.obs")
    })
    cors_to_Vi[is.na(cors_to_Vi)] <- 0  # guard: zero-variance confounder -> NA correlation
    conf.set.trio <- conf.set.trio[abs(cors_to_Vi) < V.filter.cutoff]
  }
  ### ---- END NEW STEP (2026-09-17) ----

  # All indices of variables in the trio (including counfounders)
  all_indices <- unique(c(col.indices, conf.set.trio))

  # Call infer.trio (drop rows with any NA in this trio's columns)
  trio.res <- tryCatch({
    infer.trio(na.omit(data[, all_indices]),
               use.perm = use.perm, gamma = gamma, is.CNA = is.CNA,
               alpha = alpha, nperms = nperms, verbose = verbose)
  }, error = function(e) {
    out <- data.frame(matrix(NA_real_, nrow = 1L, ncol = 13L))
    names(out) <- c("b11", "b12", "b21", "b22", "V1:T1", "V1:T2",
                    "pb11", "pb12", "pb21", "pb22", "pV1:T1", "pV1:T2", "Minor.freq")
    out$Inferred.Model <- NA_character_
    out
  })

  # If half1 = TRUE, only return the 8 columns required to correct p-values (and re-infer trio model)
  ## FIX (2026-10-07): select by name, not position 7:14. MRGN >= 1.3.9's
  ## infer.trio() inserts coef11/coef12/coef21/coef22 between Minor.freq and
  ## Inferred.Model, so trio.res now has 18 columns instead of 14 -- a fixed
  ## 7:14 slice silently grabbed pb11..Minor.freq+coef11 instead of
  ## pb11..Minor.freq+Inferred.Model, dropping Inferred.Model entirely from
  ## the qvalue/bonferroni branches' pooled correction step.
  if (half1)
    trio.res[, c("pb11", "pb12", "pb21", "pb22", "pV1:T1", "pV1:T2",
                 "Minor.freq", "Inferred.Model")]
  else
    trio.res
}

### Distinguish 'Other.1' and 'Other.2' trio structures
class.other <- function (vec) {
  # Look at pb12 and pb22
  if (any(is.na(vec[c(2,4)])) || any(vec[c(2,4)] == 1)) {
    return("Other.1") # We keep an edge ("Other.1") if the partial correlation is significant or undetermined (NA)
  }
  else {
    return("Other.2")
  }
}
