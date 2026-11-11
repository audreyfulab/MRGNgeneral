# Wrap 'analyse.triplet.i' over a set of triplets
analyse.triplet.set <- function (triplet.set, nb.triplets = NROW(triplet.set), Qlabels = NULL,
                                 alpha, FDRcontrol, fdr, lambda, lambda.step, pi0.meth = "bootstrap",
                                 data, confounders, n_t, n_v, verbose,
                                 cl = NULL, chunk.size = NULL) {
  switch(FDRcontrol,
         none = {
           all.stats <- matteApply(triplet.set,
                                          MARGIN = 1,
                                          FUN = analyse.triplet.i,
                                          Qlabels = Qlabels,
                                          data, n_t = n_t, n_v = n_v,
                                          confounders = confounders,
                                          alpha = alpha, half1 = FALSE,
                                          verbose = verbose > 1,
                                          cl = cl, chunk.size  = chunk.size)
           all.stats <- do.call('rbind', all.stats)
         },
         qvalue = { # qvalue method
           # Half trio analysis
           all.stats <- matteApply(triplet.set,
                                   MARGIN = 1,
                                   FUN = analyse.triplet.i,
                                   Qlabels = Qlabels,
                                   data, n_t = n_t, n_v = n_v,
                                   confounders = confounders,
                                   alpha = alpha, half1 = TRUE, verbose = verbose > 1,
                                   cl = cl, chunk.size  = chunk.size)

           # FDR control: apply qvalue correction
           if (is.null(lambda)) {
             lambda <- seq(min(c(all.stats, lambda.step), na.rm = T),
                           max(c(all.stats, lambda.step), na.rm = T), lambda.step)
           }

           # There can be PROBLEM in 'pi0est': there is bug at lines 28-29 (debugging line numbers)
           Qvalues <- catch.conditions (adjust.q (c(all.stats), fdr = fdr, lambda = lambda,
                                                  pi0.meth = pi0.meth))$value

           if (any(class(Qvalues) %in% c("simpleError", "error", "condition"))) {
             if (identical(pi0.meth, "smoother")) {
               if (verbose) {
                 cat("      #         * q value correction with 'pi0.meth = smoother' failled, trying 'pi0.meth = bootstrap'  \n")
               }
               Qvalues <- catch.conditions (adjust.q (c(all.stats), fdr = fdr,
                                                      lambda = lambda,
                                                      pi0.meth = "bootstrap"))$value
             }
             if (any(class(Qvalues) %in% c("simpleError", "error", "condition"))) {
               if ((ll <- length(lambda)) > 1) {
                 Qvalues <- catch.conditions (adjust.q (c(all.stats), fdr = fdr, lambda = lambda[ll],
                                                        pi0.meth = pi0.meth))$value
                 if (any(class(Qvalues) %in% c("simpleError", "error", "condition"))) {
                   Qvalues <- catch.conditions (adjust.q (c(all.stats), fdr = fdr, lambda = lambda[1],
                                                          pi0.meth = pi0.meth))$value
                 }
               }
               if (any(class(Qvalues) %in% c("simpleError", "error", "condition"))) {
                 if (verbose) {
                   cat("            -- 'qvalue::qvalue' failed. \n")
                 }
                 Qvalues <- list(qvalue = c(all.stats),
                                 significant = c(all.stats) <= alpha)
               }
             }
           }

           ## CHANGE (2026-11-11): NA is no longer treated as "not rejected" (which
           ## gave type-2 triplets M1.1, i.e. a direction). NA -> "Other" (no direction).
           model <- triplet.model(Qvalues$significant, triplet.set[,4])
           all.stats <- cbind(Qvalues$significant,
                              Qvalues$qvalue,
                              model)
           all.stats = as.data.frame(all.stats)
           colnames(all.stats) <- c("rej", "pval", "Inferred.Model")

         },
         { # BH (default), Bonferroni and all other stats::p.adjust methods
           # Half trio analysis
           all.stats <- matteApply(triplet.set,
                                   MARGIN = 1,
                                   FUN = analyse.triplet.i,
                                   Qlabels = Qlabels,
                                   data, n_t = n_t, n_v = n_v,
                                   confounders = confounders, alpha = alpha,
                                   half1 = TRUE, verbose = verbose > 1,
                                   cl = cl, chunk.size  = chunk.size)

           # p.adjust correction (BH, Bonferroni, ...)
           p.adj <- stats::p.adjust(c(all.stats), method = FDRcontrol)
           rej <- p.adj <= fdr
           ## CHANGE (2026-11-11): NA kept; NA -> "Other" (no direction), see triplet.model()
           model <- triplet.model(rej, triplet.set[,4])
           all.stats <- cbind(rej,
                              p.adj,
                              model)

           all.stats <- as.data.frame(all.stats)
           colnames(all.stats) <- c("rej", "pval", "Inferred.Model")

         })

  return(all.stats)

}

# Triplet model from the (adjusted) rejection decision (added 2026-11-11)
#   rejected (TRUE)            -> "M2.1" (Tk -> Tj)
#   not rejected (FALSE), type 2 -> "M1.1" (Tj -> Tk); other types -> "Other"
#   undetermined (NA)          -> "Other" (no direction; the edge is left as it is)
triplet.model <- function (rej, type) {
  model <- rep("Other", length(rej))
  model[rej %in% TRUE] <- "M2.1"
  model[rej %in% FALSE & type == 2] <- "M1.1"
  model
}

# Call infer.triplet for triplet analysis
analyse.triplet.i <- function (col.indices,
                               Qlabels = NULL,
                               data, n_v, n_t, confounders,
                               alpha = 0.01,
                               half1 = FALSE,
                               verbose = FALSE) {
  # Replace Adj column index by data column index if any of the three nodes is a Q-node
  n_vt <- n_v + n_t
  if (!is.null(Qlabels)) {
    if (col.indices[1] > n_vt) {
      col.indices[1] <- Qlabels[col.indices[1] - n_vt]
    }
    if (col.indices[2] > n_vt) {
      col.indices[2] <- Qlabels[col.indices[2] - n_vt]
    }
    if (col.indices[3] > n_vt) {
      col.indices[3] <- Qlabels[col.indices[3] - n_vt]
    }
  }

  # For each triplet, take the union of the confounders for the three T nodes
  conf.set.triplet <- c(if (col.indices[1] <= n_vt) confounders[[col.indices[1] - n_v]], # Take confounder index if T-node
                        if (col.indices[2] <= n_vt) confounders[[col.indices[2] - n_v]], # Result is NULL if I&C node
                        if (col.indices[3] <= n_vt) confounders[[col.indices[3] - n_v]])
  #conf.set.triplet <- unique(conf.set.triplet)

  # All indices of variables in the trio (including counfounders)
  all_indices <- unique(c(col.indices[-4], conf.set.triplet))

  infer.triplet(data[, all_indices], # c(col.indices[-4], conf.set.triplet)
                type = col.indices[4],
                alpha = alpha,
                half1 = half1,
                verbose = verbose)
}
