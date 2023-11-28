#'
#' Infer the causal networks for a set trios
#'
#' This function takes a set of trios including each a genetic variant,
#' phenotypes, and optionally confounding variables, and performs the MRGN inference.
#' This a wrapper for the function \link{infer.trio}{MRGN} of the package \code{MRGN}.
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
#' for each Wald tests. Passed to \link{infer.trio}{MRGN}
#'
#' @param FDRcontrol character, method for p-value adjustment for multiple
#' comparisons. Either \code{'qvalue'}, or any method accepted by \link{p.adjust}{stats}.
#'
#' @param fdr,lambda,pi0.meth arguments passed to \link{adjust.q}{MRGN}.
#'
#' @param lambda.step used to define the default value of \code{lambda} when the latter is \code{NULL}.
#'
#' @param use.perm,gamma,is.CNA,nperms,verbose arguments passed to \link{infer.trio}{MRGN}.
#'
#' @param data \code{data.frame} object where columns indexed in \code{trio.set},
#' \code{Qlabels}, \code{confounders} are taken from for trio analysis.
#' The first \code{q} columns are \code{V}-nodes (variants), the next \code{p}
#' columns are \code{T}-nodes (genes), and the remaining columns are confounders,
#' common children, and intermediate variables.
#'
#' @param confounders \code{p}-list of column numbers indexing confounders for each \code{T}-node in \code{data}.
#'
#' @param p,q integers, numbers of respectively \code{T}-nodes and \code{V}-nodes
#' in \code{data}.
#'
#' @param cl,chunk.size optional arguments, passed to \link{parLapply}{parallel} (when supplied) for parallel computing.
#'
#' @importFrom MRGN infer.trio
#'
#' @export analyse.trio.set
#'
#'
# Wrap 'analyse.trio.i' over a set of trios
analyse.trio.set <- function(trio.set, nb.trios = NROW(trio.set), Qlabels = NULL,
                             alpha, FDRcontrol, fdr, lambda, lambda.step,
                             pi0.meth = "bootstrap",
                             data, confounders, p, q,
                             use.perm, gamma, is.CNA,
                             nperms, verbose, cl = NULL, chunk.size = NULL) {
  switch(FDRcontrol,
         none = { # No correction
           all.stats <- matteApply(trio.set,
                                   MARGIN = 1,
                                   FUN = analyse.trio.i,
                                   Qlabels = Qlabels,
                                   data = data, confounders = confounders, p = p, q = q,
                                   use.perm = use.perm, gamma = gamma, is.CNA = is.CNA,
                                   alpha = alpha, nperms = nperms, half1 = FALSE, verbose = verbose > 1,
                                   cl = cl, chunk.size = chunk.size)
           all.stats <- do.call('rbind', all.stats)
         },
         qvalue = { # qvalue method
           # Half trio analysis
           all.stats <- matteApply(trio.set,
                                   MARGIN = 1,
                                   FUN = analyse.trio.i,
                                   Qlabels = Qlabels,
                                   data = data, confounders = confounders, p = p, q = q,
                                   use.perm = use.perm, gamma = gamma, is.CNA = is.CNA,
                                   nperms = nperms, half1 = TRUE, verbose = verbose > 1,
                                   cl = cl, chunk.size = chunk.size)
           all.stats <- do.call('rbind', all.stats)
           Inferred.Model0 <- all.stats[,8]

           # FDR control: apply qvalue correction
           if (is.null(lambda)) {
             lambda <- seq(min(c(unlist(all.stats[,1:6]), lambda.step), na.rm = T),
                           max(c(unlist(all.stats[,1:6]), lambda.step), na.rm = T), lambda.step)
           }

           Qvalues <- catch.conditions (adjust.q (unlist(all.stats[,1:6]), fdr = fdr,
                                                  lambda = lambda, pi0.meth = pi0.meth))$value
           if (any(class(Qvalues) %in% c("simpleError", "error", "condition"))) {
             if (identical(pi0.meth, "smoother")) {
               if (verbose) {
                 cat("      #         * q value correction with 'pi0.meth = smoother' failled, trying 'pi0.meth = bootstrap'  \n")
               }
               Qvalues <- catch.conditions (adjust.q (unlist(all.stats[,1:6]), fdr = fdr,
                                                      lambda = lambda, pi0.meth = "bootstrap"))$value
             }
             if (any(class(Qvalues) %in% c("simpleError", "error", "condition"))) {
               if ((ll <- length(lambda)) > 1) {
                 Qvalues <- catch.conditions (adjust.q (unlist(all.stats[,1:6]),
                                                        fdr = fdr, lambda = lambda[ll],
                                                        pi0.meth = pi0.meth))$value
                 if (any(class(Qvalues) %in% c("simpleError", "error", "condition"))) {
                   Qvalues <- catch.conditions (adjust.q (unlist(all.stats[,1:6]),
                                                          fdr = fdr, lambda = lambda[1],
                                                          pi0.meth = pi0.meth))$value
                 }
               }

               if (any(class(Qvalues) %in% c("simpleError", "error", "condition"))) {
                 Qvalues <- catch.conditions (adjust.q (unlist(all.stats[,1:6]),
                                                        fdr = fdr, lambda = 0.05,
                                                        pi0.meth = pi0.meth))$value
               }

               if (any(class(Qvalues) %in% c("simpleError", "error", "condition"))) {
                 if (verbose) {
                   cat("            -- 'qvalue::qvalue' failed. q value correction skiped. \n")
                 }
                 Qvalues <- list(qvalue = c(all.stats[,1:6]),
                                 significant = c(all.stats[,1:6]) <= alpha)
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
           all.stats$Inferred.Model <- c(matteApply(all.stats, MARGIN = 1, FUN = class.vec,
                                                    cl = cl, chunk.size = chunk.size))

           # Also return the inferred model with no correction
           all.stats$Inferred.Model0 <- Inferred.Model0

         },
         { # Bonferroni and all other Multiple Comparisons adjustment methods
           # Half trio analysis
           all.stats <- matteApply(trio.set,
                                   MARGIN = 1,
                                   FUN = analyse.trio.i,
                                   Qlabels = Qlabels,
                                   data = data, confounders = confounders, p = p, q = q,
                                   use.perm = use.perm, gamma = gamma, is.CNA = is.CNA,
                                   nperms = nperms, half1 = TRUE, verbose = verbose > 1,
                                   cl = cl, chunk.size = chunk.size)
           all.stats <- do.call('rbind', all.stats)
           Inferred.Model0 <- all.stats[,8]

           # Bonferroni correction
           p.adj <- stats::p.adjust(unlist(all.stats[,1:6]), method=FDRcontrol)
           p.adj <- matrix(p.adj, nrow = nb.trios, byrow = FALSE)
           all.stats <- cbind(p.adj <= fdr, p.adj,
                              Minor.freq = all.stats[,7])
           all.stats = as.data.frame(all.stats)
           colnames(all.stats) <- c("b11","b12", "b21","b22", "V1:T1", "V1:T2", "pb11",
                                    "pb12", "pb21","pb22","pV1:T1","pV1:T2", "Minor.freq")

           # Infer model structure (Second half trio analysis)
           all.stats$Inferred.Model <- c(matteApply(all.stats, MARGIN = 1, FUN = class.vec,
                                                    cl = cl, chunk.size = chunk.size))

           # Also return the inferred model with no correction
           all.stats$Inferred.Model0 <- Inferred.Model0

         })

  ### Distinguish 'Other.1' and 'Other.2' trio structures
  Others <- all.stats[, 14] == "Other"
  if (any(Others)) {
    all.stats[Others, 14] <- c(matteApply(all.stats[Others, 1:13, drop = FALSE],
                                          MARGIN = 1, FUN = class.other,
                                          cl = cl, chunk.size = chunk.size))
  }

  return(all.stats)

}

# Extract the correct confounder set for a trio and Call infer.trio for trio analysis
analyse.trio.i <- function (col.indices, # indicate column number in the adjacency matrix
                            data, p, q,
                            confounders, # indicate column number in the data fame
                            Qlabels = NULL, # NULL means that 'col.indices' also indicate positions in data (not just in the Adj matrix, i.e. columns of data are ordered so that they match those of Adj)
                            use.perm = TRUE, gamma = 0.05, is.CNA = FALSE,
                            alpha = 0.01, nperms = 10000, half1 = FALSE, verbose = FALSE) {

  # Replace Adj column index by data column index if any of the two T-nodes is a Q-node
  if (!is.null(Qlabels)) {
    if (col.indices[2] > (q + p)) {
      col.indices[2] <- Qlabels[col.indices[2] - (p + q)]
    }
    if (col.indices[3] > (q + p)) {
      col.indices[3] <- Qlabels[col.indices[3] - (p + q)]
    }
  }

  # For each trio, take the union of the confounders for the two T nodes
  conf.set.trio <- c(if (col.indices[2] <= q + p) confounders[[col.indices[2] - q]], # Take confounder index for T-node
                     if (col.indices[3] <= q + p) confounders[[col.indices[3] - q]]) # Result is NULL if Q-node

  #conf.set.trio <- unique(conf.set.trio)

  # Find and remove duplicates (V-nodes, T-nodes), if any
  #VTduplicates <- conf.set.trio %in% col.indices
  #if (any(VTduplicates))
  #  conf.set.trio <- conf.set.trio[!VTduplicates]

  # All indices of variables in the trio (including counfounders)
  all_indices <- unique(c(col.indices, conf.set.trio))

  # Call infer.trio
  trio.res <- infer.trio(data[, all_indices],
                         use.perm = use.perm, gamma = gamma, is.CNA = is.CNA,
                         alpha = alpha, nperms = nperms, verbose = verbose)

  # If half1 = TRUE, only return the last 8 columns required to correct p-values  (and re-infer trio model)
  if (half1)
    trio.res[,7:14]
  else
    trio.res
}

### Distinguish 'Other.1' and 'Other.2' trio structures
class.other <- function (vec) {
  # Look at pb12 and pb22
  if (any(vec[c(2,4)] == 1)) {
    return("Other.1") # We keep an edge ("Other.1") if the partial correlation is significant
  }
  else {
    return("Other.2")
  }
}
