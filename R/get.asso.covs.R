## Perform q value correction, or p-value adjustment to select confounding variables
# if required (code adapted from Jarred's 'get.conf.matrix' function)

get.asso.covs <- function(n_t,     # Number of focal variables we are selection confounding variables for (T-nodes)
                          p.mat,   # matrix of p-values
                          symmetric = FALSE,
                          FDRcontrol, adjust_by,
                          fdr, alpha,
                          lambda, pi0.method) {
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
                    qsig.mat <- adjust.q (c(p.mat),
                                          fdr = fdr,
                                          lambda = lambda,
                                          pi0.meth = pi0.method)

                    #restructure output into significance matrix of snps X covariates
                    sig.mat <- diag(n_t)
                    sig.mat[lower.tri(sig.mat, diag = FALSE)] <- qsig.mat$significant # (lower triangle)
                    sig.mat <- (sig.mat + t(sig.mat)) > 0
                    diag(sig.mat) <- FALSE
                    #restructure output into qvalue matrix of snps X covariates
                    q.mat <- diag(n_t)
                    q.mat[lower.tri(q.mat, diag = FALSE)] <- qsig.mat$qvalue # (lower triangle)
                    q.mat <- q.mat + t(q.mat)
                    diag(q.mat) <- 1
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

                  },
                  none = {
                    p.adj.mat <- q.mat <- matrix(NA, nrow = n_t, ncol = n_t)
                    sig.mat <- p.mat <= alpha
                  },
                  stop(paste0('Input to argument \'adjust_by \' = ', adjust_by, ' is not recognized. use one of \'individual\', \' all \', or \'none\' '))
           )
         })


  #-------------------obtain-final-list-of-significant-covs----------------
  #find the covs that correlated with every column of the trio matrix
  if (verbose) {
    cat("            selecting significant covariates... \n")
  }
  #extract significant covariates:
  sig.asso.covs <- apply(sig.mat,
                         MARGIN = 2,
                         FUN = function(x){which(x)}, simplify = FALSE)


  out <- list(sig.asso.covs = sig.asso.covs,
              pvalues = p.mat,
              qvalues = q.mat,
              cors = r.mat,
              sig = sig.mat,
              adj.p = p.adj.mat)

  return(out)
}
