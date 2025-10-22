
#' @importFrom propagate bigcor

## All functions here are taken from the package MRGN with slight modifications for parallel computing and time efficiency
## And controlling message display (verbose)
get.conf.matrix <- function (data = NULL,
                             cov.pool = NULL,
                             measure = c("correlation", "partial_corr"),
                             conditional.vars = NULL,
                             blocksize = 2000,
                             FDRcontrol = c("qvalue", "bonferroni", "none"),
                             apply.qval = identical(FDRcontrol[1], 'qvalue'),
                             selection_fdr = 0.05,
                             adjust_by = c("individual", "all", "none"),
                             lambda = NULL, pi0.method = "smoother",
                             alpha = 0.01,
                             verbose = TRUE,
                             cl = NULL, chunk.size = NULL,
                             save.list = FALSE,
                             save.path = "/path/to/location") {
  data = as.data.frame(data)
  n1 = dim(data)[1]
  p = dim(data)[2]
  measure <- measure[1]
  adjust_by <- adjust_by[1]
  if (!missing(FDRcontrol)) {
    FDRcontrol <- FDRcontrol[1]
    if (identical(FDRcontrol, 'qvalue')) {
      apply.qval <- TRUE
    }
    else if (identical(FDRcontrol, 'none')) {
      adjust_by <- "none"
    }
    else {
      stopifnot(FDRcontrol %in% p.adjust.methods)
    }
  }

  if (is.null(conditional.vars)) {
    if (identical(measure, "partial_corr")) {
      if (verbose) {
        cat("       Using '' measure = 'correlation' '' because 'conditional.vars' is NULL")
      }
      measure <- 'correlation'
    }
  }
  else {
    cond.p = dim(conditional.vars)[2]
  }
  cn.data = colnames(data)
  cov.pool = as.data.frame(cov.pool)
  n2 = dim(cov.pool)[1]
  q = dim(cov.pool)[2]
  cn.cov.pool = colnames(cov.pool)
  if (n1 != n2) {
    stop(paste0("cov.pool and data are not the same size! data has ",
                n1, " samples and cov.pool has ", n2, " samples!"))
  }
  non.na.vals = apply(data, 2, function(x) length(stats::na.omit(x)))
  if (any(non.na.vals <= 2)) {
    stop(paste0("some columns of \"trios\" contain <= 2 non-NA values: The column(s) is/are ",
                paste0(which(non.na.vals <= 2))))
  }
  colnames(data) = make.unique(colnames(data))
  sample.sizes = apply(data, 2, function(x) length(na.omit(x)))

  switch(measure,
         correlation = {
           if (verbose) {
           message("selected measure = correlation...")
           message(paste0("Calculating correlation matrix of size ",
                          ncol(data), " x ", ncol(cov.pool), " using ", ceiling(dim(data)[2]/blocksize),
                          " blocks"))
           }
           r.mat = propagate::bigcor(data, cov.pool, verbose = verbose,
                                      use = "pairwise.complete.obs", size = blocksize)
           r.mat = t(as.matrix(r.mat[, ]))
           p.mat = p.from.cor.mod (r.mat, n = sample.sizes)
         },
         partial_corr = {
           if (verbose) {
           message("selected measure = partial correlation")
           }
           if (is.null(conditional.vars)) {
             stop("conditional.vars is empty!")
           }
           if (n1 < (cond.p)) {
             stop(paste0("Cannot compute partial correlation test!: samples size = ",
                         n1, "< ", ((cond.p) + 3), " Consider using measure = 'correlation' instead "))
           }
           contains.nas = c(any(is.na(data)), any(is.na(cov.pool)),
                            any(is.na(conditional.vars)))
           if (any(contains.nas)) {
             object.nm = c("data", "cov.pool", "conditional.vars")
             which.contain.na = object.nm[which(contains.nas)]
             if (length(which.contain.na) > 1) {
               stop(paste0(paste(which.contain.na, collapse = " and "),
                           " contains NAs!...stopping"))
             }
             if (verbose) {
             message(paste0(which.contain.na, " contains NAs...omitting rows with missing values..."))
             }
             if (which(contains.nas) == 1) {
               omit = which(apply(data, 1, function(x) any(is.na(x))))
             }
             else if (which(contains.nas) == 2) {
               omit = which(apply(cov.pool, 1, function(x) any(is.na(x))))
             }
             else {
               omit = which(apply(conditional.vars, 1, function(x) any(is.na(x))))
             }
             if ((n1 - length(omit)) < (cond.p + 3)) {
               stop("After omitting rows with missing values, sample size is too small to compute partial correlation test!...stopping")
             }
             else {
               if (verbose) {
               message(paste0("Identified ", length(omit), " rows corresponding to missing values. Omitting rows from data, cov.pool, and conditional vars..."))
               }
             }
             data.omit = data[-omit, ]
             cov.pool.omit = cov.pool[-omit, ]
             conditional.vars.omit = conditional.vars[-omit, ]
             if (verbose) {
             message("Computing pairwise partial correlations ...")
             }
             r.mat = compute.pairwise.pcors(data = data.omit,
                                            confs = cov.pool.omit,
                                            cond.vars = conditional.vars.omit,
                                            cl = cl, chunk.size = chunk.size)
           }
           else {
             if (verbose) {
             message("Computing pairwise partial correlations ...")
             }
             r.mat = compute.pairwise.pcors(data = data,
                                            confs = cov.pool,
                                            cond.vars = conditional.vars,
                                            cl = cl, chunk.size = chunk.size)
           }
           p.mat = p.from.parcor.mod (r.mat,
                                      n = sample.sizes,
                                      S = cond.p)
         },
         stop("Argument 'measure' not specified"))

  if (apply.qval == TRUE) {
    if (verbose) {
    message(paste0("Applying qvalue correction to control the FDR at ",
                   selection_fdr))
    }
    switch(adjust_by,
           individual = {
             qsig.mat = matteLapply(1:NCOL(p.mat),
                                    FUN = function(j) {
                                      adjust.q (p.mat[,j],
                                                fdr = selection_fdr,
                                                lambda = lambda,
                                                pi0.meth = pi0.method)
                                    },
                                    cl = cl, chunk.size = chunk.size)
      sig.mat = sapply(qsig.mat, function(x) x$significant)
      q.mat = sapply(qsig.mat, function(x) x$qvalue)
    },
    all = {
      qsig.mat = adjust.q(as.vector(as.matrix(p.mat)),
                          fdr = selection_fdr, lambda = lambda, pi0.meth = pi0.method)
      sig.mat = as.data.frame(matrix(qsig.mat$significant,
                                     nrow = nrow(p.mat), ncol = ncol(p.mat), byrow = F))
      q.mat = as.data.frame(matrix(qsig.mat$qval, nrow = nrow(p.mat),
                                   ncol = ncol(p.mat), byrow = F))
    },
    none = {
      qsig.mat = q.mat = matrix(NA, nrow = nrow(p.mat),
                                ncol = ncol(p.mat))
      sig.mat = apply(p.mat, 2, function(x) x < alpha)
    }, stop(paste0("Input to argument 'adjust_by ' = ", adjust_by,
                   " is not recognized. use one of 'individual', ' all ', or 'none' ")))
    p.adj.mat = matrix(NA, nrow = nrow(p.mat), ncol = ncol(p.mat))
  }
  else {
    if (verbose) {
    message(paste0("Applying", FDRcontrol[1]," correct with threshold ",
                   alpha))
    }
    switch(adjust_by,
           individual = {
             p.adj.mat = matteApply(p.mat, MARGIN = 2,
                                    FUN = stats::p.adjust, method = FDRcontrol,
                                    cl = cl, chunk.size = chunk.size)
             sig.mat = apply(p.adj.mat, MARGIN = 2,
                             FUN = function(x) x <= alpha)
             q.mat = matrix(NA, nrow = nrow(p.mat), ncol = ncol(p.mat))
           },
           all = {
             p.adj.mat = as.data.frame(matrix(stats::p.adjust(as.vector(as.matrix(p.mat)),
                                                              method = FDRcontrol), nrow = nrow(p.mat), ncol = ncol(p.mat),
                                              byrow = F))
             sig.mat = apply(p.adj.mat, 2, function(x) x <= alpha)
             q.mat = matrix(NA, nrow = nrow(p.mat), ncol = ncol(p.mat))
           },
           none = {
             p.adj.mat = q.mat = matrix(NA, nrow = nrow(p.mat),
                                        ncol = ncol(p.mat))
             sig.mat = apply(p.mat, 2, function(x) x < alpha)
           },
           stop(paste0("Input to argument 'adjust_by ' = ", adjust_by,
                       " is not recognized. use one of 'individual', ' all ', or 'none' ")))
  }

  if (verbose) {
  message("Selecting significant covariates...")
  }

  sig.asso.covs = lapply(1:dim(sig.mat)[2], function(x, y) {
    which(sig.mat[, x])
  }, y = sig.mat)
  colnames(sig.mat) = colnames(q.mat) = colnames(r.mat) = colnames(p.mat) = colnames(p.adj.mat) = cn.data
  row.names(sig.mat) = row.names(q.mat) = row.names(p.mat) = row.names(p.adj.mat) = cn.cov.pool
  out.list = list(sig.asso.covs = sig.asso.covs, pvalues = p.mat,
                  qvalues = q.mat, cors = r.mat, sig = sig.mat, adj.p = p.adj.mat)
  if (save.list == TRUE) {
    if (verbose) {
    message(paste0("saving output at ", save.list, ".RData"))
    }
    save(out.list, file = paste0(save.path, ".RData"))
  }
  return(out.list)
}

# Only adding arguments 'cl' and 'chunk.size' for parallel computing
compute.pairwise.pcors <-
  function (data, confs, cond.vars, cl = NULL, chunk.size = NULL) {
    p = ncol(data)
    u = ncol(confs)
    z = ncol(cond.vars)
    iterable = expand.grid(Uvars = c(1:u), Pvars = c(1:p))

    pcor.mat = #catch.conditions({
      matteApply(iterable, 1,
                 FUN = function(x, p, u, z) {
                   ppcor::pcor(cbind(u[, x[1]], p[, x[2]], z))$estimate[1, 2]
                 },
                 p = data, u = confs, z = cond.vars,
                 cl = cl, chunk.size = chunk.size)
    #})$value

    pcor.final = as.data.frame(matrix(pcor.mat, nrow = p,
                                      ncol = u, byrow = T))

    return(t(pcor.final))


  error = function(e) {
    print(e)
  }
  warning = function(w) {
    mess = names(w)
    pat = "The inverse of variance-covariance matrix"
    if (grepl(pattern = pat, mess)) {
      return(t(pcor.final))
    }
    else {
      print(w)
    }
  }
}

p.from.parcor.mod <- function (r, n, S) {
  z = (sqrt(n - S - 3)/2) * log((1 + r)/(1 - r))
  p = 2 * (1 - pnorm(abs(z), lower.tail = T))
  return(p) # returning just p-values (not r as in the original function)
}

p.from.cor.mod <-
  function (r, n) {
  t = r * sqrt((n - 2)/(1 - r^2))
  p = 2 * (1 - stats::pt(q = abs(t), df = n - 2))
  return(p) # returning just p-values (not r as in the original MRGN function)
}
