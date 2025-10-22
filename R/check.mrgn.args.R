# An internal function to generate the code to check MRGN arguments
# Also saves inputs and format some arguments for the MRGN algorithm
check.mrgn.args <- function() {
  expression({
    # ============================================================================
    # Saving inputs & Checking arguments
    # ----------------------------------------------------------------------------
    # Data
    data <- as.data.frame(data)
    n <- NROW(data)
    m <- NCOL(data)

    # Check validity and consistency of block sizes
    stopifnot(is.numeric(n_v), is.numeric(n_t), is.numeric(n_q), is.numeric(n_u))
    stopifnot(n_v >= 0, n_v <= m, n_t >= 0, n_t <= m, n_q >= 0, n_q <= m, n_u >= 0, n_u < m, n_t + n_v + n_u + n_q <= m)
    if (m > n_t + n_v + n_q + n_u)
      data <- data[, 1:(n_t + n_v + n_q + n_u)]
    if (scale.data) {
      data[,-c(1:n_v)] <- scale(data[,-c(1:n_v)], center = TRUE, scale = TRUE)
    }
    m <- n_t + n_v + n_q + n_u
    nb.nodes <- n_t + n_v + n_q

    # Save original labels (if any) or label columns
    Labels <- colnames(data)

    # Labels for Q-nodes
    if(is.null(Qlabels) & n_q > 0) {
      Qlabels <- (n_t + n_v + 1):(n_t + n_v + n_q)
    }

    if (is.null(Labels)) {
      Vlabels <- if (n_v > 0) paste0('V', 1:n_v)
      Tlabels <- if (n_t > 0)  paste0('T', 1:n_t)
      WZlabels <- if (n_q > 0)  paste0('Q', 1:n_q)
      Clabels <- if (n_u > 0) paste0('C', 1:n_u)
      Labels <- c(Vlabels, Tlabels, WZlabels, Clabels)
      colnames(data) <- Labels

      # Labels for Q-nodes
      if (n_q > 0) {
        if(is.null(Qlabels)) {
          Qlabels <- (n_t + n_v + 1):(n_t + n_v + n_q)
        }
        else if (!is.numeric(Qlabels)) {
          stop("'Qlabels' must be numeric vector if 'data' columns are unnamed")
        }
      }

    }
    else {
      Labels <- make.unique(Labels)
      colnames(data) <- Labels
      if(is.null(Qlabels)) {
        Qlabels <- (n_t + n_v + 1):(n_t + n_v + n_q)
      }
      else if (!is.numeric(Qlabels)) {
        if (is.character(Qlabels)) {
          # Change to numeric
          Qlabels <- sapply(Qlabels,
                            FUN = grep,
                            x = Labels)
          if (any(Qlabels <= (n_t + n_v))) {
            stop("'Qlabels' must index columns other than V and T-nodes (first n_t+n_v columns in data)")
          }
        }
        else {
          stop("'Qlabels' must be numeric or character vector")
        }
      }
    }

    # Check the consistency of the argument 'confounders' with data and its blocks
    if (!missing(confounders) && !is.null(confounders)) {
      if (!is.list(confounders))
        stop ("The argument 'confounders' must be a list")
      if (!all(sapply(confounders, FUN = function(x) {is.numeric(x) || is.null(x)})))
        stop ("Each element of 'confounders' must be a numeric vector")
      n.con <- length(confounders)
      if (n.con == 0) {
        if (verbose & n_u > 0)
          warning(paste0("      # ", n_u, " confounders supplied, but none is associated with a node. \n"))
        confounders <- NULL
      }
      else if (n.con == 1) {
        stopifnot(all(!is.na(confounders[[1]])), all(floor(confounders[[1]]) == confounders[[1]]))
        # if (any(confounders[[1]] > n_v & confounders[[1]] <= nb.nodes) | any(confounders[[1]] > m))
        if (any(confounders[[1]] > m))
          stop ("Column index out ot bound in 'confounders'")
        confounders[1:n_t] <- confounders[1]
      }
      else if (n.con == n_t) {
        confs.range <- unlist(confounders, recursive = TRUE, use.names = FALSE)
        if (length(confs.range)) {
        stopifnot(all(!is.na(confs.range)), all(floor(confs.range) == confs.range))
        confs.range <- range(confs.range)
        # if ((confs.range[1] > n_v & confs.range[1] <= nb.nodes) | confs.range[2] > m)
        if(confs.range[2] > m)
          stop ("Column index out of bound in 'confounders'")
        }
        rm(confs.range)
      }
      else
        stop (paste0("The argument 'confounders' must be a list of length one or ",
                     n_t, ". \n"))
    }
    else {
      if (verbose & n_u > 0)
        warning(paste0("      # ", n_u, " confounders supplied, but none is associated with a node. \n"))
      confounders <- NULL
    }
    adjacency <- as.matrix(adjacency)

    # Check the supplied adjacency matrix
    stopifnot(is.adjacency.matrix(adjacency))
    if (NCOL(adjacency) %in% c(nb.nodes, m)) {
      if (NCOL(adjacency) == m)
        adjacency <- adjacency[1:nb.nodes, 1:nb.nodes]
    }
    else
      stop("The size of argument 'adjacency' is not consistent with arguments 'data', 'n_t', 'n_v', 'n_q', and 'n_u'.")

    # Orient all V -- T edges as V --> T
    adjacency[(1+n_v):nb.nodes, 1:n_v] <- 0
    dimnames(adjacency) <- list(Labels[1:nb.nodes], Labels[1:nb.nodes])

    # Save the input adjacency matrix
    Adj0 <- adjacency

    # Check FDR control arguments if required
    stopifnot(is.character(FDRcontrol))
    FDRcontrol <- FDRcontrol[1]
    stopifnot(FDRcontrol %in% c(p.adjust.methods, "qvalue"))
    stopifnot (is.numeric(alpha))
    alpha <- alpha[1]
    stopifnot (alpha > 0, alpha < .5)
    if (identical(FDRcontrol, "qvalue")) { # qvalue method
      stopifnot (is.numeric(fdr))
      fdr <- fdr[1]
      stopifnot (fdr > 0, fdr < .5)
      stopifnot(all(pi0.meth %in% c("bootstrap", "smoother")))
      stopifnot (lambda.step > 0, lambda.step < .1)
    }
    pi0.meth <- pi0.meth[1]

    stopifnot(is.logical(stringent))
    stringent <- stringent[1]

    # Set default maximum number of iteration
    if (is.null(maxiter))
      maxiter <- choose(nb.nodes, 3)

    Tlabels <- (n_v+1):(n_t+n_v)
  })
}
