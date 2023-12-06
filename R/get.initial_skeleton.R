
#' @export get.initial_skeleton

# I itialize the input graph (adjacency matrix) for a call to MRGN
get.initial_skeleton <- function (data,
                                  p, q,
                                  threshold_v = 0.2,
                                  threshold_m = 0.05,
                                  conf.sets = NULL) {

  ### Index of sets
  number.of.VT <- p + q
  indexV <- 1:q
  indexT <- (q + 1):number.of.VT

  ### Partial and marginal correlations
  MVTcor <- PTTcorGivenAllV <- NULL
  if (!is.null(conf.sets)) {
    ## V-T correlations
    MVTcor <- conf.sets$raw$Vraw$cors

    ## T-T partial correlations given all V-nodes
    PTTcorGivenAllV <- conf.sets$raw$Traw$cors
  }

  if (is.null(MVTcor))
    ## V-T correlations
    MVTcor <- cor(x = data[, indexV], y = data[, indexT])

  if (is.null(PTTcorGivenAllV))
    ## Partial correlations given all V-nodes NOT USED (too expensive to compute)
    PTTcorGivenAllV <- matrix(0, nrow = p, ncol = p)

  ## Partial T-T correlations given V,T,C
  Pcor <- partial.cor(data)
  PTTcorGivenVTC <- Pcor[indexT, indexT]

  ## Marginal T-T correlations
  MTTcor <- cor(x = data[, indexT])

  ### Absolute correlations
  MVTcor <- abs(MVTcor)
  TTMPcor <- pmax(abs(PTTcorGivenAllV),
                  abs(PTTcorGivenVTC),
                  abs(MTTcor))

  ### Number of W,Z nodes
  r <- if (!is.null(conf.sets)) length(conf.sets$WZindices) else 0

  ### Initialize result to matrix of zeros
  Adj0 <- matrix(0, nrow = number.of.VT + r, ncol = number.of.VT + r)

  ### V-T edges
  Adj0[indexV, indexT] <- (MVTcor >= threshold_v) + 0

  ### T-T edges
  Adj0[indexT, indexT] = (TTMPcor >= threshold_m) + 0

  ### Initialize edges involving common children or intermediate variables, if any
  if (r > 0) {
    ## Set of indices for Q-nodes
    indexQ <- (p + q + 1):(number.of.VT + r)

    ## Find T-Q edges
    T_WZ <- NULL
    if (!is.null(conf.sets)) {
      #* From Q-node selection results in 'conf.sets'
      T_WZ <- if (length(conf.sets$UWZconfounders) == p)
        t(sapply(conf.sets$UWZconfounders,
                 FUN = T_WZconf2adj, r = r,
                 WZindices = conf.sets$WZindices))
    }

    #* From partial T-Q correlations given all V,T,Q, if no info in 'conf.sets'
    if (is.null(T_WZ)) {
      T_WZ <- (abs (Pcor[indexT, indexQ]) > threshold_m) + 0
    }

    # Fill the input graph
    Adj0[indexT, indexQ] <- T_WZ
    Adj0[indexQ, indexT] <- t(T_WZ)
  }

  ### Ensure 'Adj0' has only zeros on its main diagonal
  diag(Adj0) <- rep(0, number.of.VT + r)

  return(structure(Adj0, class = "adjacency.matrix"))

}


# Turn a vector of selected U,W,Z-nodes into a binary vector:
# row of an T-Q adjacency matrix
T_WZconf2adj <- function (x,         # indices of selected confounding variables (U,W,Z)
                          WZindices, # indices of all W,Z-nodes
                          r = length(WZindices)) {
  TiRow <- numeric(r)
  if (length(x)) {
    # Pick 'x' elements also present in 'WZindices'
    wz <- x %in% WZindices
    y <- if(sum(wz)) x[wz]

    if (length(y)) {
      # Take the position of each y element in WZindices
      positions <- sapply(y, FUN = grep, x = WZindices)
      TiRow[positions] <- 1
    }
  }
  return (TiRow)
}
