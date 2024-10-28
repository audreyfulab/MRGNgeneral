
#' Build the input graph for MRGN inference
#'
#' Build an undirected adjacency matrix for a genomic network using a dataset
#' where rows are individuals and columns are nodes including \code{n_v} variants
#' (\code{V}-nodes), \code{n_t} phenotypes (\code{T}-nodes), and \code{n_c}
#' candidate confounding variables (\code{C}-nodes).
#'
#' @param conf.sets an optional object of class \code{conf.sets} as returned by
#' \link{get.conf.sets}.
#'
#' @details
#' This function implements \code{Step 2} of the MRGN general algorithm.
#'
#' @return an object of class \code{adjacency.matrix}, i.e. a binary square matrix
#' with \code{n_v + n_t + n_q} columns, with only zeros on the main diagonal.
#'
#' The number \code{n_q} of intermediate variables and common children is taken
#' from the argument \code{conf.sets} as \code{n_q = length(conf.sets$WZindices)}.
#' If \code{conf.sets = NULL}, then \code{n_q = 0}.
#'
#' @export get.initial.skeleton
#'
#' @examples
#' #Raw data (networkA11) & confounding variable selection result (confsetsA11)
#' library(MRGNgeneral)
#' data(networkA11)
#' data(confsetsA11)
#'
#' # Generate the input graph
#' Adjacency0 <- get.initial.skeleton (data = networkA11$data,
#'                                     n_v = networkA11$dims$n_v,
#'                                     n_t = networkA11$dims$n_t,
#'                                     threshold_v = 0.2,
#'                                     threshold_m = 0.05,
#'                                     conf.sets = confsetsA11)
#'
#' Adj0.mtx <- Adjacency0
#' class (Adj0.mtx) <- "matrix"
#' dim (Adj0.mtx)
#' # 333 333
#' # Nodes here are V and T nodes, plus the inferred Q nodes.
#' # They are a subset of all the nodes.
#' 
#' # The true adjacency matrix does not include the independent variables
#' dim (networkA11$adjacency)
#' # 500 500
#' 
#' # In order to extract the corresponding nodes from the true adjacency matrix,
#' # need to create the adj mtx for all nodes,
#' # and then extract the nodes in query mtx
#' n_all <- sum (unlist(networkA11$dims))
#' true_adj_all <- matrix(0, nrow = n_all, ncol = n_all)
#' n_vtwzu <- nrow(networkA11$adjacency)
#' true_adj_all[1:n_vtwzu, 1:n_vtwzu] <- networkA11$adjacency
#' dim (true_adj_all)
#' # 600 600
#' 
#' # The nodes in Adjacency0
#' n_v <- networkA11$dims$n_v
#' n_t <- networkA11$dims$n_t
#' index_all <- c(1:(n_v+n_t), confsetsA11$WZindices)
#' # Extract these nodes from the entire true adjacency matrix
#' true_adj <- true_adj_all[index_all, index_all]
#' 
#' # Recall and precision of the skeleton
#' Recall_Precision <- RecallPrecision(as (true_adj, 'graphNEL'),
#'                                     as (Adj0.mtx, 'graphNEL'),
#'                                     GV = networkA11$dims$n_v,
#'                                     includeGV = TRUE,
#'                                     edge.presence = 1.0,
#'                                     edge.direction = 0.5)
#'

# Initialize the input graph (adjacency matrix) for a call to MRGN
get.initial.skeleton <- function (data,
                                  n_v, n_t,
                                  threshold_v = 0.2,
                                  threshold_m = 0.05,
                                  conf.sets = NULL) {
  stopifnot(is.matrix(data) | is.data.frame(data))
  stopifnot(is.numeric(n_v), is.numeric(n_t),
            is.numeric(threshold_v), is.numeric(threshold_m))

  ### Index of sets
  n_vt <- n_t + n_v
  indexV <- 1:n_v
  indexT <- (n_v + 1):n_vt

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
    PTTcorGivenAllV <- matrix(0, nrow = n_t, ncol = n_t)

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

  ### Number of selected W,Z nodes
  n_q <- if (!is.null(conf.sets)) length(conf.sets$WZindices) else 0

  ### Initialize result to matrix of zeros
  Adj0 <- matrix(0, nrow = n_vt + n_q, ncol = n_vt + n_q)

  ### V-T edges
  Adj0[indexV, indexT] <- (MVTcor >= threshold_v) + 0

  ### T-T edges
  Adj0[indexT, indexT] = (TTMPcor >= threshold_m) + 0

  ### Initialize edges involving common children or intermediate variables, if any
  if (n_q > 0) {
    ## Set of indices for Q-nodes
    indexQ <- (n_vt + 1):(n_vt + n_q)

    ## Find T-Q edges
    T_WZ <- NULL
    if (!is.null(conf.sets)) {
      #* From Q-node selection results in 'conf.sets'
      T_WZ <- if (length(conf.sets$UWZconfounders) == n_t)
        t(sapply(conf.sets$UWZconfounders,
                 FUN = T_WZconf2adj, n_q = n_q,
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
  diag(Adj0) <- rep(0, n_vt + n_q)

  return(structure(Adj0, class = "adjacency.matrix"))

}

# Turn a vector of selected U,W,Z-nodes into a binary vector: row of an adjacency matrix for T-Q edges
T_WZconf2adj <- function (x,         # indices of selected confounding variables (U,W,Z)
                          WZindices, # indices of all W,Z-nodes
                          n_q = length(WZindices)) {
  TiRow <- numeric(n_q)
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
