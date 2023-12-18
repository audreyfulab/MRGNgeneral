#'
#' Performance of confounding variable selection
#'
#' Compute recall and precision for confounding variables
#' selected using \link{get.conf.sets}.
#'
#' @param conf.sets an object of class \code{conf.sets} as returned by \link{get.conf.sets}.
#'
#' @param adjacency numeric, adjacency matrix of the true network, i.e., a binary square matrix
#' with zeros on the main diagonal. The columns and rows of \code{adjacency} must be such
#' that the first \code{n_v} correspond to genetic variants (\code{V}-nodes),
#' the next \code{n_t} to target phenotypes/expressions values (\code{T}-nodes),
#' the next \code{n_w} to intermediate variables (\code{W}-nodes),
#' the next \code{n_z} to common children (\code{Z}-nodes), and
#' the last \code{n_u} to confounders (\code{U}-nodes).
#'
#' @param n_v,n_t,n_w,n_z,n_u integers, the numbers
#' \code{n_v} of \code{V}-nodes,
#' \code{n_t} of \code{T}-nodes,
#' \code{n_w} of \code{W}-nodes,
#' \code{n_z} of \code{Z}-node, and
#' \code{n_u} of \code{U}-nodes in the true network.
#' Note that the default values inferred from \code{conf.sets} assume
#' \code{n_w = n_z = 0}, i.e. there is no common child or intermediate variable
#' in the true network.
#'
#' @export assess.conf.selection
#'
#' @details
#' The function evaluates the performance of \link{get.conf.sets} on a single
#' genomic dataset, given the adjacency matrix of the true network.
#'
#' @return a named list with the following elements:
#' \describe{
#' \item{\code{recall}}{a 6-vector giving the recalls for the pool of
#' \code{W,Z}-nodes, and the pool of \code{U}-nodes, and the average recalls
#' for \code{V}-nodes, \code{T}-nodes, \code{Q}-nodes (mixture of \code{W} and
#' \code{Z}-nodes), and \code{U}-nodes. The raw recalls the averages are computed
#' from are in the returned element \code{raw}.}
#' \item{\code{precision}}{a 6-vector giving the precisions for the pool of
#' \code{W,Z}-nodes, and the pool of \code{U}-nodes, and the average precisions
#' for \code{V}-nodes, \code{T}-nodes, \code{Q}-nodes, and \code{U}-nodes. The
#' raw precisions the averages are computed from are in the returned element \code{raw}.}
#' \item{\code{raw}}{a list with elements \code{V}, \code{T}, \code{Q}, and
#' \code{U}, each element being a named list of two vectors (
#' \code{recall}: \code{n_t}-vector giving the recall of a particular group of
#' nodes (e.g. \code{V}, \code{T}, \code{Q}, \code{U}) for each focal \code{T}-node,
#' \code{precision}: \code{n_t}-vector giving the precisions of a particular
#' group of nodes for each focal \code{T}-node).}
#' }
#'
#' @examples
#' ## Load the 'conf.sets' object 'confsetsA11' and the raw data 'networkA11' it was obtained from
#' library(MRGNgeneral)
#' data(confsetsA11)
#' data(networkA11)
#'
#' ## Recall and precision of the selection procedure
#' Perf <- assess.conf.selection (confsetsA11,
#'                                adjacency = networkA11$adjacency,
#'                                n_v = networkA11$dims$n_v,
#'                                n_t = networkA11$dims$n_t,
#'                                n_w = networkA11$dims$n_w,
#'                                n_z = networkA11$dims$n_z,
#'                                n_u = networkA11$dims$n_u)
#'
#' Perf$recall
#' Perf$precision

assess.conf.selection <- function (conf.sets,
                                   adjacency,
                                   n_v = length(conf.sets$WZconfounders),
                                   n_t = length(conf.sets$confounders),
                                   n_w = 0, n_z = 0,
                                   n_u = NCOL(adjacency) - n_v - n_t - n_w - n_z) {
  # Check arguments
  stopifnot(is.conf.sets (conf.sets))
  n_wz <- n_w + n_z

  #### True parent sets
  ## V-nodes
  TrueVparentset <- get.true.variant.set(Adj = adjacency,
                                         n_t = n_t,
                                         n_v = n_v)
  if (!length(TrueVparentset)) {
    TrueVparentset <- vector(mode = "list", length = n_t)
  }

  ## T-nodes
  TrueTparentset <- get.true.parent.genes(Adj = adjacency,
                                          n_t = n_t,
                                          n_v = n_v)
  if (!length(TrueTparentset)) {
    TrueTparentset <- vector(mode = "list", length = n_t)
  }

  ## Q-nodes
  TrueQset <- get.true.conf.set(Adj = adjacency,
                                n_v = n_v,
                                n_t = n_t,
                                n_wz = n_wz,
                                Upool = FALSE,
                                offset = n_v + n_t)
  if (!length(TrueQset)) {
    TrueQset <- vector(mode = "list", length = n_t)
  }

  ## U-nodes
  Trueconfset <- get.true.conf.set(Adj = adjacency,
                                   n_v = n_v,
                                   n_t = n_t,
                                   n_wz = n_wz,
                                   n_u = n_u,
                                   Upool = TRUE,
                                   offset = n_v + n_t + n_wz)
  if (!length(Trueconfset)) {
    Trueconfset <- vector(mode = "list", length = n_t)
  }

  ## Pools of selected W,Z-nodes, and selected U-nodes
  Qindices <- conf.sets$WZindices
  Uindices <- setdiff(conf.sets$UWZindices, Qindices)

  ### Compute recall and precision for the W,Z pool and U pool
  ConfSelRecall <- c(Qpool = MRGNgeneral:::RecallUselection (Uset = list(Qindices),
                                                       RefSet = list((n_t+n_v+1):(n_t+n_v+n_wz))),
                     Upool = MRGNgeneral:::RecallUselection (Uset = list(Uindices),
                                                      RefSet = list((n_t+n_v+n_wz+1):(n_t+n_v+n_wz+n_u))))
  ConfSelPrecision <- c(Qpool = MRGNgeneral:::PrecisionUselection (Uset = list(Qindices),
                                                             RefSet = list((n_t+n_v+1):(n_t+n_v+n_wz))),
                        Upool = MRGNgeneral:::PrecisionUselection (Uset = list(Uindices),
                                                            RefSet = list((n_t+n_v+n_wz+1):(n_t+n_v+n_wz+n_u))))


  ### Compute recall and precision for each set of nodes
  ## V-nodes
  Vrecalls = RecallUselection (Uset = conf.sets$Vconfounders, RefSet = TrueVparentset)
  Vprecisions = PrecisionUselection (Uset = conf.sets$Vconfounders, RefSet = TrueVparentset)

  ## T-nodes
  Trecalls = RecallUselection (Uset = conf.sets$Tconfounders, RefSet = TrueTparentset)
  Tprecisions = PrecisionUselection (Uset = conf.sets$Tconfounders, RefSet = TrueTparentset)

  ## Q-nodes
  conf.sets$Qconfounders <- lapply(conf.sets$UWZconfounders,
                                   FUN = function(x) setdiff(x, Uindices))

  Qrecalls = RecallUselection (Uset = conf.sets$Qconfounders, RefSet = TrueQset)
  Qprecisions = PrecisionUselection (Uset = conf.sets$Qconfounders, RefSet = TrueQset)

  ## U-nodes
  Urecalls = RecallUselection (Uset = conf.sets$Uconfounders, RefSet = Trueconfset)
  Uprecisions = PrecisionUselection (Uset = conf.sets$Uconfounders, RefSet = Trueconfset)

  # Mean recall and precision
  ConfSelRecall = c(ConfSelRecall,
                    Vnodes = mean(Vrecalls), Tnodes = mean(Trecalls),
                    Qnodes = mean(Qrecalls), Unodes = mean(Urecalls))
  ConfSelPrecision = c(ConfSelPrecision,
                       Vnodes = mean(Vprecisions), Tnodes = mean(Tprecisions),
                       Qnodes = mean(Qprecisions), Unodes = mean(Uprecisions))

  return(list(recall = ConfSelRecall,
              precision = ConfSelPrecision,
              raw = list(V = list(recall = Vrecalls,
                                  precision = Vprecisions),
                         T = list(recall = Trecalls,
                                  precision = Tprecisions),
                         Q = list(recall = Qrecalls,
                                  precision = Qprecisions),
                         U = list(recall = Urecalls,
                                  precision = Uprecisions))))

}

# A function to compute recall for confounder selection
RecallUselection <- function (Uset, RefSet, na.rm = TRUE) {

  n_t <- length(RefSet)
  stopifnot(n_t > 0)
  stopifnot(length(Uset) == n_t)

  out <- sapply (1:n_t, FUN = RecallUselectionChild, Uset = Uset, RefSet = RefSet)

  names(Uset) -> names(out)

  if (na.rm) {
    out <- na.omit(out)

    out <- out[1:length(out)]
  }

  return(out)

}

# Child of RecallUselection
RecallUselectionChild <- function(j, Uset, RefSet) {
  if (length(RefSet[[j]]) & length(Uset[[j]])) {
    return(mean(RefSet[[j]] %in% Uset[[j]], na.rm = TRUE))
  }
  else if (length(RefSet[[j]])) {
    return(0)
  }
  else {
    return(NA)
  }
}

# A function to compute precision for confounder selection
PrecisionUselection <- function (Uset, RefSet, na.rm = TRUE) {
  n_t <- length(RefSet)
  stopifnot(n_t > 0)
  stopifnot(length(Uset) == n_t)

  out <- sapply (1:n_t, FUN = PrecisionUselectionChild, Uset = Uset, RefSet = RefSet)

  names(Uset) -> names(out)

  if (na.rm) {
    out <- na.omit(out)

    out <- out[1:length(out)]
  }

  return(out)
}

# Child of PrecisionUselection
PrecisionUselectionChild <- function(j, Uset, RefSet) {
  if (length(Uset [[j]]) & length(RefSet[[j]])) {
    return(mean(Uset[[j]] %in% RefSet[[j]], na.rm = TRUE))
  }
  else if (!length(RefSet[[j]]) & !length(RefSet[[j]])) {
    return(1)
  }
  else if (length(Uset[[j]]) & !length(RefSet[[j]])) {
    return(0)
  }
  else { # if (!length(Uset[[j]]) & length(RefSet[[j]]))
    return(NA)
  }
}

