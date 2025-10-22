
#' @title Calculate recall and precision for an inferred graph given the true graph
#'
#' @description
#' This function counts the number of true and false positives, and calculates
#' recall (i.e., power) and precision (i.e., 1-FDR), which are defined as follows:
#'
#' Recall = (# edges correctly identified in inferred graph) / (# edges in true graph).
#'
#' Precision = (# edges correctly identified in inferred graph) / (# edges in inferred graph).
#'
#' @export RecallPrecision
#'
#' @param g1 First graph object, the true graph.
#'
#' @param g2 Second graph object, the inferred graph.
#'
#' @param GV The number of genetic variants (SNPs/indels/CNV/eQTLs) in the input
#' data matrix. For example, if the data has one variant, which is in the first
#' column, then GV = 1. If there are two variants, which are in the first and
#' second Columns, then GV = 2. If there are no variants, then GV = 0.
#'
#' @param includeGV If \code{TRUE}, include edges involving genetic variants (GV)
#' when calculating recall and precision. If FALSE, exclude edges involving
#' genetic variants (GV) when calculating recall and precision.
#'
#' @param edge.presence The weight for an edge being present.
#'
#' @param edge.direction The weight for the edge direction.
#'
#' @details
#'  We consider it more important to be able to identify the presence of an edge
#'  than to also get the direct correct. Therefore, we assign 1 as the default
#'  to an edge with the correct direction and 0.5 to an edge with the wrong
#'  direction or no direction \insertCite{badsha2019learning,badsha2021mrpc}{MRGNgeneral}.
#'
#'  Column names must match between the two graphs (the order may differ).
#'
#' @return A \code{list} object that contains the following:
#' \itemize{
#' \item \code{Matrix}: Results store for TP and FP
#' \item \code{TP}: Total found edges in the inferred graph and edge exists in the true graph.
#' \item \code{FP}: Total found edges in the inferred graph but no edge exists in the true graph.
#' \item \code{NTE}: Total number of edges in the true graph.
#' \item \code{NIE}: Total number of edges in the inferred graph.
#' \item \code{Recall}: Power, or sensitivity measures how many edges from the true graph a method can recover.
#' \item \code{Precision}: Measures how many correct edges are recovered in the inferred graph.
#' }
#'
#'
#' @references
#' \insertAllCited{}
#' @author Md Bahadur Badsha (mbbadshar@gmail.com)
#'
#\references{
#  1. Badsha MB and Fu AQ (2019). Learning causal biological networks with the principle of Mendelian randomization. Frontiers in Genetics, 10:460.
#
#  2. Badsha MB, Martin EA and Fu AQ (2021). MRPC: An R package for inference of causal graphs. Frontiers in Genetics, 10:651812.
#}

#' @seealso
#' \link[MRPC]{aSHD}: adjusted Structural Hamming Distance (aSHD)
#'
#' @examples
#' # True model
#' # True graph (V1 --> T1 --> T2 --> T3)
#' # Where V1 is a genetic variant (GV) and T1, T2, and T3 are phenotypes
#' tarmat_s1 <- matrix(0,
#'                    nrow = 4,
#'                    ncol = 4)
#'
#' colnames(tarmat_s1) <- c("V1", "T1", "T2", "T3")
#' rownames(tarmat_s1) <- colnames(tarmat_s1)
#'
#' # Create an adjacency matrix for the true graph
#' tarmat_s1[1, 2] <- 1
#' tarmat_s1[2, 3] <- 1
#' tarmat_s1[3, 4] <- 1
#'
#' # Graph object of the true graph
#' Truth <- as(tarmat_s1,
#'            "graphNEL")
#'
#' # Inferred graph (V1 --> T1 <-- T2 --> T3)
#' # Where V1 is a genetic variant (GV) and T1, T2, and T3 are phenotypes
#' tarmat_s2 <- matrix(0,
#'                    nrow = 4,
#'                    ncol = 4)
#'
#' colnames(tarmat_s2) <- c("V1", "T1", "T2", "T3")
#' rownames(tarmat_s2) <- colnames(tarmat_s2)
#'
#' # Create an adjacency matrix for the inferred graph
#' tarmat_s2[1, 2] <- 1
#' tarmat_s2[3, 2] <- 1
#' tarmat_s2[3, 4] <- 1
#'
#' # Graph objects for the inferred graph
#' Inferred <- as(tarmat_s2,
#'                "graphNEL")
#'
#' # Recall and Precision
#' Recall_Precision <- RecallPrecision(Truth,
#'                                    Inferred,
#'                                    GV = 1,
#'                                    includeGV = TRUE,
#'                                    edge.presence = 1.0,
#'                                    edge.direction = 0.5)
#'
#' Recall_Precision

RecallPrecision <- function (g1, g2, GV, includeGV,
                             edge.presence = 1,
                             edge.direction = 0.5) {
  if (is(g1, "pcAlgo"))
    g1 <- g1@graph
  if (is(g2, "pcAlgo"))
    g2 <- g2@graph
  if (is(g1, "graphNEL")) {
    m1 <- wgtMatrix(g1, transpose = FALSE) # 'wgtMatrix' is from 'pcalg'
    m1[m1 != 0] <- 1
  }
  if (is.adjacency.matrix(g1)) # ADDED if statement
    m1 <- g1
  if (is(g2, "graphNEL")) {
    m2 <- wgtMatrix(g2, transpose = FALSE) # 'wgtMatrix' is from 'pcalg'
    znm2 <- m2 != 0
    if (any(znm2))  # ADDED if statement
      m2[znm2] <- 1
  }
  if (is.adjacency.matrix(g2)) # ADDED if statement
    m2 <- g2

  if (is.null(colnames(m1)))
    dimnames(m1) <- list(paste0('X', 1:NCOL(m1)), paste0('X', 1:NCOL(m1)))

  if (is.null(colnames(m2)))
    dimnames(m2) <- list(paste0('X', 1:NCOL(m2)), paste0('X', 1:NCOL(m2)))

  if (any(colnames(m1) != colnames(m2))) {
    Order_node <- match(colnames(m1), colnames(m2))
    m2 <- m2[Order_node, Order_node]
  }
  Evaluation_matrix <- matrix(0, nrow = 1, ncol = 2)
  colnames(Evaluation_matrix) <- c("TP", "FP")
  TP <- 1
  FP <- 2
  if (includeGV) {
    m11 <- m1
    for (i in 1:nrow(m11)) {
      for (j in 1:ncol(m11)) {
        if (m11[i, j] == m11[j, i]) {
          m11[i, j] <- 0
        }
      }
    }
    NTE <- length(which(m11 == 1))
    m22 <- m2
    for (i in 1:nrow(m22)) {
      for (j in 1:ncol(m22)) {
        if (m22[i, j] == m22[j, i]) {
          m22[i, j] <- 0
        }
      }
    }
    NIE <- length(which(m22 == 1))
    ind1 <- t(combn(ncol(m1), 2))
    for (i in seq_len(nrow(ind1))) {
      x <- ind1[i, 1]
      y <- ind1[i, 2]
      if ((m1[x, y] == 1 & m1[y, x] != 1) & (m2[x, y] ==
                                             1 & m2[y, x] != 1)) {
        Evaluation_matrix[TP] <- Evaluation_matrix[TP] +
          edge.presence
      }
      if ((m1[y, x] == 1 & m1[x, y] != 1) & (m2[y, x] ==
                                             1 & m2[x, y] != 1)) {
        Evaluation_matrix[TP] <- Evaluation_matrix[TP] +
          edge.presence
      }
      if ((m1[x, y] == 1 & m1[y, x] == 1) & (m2[x, y] ==
                                             1 & m2[y, x] == 1)) {
        Evaluation_matrix[TP] <- Evaluation_matrix[TP] +
          edge.presence
      }
      if ((m1[x, y] == 1 & m1[y, x] != 1) & (m2[x, y] !=
                                             1 & m2[y, x] == 1)) {
        Evaluation_matrix[TP] <- Evaluation_matrix[TP] +
          edge.direction
      }
      if ((m1[x, y] != 1 & m1[y, x] == 1) & (m2[x, y] ==
                                             1 & m2[y, x] != 1)) {
        Evaluation_matrix[TP] <- Evaluation_matrix[TP] +
          edge.direction
      }
      if ((m1[x, y] == 1 & m1[y, x] == 1) & (m2[x, y] ==
                                             1 & m2[y, x] != 1)) {
        Evaluation_matrix[TP] <- Evaluation_matrix[TP] +
          edge.direction
      }
      if ((m1[x, y] == 1 & m1[y, x] == 1) & (m2[x, y] !=
                                             1 & m2[y, x] == 1)) {
        Evaluation_matrix[TP] <- Evaluation_matrix[TP] +
          edge.direction
      }
      if ((m1[x, y] == 1 & m1[y, x] != 1) & (m2[x, y] ==
                                             1 & m2[y, x] == 1)) {
        Evaluation_matrix[TP] <- Evaluation_matrix[TP] +
          edge.direction
      }
      if ((m1[x, y] != 1 & m1[y, x] == 1) & (m2[x, y] ==
                                             1 & m2[y, x] == 1)) {
        Evaluation_matrix[TP] <- Evaluation_matrix[TP] +
          edge.direction
      }
      if ((m1[x, y] == 0 & m1[y, x] == 0) & (m2[x, y] ==
                                             1 || m2[y, x] == 1)) {
        Evaluation_matrix[FP] <- Evaluation_matrix[FP] +
          1
      }
      if (NTE != 0) {
        Recall <- Evaluation_matrix[TP]/NTE
      }
      else {
        Recall <- NA
      }
      if (NIE != 0) {
        Precision <- Evaluation_matrix[TP]/NIE
      }
      else {
        Precision <- NA
      }
    }
  }
  else {
    if (GV == 0) {
      m11 <- m1
      m22 <- m2
    }
    else {
      m11 <- m1[-c(1:GV), -c(1:GV)]
      m22 <- m2[-c(1:GV), -c(1:GV)]
    }
    for (i in 1:nrow(m11)) {
      for (j in 1:ncol(m11)) {
        if (m11[i, j] == m11[j, i]) {
          m11[i, j] <- 0
        }
      }
    }
    NTE <- length(which(m11 == 1))
    for (i in 1:nrow(m22)) {
      for (j in 1:ncol(m22)) {
        if (m22[i, j] == m22[j, i]) {
          m22[i, j] <- 0
        }
      }
    }
    NIE <- length(which(m22 == 1))
    if (GV == 0) {
      m1 <- m1
      m2 <- m2
    }
    else {
      m1 <- m1[-c(1:GV), -c(1:GV)]
      m2 <- m2[-c(1:GV), -c(1:GV)]
    }
    ind1 <- t(combn(ncol(m1), 2))
    for (i in seq_len(nrow(ind1))) {
      x <- ind1[i, 1]
      y <- ind1[i, 2]
      if ((m1[x, y] == 1 & m1[y, x] != 1) & (m2[x, y] ==
                                             1 & m2[y, x] != 1)) {
        Evaluation_matrix[TP] <- Evaluation_matrix[TP] +
          edge.presence
      }
      if ((m1[y, x] == 1 & m1[x, y] != 1) & (m2[y, x] ==
                                             1 & m2[x, y] != 1)) {
        Evaluation_matrix[TP] <- Evaluation_matrix[TP] +
          edge.presence
      }
      if ((m1[x, y] == 1 & m1[y, x] == 1) & (m2[x, y] ==
                                             1 & m2[y, x] == 1)) {
        Evaluation_matrix[TP] <- Evaluation_matrix[TP] +
          edge.presence
      }
      if ((m1[x, y] == 1 & m1[y, x] != 1) & (m2[x, y] !=
                                             1 & m2[y, x] == 1)) {
        Evaluation_matrix[TP] <- Evaluation_matrix[TP] +
          edge.direction
      }
      if ((m1[x, y] != 1 & m1[y, x] == 1) & (m2[x, y] ==
                                             1 & m2[y, x] != 1)) {
        Evaluation_matrix[TP] <- Evaluation_matrix[TP] +
          edge.direction
      }
      if ((m1[x, y] == 1 & m1[y, x] == 1) & (m2[x, y] ==
                                             1 & m2[y, x] != 1)) {
        Evaluation_matrix[TP] <- Evaluation_matrix[TP] +
          edge.direction
      }
      if ((m1[x, y] == 1 & m1[y, x] == 1) & (m2[x, y] !=
                                             1 & m2[y, x] == 1)) {
        Evaluation_matrix[TP] <- Evaluation_matrix[TP] +
          edge.direction
      }
      if ((m1[x, y] == 1 & m1[y, x] != 1) & (m2[x, y] ==
                                             1 & m2[y, x] == 1)) {
        Evaluation_matrix[TP] <- Evaluation_matrix[TP] +
          edge.direction
      }
      if ((m1[x, y] != 1 & m1[y, x] == 1) & (m2[x, y] ==
                                             1 & m2[y, x] == 1)) {
        Evaluation_matrix[TP] <- Evaluation_matrix[TP] +
          edge.direction
      }
      if ((m1[x, y] == 0 & m1[y, x] == 0) & (m2[x, y] ==
                                             1 || m2[y, x] == 1)) {
        Evaluation_matrix[FP] <- Evaluation_matrix[FP] +
          1
      }
      if (NTE != 0) {
        Recall <- Evaluation_matrix[TP]/NTE
      }
      else {
        Recall <- NA
      }
      if (NIE != 0) {
        Precision <- Evaluation_matrix[TP]/NIE
      }
      else {
        Precision <- NA
      }
    }
  }
  return(list(Matrix = Evaluation_matrix, TP = Evaluation_matrix[TP],
              FP = Evaluation_matrix[FP], NTE = NTE, NIE = NIE, Recall = Recall,
              Precision = Precision))
}
