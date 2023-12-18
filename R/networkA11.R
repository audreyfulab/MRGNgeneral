#'
#' A simulated genomic data
#'
#' An \code{.rda} structure as returned by \link{sample.graph.data}, i.e.
#' a simulated Direct Acyclic Graph (DAG), a simulated genomic data from the
#' DAG, and descriptors of the number of variants, phenotypes,
#' and confounding variables in the simulated network.
#'
#' @format  A named list with elements:
#' \itemize{
#'   \item \code{data}: a \code{data.frame} with 500 rows and 600 columns
#'   including 100 variants (named \code{V1}-\code{V100}),
#'   followed by 100 phenotypes (\code{T1}-\code{T100}),
#'   400 confounding variables. The latter includes
#'   50 intermediate variables (\code{W1}-\code{W100}),
#'   50 common children (\code{Z1}-\code{Z100}), and
#'   200 confounders (\code{U1}-\code{U100}),
#'   and 100 independent variables (\code{I1}-\code{I100}, mutually of independent all others variables).
#'   \item dims: a list with named elements \code{n_v = 100} (number of variants),
#'   \code{n_t = 100} (number of phenotypes),
#'   \code{n_w = 100} (number of intermediate variables),
#'   \code{n_z = 50} (number of common children),
#'   \code{n_u = 50} (number of unknown confounders),
#'   \code{n_k = 0} (number of known confounders),
#'   \code{n_i = 100} (number of independent variables).
#'   \item \code{b0}: conditional mean of any node given all its parent nodes.
#'   \item \code{sigma}: conditional standard deviation of any node given all its parent nodes.
#'   \item \code{adjacency}: \eqn{500 \times 500} adjacency matrix of the true network,
#'   excluding the \code{n_i} independent variables.
#'   \item \code{effects}: \eqn{500 \times 500} matrix of linear effect sizes
#'   (coefficients in the linear model used to simulate the data).
#'   \item \code{igraph}: \code{igraph} object (from the package \code{igraph})
#'   corresponding to the true network.
#' }
#'
#' @docType data
#' @keywords datasets
#'
#' @source Simulated in the package \code{MRGNgeneral}.
#'
#' @examples
#' ## Load the data
#' library(MRGNgeneral)
#' data(networkA11)
#'
#' ## Adjacency matrix of a subset of the network
#' adjacency <- structure(
#'   networkA11$adjacency[c('V39', 'T39', 'T43', 'T52', 'W11', 'Z2', 'U28'),
#'                        c('V39', 'T39', 'T43', 'T52', 'W11', 'Z2', 'U28')],
#'   class = 'adjacency.matrix')
#' adjacency
#'
#' ## Plot the graph of the subset
#' plot (adjacency)
#'
#' ## Scatterplot of the data for the subset
#' pairs (networkA11$data[, c('V39', 'T39', 'T43', 'T52', 'W11', 'Z2', 'U28')])
#'
#'
"networkA11"

#library(MRGNgeneral)
#set.seed(11000)
#networkA11 <- sample.graph.data (n_t = 100,
#                              n_v.t = 1,
#                              family.n_v = NULL,
#                              conf.num.vec = c(W = 50, Z = 50,
#                                               U = 200, K = 0, I = 100),
#                              graph_type = "scale-free",
#                              degree = 3,
#                              theta = .4,
#                              b0 = 0,
#                              b.snp = c(-0.5, 0.5),
#                              b.med = c(-0.8, 0.8),
#                              sigma = 0.1,
#                              neg.freq = 0.5,
#                              conf.coef.ranges = list(W = c(0.4, 0.5),
#                                                      Z = c(1, 1.5),
#                                                      U = c(0.4, 0.5),
#                                                      K = c(0.01, 0.1)),
#                              scale = FALSE,
#                              sample.size = 500)
#
#save(networkA11, file="data/networkA11.rda")
#
