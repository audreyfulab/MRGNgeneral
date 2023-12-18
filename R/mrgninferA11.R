#'
#' An object of class \code{'MRGN'}
#'
#' An \code{.rda} structure as returned by \link{MRGN}, i.e.
#' an object of class \code{'MRGN'}.
#'
#' @format See the section "Value" of \link{MRGN} for a description.
#'
#' @docType data
#' @keywords datasets
#'
#' @source Generated in the package \code{MRGNgeneral}.
#'
#' @examples
#' ## Load data (mrgninferA11)
#' library(MRGNgeneral)
#' data(mrgninferA11)
#'
#' ## Print a summary
#' mrgninferA11
#'
"mrgninferA11"

#Adjacency0 <- get.initial.skeleton (data = networkA11$data,
#                                    n_v = networkA11$dims$n_v,
#                                    n_t = networkA11$dims$n_t,
#                                    threshold_v = 0.2,
#                                    threshold_m = 0.05,
#                                   conf.sets = confsetsA11)
#
#
#mrgninferA11 <- MRGN(data = networkA11$data,
#                     n_v = networkA11$dims$n_v,
#                     n_t = networkA11$dims$n_t,
#                     Qlabels = confsetsA11$WZindices,
#                     n_q = length(confsetsA11$WZindices),
#                     n_u = networkA11$dims$n_w + networkA11$dims$n_z +
#                       networkA11$dims$n_u + networkA11$dims$n_k +
#                       networkA11$dims$n_i - length(confsetsA11$WZindices),
#                     adjacency = Adjacency0,
#                     confounders = confsetsA11$confounders,
#                     alpha = 0.01,
#                     FDRcontrol = 'bonferroni',
#                     fdr = 0.05,
#                     verbose = TRUE)
#save(mrgninferA11, file="data/mrgninferA11.rda")
#
