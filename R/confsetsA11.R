#'
#' An object of class \code{'conf.sets'}
#'
#' An \code{.rda} structure as returned by \link{get.conf.sets}, i.e.
#' an object of class \code{'conf.sets'}.
#'
#' @format  A named list with elements \code{confounder}, \code{Vconfounder},
#' \code{Tconfounder}, \code{Uconfounder}, \code{confounder},
#' \code{UWZconfounders}, \code{WZconfounders}, \code{UWZindices},
#' \code{WZindices}, \code{raw}, \code{time}, and \code{call}.
#' See the section "Value" of \link{get.conf.sets} for a description of these elements.
#'
#' @docType data
#' @keywords datasets
#'
#' @source Generated in the package \code{MRGNgeneral}.
#'
#' @examples
#' ## Load the raw data (networkA11) used to obtain 'confsetsA11'
#' library(MRGNgeneral)
#' data(networkA11)
#' data(confsetsA11)
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
#'
"confsetsA11"
#confsetsA11 <- get.conf.sets(data = networkA11$data,
#                             n_v = networkA11$dims$n_v,
#                             n_t = networkA11$dims$n_t,
#                             n_c = NCOL(networkA11$data) - networkA11$dims$n_v - networkA11$dims$n_t,
#                             blocksize = min(networkA11$dims$n_v, networkA11$dims$n_t, 100),
#                             T.measure = 'partial',
#                             C.measure = 'partial',
#                             FDRcontrol = 'qvalue',
#                             adjust_by = 'individual',
#                             alpha = 0.01,
#                             fdr = 0.05,
#                             lambda = 0.05,
#                             pi0.method = 'smoother')
#save(confsetsA11, file="data/confsetsA11.rda")
#
