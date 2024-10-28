#'
#' Auxiliaries for \code{'conf.sets'} class object manipulations
#'
#' \code{reorder.conf.sets} updates a \code{'conf.sets'} class object after
#' the ordering of of the columns dataset used to obtain the \code{'conf.sets'}
#' is modified. This is achieved by replace column labels (numbers) of selected
#' confounding variables by new labels (positions) in the column-reordered dataset
#'
#' @param conf.sets an object of class 'conf.sets'
#'
#' @param new.order a numeric vector of the positions of variables in the original dataset,
#' i.e., if e.g. the first element of \code{new.order} is 10, then the 10th variable
#' in the original dataset is now in the first position in the reordered dataset.
#' Note that the elements of \code{new.order} must be unique, failure of this condition
#' produces an error.
#'
#' @usage
#' is.conf.sets (x)
#'
#' reorder.conf.sets (conf.sets, new.order)
#'
#' @details
#' The function \code{reorder.conf.sets} is useful to test the stability of
#' network inference methods.
#'
#' @export reorder.conf.sets
#' @export reorder.set
#'
#' @return
#' \code{reorder.conf.sets} returns an object of class 'conf.sets' with column
#' identifiers updated based on the argument \code{new.order}. However, the slot \code{raw} of
#' the input 'conf.sets' is not updated. An additional slot named \code{new.order} is added
#' to the returned object (the presence of this optional slot will indicate that
#' all other slots have been altered using \code{new.order}).
#'
#' @seealso \link{assess.conf.selection} to assess the performance of the selection
#' procedure given the adjacency matrix of the true network.
#'
#' @examples
#' ## Load data
#' library(MRGNgeneral)
#' data(confsetsA11)
#'
#' ## Test if 'confsetsA11' is a 'conf.sets' object
#' is.conf.sets (confsetsA11)
#'
#' ## Test if a list of one element 'a = 0' is a 'conf.sets' object
#' is.conf.sets (list(a = 0))
#'
#' ## Create a vector of new column order
#' set.seed(167)
#' Vorder <- sample(100, size = 100, replace = FALSE) # V-nodes
#' Torder <- sample(100, size = 100, replace = FALSE) + 100 # T-nodes
#' Corder <- sample(400, size = 400, replace = FALSE) + 200 # C-nodes
#' VTCorder <- c(Vorder, Torder, Corder)
#'
#' ## Update 'confsetsA11'
#' new.confsetsA11 <- reorder.conf.sets (confsetsA11, new.order = VTCorder)
#'
#' ## Compare the result with the expectation for V-nodes selected for individual 1
#' # old = indices of selected V-nodes in the original dataset
#' # new = indices of the same selected V-nodes 
#' #       in the dataset with columns re-ordered using VTCorder
#' # expected.new = what new is supposed to be
#' cbind(old = confsetsA11$Vconfounders[[1]],
#'       new = new.confsetsA11$Vconfounders[[1]],
#'       expected.new = sapply(confsetsA11$Vconfounders[[1]],
#'                      FUN = function(x) which(VTCorder == x)))
#

reorder.conf.sets <- function (conf.sets, new.order) {
  # Check object class
  stopifnot(is.conf.sets (conf.sets))
  stopifnot(length(new.order) == length(unique(new.order)))

  # Re-order column numbers in each component of conf.sets
  new.conf.sets <- conf.sets
  new.conf.sets$Vconfounders <- lapply(conf.sets$Vconfounders, # For each Tj
                                       FUN = reorder.set,
                                       new.order = new.order)
  new.conf.sets$Tconfounders <- lapply(conf.sets$Tconfounders, # For each Tj
                                       FUN = reorder.set,
                                       new.order = new.order)
  new.conf.sets$Uconfounders <- lapply(conf.sets$Uconfounders, # For each Tj
                                       FUN = reorder.set,
                                       new.order = new.order)
  new.conf.sets$confounders <- mapply (FUN = function(x,y,z) {c(x, y, z)},
                                       new.conf.sets$Vconfounders,
                                       new.conf.sets$Tconfounders,
                                       new.conf.sets$Uconfounders,
                                       SIMPLIFY = FALSE)
  new.conf.sets$UWZconfounders <- lapply(conf.sets$UWZconfounders, # For each Tj
                                         FUN = reorder.set,
                                         new.order = new.order)
  new.conf.sets$WZconfounders <- lapply(conf.sets$WZconfounders, # For each Vj
                                        FUN = reorder.set,
                                        new.order = new.order)
  new.conf.sets$UWZindices <- reorder.set (conf.sets$UWZindices,
                                           new.order = new.order)
  new.conf.sets$WZindices <- reorder.set (conf.sets$WZindices,
                                          new.order = new.order)
  new.conf.sets$new.order <- new.order

  return(new.conf.sets)
}

# Replace column numbers by new numbers in a reordered dataset
# x         = numeric vector of some of the original column numbers
# new.order = numeric vector where the ith element gives the original (old) position of the new ith column
# For instance, x = c(2,5,7); and new.order = c(5,10,2,9,7,1,8,4,6,3)
# is transformed into y = c(3, 1, 5) because: 2 in x is in the 3rd position in new.order,
#                                             5 in x is in the 1st position in new.order, and
#                                             7 in x is in the 5th position in new.order
reorder.set <- function(x, new.order) {
  if (is.null(x)) {
    return(NULL)
  }
  y <- sapply(x, # For each Vi selected for Tj
              FUN = function(xj) {
                which(new.order == xj)
              })
  return(y)
}