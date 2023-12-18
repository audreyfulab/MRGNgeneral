#'
#' Auxiliaries for \code{'conf.sets'} class object manipulations
#'
#' \code{is.conf.sets} tests if an object is of class \code{'conf.sets'}.
#' \code{reorder.conf.sets} updates a \code{'conf.sets'} class object after
#' the ordering of of the columns dataset used to obtain the \code{'conf.sets'}
#' is modified. This is achieved by replace column labels (numbers) of selected
#' confounding variables by new labels (positions) in the column-reordered dataset
#'
#' @param conf.sets,x an object of class 'conf.sets'
#'
#' @param new.order a numeric vector of the positions of variables in the original dataset,
#' i.e., if e.g. the first element of \code{new.order} is 10, then the 10th variable
#' in the original dataset is now in the first position in the reordered dataset.
#' Note that the elements of \code{new.order} must be unique, failure of this condition
#' produces an error.
#'
#' @param ... not used.
#'
#' @usage
#' is.conf.sets (x)
#'
#' ## S3 method for class 'conf.sets'
#' reorder (conf.sets, new.order, ...)
#'
#' @details
#' The function \code{is.conf.sets} helps stop computations early if a user specified
#' \code{'conf.sets'} object is not licit.
#'
#' The function \code{reorder.conf.sets} is useful to test the stability of
#' network inference methods.
#'
#' @export is.conf.sets
#' @export reorder.conf.sets
#'
#' @return
#' \code{is.conf.sets} returns a logical, \code{TRUE} is \code{x} is \code{'conf.sets'} class object.
#'
#' \code{reorder.conf.sets} returns an object of class 'conf.sets' with column
#' identifiers updated based on the argument \code{new.order}. However, the slot code{raw} of
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
#' data(confsetsa11)
#'
#' ## Test if 'confsetsa11' is a 'conf.sets' object
#' is.conf.sets (confsetsa11)
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
#' ## Update 'confsetsa11'
#' new.confsetsa11 <- reorder.conf.sets (confsetsa11, new.order = VTCorder)
#'
#' ## Compare the result with the expectation for V-nodes selected for individual 1
#' cbind(# old = indices of selected V-nodes in the original dataset,
#'       old = confsetsa11$Vconfounders[[1]],
#'       # new = indices of the same selected V-nodes in the dataset with columns re-ordered using VTCorder,
#'       new = new.confsetsa11$Vconfounders[[1]],
#'       # expected.new = what new is supposed to be.
#'       expected.new = sapply(confsetsa11$Vconfounders[[1]],
#'                      FUN = function(x) which(VTCorder == x)))
#

reorder.conf.sets <- function (conf.sets, new.order, ...) {
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

# Test is an object is of class \code{'conf.sets'}
#
# Tests if its argument is a \code{'conf.sets'} class object
#
is.conf.sets <- function(x) {
  # First check that 'x' is a list
  if (!is.list(x)) {
    return(FALSE)
  }

  # Check that all required fields/slots are presents
  required_slots <- c('Vconfounders', 'Tconfounders',
                      'Uconfounders', 'WZconfounders')

  object_slots <- names(x)
  yes <- all(required_slots %in% object_slots)
  if (!yes) {
    return(FALSE)
  }

  # Check that these slots are NULL or lists
  NL <- sapply(x[c('Vconfounders', 'Tconfounders',
                        'Uconfounders', 'WZconfounders')],
               FUN = function(xj) is.null(xj) | is.list(xj))
  yes <- all(NL)
  if (!yes) {
    return(FALSE)
  }

  # Check that these slots are NULL or lists
  nn <- sapply(x[c('Vconfounders', 'Tconfounders',
                        'Uconfounders', 'WZconfounders')],
               FUN = length)
  n_t <- max(nn)
  yes <- all(nn[1:3] %in% c(0, n_t))
  if (!yes) {
    return(FALSE)
  }

  n_v <- nn[4]
  if (n_v) {
    if (nn[1] & NL[1]) {
      yes <- all(unlist(x$Vconfounders) %in% 1:n_v)
      if (!yes) {
        return(FALSE)
      }
    }
  }

  if (n_t) {
    if (nn[2] & NL[2]) {
      yes <- all(unlist(x$Tconfounders) %in% (1+n_t):(n_v+n_t))
      if (!yes) {
        return(FALSE)
      }
    }

    if (nn[3] & NL[3]) {
      yes <- all(unlist(x$Uconfounders) > (n_v+n_t))
      if (!yes) {
        return(FALSE)
      }
    }
  }

  return(TRUE)
}
