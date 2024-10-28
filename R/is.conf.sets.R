#' \code{is.conf.sets} tests if an object is of class \code{'conf.sets'}.
#' 
#' @param x an object of class 'conf.sets'
#' 
#' @details
#' The function \code{is.conf.sets} helps stop computations early if a user specified
#' \code{'conf.sets'} object is not licit.
#'
#' 
#' @usage
#' is.conf.sets (x)
#' @return
#' \code{is.conf.sets} returns a logical, \code{TRUE} is \code{x} is \code{'conf.sets'} class object.
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
#' @export
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
