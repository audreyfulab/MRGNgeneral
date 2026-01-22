### A function to convert a list of selected confounders (C-nodes)
### into an adjacency matrix
# n_c = Total number of candidate confounding variables
# n_v = Number of variants (V-nodes)
# Return a n_c X n_t matrix of 0/1 where
# n_t in the number of genes (T-nodes); n_t = length(conf.list)

#' @keywords internal
#' @noRd

get.adj.from.conf.list <- function (conf.list,
                                    n_c,
                                    offset = n_t + n_v,
                                    n_v = 0) {
  # Number of genes (T-nodes)
  n_t <- length(conf.list)

  # Convert each element of 'conf.list' (vector of indices of U-nodes selected for a T-node)
  # into a matrix column using 'confounder.vec2adj'
  Adj <- sapply(conf.list,
                FUN = confounder.vec2adj,
                n_c = n_c,
                offset = offset)
  dimnames(Adj) <- list(paste0('C', 1:n_c), paste0('T', 1:n_t))

  # Return an 'adjacency.matrix' object
  return(structure(Adj, class = "adjacency.matrix"))
}

# Child of 'get.adj_from_conf_list' to convert vector of selected
# confounders into a column in an adjacency matrix
confounder.vec2adj <- function (confounder.vec, n_c, offset) {
  # Initialize the column
  adj_col <- numeric(n_c)

  # Set ones for selected confounders
  if(!is.null(confounder.vec))
    adj_col[confounder.vec - offset] <- 1

  return(adj_col)
}
