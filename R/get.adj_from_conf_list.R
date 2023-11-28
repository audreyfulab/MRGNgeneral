### A function to convert a list of selected confounders (U-nodes)
### into an adjacency matrix
# u = Total number of candidate confounders
# q = Number of variants (V-nodes)
# Return a u X p matrix of 0/1 where p in the number of genes (T-nodes)

#' @export get.adj_from_conf_list

get.adj_from_conf_list <- function (confounder.list,
                                    n_c,
                                    offset = p + q,
                                    q = 0) {
  # Number of genes (T-nodes)
  p <- length(confounder.list)

  # Convert each element of 'confounder.list' (vector of indices of U-nodes selected for a T-node)
  # into a matrix column using 'confounder.vec2adj'
  Adj <- sapply(confounder.list,
                FUN = confounder.vec2adj,
                n_c = n_c,
                offset = offset)
  dimnames(Adj) <- list(paste0('C', 1:n_c), paste0('T', 1:p))

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
