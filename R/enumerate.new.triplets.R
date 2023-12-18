# triplet.set = trios that were inferred "Other" from a previous analysis

get.new.triplet.type <- function (triplet.set,
                                  adjacency,
                                  cl, chunk.size = NULL) {
  # Get update-able triplets (and types of triplet) if any
  triplets <- matteApply(triplet.set[,1:3, drop = FALSE],
                              MARGIN = 1,
                              FUN = get.new.triplet.type.i,
                              adjacency = adjacency,
                              cl = cl, chunk.size = chunk.size)
  if (NCOL(triplets) > 1) {
    triplets <- t(triplets)
  }
  else {
    triplets <- matrix(c(triplets), nrow = 1)
  }

  colnames(triplets) <- c("Ti", "Tj", "Tk", "types")

  return(triplets)
}

# x is a vector of three elements such that we have T1-T2-T3
get.new.triplet.type.i <- function(x, adjacency) {
  Aijk <- adjacency[x, x]

  if (all(Aijk == matrix(c(0, 1, 0,
                           1, 0, 1,
                           0, 1, 0), nrow = 3))) { # Type 1
    return(c(x, 1))
  }
  else if (all(Aijk == matrix(c(0, 1, 0,
                                0, 0, 1,
                                0, 1, 0), nrow = 3, byrow = TRUE))) { # Type 2
    return(c(x, 2))
  }
  else { # Type 3
    return(c(x, 3))
  }
}
