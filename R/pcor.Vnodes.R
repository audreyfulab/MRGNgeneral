
# Compute T-T node partial correlations, possibly conditioning only on V-nodes (all)
partial.cor.given.Vnodes <- function(data, p, q,
                                     onlyVnodes = TRUE,
                                     Mcov = NULL,
                                     parallel = FALSE,
                                     cl =  parallel::getDefaultCluster(),
                                     chunk.size = NULL) {
  if (is.null(Mcov)) {
    Mcov <- cov(data)
  }

  if (onlyVnodes) {
    ### Take each pair of T nodes
    resultat <- matrix(1, nrow = p, ncol = p)
    indices <- which(lower.tri(resultat, diag = FALSE), arr.ind = TRUE) + q
    if (parallel) {
      # Revise to compute the var-cov matrix once, and take the appropriate subset to get the partial cor for each pair
      corvalues <- matteApply (X = indices, MARGIN = 1,
                                      FUN = childpartial.cor.given.Vnodes,
                                      Mcov = Mcov,
                                      q = q,
                                      cl = cl,
                                      chunk.size = chunk.size, simplify = TRUE)
    }
    else {
      corvalues <- apply(indices,
                         MARGIN = 1,
                         FUN = childpartial.cor.given.Vnodes,
                         Mcov = Mcov,
                         q = q)
    }
    indices <- indices - q
    resultat[indices] <- corvalues
    resultat[indices[,2:1]] <- corvalues
  }
  else {
    resultat <- - cov2cor(mpinv(Mcov))
    diag(resultat) <- 1
    dimnames(resultat) <- dimnames(Mcov)

    resultat <- resultat[(q+1):(q+p), (q+1):(q+p)]
  }
  return(resultat)
}

childpartial.cor.given.Vnodes <- function(index, q, Mcov) {
  H <- - cov2cor(mpinv(Mcov[c(index, 1:q), c(index, 1:q)]))
  return(H[1,2])
}
