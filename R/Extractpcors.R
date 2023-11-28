# Compute partial correlations between a V-node indexed 'v' and each T-node,
# given all other T-nodes, starting from a covariance matrix
extractpcorVT <- function(v, cov.mat, p, q) {
  vTcov.mat <- cov.mat[c(v, (q+1):(p+q)), c(v, (q+1):(p+q))]
  H <- - cov2cor(mpinv(vTcov.mat))
  H[1, -1]
}

# Compute partial correlations between each U-node indexed 'u' (from a pool 'U.pool')
# and each T-node indexed 't', given all V-nodes, starting from a covariance matrix
extractpcorUT <- function(t, cov.mat, U.pool, p, q) {
  sapply(U.pool,
         FUN = extractpcorUTi, t = t,
         cov.mat = cov.mat, U.pool = U.pool, p = p, q = q)
}

extractpcorUTi <- function(u, t, cov.mat, U.pool, p, q) {
  uTcov.mat <- cov.mat[c(t, u, 1:q), c(t, u, 1:q)]
  H <- - cov2cor(mpinv(uTcov.mat))
  H[1, 2]
}
