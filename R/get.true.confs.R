# Extract true confounding variables using an adjacency matrix

get.true.confs <- function (adjacency,
                               n_v = length(conf.sets$WZconfounders),
                               n_t = length(conf.sets$confounders),
                               n_w = 0, n_z = 0,
                               n_u = NCOL(adjacency) - n_v - n_t - n_w - n_z) {
  ### True parent sets
  ## V-nodes
  TrueVparentset <- get.true.variant.set(Adj = adjacency,
                                         n_t = n_t,
                                         n_v = n_v)
  if (!length(TrueVparentset)) {
    TrueVparentset <- vector(mode = "list", length = n_t)
  }

  ## T-nodes
  TrueTparentset <- get.true.parent.genes(Adj = adjacency,
                                          n_t = n_t,
                                          n_v = n_v)
  if (!length(TrueTparentset)) {
    TrueTparentset <- vector(mode = "list", length = n_t)
  }

  ## Q-nodes
  n_wz <- n_w + n_z
  TrueQset <- get.true.conf.set(Adj = adjacency,
                                n_v = n_v,
                                n_t = n_t,
                                n_wz = n_wz,
                                Upool = FALSE,
                                offset = n_v + n_t)
  if (!length(TrueQset)) {
    TrueQset <- vector(mode = "list", length = n_t)
  }

  ## U-nodes
  Trueconfset <- get.true.conf.set(Adj = adjacency,
                                   n_v = n_v,
                                   n_t = n_t,
                                   n_wz = n_wz,
                                   n_u = n_u,
                                   Upool = TRUE,
                                   offset = n_v + n_t + n_wz)
  if (!length(Trueconfset)) {
    Trueconfset <- vector(mode = "list", length = n_t)
  }

  return(list(Vset = TrueVparentset,
              Tset = TrueTparentset,
              Qset = TrueQset,
              Uset = Trueconfset))

}

