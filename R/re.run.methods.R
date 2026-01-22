
# Function to test the stability of network inference method
# Re-order data columns and re-run network inference methods

#' @keywords internal
#' @noRd
reorder_nre_run_methods <- function (seed = NULL,
                                     simdata, conf.sets, Adjacency0,
                                     alpha, FDRcontrol, fdr, lambda,
                                     analyse.triplets, stringent,
                                     nb.cl, indexT,
                                     bn.methods = "none",
                                     whitelist0,
                                     blacklist0,
                                     restart.hc,
                                     verbose) {
  n_v <- simdata$dims$n_v
  n_t <- simdata$dims$n_t
  n_vt <- n_v + n_t
  n_w <- simdata$dims$n_w
  n_z <- simdata$dims$n_z
  n_wz <- n_w + n_z
  n_u <- simdata$dims$n_u + simdata$dims$n_k
  n_c <- n_w + n_z + n_u + simdata$dims$n_i
  n_vtc <- n_vt + n_c
  n_q <- length(conf.sets$WZindices)
  n_vtq <- n_vt + n_q

  ### Save random generator state to restitute it later
  saved.seed <- .Random.seed
  if(!is.null(seed)) {
    set.seed(seed)
  }

  # Generate the new ordering for each set of variables (V, T, and C-nodes)
  Vorder <- sample(n_v, size = n_v, replace = FALSE)
  Torder <- sample(n_t, size = n_t, replace = FALSE) + n_v
  Corder <- sample(n_c, size = n_c, replace = FALSE) + n_vt
  VTCorder <- c(Vorder, Torder, Corder)

  # Re-order columns in 'simdata'
  simdata$data <- simdata$data[, VTCorder]
  simdata$adjacency <- simdata$adjacency[VTCorder, VTCorder]

  # Re-order column indices in 'conf.sets'
  conf.sets <- reorder_conf_sets (conf.sets, new.order = VTCorder)
  Adjacency0 <- Adjacency0[VTCorder, VTCorder]

  # Re-order column indices in 'whitelist0' and 'blacklist0'
  whitelist0[,1] <- reorder_set (whitelist0[,1], new.order = VTCorder)
  whitelist0[,2] <- reorder_set (whitelist0[,2], new.order = VTCorder)
  blacklist0[,1] <- reorder_set (blacklist0[,1], new.order = VTCorder)
  blacklist0[,2] <- reorder_set (blacklist0[,2], new.order = VTCorder)

  # Data for bnlearn calls
  bndata <- simdata$data
  colnames(bndata) <- paste0('G', 1:dim(simdata$data)[2])
  dnames <- colnames(simdata$data)



  ### Restitute random generator state
  .Random.seed <- saved.seed


}
