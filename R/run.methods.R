#' @name run.methods
#' @title Run MRGN and bnlearn methods
#' @description Run MRGN and selected \code{bnlearn} methods on one simulated genomic dataset.
#' Depends on package \code{bnlearn} which is not imported and loaded with MRGNgeneral.
#'
#' @param simdata Simulated data object from \code{sample.graph.data}, or NULL to generate new data
#' @param n_t Integer, number of T-nodes (phenotypes)
#' @param n_v.t Numeric, number of V-nodes per T-node
#' @param conf.num.vec Named numeric vector of confounder counts
#' @param graph_type Character, type of graph structure
#' @param degree Integer, average node degree
#' @param theta Numeric, minor allele frequency
#' @param b0 Numeric, baseline value
#' @param b.snp Numeric vector, range for SNP effects
#' @param b.med Numeric vector, range for mediator effects
#' @param sigma Numeric, noise standard deviation
#' @param neg.freq Numeric, frequency of negative effects
#' @param conf.coef.ranges List of coefficient ranges for confounders
#' @param sample.size Integer, number of samples
#' @param seed Integer, random seed
#' @param conf.sets Confounder sets from \code{get.conf.sets}, or NULL to compute
#' @param blocksize Integer, block size for parallel processing
#' @param T.measure Character, measure type for T-nodes
#' @param C.measure Character, measure type for C-nodes
#' @param selection_FDRcontrol Character, FDR control method for selection
#' @param T.FDRcontrol Character, FDR control for T-nodes
#' @param C.FDRcontrol Character, FDR control for C-nodes
#' @param V.FDRcontrol Character, FDR control for V-nodes
#' @param adjust_by Character, adjustment method
#' @param T.adjust_by Character, adjustment method for T-nodes
#' @param C.adjust_by Character, adjustment method for C-nodes
#' @param V.adjust_by Character, adjustment method for V-nodes
#' @param selection_alpha Numeric, alpha level for selection
#' @param selection_fdr Numeric, FDR level for selection
#' @param lambda Numeric, tuning parameter for q-value estimation
#' @param pi0.method Character, method for pi0 estimation
#' @param bn.methods Character vector of bnlearn methods to run
#' @param MRGNfit Pre-computed MRGN fit object, or NULL
#' @param threshold_v Numeric, threshold for V-node associations
#' @param threshold_m Numeric, threshold for marginal associations
#' @param alpha Numeric, significance level
#' @param FDRcontrol Character, FDR control method for MRGN
#' @param fdr Numeric, FDR level for MRGN
#' @param analyse.triplets Logical, whether to analyze triplets
#' @param stringent Logical, use stringent edge criteria
#' @param restart.hc Integer, number of restarts for hill-climbing
#' @param n.reorder Integer, number of reordering iterations
#' @param verbose Integer, verbosity level
#' @param nb.cl Integer, number of parallel clusters
#' @param savetopath Logical, whether to save results to path
#' @param path Character, base path for saving
#' @param setpath Character, subdirectory for saving
#'
#' @return A list containing performance metrics, simulated data, confounder sets, and model fits
#'
#' @export run.methods
#'
run.methods <- function (simdata = NULL,
                         n_t = 100,
                         n_v.t = 1,
                         conf.num.vec = c(W = 0, Z = 0, U = 200, K = 0, I = 100),
                         graph_type = "scale-free", degree = 3,
                         theta = 0.4, b0 = 0, b.snp = c(-.5, .5),
                         b.med = c(-.8, .8), sigma = 0.1, neg.freq = 0.5,
                         conf.coef.ranges = list(W = c(0.4, 0.5),
                                                 Z = c(1, 1.5),
                                                 U = c(0.4, 0.5),
                                                 K = c(0.4, 0.5)),
                         sample.size = 500,
                         seed = NULL,
                         conf.sets = NULL,
                         blocksize = min(n_v, n_t, 100),
                         T.measure = 'partial',
                         C.measure = 'partial',
                         selection_FDRcontrol = 'qvalue',
                         T.FDRcontrol = selection_FDRcontrol,
                         C.FDRcontrol = selection_FDRcontrol,
                         V.FDRcontrol = selection_FDRcontrol,
                         adjust_by = 'individual',
                         T.adjust_by = adjust_by,
                         C.adjust_by = adjust_by,
                         V.adjust_by = adjust_by,
                         selection_alpha = 0.01,
                         selection_fdr = 0.05,
                         lambda = 0.05,
                         pi0.method = 'smoother',
                         bn.methods = c('none', 'tabu', 'hc', 'pc.stable', 'mmhc'),
                         MRGNfit = NULL,
                         threshold_v = 0.2,
                         threshold_m = 0.05,
                         alpha = 0.05,
                         FDRcontrol = "bonferroni", fdr = 0.05,
                         analyse.triplets = TRUE, stringent = TRUE,
                         restart.hc = 10,
                         n.reorder = 0,
                         verbose = 2L, nb.cl = 4,
                         savetopath = FALSE, path = '', setpath = '') {
  stopifnot(is.character(bn.methods))
  bn.methods <- tolower(bn.methods)
  stopifnot(all(bn.methods %in% c('none', 'tabu', 'hc', 'pc.stable', 'mmhc')))
  if (!all(bn.methods == 'none')) {
    stopifnot('bnlearn' %in% installed.packages()[,1])
  }

  ### Create the directory path to save results (if non existent)
  if (savetopath) {
    savepath <- paste0(path, setpath, "/")
    if (!missing(path)) {
      if (!dir.exists(savepath))
        dir.create(savepath)
    }
  }

  ### Save random generator state to restitute it later
  if(!is.null(seed)) {
    saved.seed <- .Random.seed
    set.seed(seed)
  }

  ### Generate data if not provided
  if (is.null(simdata)) {
    # No known confounder is allowed in these simulation: all C-nodes are candidate for selection
    if (length(conf.num.vec) >= 4) {
      conf.num.vec[3] <- conf.num.vec[3] + conf.num.vec[4]
      conf.num.vec[4] <- 0
    }

    simdata <- sample.graph.data (n_t = n_t,
                                  n_v.t = n_v.t,
                                  conf.num.vec = conf.num.vec,
                                  graph_type = graph_type,
                                  degree = degree,
                                  theta = theta,
                                  b0 = b0,
                                  b.snp = b.snp,
                                  b.med = b.med,
                                  sigma = sigma,
                                  neg.freq = neg.freq,
                                  conf.coef.ranges = conf.coef.ranges,
                                  scale = FALSE,
                                  sample.size = sample.size,
                                  seed = NULL)
  }
  if (savetopath){
    save(simdata, file = paste0(path, setpath, "/", "simdata", ".RData"))
  }
  n_v <- simdata$dims$n_v
  n_t <- simdata$dims$n_t
  n_vt <- n_v + n_t
  n_w <- simdata$dims$n_w
  n_z <- simdata$dims$n_z
  n_wz <- n_w + n_z
  n_u <- simdata$dims$n_u + simdata$dims$n_k # (U and K-nodes are considered U-nodes)
  n_c <- n_w + n_z + n_u + simdata$dims$n_i
  n_vtc <- n_vt + n_c

  ### Select U, V, T
  if (is.null(conf.sets)) {
    if (!missing(blocksize))
      blocksize <- min(blocksize, n_t, n_v, 100)

    mycl <- if (nb.cl > 0) parallel::makeCluster(nb.cl) else NULL
    if (nb.cl > 0) {
      setup_cluster(mycl, packages = c("MRGN", "ppcor", "propagate"), seed = seed)
    }
    conf.sets <- get.conf.sets (data = simdata$data[, 1:n_vtc], scale.data = TRUE,
                               n_v = n_v, n_t = n_t, n_c = n_c,
                               T.measure = T.measure,
                               C.measure = C.measure,
                               blocksize = blocksize,
                               FDRcontrol = selection_FDRcontrol,
                               V.FDRcontrol = V.FDRcontrol,
                               T.FDRcontrol = T.FDRcontrol,
                               C.FDRcontrol = C.FDRcontrol,
                               adjust_by = adjust_by,
                               V.adjust_by = V.adjust_by,
                               T.adjust_by = T.adjust_by,
                               C.adjust_by = C.adjust_by,
                               fdr = selection_fdr,
                               lambda = lambda,
                               pi0.method = pi0.method,
                               alpha = selection_alpha,
                               parallel = nb.cl > 0, cl = mycl,
                               verbose = verbose,
                               save.list = FALSE)
    MRGNgeneral:::catch.conditions({
      if (nb.cl > 0) parallel::stopCluster(cl = mycl)
    })
  }

  if (savetopath){
    save(conf.sets, file = paste0(path, setpath, "/", "conf.sets", ".RData"))
  }

  ConfSelection <- MRGNgeneral::assess.conf.selection (conf.sets,
                                                       adjacency = simdata$adjacency,
                                                       n_v = n_v,
                                                       n_t = n_t,
                                                       n_w = n_w, n_z = n_z,
                                                       n_u = n_u)

  cat("\n       Recall and Precision of selected confounders:       \n")
  print(rbind(`Recall` = ConfSelection$recall,
              `Precision` = ConfSelection$precision))

  ## Number of selected W,Z-nodes (in place of n_wz, the truth)
  n_q <- length(conf.sets$WZindices)
  n_vtq <- n_vt + n_q

  ### Input adjacency
  Adjacency0 <- MRGNgeneral::get.initial.skeleton (data = simdata$data,
                                             n_t = n_t, n_v = n_v,
                                             threshold_v = threshold_v,
                                             threshold_m = threshold_m,
                                             conf.sets = conf.sets)
  ### Index for sets of variables
  indexV <- 1:n_v
  indexT <- (n_v + 1):n_vt
  #indexVT <- 1:n_vt

  if (is.null(MRGNfit)) {
    ## Run MRGN: including selected T, U and V-nodes as confounders
    cat(paste0("\n ************ MRGN general ... \n"))
    mycl <- if (nb.cl > 0) parallel::makeCluster(nb.cl) else NULL
    if (nb.cl > 0) {
      setup_cluster(mycl, packages = c("MRGN", "ppcor"), seed = if (!is.null(seed)) seed + 1000 else NULL)
    }
    MRGNTime <- system.time({
      MRGNfit <- MRGN(data = simdata$data[, 1:n_vtc, drop = FALSE],
                      Qlabels = conf.sets$WZindices,
                      n_t = n_t,
                      n_v = n_v,
                      n_q = n_q,
                      n_u = n_c - n_q,
                      adjacency = Adjacency0,
                      confounders = conf.sets$confounders,
                      alpha = alpha,
                      FDRcontrol = FDRcontrol,
                      fdr = fdr,
                      lambda = lambda,
                      analyse.triplets = analyse.triplets,
                      stringent = stringent,
                      verbose = verbose,
                      parallel = nb.cl > 0, cl = mycl)
    })

    MRGNgeneral:::catch.conditions({
      parallel::stopCluster(cl = mycl)
    })

    MRGNfit$time <- MRGNTime

    cat("\n ............ DONE. \n")
  }
  if (is.null(MRGNfit)) {
    MRGNfit <- list(time = rep(NA, 3))
  }
  else {
    if (is.null(MRGNfit$time)) {
      MRGNfit$time <- rep(NA, 3)
    }

    MRGNfit$Performance <- unlist(MRGNgeneral::RecallPrecision(g1 = as (simdata$adjacency[indexT, indexT], "graphNEL"),
                                                               g2 = as (MRGNfit$adjacency[indexT, indexT], "graphNEL"),
                                                               GV = 0,
                                                               includeGV = FALSE)[-1])[-2]

    if (savetopath){
      save(MRGNfit, file = paste0(path, setpath, "/", "MRGN", ".RData"))
    }
  }

  print(rbind(c(MRGNfit$Performance, MRGNfit$time[3])))

  ##################### bnlearn methods##################
  if(any(c('tabu', 'hc', 'pc.stable', 'mmhc') %in% bn.methods)) {
    ## Build matrix of Q,U -> T edges from confounding variable selection results
    AdjQU_T <- get.adj.from.conf.list (conf.sets$UWZconfounders,
                                                     n_c = n_c,
                                                     offset = n_vt)

    ### matrix indicating fixed and absent edges
    # Initialize a matrix of zeros
    blackmatrix0 = matrix(0,
                          nrow = n_vtc,
                          ncol = n_vtc)

    ### Index for sets of variables
    indexQU <- (n_vt + 1):n_vtc
    indexQ <- conf.sets$WZindices
    indexU <- conf.sets$UWZindices[!(conf.sets$UWZindices %in% indexQ)]

    #indexQ <- (n_vt + 1):n_vtq
    indexVTQU <- 1:n_vtc

    #  indexVTQ <- 1:n_vtq
    #  indexUK <- (n_vtq + 1):(n_vtq + n_u)

    # V - V edges must be absent
    blackmatrix0[indexV, indexV] = 1

    # Q,U - Q,U edges must be absent
    blackmatrix0[indexQU, indexQU] = 1

    # V -> Q,U and Q,U -> V edges must be absent
    blackmatrix0[indexV, indexQU] = 1
    blackmatrix0[indexQU, indexV] = 1

    # T -> V edges must be absent
    blackmatrix0[indexT, indexV] = 1

    # T -> U edges must be absent
    blackmatrix0[indexT, indexU] = 1

    # U -> T edges for non-selected confounders must be absent
    blackmatrix0[indexU, indexT] = 1 - AdjQU_T[indexU - n_vt, ]

    ### matrix indicating fixed and present edges
    # Initialize a matrix of zeros
    whitematrix0 = matrix(0, nrow = dim(simdata$data)[2],
                          ncol = dim(simdata$data)[2])

    # U -> T edges for selected Q and U must be present
    whitematrix0[indexU, indexT] <- AdjQU_T[indexU - n_vt, ] # replace indexQU by the index for just U

    # Check that No duplicate in whitematrix0 and blackmatrix0
    # sum(whitematrix0 * blackmatrix0 == 1)

    # White and black lists
    whitelist0 = which(whitematrix0 == 1, arr.ind = TRUE)
    blacklist0 = which(blackmatrix0 == 1, arr.ind = TRUE)

    # Data for bnlearn calls
    bndata <- simdata$data
    colnames(bndata) <- paste0('G', 1:dim(simdata$data)[2])
    dnames <- colnames(simdata$data)

    bnfits = bnlearning (data = bndata, dnames = dnames,
                         bn.methods = bn.methods,
                         whitelist0 = whitelist0,
                         blacklist0 = blacklist0,
                         restart = restart.hc,
                         adjacency = simdata$adjacency,
                         indexT = indexT,
                         verbose = verbose, nb.cl = nb.cl,
                         savetopath = savetopath,
                         path = path, setpath = setpath)
  }
  else {
    bnfits <- list(TABU = list(time = rep(NA, 3),
                                Performance = rep(NA, 5)),
                    HC = list(time = rep(NA, 3),
                              Performance = rep(NA, 5)),
                    PCSTABLE = list(time = rep(NA, 3),
                                    Performance = rep(NA, 5)),
                    MMHC = list(time = rep(NA, 3),
                                Performance = rep(NA, 5)))
  }

  # Output list
  out <- list(Performance = c(ConfSelection$recall,
                              ConfSelection$precision,
                              MRGNfit$Performance, MRGNfit$time[3],
                              bnfits$TABU$Performance, bnfits$TABU$time[3],
                              bnfits$HC$Performance, bnfits$HC$time[3],
                              bnfits$PCSTABLE$Performance, bnfits$PCSTABLE$time[3],
                              bnfits$MMHC$Performance, bnfits$MMHC$time[3]),
              simdata = simdata,
              conf.sets = conf.sets,
              conf.perf = ConfSelection,
              fits = list(MRGN = MRGNfit,
                          bn = bnfits))



  ### Reorder columns of data and re-run the all inferences
  out$reorder <- list()
  if (n.reorder) {
    if(any(c('tabu', 'hc', 'pc.stable', 'mmhc') %in% bn.methods)) {
      for (k in 1:n.reorder) { # Reorder and re-run
        out$reorder[[k]] <- reorder_nre_run_methods (seed = NULL,
                                                     simdata = simdata,
                                                     conf.sets = conf.sets, Adjacency0 = Adjacency0,
                                                    alpha = alpha, FDRcontrol = FDRcontrol,
                                                    fdr = fdr, lambda = lambda,
                                                    analyse.triplets = analyse.triplets,
                                                    stringent = stringent, verbose = verbose,
                                                    nb.cl = nb.cl, indexT = indexT,
                                                    bn.methods = bn.methods, restart.hc = restart.hc,
                                                    whitelist0 = whitelist0, blacklist0 = blacklist0)
      }
    }
    else {
      for (k in 1:n.reorder) { # Reorder and re-run
        out$reorder[[k]] <- reorder_nre_run_methods (seed = NULL,
                                                     simdata = simdata,
                                                     conf.sets = conf.sets, Adjacency0 = Adjacency0,
                                                     alpha = alpha, FDRcontrol = FDRcontrol,
                                                     fdr = fdr, lambda = lambda,
                                                     analyse.triplets = analyse.triplets,
                                                     stringent = stringent, verbose = verbose,
                                                     nb.cl = nb.cl, indexT = indexT)
      }
    }

  }

  if (savetopath) {
    save(out, file = paste0(savepath, "out", ".RData"))
  }

  ### Restitute random generator state
  if(!is.null(seed)) {
    .Random.seed <- saved.seed
  }

  return(out)

}
