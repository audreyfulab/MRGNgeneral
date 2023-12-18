# Run a few bnlearn network inference algorithms
# Internal routine for 'run.methods'
bnlearning <- function (data, # data.frame/matrix, typically the \code{data} element of an output from \link{sample.graph.data}.
                        dnames = colnames(data), # character vector, names for nodes corresponding to colnames of 'data'
                        bn.methods, # list, bnlearn network inference algorithms, one of: 'none', 'tabu', 'hc', 'pc.stable' or 'mmhc'
                        whitelist0, # two-column matrix, each row indicates an edge direction that must be present in the final inferred network
                        blacklist0, # two-column matrix, each row indicates an edge direction that must be absent in the final inferred network
                        restart,    # an integer, the number of random restarts.
                        adjacency = NULL,  # adjacency matrix of the true network from which \code{data} was generated
                        indexT = 1:NCOL(data),     # integer vector, indices of the columns of \code{data} corresponding to genes/phenotypes of interest and on which performance measures are computed.
                        verbose, nb.cl, #
                        savetopath,     # logical, save the output of each method, and the final output?
                        path, setpath) {

  # Initialize outputs
  Tabu <- HC <- PCstable <- MMHC <- list(time = rep(NA, 3),
                                         Performance = rep(NA, 5))

  # Prepare white/black lists
  whitelist1 = cbind(paste0('G', whitelist0[,1]), paste0('G', whitelist0[,2]))
  blacklist1 = cbind(paste0('G', blacklist0[,1]), paste0('G', blacklist0[,2]))

  if(any(c('tabu', 'hc', 'pc.stable', 'mmhc') %in% bn.methods)) {
    require(bnlearn)
    ## Tabu Search -- Score-based (modified Hill-Climbing to escape local optima)
    if('tabu' %in% bn.methods) {
      if (verbose)
        cat(paste0("\n ************ Tabu Search ... \n"))
      TimeTabu = system.time({
        Tabu =  tabu (x = data,
                      whitelist = whitelist1,
                      blacklist = blacklist1)
      })
      # Error in build.whitelist(whitelist, nodes = names(x), data = x, algo = heuristic, :
      # score-based algorithms do not support whitelisting both directions of an arc.

      Tabu$time = TimeTabu

      if (verbose)
        cat("\n ............ DONE. \n")

      Tabu$adjacency = amat(Tabu)
      dimnames(Tabu$adjacency) <- list(dnames, dnames)

      if (!is.null(adjacency)) {
      Tabu$Performance <- unlist(MRGNgeneral::RecallPrecision(g1 = as (adjacency[indexT, indexT], "graphNEL"),
                                                              g2 = as (Tabu$adjacency[indexT, indexT], "graphNEL"),
                                                              GV = 0,
                                                              includeGV = FALSE)[-1])[-2]

      if (verbose)
        print(rbind(c(Tabu$Performance, Tabu$time[3])))
      }
      else {
        Tabu$Performance <- Performance = rep(NA, 5)
      }

      if (savetopath){
        save(Tabu, file = paste0(path, setpath, "Tabu", ".RData"))
      }
    }

    ## Run HC -- Score-based
    if('hc' %in% bn.methods) {
      if (verbose)
        cat(paste0("\n ************ HC ... \n"))
      TimeHC = system.time({
        HC =  hc (x = data,
                  whitelist = whitelist1,
                  blacklist = blacklist1,
                  restart = restart,
                  perturb = 1)
      })

      HC$time = TimeHC

      if (verbose)
        cat("\n ............ DONE. \n")

      HC$adjacency = amat(HC)
      dimnames(HC$adjacency ) <- list(dnames, dnames)

      if (!is.null(adjacency)) {
      HC$Performance <- unlist(MRGNgeneral::RecallPrecision(g1 = as (adjacency[indexT, indexT], "graphNEL"),
                                                            g2 = as (HC$adjacency[indexT, indexT], "graphNEL"),
                                                            GV = 0,
                                                            includeGV = FALSE)[-1])[-2]
      if (verbose)
        print(rbind(c(HC$Performance, HC$time[3])))
      }
      else {
        HC$Performance <- Performance = rep(NA, 5)
      }

      if (savetopath){
        save(HC, file = paste0(path, setpath, "HC", restart, ".RData"))
      }


    }

    ## Run PC (stable) -- Constraint-based
    if('pc.stable' %in% bn.methods) {
      if (verbose)
        cat(paste0("\n ************ PC (stable) ... \n"))
      mycl <- if (nb.cl > 0) parallel::makeCluster(nb.cl) else NULL
      TimePCstable = system.time({
        PCstable =  pc.stable (x = data,
                               cluster = mycl,
                               whitelist = whitelist1,
                               blacklist = blacklist1)
      })

      PCstable$time <- TimePCstable

      if (verbose)
        cat("\n ............ DONE. \n")

      MRGNgeneral:::catch.conditions({
        parallel::stopCluster(cl = mycl)
      })

      PCstable$adjacency <- amat(PCstable)
      dimnames(PCstable$adjacency) <- list(dnames, dnames)

      if (!is.null(adjacency)) {
      PCstable$Performance <- unlist(MRGNgeneral::RecallPrecision(g1 = as (adjacency[indexT, indexT], "graphNEL"),
                                                                  g2 = as (PCstable$adjacency[indexT, indexT], "graphNEL"),
                                                                  GV = 0,
                                                                  includeGV = FALSE)[-1])[-2]
      if (verbose)
        print(rbind(c(PCstable$Performance, PCstable$time[3])))
      }
      else {
        PCstable$Performance <- Performance = rep(NA, 5)
      }

      if (savetopath){
        save(PCstable, file = paste0(path, setpath, "PCstable", ".RData"))
      }

    }

    ## Max-Min Hill-Climbing -- Hybrid
    if (verbose)
      cat(paste0("\n ************ Max-Min Hill-Climbing ... \n"))
    if('mmhc' %in% bn.methods) {
      TimeMMHC = system.time({
        MMHC =  mmhc (x = data,
                      whitelist = whitelist1,
                      blacklist = blacklist1,
                      restrict.args = list(alpha = 0.01))
      })
      # score-based algorithms do not support whitelisting both directions of an arc
      MMHC$time = TimeMMHC

      if (verbose)
        cat("\n ............ DONE. \n")

      MMHC$adjacency <- amat(MMHC)
      dimnames(MMHC$adjacency) <- list(dnames, dnames)

      if (!is.null(adjacency)) {
      MMHC$Performance <- unlist(MRGNgeneral::RecallPrecision(g1 = as (adjacency[indexT, indexT], "graphNEL"),
                                                              g2 = as (MMHC$adjacency[indexT, indexT], "graphNEL"),
                                                              GV = 0,
                                                              includeGV = FALSE)[-1])[-2]

      if (verbose)
        print(rbind(c(MMHC$Performance, MMHC$time[3])))
      }
      else {
        MMHC$Performance <- Performance = rep(NA, 5)
      }

      if (savetopath){
        save(MMHC, file = paste0(path, setpath, "MMHC", ".RData"))
      }

    }
  }

  return(list(TABU = Tabu,
              HC = HC,
              PCSTABLE = PCstable,
              MMHC = MMHC,
              whitelist = whitelist1,
              blacklist = blacklist1))
}
