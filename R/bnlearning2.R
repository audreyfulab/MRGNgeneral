# Run a few bnlearn network inference algorithms
# Internal routine for 'run.methods'
bnlearning2 <- function (data, # data.frame/matrix, typically the \code{data} element of an output from \link{sample.graph.data}.
                        dnames = colnames(data), # character vector, names for nodes corresponding to colnames of 'data'
                        bn.methods, # list, bnlearn network inference algorithms, one of: 'none', 'tabu', 'hc', 'pc.stable', 'mmhc', 'fast.iamb', 'mmpc', 'hpc'
                        whitelist0, # two-column matrix, each row indicates an edge direction that must be present in the final inferred network
                        blacklist0, # two-column matrix, each row indicates an edge direction that must be absent in the final inferred network
                        restart,    # an integer, the number of random restarts.
                        adjacency = NULL,  # adjacency matrix of the true network from which \code{data} was generated
                        indexT = 1:NCOL(data),     # integer vector, indices of the columns of \code{data} corresponding to genes/phenotypes of interest and on which performance measures are computed.
                        verbose,
                        nb.cl, # Only used by pc.stable and other constraint-based methods
                        savetopath,     # logical, save the output of each method, and the final output?
                        path, setpath) {
  
  # Initialize outputs
  Tabu <- HC <- PCstable <- MMHC <- FastIAMB <- MMPC <- HPC <- list(time = rep(NA, 3),
                                                                    Performance = rep(NA, 5))
  
  # Prepare white/black lists
  whitelist1 = cbind(paste0('G', whitelist0[,1]), paste0('G', whitelist0[,2]))
  blacklist1 = cbind(paste0('G', blacklist0[,1]), paste0('G', blacklist0[,2]))
  
  if(any(c('tabu', 'hc', 'pc.stable', 'mmhc', 'fast.iamb', 'mmpc', 'hpc') %in% bn.methods)) {
    if (!requireNamespace("bnlearn", quietly = TRUE)) {
      stop("Package 'bnlearn' is needed. Please install it.")
    }
    
    ## Tabu Search -- Score-based (modified Hill-Climbing to escape local optima)
    if('tabu' %in% bn.methods) {
      if (verbose)
        cat(paste0("\n ************ Tabu Search ... \n"))
      TimeTabu = system.time({
        Tabu =  tabu (x = data,
                      whitelist = whitelist1,
                      blacklist = blacklist1)
      })
      
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
      if (nb.cl > 0) {
        setup_cluster(mycl, packages = c("bnlearn"), seed = NULL)
      }
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
    if('mmhc' %in% bn.methods) {
      if (verbose)
        cat(paste0("\n ************ Max-Min Hill-Climbing ... \n"))
      TimeMMHC = system.time({
        MMHC =  mmhc (x = data,
                      whitelist = whitelist1,
                      blacklist = blacklist1,
                      restrict.args = list(alpha = 0.01))
      })
      
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
    
    ## Fast IAMB -- Constraint-based
    if('fast.iamb' %in% bn.methods) {
      if (verbose)
        cat(paste0("\n ************ Fast IAMB ... \n"))
      mycl <- if (nb.cl > 0) parallel::makeCluster(nb.cl) else NULL
      if (nb.cl > 0) {
        setup_cluster(mycl, packages = c("bnlearn"), seed = NULL)
      }
      TimeFastIAMB = system.time({
        FastIAMB =  fast.iamb (x = data,
                               cluster = mycl,
                               whitelist = whitelist1,
                               blacklist = blacklist1,
                               alpha = 0.05)
      })
      
      FastIAMB$time <- TimeFastIAMB
      
      if (verbose)
        cat("\n ............ DONE. \n")
      
      MRGNgeneral:::catch.conditions({
        if (nb.cl > 0) parallel::stopCluster(cl = mycl)
      })
      
      FastIAMB$adjacency <- amat(FastIAMB)
      dimnames(FastIAMB$adjacency) <- list(dnames, dnames)
      
      if (!is.null(adjacency)) {
        FastIAMB$Performance <- unlist(MRGNgeneral::RecallPrecision(g1 = as (adjacency[indexT, indexT], "graphNEL"),
                                                                    g2 = as (FastIAMB$adjacency[indexT, indexT], "graphNEL"),
                                                                    GV = 0,
                                                                    includeGV = FALSE)[-1])[-2]
        if (verbose)
          print(rbind(c(FastIAMB$Performance, FastIAMB$time[3])))
      }
      else {
        FastIAMB$Performance <- Performance = rep(NA, 5)
      }
      
      if (savetopath){
        save(FastIAMB, file = paste0(path, setpath, "FastIAMB", ".RData"))
      }
    }
    
    ## MMPC -- Constraint-based
    if('mmpc' %in% bn.methods) {
      if (verbose)
        cat(paste0("\n ************ MMPC ... \n"))
      mycl <- if (nb.cl > 0) parallel::makeCluster(nb.cl) else NULL
      if (nb.cl > 0) {
        setup_cluster(mycl, packages = c("bnlearn"), seed = NULL)
      }
      TimeMMPC = system.time({
        MMPC =  mmpc (x = data,
                      cluster = mycl,
                      whitelist = whitelist1,
                      blacklist = blacklist1,
                      alpha = 0.05,
                      undirected = TRUE)
      })
      
      MMPC$time <- TimeMMPC
      
      if (verbose)
        cat("\n ............ DONE. \n")
      
      MRGNgeneral:::catch.conditions({
        if (nb.cl > 0) parallel::stopCluster(cl = mycl)
      })
      
      MMPC$adjacency <- amat(MMPC)
      dimnames(MMPC$adjacency) <- list(dnames, dnames)
      
      if (!is.null(adjacency)) {
        MMPC$Performance <- unlist(MRGNgeneral::RecallPrecision(g1 = as (adjacency[indexT, indexT], "graphNEL"),
                                                                g2 = as (MMPC$adjacency[indexT, indexT], "graphNEL"),
                                                                GV = 0,
                                                                includeGV = FALSE)[-1])[-2]
        if (verbose)
          print(rbind(c(MMPC$Performance, MMPC$time[3])))
      }
      else {
        MMPC$Performance <- Performance = rep(NA, 5)
      }
      
      if (savetopath){
        save(MMPC, file = paste0(path, setpath, "MMPC", ".RData"))
      }
    }
    
    ## HPC -- Constraint-based
    if('hpc' %in% bn.methods) {
      if (verbose)
        cat(paste0("\n ************ HPC ... \n"))
      mycl <- if (nb.cl > 0) parallel::makeCluster(nb.cl) else NULL
      if (nb.cl > 0) {
        setup_cluster(mycl, packages = c("bnlearn"), seed = NULL)
      }
      TimeHPC = system.time({
        HPC =  hpc (x = data,
                    cluster = mycl,
                    whitelist = whitelist1,
                    blacklist = blacklist1,
                    alpha = 0.05,
                    undirected = TRUE)
      })
      
      HPC$time <- TimeHPC
      
      if (verbose)
        cat("\n ............ DONE. \n")
      
      MRGNgeneral:::catch.conditions({
        if (nb.cl > 0) parallel::stopCluster(cl = mycl)
      })
      
      HPC$adjacency <- amat(HPC)
      dimnames(HPC$adjacency) <- list(dnames, dnames)
      
      if (!is.null(adjacency)) {
        HPC$Performance <- unlist(MRGNgeneral::RecallPrecision(g1 = as (adjacency[indexT, indexT], "graphNEL"),
                                                               g2 = as (HPC$adjacency[indexT, indexT], "graphNEL"),
                                                               GV = 0,
                                                               includeGV = FALSE)[-1])[-2]
        if (verbose)
          print(rbind(c(HPC$Performance, HPC$time[3])))
      }
      else {
        HPC$Performance <- Performance = rep(NA, 5)
      }
      
      if (savetopath){
        save(HPC, file = paste0(path, setpath, "HPC", ".RData"))
      }
    }
  }
  
  return(list(TABU = Tabu,
              HC = HC,
              PCSTABLE = PCstable,
              MMHC = MMHC,
              FASTIAMB = FastIAMB,
              MMPC = MMPC,
              HPC = HPC,
              whitelist = whitelist1,
              blacklist = blacklist1))
}