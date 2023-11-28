
# Keep the function simple for now
# Assume we know (as input) the set of confounding variables for each T node
# Know the type of each confounding variable (actual confounder, or one of intermediate or common child variable)
# Dealing with missing values ???

#' @name MRGN
#'
#' @title Infer a causal network using the MRGN algorithm
#'
#' @description Infer a causal network/graph with directed and undirected edges from observational data.
#'
#' @param data \code{data.frame}
#'
#' @param Qlabels labels indicating columns of \code{data} which are
#' \code{Q}-nodes (common children or intermediate variables).
#'
#' @param confounders a list of length \code{p} giving the set of known confounders
#' for each phenotype, as returned by \link{get.conf.sets} with the option
#' \code{return.for.trios = FALSE}.
#'
#' @param Adj input adjacency matrix (binary)
#'
#' @param stringent logical, should edges not appearing in \code{'M1'} and
#' \code{'M2'} trio structures be dropped during adjacency matrix update?
#' Passed to \link{update.adjacency.matrix}.
#'
#'
#' @seealso \link{infer.trio} for inferring edges in small networks of one
#'  genetic variant and only two genes.

#' @export MRGNgeneral
#' @import parallel

# We have four blocks of columns in data: V, T, (I&C) and confounders

# Post filtering based on marginal tests ???


# Keep track of all added edges (directions) so that
# at each iteration, I can check for conflicts with the current adds
# but also previous adds

# Use one 'control' argument, a list of (all?) control arguments (which ones?)

# To do: make dropped.edges and added.edges lists, so as to save the edges dropped/added at each iteration

MRGN <- function (data, # input n-by-m data matrix: 'q' Variants, 'p' Phenotypes, r Intermediate&Common_Child and 'u' confounders
                         scale.data = TRUE,
                         Qlabels = NULL, # column number or names of W,Z-nodes in data
                         p, q, r = length(Qlabels),
                         u = m - p - q - r,
                         confounders = list(), # a list of length p giving the set of confounders for each phenotype
                         Adj, # a (p+q+r)-by-(p+q+r) adjacency matrix of an undirected acyclic graph (zeros on diagonal).
                         is.CNA = FALSE, # logical indicating if a genetic variant is a copy number alteration.
                         alpha = 0.01, # in (0, .5) type I error rate for individual tests, in the open (0, .5)
                         use.perm = TRUE, # See infer.trio
                         gamma = 0.05, #
                         nperms = 10000, #
                         FDRcontrol = c("bonferroni", "qvalue", "none"),
                         fdr = 0.05, # in (0, .5) False discovery rate for all trio analyses at one iteration, in the open (0, .5)
                         lambda = NULL, # The value of the tuning parameter for q-value estimation
                         lambda.step = 0.05, # in (0, .1); Used to initialize "lambda" when lambda = NULL
                         pi0.meth = c("bootstrap", "smoother"),
                         stringent = TRUE,
                         analyse.triplets = TRUE,
                         solve.conflicts = TRUE, # Just used for testing, to be removed?
                         method = "naive",
                         maxiter = NULL,
                         parallel = FALSE,
                         cl = parallel:::getDefaultCluster(),
                         chunk.size = NULL, # scalar number; number of invocations of fun or FUN in one chunk; a chunk is a unit for scheduling.
                         verbose = 0L) {
  # ============================================================================
  # Save the call to MRGN
  # ----------------------------------------------------------------------------
  mrgncl <- match.call()

  # ============================================================================
  # Parallelize computations or not?
  # ----------------------------------------------------------------------------
  stopifnot(is.logical(parallel))
  if (parallel[1]) {
    if (is.null(cl)) {
      # Use option cl.cores to choose an appropriate cluster size
      cl <- parallel::makeCluster(getOption("cl.cores", 2L))
    }
    if (!is.null(chunk.size)) {
      stopifnot(is.numeric(chunk.size))
      chunk.size <- chunk.size[1]
      stopifnot(chunk.size > 0)
      chunk.size <- ceiling(chunk.size)
    }
  }
  else {
    cl <- NULL
  }

  # ============================================================================
  # Check arguments & Pre-process some inputs
  # ----------------------------------------------------------------------------
  eval(check.mrgn.args())

  # ============================================================================
  ### Save the current Adjacency matrix
  old.Adj <- Adj

  # Initialize some outputs as empty vector/lists
  trio.analysis <- NULL
  dropped.edges <- added.edges <- NULL # list()
  Qtrio.set <- Qtrio.analysis <- NULL
  triplet.set <- triplet.analysis <- NULL
  Qtriplet.set <- Qtriplet.analysis <- NULL
  iter.nb.triplets <- iter.nb.Qtriplets <- NULL

  # Initialize iteration-wise counters/indicators
  new.edges <- NULL # indicator of trios with new added edges at an iteration (length = NROW(trio.set))
  new_edge_added <- FALSE # scalar indicator of at least one new edge added # new_edge_added <- any(new.edges)
  nb.new_edges <- 0 # total number of new edges added into the network at an iteration
  iter <- 1

  ### Step 3: List trios involving genetic variants and direct related edges
  # ----------------------------------------------------------------------------
  ## Step 3.1: Identify target trios, i.e. generate trio list & Bind trios identifiers as rows
  if (verbose) {
    cat("\n      # Generating an initial list of trios with genetic variants. \n")
  }
  trio.set <- enumerate.trios(i = 1:q, Adj = Adj,
                              p = p, q = q,
                              r = r,
                              VTT = TRUE,
                              cl = cl,
                              chunk.size = chunk.size)
  trio.set <- do.call('rbind', trio.set)

  # Number of initial trios
  nb.trios <- NROW(trio.set)

  if (verbose) {
    cat(paste0("          * ",
               if (nb.trios == 0) "zero trio"
               else if (nb.trios == 1) "one [Vi,Tj,Tk] trio"
               else paste0(nb.trios, " [Vi,Tj,Tk] trios"),
               " generated. \n"))
  }

  ## Step 3.2-3.3: Determine the directed structure for each trio and update Adj
  if (nb.trios > 0) {
    ## Step 3.2: Determine the directed structure for each trio
    if (verbose)
      cat("      # Inferring the causal network for each [Vi,Tj,Tk] trio. \n")
    trio.analysis <- analyse.trio.set (trio.set, nb.trios = nb.trios, Qlabels = Qlabels,
                                       alpha = alpha, FDRcontrol = FDRcontrol,
                                       fdr = fdr, lambda = lambda, lambda.step = lambda.step, pi0.meth = pi0.meth,
                                       data = data, confounders = confounders, p = p, q = q,
                                       use.perm = use.perm, gamma = gamma, is.CNA = is.CNA,
                                       nperms = nperms, verbose = verbose,
                                       cl = cl, chunk.size = chunk.size)

    ## Step 3.3: Update the adjacency matrix
    if (verbose) {
      cat(paste0("      # Updating the adjacency matrix: adding edge ",
                 if (nb.trios == 1) "direction for one [Vi,Tj,Tk] trio. \n"
                 else "directions for ", nb.trios, " [Vi,Tj,Tk] trios. \n"))
    }
    update.object <- update.adjacency.matrix (Adj = Adj, p = p, q = q, r = r,
                                              trio.set = trio.set,
                                              inferred.models = trio.analysis$Inferred.Model,
                                              stringent = stringent, add.edges = TRUE,
                                              solve.conflicts = solve.conflicts,
                                              method = method, cl = cl, chunk.size = chunk.size)
    # Update dropped.edges and added.edges
    added.edges <- update.object$added.edges
    dropped.edges <- update.object$dropped.edges

    # Report conflicts, if any
    if (verbose & solve.conflicts) {
      if (update.object$nb.conflicts)
        cat(paste0("          * ",
                   if (update.object$nb.conflicts == 1) "one conflict detected"
                   else paste0(update.object$nb.conflicts, " conflicts detected"),
                   " and handled with 'method = ", method, "'. \n"))
    }

    # Logical vector indicating trios in which new edges were added
    nb.new_edges <- update.object$nb.new_edges
    new.edges <- (update.object$new.edges == 1) & (nb.new_edges > 0)

    # Check if new edges were added from [Vi,Tj,Tk] trio analysis
    new_edge_added <- any(new.edges)
    if (verbose & new_edge_added) {
      cat(paste0("          * ",
                 if (nb.new_edges == 1) "one new edge"
                 else paste0(nb.new_edges, " new edges"),
                 " added into the network. \n"))
    }

    # Save the new Adj
    Adj <- update.object$Adj
  }

  ## Step 3.4: Generating an initial list of trios involving a genetic variant and a Q-node
  if (r > 0) {
    if (verbose) {
      cat("\n      # Generating an initial list of [Vi,Tj,Qk] trios. \n")
    }
    Qtrio.set <- enumerate.trios(i = 1:q,
                                 Adj = Adj,
                                 p = p, q = q,
                                 r = r,
                                 VTT = FALSE,
                                 cl = cl, chunk.size = chunk.size)
    Qtrio.set <- do.call('rbind', Qtrio.set)

    # Number of initial trios
    nb.Qtrios <- NROW(Qtrio.set)

    if (verbose) {
      cat(paste0("          * ",
                 if (nb.Qtrios == 0) "zero trio"
                 else if (nb.Qtrios == 1) "one [Vi,Tj,Qk] trio"
                 else paste0(nb.Qtrios, " [Vi,Tj,Qk] trios"),
                 " generated. \n"))
    }

    if (nb.Qtrios > 0) {
      # Step 3.5: Determine the directed structure for each trio
      if (verbose)
        cat("      # Inferring the causal network for each trio. \n")
      Qtrio.analysis <- analyse.trio.set (Qtrio.set, nb.trios = nb.Qtrios, Qlabels = Qlabels,
                                          alpha = alpha, FDRcontrol = FDRcontrol,
                                          fdr = fdr, lambda = lambda, lambda.step = lambda.step, pi0.meth = pi0.meth,
                                          data = data, confounders = confounders, p = p, q = q,
                                          use.perm = use.perm, gamma = gamma, is.CNA = is.CNA,
                                          nperms = nperms, verbose = verbose,
                                          cl = cl, chunk.size = chunk.size)

      ## Step 3.6: Update the adjacency matrix
      if (verbose) {
        cat(paste0("      # Updating the adjacency matrix: adding edge ",
                   if (nb.Qtrios == 1) "direction for one [Vi,Tj,Qk] trio. \n"
                   else "directions for ", nb.Qtrios, " [Vi,Tj,Qk] trios. \n"))
      }
      Qupdate.object <- update.adjacency.matrix (Adj = Adj,
                                                 p = p, q = q, r = r,
                                                 trio.set = Qtrio.set,
                                                 inferred.models = Qtrio.analysis$Inferred.Model,
                                                 stringent = stringent, add.edges = FALSE,
                                                 solve.conflicts = solve.conflicts, method = method,
                                                 added.edges = added.edges, dropped.edges = dropped.edges,
                                                 cl = cl, chunk.size = chunk.size)
      dropped.edges <- Qupdate.object$dropped.edges

      # Update Adj
      Adj <- Qupdate.object$Adj
    }
  }
  # ----------------------------------------------------------------------------

  ## Iterate Step 3 each time new edges are added (and the network is not fully directed)
  # Indicator of the presence of undirected edges in the updated network
  Not.fully.directed <- any(rowSums((Adj + t(Adj)) == 2) > 0)

  # Number of [Vi,Tj,Tk] trios generated at each iteration
  iter.nb.trios <- nb.trios

  # Number of [Vi,Tj,Qk] trios generated at each iteration
  iter.nb.Qtrios <- if (r > 0) nb.Qtrios

  # Initialize the number of new trios to nb.trios
  nb.new.trios <- nb.trios

  while ((new_edge_added & Not.fully.directed) & iter < maxiter) {
    ## Step 3.1: Form additional trios
    if (verbose) {
      cat("\n      # Generating additional [Vi,Tj,Tk] trios. \n")
    }
    new.trio.set <- enumerate.new.trios (trio.set = trio.set,
                                         new.edges = c(rep(FALSE, nb.trios - nb.new.trios),
                                                       new.edges),
                                         old.Adj = old.Adj, Adj = Adj,
                                         p = p, q = q, r = r,
                                         VTT = TRUE,
                                         cl = cl, chunk.size = chunk.size)

    # Stop if no new trio generated
    nb.new.trios <- NROW(new.trio.set)
    if (nb.new.trios == 0) {
      new_edge_added <- FALSE
      if (verbose) {
        cat("          * no additional [Vi,Tj,Tk] trio found.")
      }

      break
    }

    if (verbose) {
      cat(paste0("          * ",
                 if (nb.new.trios == 1) "one new [Vi,Tj,Tk] trio"
                 else paste0(nb.new.trios, " new [Vi,Tj,Tk] trios"),
                 " generated. \n"))
    }

    # Save the current Adjacency matrix (at iter)
    old.Adj <- Adj

    # iteration counter
    iter <- iter + 1

    ## Step 3.2: Determine the directed structure for each trio
    if (verbose) {
      cat(paste0("      # Inferring the causal network for each new [Vi,Tj,Tk] trio (iteration ",
                 iter, "). \n"))
    }
    new.trio.analysis <- analyse.trio.set (new.trio.set, nb.trios = nb.new.trios, Qlabels = Qlabels,
                                           alpha = alpha, FDRcontrol = FDRcontrol,
                                           fdr = fdr, lambda = lambda,
                                           lambda.step = lambda.step, pi0.meth = pi0.meth,
                                           data = data, confounders = confounders, p = p, q = q,
                                           use.perm = use.perm, gamma = gamma, is.CNA = is.CNA,
                                           nperms = nperms, verbose = verbose,
                                           cl = cl, chunk.size = chunk.size)

    ## Step 3.3: Update the adjacency matrix
    if (verbose) {
      cat(paste0("      # Updating the adjacency matrix: adding edge ",
                 if (nb.new.trios == 1) "direction for one [Vi,Tj,Tk] trio. \n"
                 else "directions for ", nb.new.trios, " [Vi,Tj,Tk] trios. \n"))
    }
    update.object <- update.adjacency.matrix (Adj = Adj,
                                              p = p, q = q, r = r,
                                              trio.set = new.trio.set,
                                              inferred.models = new.trio.analysis$Inferred.Model,
                                              stringent = stringent, add.edges = TRUE,
                                              solve.conflicts = solve.conflicts, method = method,
                                              added.edges = added.edges,
                                              dropped.edges = dropped.edges,
                                              cl = cl, chunk.size = chunk.size)

    # Update dropped.edges and added.edges
    added.edges <- update.object$added.edges
    dropped.edges <- update.object$dropped.edges

    # Report conflicts, if any
    if (verbose & solve.conflicts) {
      if (update.object$nb.conflicts)
        cat(paste0("          * ",
                   if (update.object$nb.conflicts == 1) "one conflict detected"
                   else paste0(update.object$nb.conflicts, " conflicts detected"),
                   " and handled with 'method = ", method, "'. \n"))
    }

    # Logical vector indicating trios in which new edges were added
    nb.new_edges <- update.object$nb.new_edges
    new.edges <- (update.object$new.edges == 1) & (nb.new_edges > 0)

    # Check if new edges were added from [V,T1,T2] trio analysis
    new_edge_added <- any(new.edges)
    if (verbose & new_edge_added) {
      cat(paste0("          * ",
                 if (nb.new_edges == 1) "one new edge"
                 else paste0(nb.new_edges, " new edges"),
                 " added into the network. \n"))
    }

    # Save the new Adj
    Adj <- update.object$Adj

    ## Step 3.4: Generating a list of trios involving a genetic variant and a Q-node
    if (r > 0) {
      if (verbose) {
        cat("\n      # Generating additional [Vi,Tj,Qk] trios. \n")
      }

      new.Qtrio.set <- enumerate.new.trios (trio.set = new.trio.set,
                                            new.edges = new.edges,
                                            old.Adj = old.Adj, Adj = Adj,
                                            p = p, q = q, r = r,
                                            VTT = FALSE,
                                            Qtrio.set = Qtrio.set,
                                            cl = cl, chunk.size = chunk.size)

      # Number of initial trios
      nb.new.Qtrios <- NROW(new.Qtrio.set)

      if (verbose) {
        cat(paste0("          * ",
                   if (nb.new.Qtrios == 0) "zero new trio"
                   else if (nb.new.Qtrios == 1) "one new [Vi,Tj,Qk] trio"
                   else paste0(nb.new.Qtrios, " new [Vi,Tj,Qk] trios"),
                   " generated. \n"))
      }

      if (nb.new.Qtrios > 0) {
        # Step 3.5: Determine the directed structure for each trio
        if (verbose)
          cat("      # Inferring the causal network for each trio. \n")
        new.Qtrio.analysis <- analyse.trio.set (new.Qtrio.set, nb.trios = nb.new.Qtrios, Qlabels = Qlabels,
                                                alpha = alpha, FDRcontrol = FDRcontrol,
                                                fdr = fdr, lambda = lambda, lambda.step = lambda.step, pi0.meth = pi0.meth,
                                                data = data, confounders = confounders, p = p, q = q,
                                                use.perm = use.perm, gamma = gamma, is.CNA = is.CNA,
                                                nperms = nperms, verbose = verbose,
                                                cl = cl, chunk.size = chunk.size)

        ## Step 3.6: Update the adjacency matrix
        if (verbose) {
          cat(paste0("      # Updating the adjacency matrix: adding edge ",
                     if (nb.new.Qtrios == 1) "direction for one [Vi,Tj,Qk] trio. \n"
                     else "directions for ", nb.new.Qtrios, " [Vi,Tj,Qk] trios. \n"))
        }
        Qupdate.object <- update.adjacency.matrix (Adj = Adj,
                                                   p = p, q = q, r = r,
                                                   trio.set = new.Qtrio.set,
                                                   inferred.models = new.Qtrio.analysis$Inferred.Model,
                                                   stringent = stringent, add.edges = FALSE,
                                                   solve.conflicts = solve.conflicts,
                                                   method = method, dropped.edges = dropped.edges,
                                                   cl = cl, chunk.size = chunk.size)
        dropped.edges <- Qupdate.object$dropped.edges

        # Update Adj
        Adj <- Qupdate.object$Adj

        # Save the additional trios in 'Qtrio.set'
        Qtrio.set <- rbind(Qtrio.set, new.Qtrio.set)

        # Save the additional trios analysis result in 'Qtrio.analysis'
        Qtrio.analysis <- rbind(Qtrio.analysis, new.Qtrio.analysis) # (not useful at this point)
      }

      # Save the number of additional trios generated
      iter.nb.Qtrios <- c(iter.nb.Qtrios, nb.new.Qtrios)
    }

    # Save the number of additional trios generated
    iter.nb.trios <- c(iter.nb.trios, nb.new.trios)

    # Save the additional trios in 'trio.set'
    trio.set <- rbind(trio.set, new.trio.set)

    # Save the additional trios analysis result in 'trio.analysis'
    trio.analysis <- rbind(trio.analysis, new.trio.analysis) # (not useful at this point)

    # Update the total number of trios analysed so far
    nb.trios <- nb.trios + nb.new.trios

    # Indicator of the presence of undirected edges in the updated network
    Not.fully.directed <- any(rowSums((Adj + t(Adj)) == 2) > 0)
  }

  ### Step 4: List trios involving no genetic variant and direct related edges
  # ----------------------------------------------------------------------------

  if (all(c(Not.fully.directed, analyse.triplets, p + r >= 3))) { # p >= 3 means at least one triplet can be formed
    ## Step 4.1: Form [Ti,Tj,Tk]  trios
    if (verbose) {
      cat("\n      # Generating an initial list of triplets involving only T-nodes. \n")
    }
    triplet.set <- enumerate.triplets (Adj,
                                       p = p,
                                       q = q,
                                       r = r,
                                       TTT = TRUE,
                                       cl = cl,
                                       chunk.size = chunk.size)

    # Number of initial triplets
    nb.triplets <- length(triplet.set$types)
    if (!nb.triplets) {
      if (verbose) {
        cat("          * no updateable [Ti,Tj,Tk] trio found. \n")
      }
      triplet.set <- triplet.analysis <- NULL
    }
    else {
      triplet.set <- cbind(triplet.set$triplets, type = triplet.set$types)

      if (verbose) {
        cat(paste0("          * ",
                   if (nb.triplets == 1) "one updateable [Ti,Tj,Tk] trio"
                   else paste0(nb.triplets, " updateable  [Ti,Tj,Tk] trios"),
                   " generated. \n"))
      }

      # Step 4.2: Determine the directed structure for each triplet
      if (verbose)
        cat("      # Inferring the causal network for each [Ti,Tj,Tk] trio. \n")
      triplet.analysis <- analyse.triplet.set (triplet.set, nb.triplets = nb.triplets,
                                               alpha = alpha, FDRcontrol = FDRcontrol, fdr = fdr,
                                               lambda = lambda, lambda.step = lambda.step,
                                               pi0.meth = pi0.meth,
                                               data = data, confounders = confounders,
                                               p = p, q = q, verbose = verbose,
                                               cl = cl, chunk.size = chunk.size)

      # Step 4.3: Update the adjacency matrix
      if (verbose) {
        cat(if (nb.triplets == 1) paste0("      # Updating the adjacency matrix: adding edge ",
                                         "direction for one [Ti,Tj,Tk] trio. \n")
            else paste0("      # Updating the adjacency matrix: adding edge directions for ",
                        nb.triplets, " [Ti,Tj,Tk] trios. \n"))
      }
      update.object <- update.adjacency.matrix (Adj = Adj,
                                                q = q,
                                                trio.set = triplet.set,
                                                inferred.models = triplet.analysis$Inferred.Model,
                                                stringent = stringent,
                                                add.edges = FALSE,
                                                solve.conflicts = solve.conflicts,
                                                method = method,
                                                added.edges = added.edges,
                                                dropped.edges = dropped.edges,
                                                cl = cl, chunk.size = chunk.size)
      dropped.edges <- update.object$dropped.edges

      # Update Adj
      Adj <- update.object$Adj
    }

    # Step 4.4: Form trio with T-nodes and Q-nodes
    if (r > 0) {
      if (verbose) {
        cat("\n      # Generating an initial list of [Ti,Tj,Qk] trios. \n")
      }
      Qtriplet.set <- enumerate.triplets(Adj = Adj,
                                         p = p, q = q,
                                         r = r,
                                         TTT = FALSE,
                                         cl = cl, chunk.size = chunk.size)

      # Number of initial trios
      nb.Qtriplets <- length(Qtriplet.set$types)
      if (!nb.Qtriplets) {
        if (verbose) {
          cat("          * no updateable [Ti,Tj,Qk] trio found. \n")
        }
        Qtriplet.set <- Qtriplet.analysis <- NULL
      }
      else {
        Qtriplet.set <- cbind(Qtriplet.set$triplets, type = Qtriplet.set$types)

        if (verbose) {
          cat(paste0("          * ",
                     if (nb.Qtriplets == 1) "one updateable [Ti,Tj,Tk] trio"
                     else paste0(nb.Qtriplets, " updateable  [Ti,Tj,Tk] trios"),
                     " generated. \n"))
        }

        # Step 4.5: Determine the directed structure for each trio
        if (verbose)
          cat("      # Inferring the causal network for each [Ti,Tj,Qk] trio. \n")
        Qtriplet.analysis <- analyse.triplet.set (Qtriplet.set, nb.triplets = nb.Qtriplets,
                                                  alpha = alpha, FDRcontrol = FDRcontrol, fdr = fdr,
                                                  lambda = lambda, lambda.step = lambda.step,
                                                  pi0.meth = pi0.meth,
                                                  data = data, confounders = confounders,
                                                  p = p, q = q, verbose = verbose,
                                                  cl = cl, chunk.size = chunk.size)

        # Step 4.3: Update the adjacency matrix
        if (verbose) {
          cat(if (nb.Qtriplets == 1) paste0("      # Updating the adjacency matrix: adding edge ",
                                            "direction for one [Ti,Tj,Qk] trio. \n")
              else paste0("      # Updating the adjacency matrix: adding edge directions for ",
                          nb.Qtriplets, " [Ti,Tj,Qk] trios. \n"))
        }
        Qupdate.object <- update.adjacency.matrix (Adj = Adj,
                                                   q = q,
                                                   trio.set = Qtriplet.set,
                                                   inferred.models = Qtriplet.analysis$Inferred.Model,
                                                   stringent = stringent,
                                                   add.edges = FALSE,
                                                   solve.conflicts = solve.conflicts,
                                                   method = method,
                                                   added.edges = added.edges,
                                                   dropped.edges = dropped.edges,
                                                   cl = cl, chunk.size = chunk.size)
        dropped.edges <- Qupdate.object$dropped.edges

        # Update Adj
        Adj <- Qupdate.object$Adj
      }
    }

    ## Iterate Step 4
    # Indicator of the presence of undirected edges in the updated network
    Not.fully.directed <- any(rowSums((Adj + t(Adj)) == 2) > 0)

    # Repeat step 4: re-evaluate the type of each triplet with two edges, and one undirected.
    # Re-interpret the result of triplet analysis, and update Adj.
    # Number of [Vi,Tj,Tk] trios generated at each iteration
    iter.nb.triplets <- nb.triplets
    new.triplet.set <- triplet.set
    new.triplet.analysis <- triplet.analysis
    other.triplets <- if (nb.triplets) new.triplet.analysis$Inferred.Model == "Other" else FALSE
    new_directions <- TRUE
    if (r > 0) {
      # Number of [Vi,Tj,Qk] trios generated at each iteration
      iter.nb.Qtriplets <- if (r > 0) nb.Qtriplets
      new.Qtriplet.set <- Qtriplet.set
      new.Qtriplet.analysis <- Qtriplet.analysis
      other.Qtriplets <- if (nb.Qtriplets) new.Qtriplet.analysis$Inferred.Model == "Other" else FALSE
    }
    else {
      other.Qtriplets <- FALSE
    }

    Not.stop <- Not.fully.directed & (any(other.Qtriplets) | any(other.triplets)) # & FALSE

    while(Not.stop & iter < maxiter) {

      if (any(other.triplets)) {
        ## Step 4.1: Extract [Ti,Tj,Tk] trios to be updated
        if (verbose) {
          cat("\n      # Identifying updateable [Ti,Tj,Tk] trios. \n")
        }

        # Save old triplet types
        old.types <- new.triplet.set[other.triplets, 4]

        # Find new triplet types based on the updated Adj
        new.triplet.set <- get.new.triplet.type (triplet.set = new.triplet.set[other.triplets,, drop = FALSE],
                                                 Adj = Adj,
                                                 cl = cl,
                                                 chunk.size = chunk.size)
        keep.new.triplets <- (new.triplet.set[,4] <= 2) & (new.triplet.set[,4] != old.types)
        new.triplet.set <- new.triplet.set[keep.new.triplets,,drop = FALSE]

        # Number of initial triplets
        nb.new.triplets <- sum(keep.new.triplets)

        # Stop if no new updateble trio found
        if (nb.new.triplets == 0) {
          new_directions <- FALSE
          if (verbose) {
            cat("          * no additional [Ti,Tj,Tk] trio found.")
          }
          other.triplets <- FALSE
        }

        if (nb.new.triplets > 0) {
          if (verbose) {
            cat(paste0("          * ",
                       if (nb.new.triplets == 1) "one updateable [Ti,Tj,Tk] trio"
                       else paste0(nb.new.triplets, " updateable  [Ti,Tj,Tk] trios"),
                       " found \n"))
          }

          # Step 4.2: Extract previous triplet analysis results for kept triplets
          new.triplet.analysis <- new.triplet.analysis[other.triplets,, drop = FALSE]
          new.triplet.analysis <- new.triplet.analysis[keep.new.triplets,, drop = FALSE]

          # Step 4.3: Update the adjacency matrix
          if (verbose) {
            cat(if (nb.new.triplets == 1) paste0("      # Updating the adjacency matrix: adding edge ",
                                                 "direction for one [Ti,Tj,Tk] trio. \n")
                else paste0("      # Updating the adjacency matrix: adding edge directions for ",
                            nb.new.triplets, " [Ti,Tj,Tk] trios. \n"))
          }
          update.object <- update.adjacency.matrix (Adj = Adj,
                                                    q = q,
                                                    trio.set = new.triplet.set,
                                                    inferred.models = new.triplet.analysis$Inferred.Model,
                                                    stringent = stringent,
                                                    add.edges = FALSE,
                                                    solve.conflicts = solve.conflicts,
                                                    method = method,
                                                    added.edges = added.edges,
                                                    dropped.edges = dropped.edges,
                                                    cl = cl, chunk.size = chunk.size)
          dropped.edges <- update.object$dropped.edges

          # Update Adj
          Adj <- update.object$Adj

          # Indicator of the presence of triplets of type "Other"
          other.triplets <- new.triplet.analysis$Inferred.Model == "Other"

          # Save the number of additional trios generated
          iter.nb.triplets <- c(iter.nb.triplets, nb.new.trios)

          # Save the additional trios in 'trio.set'
          triplet.set <- rbind(triplet.set, new.triplet.set)

          # Save the additional trios analysis result in 'trio.analysis'
          triplet.analysis <- rbind(triplet.analysis, new.triplet.analysis) # (not useful at this point)

        }

      }

      if (r > 0 & any(other.Qtriplets)) {
        ## Step 3.4: Extract [Ti,Tj,Qk] trios to be updated
        if (verbose) {
          cat("\n      # Identifying updateable [Ti,Tj,Qk] trios. \n")
        }

        # Save old triplet types
        old.Qtypes <- new.Qtriplet.set[other.Qtriplets, 4]

        # Find new triplet types based on the updated Adj
        new.Qtriplet.set <- get.new.triplet.type (triplet.set = new.Qtriplet.set[other.Qtriplets,, drop = FALSE],
                                                  Adj = Adj,
                                                  cl = cl,
                                                  chunk.size = chunk.size)
        keep.new.Qtriplets <- (new.Qtriplet.set[,4] <= 2) & (new.Qtriplet.set[,4] != old.Qtypes)
        new.Qtriplet.set <- new.Qtriplet.set[keep.new.Qtriplets,,drop = FALSE]

        # Number of initial triplets
        nb.new.Qtriplets <- sum(keep.new.Qtriplets)

        # Stop if no new updateble trio found
        if (nb.new.Qtriplets == 0) {
          new_directions <- FALSE
          if (verbose) {
            cat("          * no additional [Ti,Tj,Qk] trio found.")
          }

          other.Qtriplets <- FALSE
        }

        if (nb.new.Qtriplets > 0) {
          if (verbose) {
            cat(paste0("          * ",
                       if (nb.new.Qtriplets == 1) "one updateable [Ti,Tj,Qk] trio"
                       else paste0(nb.new.Qtriplets, " updateable  [Ti,Tj,Qk] trios"),
                       " found. \n"))
          }

          # Step 4.2: Extract previous triplet analysis results for kept triplets
          new.Qtriplet.analysis <- new.Qtriplet.analysis[other.Qtriplets,, drop = FALSE]
          new.Qtriplet.analysis <- new.Qtriplet.analysis[keep.new.Qtriplets,, drop = FALSE]

          # Step 4.3: Update the adjacency matrix
          if (verbose) {
            cat(if (nb.new.Qtriplets == 1) paste0("      # Updating the adjacency matrix: adding edge ",
                                                  "direction for one [Ti,Tj,Qk] trio. \n")
                else paste0("      # Updating the adjacency matrix: adding edge directions for ",
                            nb.new.Qtriplets, " [Ti,Tj,Qk] trios. \n"))
          }
          update.object <- update.adjacency.matrix (Adj = Adj,
                                                    q = q,
                                                    trio.set = new.Qtriplet.set,
                                                    inferred.models = new.Qtriplet.analysis$Inferred.Model,
                                                    stringent = stringent,
                                                    add.edges = FALSE,
                                                    solve.conflicts = solve.conflicts,
                                                    method = method,
                                                    added.edges = added.edges,
                                                    dropped.edges = dropped.edges,
                                                    cl = cl, chunk.size = chunk.size)
          dropped.edges <- update.object$dropped.edges

          # Update Adj
          Adj <- update.object$Adj

          # Save the additional trios in 'Qtrio.set'
          Qtriplet.set <- rbind(Qtriplet.set, new.Qtriplet.set)

          # Save the additional trios analysis result in 'Qtrio.analysis'
          Qtriplet.analysis <- rbind(Qtriplet.analysis, new.Qtriplet.analysis) # (not useful at this point)

          # Indicator of the presence of triplets of type "Other"
          other.Qtriplets <- new.Qtriplet.analysis$Inferred.Model == "Other"
        }

        # Save the number of additional trios generated
        iter.nb.Qtriplets <- c(iter.nb.Qtriplets, nb.new.Qtriplets)
      }

      # iteration counter
      iter <- iter + 1

      # Indicator of the presence of undirected edges in the updated network
      Not.fully.directed <- any(rowSums((Adj + t(Adj)) == 2) > 0)

      Not.stop <- Not.fully.directed & any(c(other.triplets, other.Qtriplets))
    }

  }

  # Report on the presence of undirected edges in the final network
  if (verbose) {
    if (Not.fully.directed)
      cat("\n      # the final network is not fully directed. \n")
    else
      cat("\n      # the final network is fully directed. \n")
  }

  # ======================================================================================
  structure(
    list(Adj = structure(Adj, class = 'adjacency.matrix'),
         Adj0 = structure(Adj0, class = 'adjacency.matrix'),
         VTTtrio.set = trio.set,
         iter.nb.VTTtrios = iter.nb.trios,
         VTTtrio.analysis = trio.analysis,
         VTQtrio.set = Qtrio.set,
         iter.nb.VTQtrios = iter.nb.Qtrios,
         VTQtrio.analysis = Qtrio.analysis,
         TTTtrio.set = triplet.set,
         iter.nb.TTTtrios = iter.nb.triplets,
         TTTtrio.analysis = triplet.analysis,
         TTQtrio.set = Qtriplet.set,
         iter.nb.TTQtrios = iter.nb.Qtriplets,
         TTQtrio.analysis = Qtriplet.analysis,
         added.edges = added.edges,
         dropped.edges = dropped.edges,
         Qlabels = Qlabels, niter = iter,
         cl = cl, chunk.size = chunk.size,
         call = mrgncl),
    class = "MRGNg")
}


# A brief print method for MRGN output (class 'MRGNg')
# Issue with the count of edges: not consistent with graphNEL result
# when there are some undirected edges
print.MRGNg <- function (x, TTonly = FALSE, format = "adjacency", ...) {
  p <- eval(x$call$p)
  q <- eval(x$call$q)
  r <- if (TTonly) 0 else eval(x$call$r)
  Adj <- as.matrix(x$Adj[1:(p+q+r), 1:(p+q+r)])

  switch(format,
         graphNEL = {
           print(as(Adj, "graphNEL"), ...)
         },
         igraph = { # igraph object
           print(igraph::graph_from_adjacency_matrix(Adj), ...)
         },
         {
           cat("\nCall:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"),
               "\n\n", sep = "")

           #AugAdj <- Adj + t(Adj)
           #AugAdj[upper.tri(AugAdj)] <- 0

           #DirEdges <- which(AugAdj == 1, arr.ind = TRUE)
           #UndirEdges <- which(AugAdj == 2, arr.ind = TRUE)

           #NROW(DirEdges)
           #NROW(UndirEdges)

           EdgeDirCount <- Adj[lower.tri(Adj)] + Adj[upper.tri(Adj)]
           unidir <- EdgeDirCount == 1
           nb.uni <- sum(unidir)
           bidir <- EdgeDirCount == 2
           nb.bi <- sum(bidir)
           nb.edges <- nb.uni + nb.bi
           cat(paste0(" MRGN Inferred a Direct Acyclic Graph with: \n"))
           cat(paste0("    ", p+q+r, " nodes \n"))
           cat(paste0("        ", q, " variants (V-nodes) \n"))
           cat(paste0("        ", p, " genes (T-nodes) \n"))
           if (r > 0)
             cat(paste0("        ", r, " intermadiate variables or common children (Q-nodes) \n"))
           if (nb.uni > 0 & nb.bi > 0)
             cat(paste0("    ", nb.edges, " edges \n"))
           if (nb.uni > 0) {
             if (nb.bi > 0)
               cat(paste0("        ", nb.uni, " directed edges \n"))
             else
               cat(paste0("    ", nb.uni, " directed edges \n"))
           }
           if (nb.bi > 0) {
             if (nb.uni > 0)
               cat(paste0("        ", nb.bi, " undirected edges \n"))
             else
               cat(paste0("    ", nb.bi, " undirected edges \n"))
           }

         })
}
