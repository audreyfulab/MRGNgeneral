# Update all edges at once (not sequentially) based on trio analysis
# Internal child function of "update.adjacency.matrix
update.adjacency.matrix_trio <-
  function (Adj, p, q, r, trio.set, inferred.models,
            stringent = FALSE,
            add.edges = TRUE,
            added.edges = NULL,
            dropped.edges = NULL,
            solve.conflicts = TRUE, method = "naive",
            cl = NULL, chunk.size = NULL,
            ...) {
    # Save the input Adjacency matrix and added/dropped edges
    Adj.old <- Adj
    dropped.edges.old <- dropped.edges

    # A vector to record the addition of new edges; Useful for generating new trios
    edges.new <- numeric(length(inferred.models))

    # Initialize two 2-columns matrices to track edge adds (or confirmation) and edge drops
    edge_add <- edge_drop <- NULL # no need of matrices; NULL is OK

    #=============================================================================
    #------------------------------- M0.1 trios ----------------------------------
    #=============================================================================
    index <- inferred.models == 'M0.1'
    if (any(index)) {
      subtrio.set <- trio.set[index,, drop = FALSE]

      # Set edge updating plan given trio analysis
      if (add.edges)
        edge_add <- rbind(edge_add,
                          cbind(subtrio.set[,1], subtrio.set[,2]))
      edge_drop <- rbind(edge_drop,
                         cbind(subtrio.set[,1], subtrio.set[,3]),
                         cbind(subtrio.set[,2], subtrio.set[,3]),
                         cbind(subtrio.set[,3], subtrio.set[,2]))

      # Check for new edge(s) between nodes (V, T)
      edges.new[index] <- !Adj.old[cbind(subtrio.set[,1], subtrio.set[,2])]
    }

    #=============================================================================
    #------------------------------- M0.2 trios ----------------------------------
    #=============================================================================
    index <- inferred.models == 'M0.2'
    if (any(index)) {
      subtrio.set <- trio.set[index,, drop = FALSE]

      if (add.edges)
        edge_add <- rbind(edge_add,
                        cbind(subtrio.set[,1], subtrio.set[,3]))
      edge_drop <- rbind(edge_drop,
                         cbind(subtrio.set[,1], subtrio.set[,2]),
                         cbind(subtrio.set[,2], subtrio.set[,3]),
                         cbind(subtrio.set[,3], subtrio.set[,2]))

      # Check for new edge(s) between V and T nodes
      edges.new[index] <- !Adj.old[cbind(subtrio.set[,1], subtrio.set[,3])]
    }

    #=============================================================================
    #------------------------------- M1.1 trios ----------------------------------
    #=============================================================================
    index <- inferred.models == 'M1.1'
    if (any(index)) {
      subtrio.set <- trio.set[index,, drop = FALSE]

      if (add.edges)
        edge_add <- rbind(edge_add,
                        cbind(subtrio.set[,1], subtrio.set[,2]),
                        cbind(subtrio.set[,2], subtrio.set[,3]))

      # Drop edges not appearing in 'M1.1'?
      if (stringent) {
       edge_drop <- rbind(edge_drop,
                         cbind(subtrio.set[,1], subtrio.set[,3]),
                         cbind(subtrio.set[,3], subtrio.set[,2]))
      }

      # Check for new edge(s) between V and T nodes
      edges.new[index] <- Adj.old[cbind(subtrio.set[,1], subtrio.set[,2])] *
        ((Adj.old[cbind(subtrio.set[,2], subtrio.set[,3])] +
            Adj.old[cbind(subtrio.set[,3], subtrio.set[,2])]) == 0)
    }

    #=============================================================================
    #------------------------------- M1.2 trios ----------------------------------
    #=============================================================================
    index <- inferred.models == 'M1.2'
    if (any(index)) {
      subtrio.set <- trio.set[index,, drop = FALSE]

      if (add.edges)
        edge_add <- rbind(edge_add,
                        cbind(subtrio.set[,1], subtrio.set[,3]),
                        cbind(subtrio.set[,3], subtrio.set[,2]))

      # Drop edges not appearing in 'M1.2'?
      if (stringent) {
      edge_drop <- rbind(edge_drop,
                         cbind(subtrio.set[,1], subtrio.set[,2]),
                         cbind(subtrio.set[,2], subtrio.set[,3]))
      }

      # Check for new edge(s) between V and T nodes
      edges.new[index] <- Adj.old[cbind(subtrio.set[,1], subtrio.set[,3])] *
        (Adj.old[cbind(subtrio.set[,2], subtrio.set[,3])] +
           Adj.old[cbind(subtrio.set[,3], subtrio.set[,2])]) == 0
    }

    #=============================================================================
    #------------------------------- M2.1 trios ----------------------------------
    #=============================================================================
    index <- inferred.models == 'M2.1'
    if (any(index)) {
      subtrio.set <- trio.set[index,, drop = FALSE]

      if (add.edges)
        edge_add <- rbind(edge_add,
                        cbind(subtrio.set[,1], subtrio.set[,2]),
                        cbind(subtrio.set[,3], subtrio.set[,2]))

      # Drop edges not appearing in 'M2.1'?
      if (stringent) {
      edge_drop <- rbind(edge_drop,
                         cbind(subtrio.set[,1], subtrio.set[,3]),
                         cbind(subtrio.set[,2], subtrio.set[,3]))
      }

      # Check for new edge(s) between V and T nodes
      edges.new[index] <- Adj.old[cbind(subtrio.set[,1], subtrio.set[,2])] *
        (Adj.old[cbind(subtrio.set[,2], subtrio.set[,3])] +
           Adj.old[cbind(subtrio.set[,3], subtrio.set[,2])]) == 0
    }

    #=============================================================================
    #------------------------------- M2.2 trios ----------------------------------
    #=============================================================================
    index <- inferred.models == 'M2.2'
    if (any(index)) {
      subtrio.set <- trio.set[index,, drop = FALSE]
      subtrio.set <- trio.set[index,, drop = FALSE]

      if (add.edges)
        edge_add <- rbind(edge_add,
                        cbind(subtrio.set[,1], subtrio.set[,3]),
                        cbind(subtrio.set[,2], subtrio.set[,3]))

      # Drop edges not appearing in 'M2.2'?
      if (stringent) {
      edge_drop <- rbind(edge_drop,
                         cbind(subtrio.set[,1], subtrio.set[,2]),
                         cbind(subtrio.set[,3], subtrio.set[,2]))
      }

      # Check for new edge(s) between V and T nodes
      edges.new[index] <- Adj.old[cbind(subtrio.set[,1], subtrio.set[,3])] *
        (Adj.old[cbind(subtrio.set[,2], subtrio.set[,3])] +
           Adj.old[cbind(subtrio.set[,3], subtrio.set[,2])]) == 0
    }

    #=============================================================================
    #-------------------------------- M3 trios -----------------------------------
    #=============================================================================
    index <- inferred.models == 'M3'
    if (any(index)) {
      subtrio.set <- trio.set[index,, drop = FALSE]

      if (add.edges)
        edge_add <- rbind(edge_add,
                        cbind(subtrio.set[,1], subtrio.set[,2]),
                        cbind(subtrio.set[,1], subtrio.set[,3]))
      edge_drop <- rbind(edge_drop,
                         cbind(subtrio.set[,2], subtrio.set[,3]),
                         cbind(subtrio.set[,3], subtrio.set[,2]))

      # Check for new edge(s) between V and T nodes
      edges.new[index] <- Adj.old[cbind(subtrio.set[,1], subtrio.set[,2])] *
        Adj.old[cbind(subtrio.set[,1], subtrio.set[,3])] == 0

    }

    #=============================================================================
    #------------------------------- M4 trios ----------------------------------
    #=============================================================================
    index <- inferred.models == 'M4'
    if (any(index)) {
      subtrio.set <- trio.set[index,, drop = FALSE]

      if (add.edges)
        edge_add <- rbind(edge_add,
                        cbind(subtrio.set[,1], subtrio.set[,2]),
                        cbind(subtrio.set[,1], subtrio.set[,3]),
                        cbind(subtrio.set[,2], subtrio.set[,3]),
                        cbind(subtrio.set[,3], subtrio.set[,2]))

      # Check for new edge(s) between V and T nodes
      edges.new[index] <- Adj.old[cbind(subtrio.set[,1], subtrio.set[,2])] *
        Adj.old[cbind(subtrio.set[,1], subtrio.set[,3])] *
        (Adj.old[cbind(subtrio.set[,2], subtrio.set[,3])] +
           Adj.old[cbind(subtrio.set[,3], subtrio.set[,2])]) == 0
    }

    #=============================================================================
    #------------------------------- Other.1 trios ----------------------------------
    #=============================================================================
    index <- inferred.models == 'Other.1'
    if (any(index)) {
      subtrio.set <- trio.set[index,, drop = FALSE]

      if (add.edges)
        edge_add <- rbind(edge_add,
                        cbind(subtrio.set[,2], subtrio.set[,3]),
                        cbind(subtrio.set[,3], subtrio.set[,2]))

      # Check for new edge(s) between V and T nodes
      edges.new[index] <- (Adj.old[cbind(subtrio.set[,2], subtrio.set[,3])] +
           Adj.old[cbind(subtrio.set[,3], subtrio.set[,2])]) == 0
    }

    #=============================================================================
    #------------------------------- Other.2 trios ----------------------------------
    #=============================================================================
    index <- inferred.models == 'Other.2'
    if (any(index)) {
      subtrio.set <- trio.set[index,, drop = FALSE]
      edge_drop <- rbind(edge_drop,
                        cbind(subtrio.set[,2], subtrio.set[,3]),
                        cbind(subtrio.set[,3], subtrio.set[,2]))

    }

    # Remove duplicates
    edge_drop <- unique(edge_drop)
    if (add.edges)
      edge_add <- unique(edge_add)

    # Update the list of dropped edge directions
    dropped.edges <- c(dropped.edges, paste0(edge_drop[,1], '.', edge_drop[,2]))

    ### Update edge directions

    ## List of edge directions to be added
    add.edges_new <- if (add.edges) paste0(edge_add[,1], '.', edge_add[,2])

    # logical vector indicating conflict or not
    if (add.edges) {
      conflicts.adds <- is.element(add.edges_new, dropped.edges)
      nb.add_conflicts <- sum(conflicts.adds)
      if (nb.add_conflicts) {# If any conflict
        edge_add <- edge_add[!conflicts.adds, , drop = FALSE] # solve conflict: remove conflict edge directions from 'edge_add'
      }
    }
    else {
      nb.add_conflicts <- 0
    }

    ## Add edges
    if (NROW(edge_add)) {
      Adj[edge_add] <- 1
    }

    ## Drop edges
    if (NROW(edge_drop)) {
      Adj[edge_drop] <- 0
    }
    nb.drop_conflicts <- 0

    # Count newly added edges
    if (!is.null(edge_add)) {
      ## Set of new edges (all directions)
      add.edges_new <- paste0(edge_add[,1], '.', edge_add[,2])
      add.edges_new.rv <- paste0(edge_add[,2], '.', edge_add[,1])

      ## Count new edges
      if (!is.null(added.edges))
        nb.new_edges <- sum((add.edges_new %in% added.edges) + (add.edges_new.rv %in% added.edges) == 0)
      else
        nb.new_edges <- length(edge_add)

      ## Update the list 'added.edges'
      added.edges <- unique(c(added.edges, add.edges_new))
    }
    else
      nb.new_edges <- 0

    # Count dropped edges
    if (!is.null(edge_drop)) {
      ## Set of dropped edges (all directions)
      drop.edges_new <- paste0(edge_drop[,1], '.', edge_drop[,2])
      drop.edges_new.rv <- paste0(edge_drop[,2], '.', edge_drop[,1])

      ## Count new drops
      nb.dropped_edges <- sum((drop.edges_new %in% dropped.edges.old) + (drop.edges_new.rv %in% dropped.edges.old) == 0)
    }
    else
      nb.dropped_edges <- 0

    return(list(Adj = Adj,
                nb.new_edges = nb.new_edges,
                nb.dropped_edges = nb.dropped_edges,
                new.edges = edges.new,
                added.edges = added.edges,
                dropped.edges = dropped.edges,
                nb.conflicts = nb.drop_conflicts + nb.add_conflicts))
  }

dropedges <- function (edge_drop, Adj, dropped.edges) {

  if (!is.null(edge_drop)) {
    # Remove duplicates
    edge_drop <- unique(edge_drop)

    # Update edge directions
    Adj[edge_drop] <- 0

    # Update the list of dropped edge directions
    drop.edges_new <- paste0(edge_drop[,1], '.', edge_drop[,2])
    dropped.edges <- c(dropped.edges, drop.edges_new)
  }

  list(Adj = Adj,
       edge_drop = edge_drop,
       dropped.edges = dropped.edges,
       nb.drop_conflicts = 0)

}

### Add edges
addedges <- function (edge_add, dropped.edges,
                      q, Adj, edge_drop,
                      solve.conflicts,
                      method, cl, chunk.size) {

  if (!is.null(edge_add)) {
    ## Remove duplicates
    edge_add <- unique(edge_add)

    ## Check for conflicts with previous drops before adding edge directions
    ## (Check that we are not adding an edge we removed earlier)

    # List of edge directions to be added
    add.edges_new <- paste0(edge_add[,1], '.', edge_add[,2])
    if (!is.null(dropped.edges)) {
      # logical vector indicating conflict or not
      conflicts.adds <- is.element(add.edges_new, dropped.edges)
      nb.add_conflicts <- sum(conflicts.adds)

      if (nb.add_conflicts) {# If any conflict

          # Does the conflicting edge involve a genetic variant or just T-nodes?
          V.conflicts.adds <-  edge_add[conflicts.adds,1] <= q


          #

          if (solve.conflicts) {
          # Use 'conflict.add.fun' to solve each conflict
          confirm.add <- matteApply(X = cbind(edge_add[conflicts.adds,, drop = FALSE],
                                              V.conflicts.adds),
                                    MARGIN = 1,
                                    FUN = conflict.add.fun,
                                    Adj = Adj,
                                    method = method,
                                    cl = cl, chunk.size = chunk.size)
        }
        else {
          # Keep any conflicting edge that involves a genetic variant; drop any T-T edge
          confirm.add <- V.conflicts.adds
        }

        # Update 'edge_add' with the result
        edge_add <- rbind(edge_add[!conflicts.adds,], # non-conflicting adds
                          edge_add[conflicts.adds,,drop = FALSE][confirm.add,]) # keep only confirmed adds

        # Update the list of newly added edges
        add.edges_new <- paste0(edge_add[,1], '.', edge_add[,2])

        # Did we re-added an edge direction we dropped?
        if (any(re.add <- is.element(dropped.edges, add.edges_new))) {
          # Update the list of dropped edge directions
          dropped.edges <- dropped.edges[!re.add]
        }

        # List of newly dropped edges
        drop.edges_new <- paste0(edge_drop[,1], '.', edge_drop[,2])
        # Did we re-added an edge direction we just dropped?
        if (any(re.add <- is.element(drop.edges_new, add.edges_new))) {
          # Update the list of dropped edge directions
          edge_drop <- edge_drop[!re.add, , drop = FALSE]
        }

      }
    }
    else {
      nb.add_conflicts <- 0
    }

    # Update edge directions
    Adj[edge_add] <- 1
  }
  else {
    nb.add_conflicts <- 0
  }

  list(Adj = Adj,
       edge_add = edge_add,
       edge_drop = edge_drop,
       dropped.edges = dropped.edges,
       nb.add_conflicts = nb.add_conflicts)
}

# Update edge directions
#Adj <- addedges (edge_add = edge_add,
#                 dropped.edges = dropped.edges,
#                 q = q,
#                 Adj = Adj,
#                 edge_drop = edge_drop,
#                 solve.conflicts = solve.conflicts,
#                 method = method,
#                 cl = cl,
#                 chunk.size = chunk.size)
#edge_add <- Adj$edge_add
#dropped.edges <- Adj$dropped.edges
#nb.add_conflicts <- Adj$nb.add_conflicts
#edge_drop <- Adj$edge_drop
#Adj <- Adj$Adj
## Drop edges
#Adj <- dropedges (edge_drop = edge_drop,
#                  dropped.edges = dropped.edges,
#                  Adj = Adj)
#edge_drop <- Adj$edge_drop
#dropped.edges <- Adj$dropped.edges
#nb.drop_conflicts <- Adj$nb.drop_conflicts
#Adj <- Adj$Adj

