# Update all edges at once (not sequentially) based on trio analysis
# Internal child function of "update.adjacency.matrix
#' @export
update.adjacency.matrix_trio <-
  function (adjacency,
            n_v, n_t, n_q, # Currently not used,
            trio.set, inferred.models,
            stringent = FALSE,
            add.edges = TRUE,
            added.edges = NULL,
            dropped.edges = NULL,
            solve.conflicts = TRUE, method = "conservative", # Currently not used: this is the only available method
            cl = NULL, chunk.size = NULL,
            ...) {
    # Save the input Adjacency matrix and added/dropped edges
    adjacency.old <- adjacency
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
      edges.new[index] <- !adjacency.old[cbind(subtrio.set[,1], subtrio.set[,2])]
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
      edges.new[index] <- !adjacency.old[cbind(subtrio.set[,1], subtrio.set[,3])]
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
      edges.new[index] <- adjacency.old[cbind(subtrio.set[,1], subtrio.set[,2])] *
        ((adjacency.old[cbind(subtrio.set[,2], subtrio.set[,3])] +
            adjacency.old[cbind(subtrio.set[,3], subtrio.set[,2])]) == 0)
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
      edges.new[index] <- adjacency.old[cbind(subtrio.set[,1], subtrio.set[,3])] *
        (adjacency.old[cbind(subtrio.set[,2], subtrio.set[,3])] +
           adjacency.old[cbind(subtrio.set[,3], subtrio.set[,2])]) == 0
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
      edges.new[index] <- adjacency.old[cbind(subtrio.set[,1], subtrio.set[,2])] *
        (adjacency.old[cbind(subtrio.set[,2], subtrio.set[,3])] +
           adjacency.old[cbind(subtrio.set[,3], subtrio.set[,2])]) == 0
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
      edges.new[index] <- adjacency.old[cbind(subtrio.set[,1], subtrio.set[,3])] *
        (adjacency.old[cbind(subtrio.set[,2], subtrio.set[,3])] +
           adjacency.old[cbind(subtrio.set[,3], subtrio.set[,2])]) == 0
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
      edges.new[index] <- adjacency.old[cbind(subtrio.set[,1], subtrio.set[,2])] *
        adjacency.old[cbind(subtrio.set[,1], subtrio.set[,3])] == 0

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
      edges.new[index] <- adjacency.old[cbind(subtrio.set[,1], subtrio.set[,2])] *
        adjacency.old[cbind(subtrio.set[,1], subtrio.set[,3])] *
        (adjacency.old[cbind(subtrio.set[,2], subtrio.set[,3])] +
           adjacency.old[cbind(subtrio.set[,3], subtrio.set[,2])]) == 0
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
      edges.new[index] <- (adjacency.old[cbind(subtrio.set[,2], subtrio.set[,3])] +
           adjacency.old[cbind(subtrio.set[,3], subtrio.set[,2])]) == 0
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
      adjacency[edge_add] <- 1
    }

    ## Drop edges
    if (NROW(edge_drop)) {
      adjacency[edge_drop] <- 0
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

    return(list(adjacency = adjacency,
                nb.new_edges = nb.new_edges,
                nb.dropped_edges = nb.dropped_edges,
                new.edges = edges.new,
                added.edges = added.edges,
                dropped.edges = dropped.edges,
                nb.conflicts = nb.drop_conflicts + nb.add_conflicts))
  }
