# Internal child function of "update.adjacency.matrix
update.adjacency.matrix_triplet <- function (adjacency,
                                             n_v, n_t, n_q, # Currently not used
                                             triplet.set,
                                             inferred.models,
                                             stringent = FALSE,
                                             add.edges = TRUE,
                                             solve.conflicts = TRUE,
                                             method = "conservative", # Currently not used: this is the only available method
                                             added.edges = NULL,
                                             dropped.edges = NULL,
                                             cl = NULL,
                                             chunk.size = NULL,
                                             ...) {

  # Save the input Adjacency matrix
  adjacency.old <- adjacency
  dropped.edges.old <- dropped.edges

  # Initialize two 2-columns matrices to track edge adds (or confirmation) and edge drops
  edge_add <- edge_drop <- NULL # no need of matrices; NULL is OK

  #=============================================================================
  #------------------------------- M1.1 trios ----------------------------------
  #=============================================================================
  index <- inferred.models == 'M1.1'
  if (any(index)) {
    subtriplet.set <- triplet.set[index,, drop = FALSE]

    # Triplets of type 2
    if (any(type2 <- subtriplet.set[,4] == 2)) {

      if (add.edges) {
      edge_add <- rbind(edge_add,
                        cbind(subtriplet.set[type2,1], subtriplet.set[type2,3])) # Yes (Ti --> Tk) # useless, Ti --> Tk is always 1 for a type 2 triplet
      }

      # Drop edges not appearing in 'M1.1'?
      if (stringent) {
        edge_drop <- rbind(edge_drop,
                           cbind(subtriplet.set[type2,3], subtriplet.set[type2,1]))
      }
    }
  }

  #=============================================================================
  #------------------------------- M2.1 trios ----------------------------------
  #=============================================================================
  index <- inferred.models == 'M2.1'
  if (any(index)) {
    subtriplet.set <- triplet.set[index,, drop = FALSE]

    if (add.edges) {
    edge_add <- rbind(edge_add,
                      cbind(subtriplet.set[,3], subtriplet.set[,1]), # Yes (Ti --> Tk) # useless, Ti --> Tk is always 1 for a type 2 triplet
                      cbind(subtriplet.set[,2], subtriplet.set[,1])) # Yes (Ti --> Tk) # useless, Ti --> Tk is always 1 for a type 2 triplet
    }

    # Drop edges not appearing in 'M2.1'?
    if (stringent) {
      edge_drop <- rbind(edge_drop,
                         cbind(subtriplet.set[,1], subtriplet.set[,3]), # No  (Tk --> Ti)
                         cbind(subtriplet.set[,1], subtriplet.set[,2])) # No  (Tk --> Ti)
    }
  }

  # Remove duplicates
  edge_drop <- unique(edge_drop)
  edge_add <- unique(edge_add)

  # Update the list of dropped edge directions
  dropped.edges <- c(dropped.edges, paste0(edge_drop[,1], '.', edge_drop[,2]))

  ## Update edge directions
  nb.drop_conflicts <- nb.add_conflicts <- 0
  ## Drop edges
  if (NROW(edge_drop)) {
    adjacency[edge_drop] <- 0
  }

  # List of edge directions to be added
  if (!is.null(edge_add)) {
    add.edges_new <- paste0(edge_add[,1], '.', edge_add[,2])
    # logical vector indicating conflict or not
    conflicts.adds <- is.element(add.edges_new, dropped.edges)
    nb.add_conflicts <- sum(conflicts.adds)
    if (nb.add_conflicts) {# If any conflict
      edge_add <- edge_add[!conflicts.adds, , drop = FALSE] # solve conflict: remove conflict edge directions from 'edge_add'
    }
    ## Add edges
    if (NROW(edge_add)) {
      adjacency[edge_add] <- 1
    }

  # Count newly added edges
    ## Set of new edges (all directions)
    add.edges_new <- paste0(edge_add[,1], '.', edge_add[,2])
    add.edges_new.rv <- paste0(edge_add[,2], '.', edge_add[,1])

    ## Count new edges
    if (!is.null(added.edges))
      nb.new_edges <- sum((add.edges_new %in% added.edges) + (add.edges_new.rv %in% added.edges) == 0)
    else
      nb.new_edges <- length(edge_add)

    ## Update the list 'added.edges'
    added.edges <- c(added.edges, edge_add)
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
              nb.dropped_edges = nb.dropped_edges,
              nb.new_edges = length(edge_add),
              added.edges = added.edges,
              dropped.edges = dropped.edges,
              nb.conflicts = nb.drop_conflicts + nb.add_conflicts))
}

