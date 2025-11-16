# Internal helper function to set up parallel clusters with proper library paths and packages
# This centralizes cluster setup logic to avoid code duplication
setup_cluster <- function(cl, packages = c("MRGN", "ppcor", "propagate"), seed = NULL) {
  if (is.null(cl)) return(NULL)

  # Set library paths on workers to match main session
  # Find where the required packages are actually installed
  dep_paths <- unique(dirname(find.package(packages, quiet = TRUE)))
  parallel::clusterCall(cl, function(paths, dep_paths) {
    .libPaths(unique(c(paths, dep_paths)))
  }, .libPaths(), dep_paths)

  # Set random seed on workers for reproducibility (if provided)
  if (!is.null(seed)) {
    parallel::clusterSetRNGStream(cl, seed)
  }

  # Load required packages on cluster workers
  # Pass the package list to workers and load them dynamically
  parallel::clusterCall(cl, function(pkgs) {
    for (pkg in pkgs) {
      suppressPackageStartupMessages(library(pkg, character.only = TRUE))
    }
  }, packages)

  invisible(cl)
}
