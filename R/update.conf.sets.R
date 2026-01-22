# A function to update an object of class 'conf.sets' as returned by 'get.conf.sets'
# Change can be made to the selection measure (for T and C-nodes), the
# FDR control method, the grouping for adjustment (adjust_by),
# the level of adjustment (fdr) or the individual p-value significance level (alpha)
#' @keywords internal
#' @noRd
update_conf_sets <- function (object,
                             T.measure = c("partial", "marginal"), # Selection measure
                             C.measure = T.measure,
                             FDRcontrol = c("qvalue", "bonferroni", "none"),
                             V.FDRcontrol = FDRcontrol,
                             T.FDRcontrol = FDRcontrol,
                             C.FDRcontrol = FDRcontrol,
                             WZ.FDRcontrol = FDRcontrol,
                             adjust_by = c("individual", "all", "none"),
                             V.adjust_by = adjust_by,
                             T.adjust_by = adjust_by,
                             C.adjust_by = adjust_by,
                             WZ.adjust_by = adjust_by,
                             fdr = 0.05,
                             lambda = 0.05,
                             pi0.method = c("smoother", "boostrap"),
                             alpha = 0.01, # Only used if FDRcontrol = 'none' or 'bonferroni'
                             verbose = 0,
                             save.list = FALSE, # Only used if 'save.path' is not missing
                             save.path = "/path/to/save/location/") {
  ### Write diverse child functions to perform the tasks in 'get.conf.sets'

  ### Determine the components of 'object' that need updates

  ### Call the appropriate one here to update 'object'
}


