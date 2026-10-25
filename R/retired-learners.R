# Retirement concerns the binary/survival learners, not similarly named engines.
sb_retired_learner_ids <- function() {
  c("gam_cox_reml", "gam_binary_reml", "slope_binary")
}

sb_reject_retired_learner <- function(learner) {
  if (is.character(learner) && any(learner %in% sb_retired_learner_ids())) {
    sb_abort(
      paste0(
        "Learner `", paste(learner[learner %in% sb_retired_learner_ids()], collapse = ", "),
        "` was removed in 0.3.38 because of excessive runtime and computational cost."
      ),
      "stormbreaker_retired_learner"
    )
  }
  invisible(NULL)
}
