sb_survival_menu_arguments <- function(row, selector_features, feature_set_id,
                                       selector_n_features) {
  args <- list(selector = row$selector, learner = row$learner)
  parameter <- row$selector_parameter
  if (identical(parameter, "selector_n_features")) {
    args$selector_n_features <- selector_n_features
  }
  if (parameter %in% c(
    "stepcox_direction", "selector_rsf_nodesize_seed_and_stepcox_direction",
    "selector_inner_folds_seed_and_stepcox_direction"
  )) {
    args$stepcox_direction <- "both"
  }
  if (parameter %in% c(
    "selector_rsf_nodesize_and_seed",
    "selector_rsf_nodesize_seed_and_stepcox_direction"
  )) {
    args$selector_rsf_nodesize <- 5L
    args$selector_seed <- 1L
  }
  if (parameter %in% c(
    "selector_inner_folds_and_seed",
    "selector_inner_folds_seed_and_stepcox_direction"
  )) {
    args$selector_inner_folds <- 10L
    args$selector_seed <- 1L
  }
  if (identical(parameter, "selector_features_and_feature_set_id")) {
    args$selector_features <- selector_features
    args$feature_set_id <- feature_set_id
  }
  if (identical(row$tuning_semantics, "lambda_selected_by_inner_cv")) {
    args$inner_folds <- 10L
    args$inner_seed <- 1L
    args$glmnet_standardize <- TRUE
    args$glmnet_lambda_rule <- "lambda.min"
    args$glmnet_tuning_measure <- "deviance"
  }
  if (identical(row$tuning_semantics, "analysis_fold_internal_cv")) {
    args$learner_tuning_folds <- if (row$learner == "plsrcox_cv") 5L else 10L
    args$learner_seed <- 1L
    if (row$learner == "gbm_cox_cv") {
      args$gbm_cv_max_trees <- 10000L
    }
  }
  args
}

#' Construct the full MIME-aligned survival pipeline menu
#'
#' Creates one specification for every registered survival selector-learner
#' pair and includes all 117 MIME reference variants. Identical complete
#' pipeline IDs are deduplicated. This constructor only creates specifications;
#' it does not read expression data, fit models, rank methods or select genes.
#'
#' Common MIME settings are saved in the returned pipeline objects: 10 inner
#' folds for glmnet, learned Lasso/CoxBoost selectors, CoxBoost, GBM and SuperPC;
#' 5 folds for PLS-Cox; a common declared seed of 1; and up to 10000 GBM trees.
#' The reference suite retains all Elastic Net alphas and StepCox directions.
#' Other registered methods keep their native package defaults, including the
#' default Boruta budget. Fixed-component and fixed-step extensions remain
#' separate from MIME's tuned learners. The low-level [sb_pipeline()] API keeps
#' its explicit parameter requirements and is not changed by this preset.
#'
#' Parameter alignment does not imply identical internal fold assignments or
#' bug-for-bug MIME implementation. It does not turn a training-selected input
#' panel into independent or unbiased validation. See [sb_mime117_pipelines()]
#' for the reference suite's implementation boundaries.
#'
#' @param selector_features Ordered, unique names for the caller-declared
#'   fixed-feature routes. Every fixed-feature pipeline uses this same vector.
#' @param feature_set_id Non-empty label for the declared feature block.
#' @param selector_n_features Positive integer, at least 2, for Top-K selectors.
#'   The default is 30. This number must be declared before evaluation; the
#'   constructor does not choose it from outcomes or clip it to observed data.
#'
#' @return A data frame with a `pipeline` list-column, unique `config_id`,
#'   `selector`, `learner`, `mime_labels` and `comparison_role` columns.
#'   The `mime117_mapping` attribute maps the 117 reference labels to menu
#'   configuration IDs. The `registry_pairs` attribute records registry coverage.
#'   Menu size can increase when the package admits additional methods; counts
#'   describe configurations, not independent algorithm families.
#' @export
sb_survival_menu_pipelines <- function(selector_features, feature_set_id,
                                       selector_n_features = 30L) {
  selector_features <- sb_validate_fixed_features(selector_features)
  feature_set_id <- sb_validate_feature_set_id(feature_set_id)
  selector_n_features <- sb_check_positive_integer(
    selector_n_features, "selector_n_features", minimum = 2L
  )
  registry <- sb_pipeline_registry(task_type = "survival")
  registry <- registry[order(
    registry$selector != "none", registry$selector, registry$learner
  ), , drop = FALSE]
  pipelines <- lapply(seq_len(nrow(registry)), function(i) {
    args <- sb_survival_menu_arguments(
      registry[i, , drop = FALSE], selector_features, feature_set_id,
      selector_n_features
    )
    do.call(sb_pipeline, args)
  })
  reference <- sb_mime117_pipelines(
    inner_folds = 10L, inner_seed = 1L,
    selector_inner_folds = 10L, selector_seed = 1L, rsf_nodesize = 5L,
    learner_tuning_folds = 10L, plsrcox_tuning_folds = 5L,
    learner_seed = 1L, gbm_cv_max_trees = 10000L
  )
  pipelines <- c(pipelines, reference$pipeline)
  ids <- vapply(pipelines, `[[`, "", "id")
  pipelines <- pipelines[!duplicated(ids)]
  ids <- vapply(pipelines, `[[`, "", "id")
  config_ids <- sprintf("aligned_%04d", seq_along(pipelines))
  mapping <- reference[setdiff(names(reference), "pipeline")]
  mapping$comparison_config_id <- config_ids[match(
    vapply(reference$pipeline, `[[`, "", "id"), ids
  )]
  menu <- data.frame(
    config_id = config_ids,
    selector = vapply(pipelines, `[[`, "", "selector"),
    learner = vapply(pipelines, `[[`, "", "learner")
  )
  menu$mime_labels <- vapply(config_ids, function(id) {
    paste(mapping$official_model_label[mapping$comparison_config_id == id],
      collapse = " | "
    )
  }, "")
  menu$comparison_role <- ifelse(
    nzchar(menu$mime_labels), "MIME_structure_matched", "STORMBREAKER_extension"
  )
  menu$pipeline <- pipelines
  attr(menu, "mime117_mapping") <- mapping
  attr(menu, "registry_pairs") <- nrow(registry)
  menu
}
