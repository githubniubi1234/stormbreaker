sb_binary_menu_arguments <- function(row, selector_features, feature_set_id,
                                     selector_n_features, inner_folds,
                                     inner_seed) {
  args <- list(selector = row$selector, learner = row$learner)
  if (identical(row$selector_parameter, "selector_n_features")) {
    args$selector_n_features <- selector_n_features
  }
  if (identical(row$selector_parameter, "selector_features_and_feature_set_id")) {
    args$selector_features <- selector_features
    args$feature_set_id <- feature_set_id
  }
  if (identical(row$tuning_semantics, "lambda_selected_by_inner_cv")) {
    args$inner_folds <- inner_folds
    args$inner_seed <- inner_seed
  }
  if (identical(row$learner, "stepglm_binary")) {
    args$stepglm_direction <- "both"
  }
  args
}

#' Construct the complete binary classification pipeline menu
#'
#' Creates one specification for every registered binary selector-learner pair.
#' The constructor does not read expression data, fit models, rank methods or
#' choose genes. All fixed-feature routes use the same caller-declared vector.
#' Learned selectors retain their existing native defaults, including their
#' own tuning folds, seeds and Boruta budgets. Top-K size and glmnet learner
#' tuning controls are explicitly recorded in the returned pipelines.
#'
#' The seven MIME classification families are identified for comparison, not
#' reproduced by this constructor. These fixed STORMBREAKER configurations do
#' not implement MIME's caret tuning grids or repeated cross-validation; run
#' the original MIME methods separately for that benchmark. The mapping is
#' deliberately labelled `algorithm_family_only`.
#'
#' The `none` and `fixed_features` routes remain separate specifications.
#' They may use identical inputs when the task's predictors equal the fixed
#' feature vector. Use the shared comparison/reuse interfaces rather than
#' treating such duplicate inputs as additional independent evidence.
#'
#' @param selector_features Ordered, unique names for fixed-feature routes.
#' @param feature_set_id Non-empty label for the declared feature block.
#' @param selector_n_features Positive integer, at least 2, for Top-K screening.
#'   Defaults to 10 and is not selected from outcomes or clipped to the data.
#' @param inner_folds Number of inner folds for the glmnet logistic learner,
#'   at least 3. Defaults to 5. It does not override selector-native tuning.
#' @param inner_seed Positive integer seed for glmnet learner tuning. Other
#'   learners and selectors retain their saved native defaults.
#'
#' @return A data frame with unique `config_id`, `selector`, `learner`,
#'   `mime_family`, `comparison_role` and a `pipeline` list-column. The
#'   `mime7_family_mapping` attribute records the seven family references and
#'   their non-equivalence to MIME tuning. `registry_pairs` records coverage.
#'   Counts describe configurations, not independent algorithm families.
#' @export
sb_binary_menu_pipelines <- function(selector_features, feature_set_id,
                                      selector_n_features = 10L,
                                      inner_folds = 5L, inner_seed = 1L) {
  selector_features <- sb_validate_fixed_features(selector_features)
  feature_set_id <- sb_validate_feature_set_id(feature_set_id)
  selector_n_features <- sb_check_positive_integer(
    selector_n_features, "selector_n_features", minimum = 2L
  )
  inner_folds <- sb_check_positive_integer(
    inner_folds, "inner_folds", minimum = 3L
  )
  inner_seed <- sb_check_positive_integer(inner_seed, "inner_seed", minimum = 1L)
  registry <- sb_pipeline_registry(task_type = "binary")
  registry <- registry[order(
    registry$selector != "none", registry$selector, registry$learner
  ), , drop = FALSE]
  pipelines <- lapply(seq_len(nrow(registry)), function(i) {
    args <- sb_binary_menu_arguments(
      registry[i, , drop = FALSE], selector_features, feature_set_id,
      selector_n_features, inner_folds, inner_seed
    )
    do.call(sb_pipeline, args)
  })
  menu <- data.frame(
    config_id = sprintf("binary_%04d", seq_along(pipelines)),
    selector = registry$selector,
    learner = registry$learner,
    stringsAsFactors = FALSE
  )
  reference <- data.frame(
    official_method = c("nb", "svmRadialWeights", "rf", "kknn",
      "adaboost", "LogitBoost", "cancerclass"),
    learner = c("kernel_naive_bayes_binary", "svm_radial_weighted_binary",
      "random_forest_binary", "kernel_weighted_knn_binary",
      "adaboost_binary", "logitboost_binary", "cancerclass_binary"),
    relationship = "algorithm_family_only",
    stringsAsFactors = FALSE
  )
  keys <- paste(menu$selector, menu$learner)
  reference$config_id <- menu$config_id[match(paste("none", reference$learner), keys)]
  menu$mime_family <- reference$official_method[match(menu$config_id, reference$config_id)]
  menu$mime_family[is.na(menu$mime_family)] <- ""
  menu$comparison_role <- ifelse(nzchar(menu$mime_family),
    "MIME_family_reference", "STORMBREAKER_extension")
  menu$pipeline <- pipelines
  attr(menu, "mime7_family_mapping") <- reference
  attr(menu, "registry_pairs") <- nrow(registry)
  menu
}
