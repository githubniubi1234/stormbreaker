make_sb_curve_view_fixture <- function(task_type = "binary", comparison = FALSE) {
  thresholds <- c(0.1, 0.2, 0.3, 0.5, 0.7)
  predictors <- c("gene_a", "gene_b")
  if (identical(task_type, "binary")) {
    d <- sb_example_binary(120, seed = 1)
    task <- sb_binary_task(d, "patient_id", "outcome", predictors,
                           decision_thresholds = thresholds)
    pipelines <- list(sb_pipeline(learner = "logistic"))
    if (comparison) {
      pipelines[[2]] <- sb_pipeline(learner = "gaussian_naive_bayes")
    }
    external <- sb_example_binary(100, seed = 2)
    external$patient_id <- paste0("EXT_", external$patient_id)
    endpoint <- sb_binary_endpoint(external, "patient_id", "outcome")
  } else {
    d <- sb_example_survival(120, seed = 1)
    task <- sb_surv_task(d, "patient_id", "time", "event", predictors, "month",
                         prediction_horizon = 6, decision_thresholds = thresholds)
    pipelines <- list(sb_pipeline())
    if (comparison) {
      pipelines[[2]] <- sb_pipeline(selector = "univ_cox_top_k",
                                    selector_n_features = 1)
    }
    external <- sb_example_survival(100, seed = 2)
    external$patient_id <- paste0("EXT_", external$patient_id)
    endpoint <- sb_surv_endpoint(external, "patient_id", "time", "event", "month",
                                 prediction_horizon = 6)
  }
  predictions <- lapply(pipelines, function(pipeline) {
    fit <- sb_fit_final(task, pipeline)
    sb_predict_final(fit, external[c("patient_id", predictors)],
                     id_col = "patient_id")
  })
  if (comparison) {
    cv <- sb_compare_cv(task, pipelines, folds = 3, seed = 1)
    evaluated <- sb_compare_external(predictions, endpoint)
  } else {
    cv <- sb_run_cv(task, pipelines[[1]], folds = 3, seed = 1)
    evaluated <- sb_evaluate_external(predictions[[1]], endpoint)
  }
  list(cv = cv, external = evaluated, task = task)
}
