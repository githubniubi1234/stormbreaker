# STORMBREAKER minimum end-to-end example
# Synthetic data only: no result below is clinical evidence.

library(stormbreaker)

# ---- Right-censored survival workflow --------------------------------------

survival_development <- surv_example(n = 100, seed = 20260812)
survival_task <- surv_task(
  data = survival_development,
  id_col = "patient_id",
  time_col = "time",
  event_col = "event",
  predictors = c("gene_a", "gene_b", "gene_c"),
  time_unit = "month",
  prediction_horizon = 6
)

# Both specifications are declared before resampling and use identical folds.
survival_pipelines <- list(
  all_features_cox = sb_pipeline(
    selector = "none",
    learner = "coxph"
  ),
  top2_cox = sb_pipeline(
    selector = "univ_cox_top_k",
    learner = "coxph",
    selector_n_features = 2L
  )
)
survival_cv <- compare_cv(
  survival_task,
  survival_pipelines,
  folds = 5L,
  seed = 5201314L
)
print(survival_cv)
survival_cv$comparison

# STORMBREAKER does not choose a winner. The study team supplies the final
# specification; this synthetic example uses the predeclared all-feature Cox.
survival_final_pipeline <- survival_pipelines[["all_features_cox"]]
survival_fit <- fit_final(survival_task, survival_final_pipeline)

survival_external <- surv_example(n = 40, seed = 20260820)
survival_external$patient_id <- paste0(
  "survival_external_",
  seq_len(nrow(survival_external))
)
survival_new_data <- survival_external[
  c("patient_id", "gene_a", "gene_b", "gene_c")
]

# Prediction is outcome-blind: time and event are not supplied here.
survival_prediction <- predict_final(
  survival_fit,
  survival_new_data,
  id_col = "patient_id"
)

# Outcomes enter only after predictions have been saved.
survival_endpoint <- surv_endpoint(
  survival_external,
  id_col = "patient_id",
  time_col = "time",
  event_col = "event",
  time_unit = "month",
  prediction_horizon = 6
)
survival_evaluation <- eval_external(
  survival_prediction,
  survival_endpoint
)
print(survival_evaluation)
survival_evaluation$metrics$harrell_c
survival_evaluation$metrics$ipcw_brier

# ---- Binary workflow -------------------------------------------------------

binary_development <- bi_example(n = 100, seed = 20260819)
binary_task <- bi_task(
  data = binary_development,
  id_col = "patient_id",
  outcome_col = "outcome",
  predictors = c("gene_a", "gene_b", "gene_c")
)

binary_pipelines <- list(
  logistic = sb_pipeline(learner = "logistic"),
  gaussian_naive_bayes = sb_pipeline(learner = "gaussian_naive_bayes")
)
binary_cv <- compare_cv(
  binary_task,
  binary_pipelines,
  folds = 5L,
  seed = 5201314L
)
print(binary_cv)
binary_cv$comparison

# Again, the package does not select from the comparison automatically.
binary_final_pipeline <- binary_pipelines[["logistic"]]
binary_fit <- fit_final(binary_task, binary_final_pipeline)

binary_external <- bi_example(n = 40, seed = 20260821)
binary_external$patient_id <- paste0(
  "binary_external_",
  seq_len(nrow(binary_external))
)
binary_new_data <- binary_external[
  c("patient_id", "gene_a", "gene_b", "gene_c")
]
binary_prediction <- predict_final(
  binary_fit,
  binary_new_data,
  id_col = "patient_id"
)
binary_endpoint <- bi_endpoint(
  binary_external,
  id_col = "patient_id",
  outcome_col = "outcome"
)
binary_evaluation <- eval_external(
  binary_prediction,
  binary_endpoint
)
print(binary_evaluation)
binary_evaluation$metrics$roc_auc
binary_evaluation$metrics$brier_score
binary_evaluation$metrics$calibration

# Minimal integration assertions. Detailed statistical reconstruction tests
# belong to tests/testthat rather than this user-facing example.
stopifnot(
  identical(survival_cv$status, "success"),
  identical(survival_fit$status, "success"),
  identical(survival_prediction$status, "success"),
  identical(survival_evaluation$status, "success"),
  !any(c("time", "event") %in% names(survival_new_data)),
  identical(binary_cv$status, "success"),
  identical(binary_fit$status, "success"),
  identical(binary_prediction$status, "success"),
  identical(binary_evaluation$status, "success"),
  !"outcome" %in% names(binary_new_data)
)
