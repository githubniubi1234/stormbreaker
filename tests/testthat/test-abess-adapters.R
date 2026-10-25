test_that("ABESS components declare one fixed GIC specification per task", {
  expect_identical(
    sb_abess_learner_ids(),
    c("abess_binary_gic", "abess_cox_gic")
  )
  expect_identical(
    sb_abess_selector_ids(),
    c("abess_binary_gic_nonzero", "abess_cox_gic_nonzero")
  )

  binary <- sb_abess_learner_params("abess_binary_gic")
  survival <- sb_abess_learner_params("abess_cox_gic")
  expect_identical(binary$family, "binomial")
  expect_identical(binary$normalize, 2L)
  expect_identical(survival$family, "cox")
  expect_identical(survival$normalize, 3L)
  expect_identical(binary$tune_path, "sequence")
  expect_identical(binary$tune_type, "gic")
  expect_identical(binary$minimum_support, 0L)
  expect_identical(binary$maximum_support, 10L)
  expect_identical(binary$num_threads, 1L)
  expect_identical(binary$seed, 1L)
  expect_identical(sb_abess_dependency()$package, "abess")
})

test_that("ABESS learners match native predictions and portable coefficients", {
  skip_if_not_installed("abess", minimum_version = "0.4.11")

  binary_data <- make_sb_binary_fixture(n = 84L, seed = 8611L)
  binary_x <- binary_data[c("gene_a", "gene_b", "gene_c")]
  binary_new <- binary_x[1:7, , drop = FALSE]
  binary_adapter <- sb_abess_learner_adapter("abess_binary_gic")

  set.seed(17L)
  caller_rng <- .Random.seed
  binary_model <- binary_adapter$fit(binary_x, binary_data$outcome)
  expect_identical(.Random.seed, caller_rng)
  expect_s3_class(binary_model$engine, "abess")
  expect_null(binary_model$engine$call)

  binary_observed <- binary_adapter$predict(binary_model, binary_new)
  binary_native <- stats::predict(
    binary_model$engine,
    as.matrix(binary_new),
    type = "response"
  )
  binary_coefficients <- binary_adapter$extract_coefficients(binary_model)
  binary_manual <- stats::plogis(
    binary_coefficients$coefficient[[1L]] +
      as.vector(as.matrix(binary_new) %*%
        binary_coefficients$coefficient[-1L])
  )
  binary_state <- binary_adapter$extract_state(binary_model, list())
  binary_replay <- binary_adapter$predict_state(binary_state, binary_new)

  expect_equal(binary_observed, unname(as.numeric(binary_native)), tolerance = 1e-12)
  expect_equal(binary_observed, unname(binary_manual), tolerance = 1e-12)
  expect_equal(binary_replay$probability, binary_observed, tolerance = 1e-12)
  expect_null(binary_state$engine)
  expect_null(binary_state$call)

  survival_data <- make_sb_fixture(n = 90L, seed = 8612L)
  survival_x <- survival_data[c("gene_a", "gene_b", "gene_c")]
  survival_new <- survival_x[1:7, , drop = FALSE]
  survival_adapter <- sb_abess_learner_adapter("abess_cox_gic")
  survival_model <- survival_adapter$fit(
    survival_x,
    survival_data$time,
    survival_data$event
  )
  survival_observed <- survival_adapter$predict(survival_model, survival_new)
  survival_native <- stats::predict(
    survival_model$engine,
    as.matrix(survival_new),
    type = "link"
  )
  horizon <- stats::median(survival_data$time)
  survival_probability <- survival_adapter$predict_survival(
    survival_model,
    survival_x,
    survival_data$time,
    survival_data$event,
    survival_new,
    horizon
  )
  survival_state <- survival_adapter$extract_state(
    survival_model,
    list(
      x = survival_x,
      time = survival_data$time,
      event = survival_data$event,
      prediction_horizon = horizon,
      time_unit = "month"
    )
  )
  survival_replay <- survival_adapter$predict_state(
    survival_state,
    survival_new
  )

  expect_equal(survival_observed, unname(as.numeric(survival_native)), tolerance = 1e-12)
  expect_equal(survival_replay$risk, survival_observed, tolerance = 1e-12)
  expect_equal(
    survival_replay$survival_probability,
    survival_probability$probability,
    tolerance = 1e-12
  )
  expect_identical(survival_state$learner, "abess_cox_gic")
  expect_null(survival_state$engine)
  expect_null(survival_state$call)
})

test_that("ABESS Cox reports its native one-input boundary exactly", {
  skip_if_not_installed("abess", minimum_version = "0.4.11")
  survival_data <- make_sb_fixture(n = 90L, seed = 8613L)
  binary_data <- make_sb_binary_fixture(n = 84L, seed = 8614L)

  expect_error(
    sb_abess_learner_adapter("abess_cox_gic")$fit(
      survival_data["gene_a"],
      survival_data$time,
      survival_data$event
    ),
    class = "stormbreaker_abess_cox_too_few_predictors"
  )
  expect_s3_class(
    sb_abess_learner_adapter("abess_binary_gic")$fit(
      binary_data["gene_a"], binary_data$outcome
    ),
    "stormbreaker_abess_model"
  )
})

test_that("ABESS selectors retain only fold-local nonzero feature state", {
  skip_if_not_installed("abess", minimum_version = "0.4.11")

  binary_data <- make_sb_binary_fixture(n = 90L, seed = 8621L)
  binary_x <- binary_data[c("gene_a", "gene_b", "gene_c")]
  binary_adapter <- sb_abess_selector_adapter(
    "abess_binary_gic_nonzero"
  )
  binary_fit <- binary_adapter$fit(binary_x, binary_data$outcome)

  expect_s3_class(binary_fit, "stormbreaker_selector_fit")
  expect_gte(length(binary_fit$selected_features), 2L)
  expect_identical(
    binary_fit$selected_features,
    binary_fit$scores$feature[binary_fit$scores$selected]
  )
  expect_identical(
    names(binary_adapter$apply(binary_fit, binary_x)),
    binary_fit$selected_features
  )
  expect_null(binary_fit$engine)
  expect_null(binary_fit$call)

  survival_data <- make_sb_fixture(n = 90L, seed = 8622L)
  survival_x <- survival_data[c("gene_a", "gene_b", "gene_c")]
  survival_adapter <- sb_abess_selector_adapter(
    "abess_cox_gic_nonzero"
  )
  survival_fit <- survival_adapter$fit(
    survival_x,
    survival_data$time,
    survival_data$event
  )

  expect_gte(length(survival_fit$selected_features), 2L)
  expect_identical(
    names(survival_adapter$apply(survival_fit, survival_x)),
    survival_fit$selected_features
  )
  expect_null(survival_fit$engine)
  expect_null(survival_fit$call)
})

test_that("ABESS routes expand the registry without redundant self chains", {
  registry <- sb_pipeline_registry()
  restricted_internal_cv <- c(
    setdiff(sb_rule_aft_learner_ids(), "bart_aft_survival"),
    sb_slope_binary_learner_ids(),
    "cancerclass_binary"
  )
  base_binary <- registry$learner[
    registry$task_type == "binary" & registry$selector == "none"
  ]
  base_survival <- registry$learner[
    registry$task_type == "survival" & registry$selector == "none"
  ]
  expect_setequal(
    registry$learner[registry$selector == "abess_binary_gic_nonzero"],
    setdiff(base_binary, c("abess_binary_gic", restricted_internal_cv))
  )
  expect_setequal(
    registry$learner[registry$selector == "abess_cox_gic_nonzero"],
    setdiff(
      base_survival,
      c(
        "abess_cox_gic",
        sb_analysis_tuned_survival_learner_ids(),
        restricted_internal_cv
      )
    )
  )
  expect_identical(
    nrow(registry[
      registry$selector == "abess_binary_gic_nonzero" &
        registry$learner %in% c(
          "abess_binary_gic",
          restricted_internal_cv
        ),
      ,
      drop = FALSE
    ]),
    0L
  )
  expect_identical(
    nrow(registry[
      registry$selector == "abess_cox_gic_nonzero" &
        registry$learner %in% c(
          "abess_cox_gic",
          sb_analysis_tuned_survival_learner_ids(),
          restricted_internal_cv
        ),
      ,
      drop = FALSE
    ]),
    0L
  )

  binary <- sb_pipeline(
    selector = "abess_binary_gic_nonzero",
    learner = "glmnet_logistic",
    inner_folds = 3L,
    inner_seed = 1L
  )
  survival <- sb_pipeline(
    selector = "abess_cox_gic_nonzero",
    learner = "glmnet_cox",
    inner_folds = 3L,
    inner_seed = 1L
  )
  expect_null(binary$selector_params)
  expect_null(survival$selector_params)
  expect_match(binary$id, "abess_binary_gic_nonzero", fixed = TRUE)
  expect_match(survival$id, "abess_cox_gic_nonzero", fixed = TRUE)
})

test_that("ABESS learner and selector roles run through final fitting", {
  skip_if_not_installed("abess", minimum_version = "0.4.11")

  binary_data <- make_sb_binary_fixture(n = 84L, seed = 8631L)
  binary_task <- sb_binary_task(
    binary_data,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b", "gene_c")
  )
  binary_fit <- sb_fit_final(
    binary_task,
    sb_pipeline(learner = "abess_binary_gic")
  )
  binary_selected_fit <- sb_fit_final(
    binary_task,
    sb_pipeline(
      selector = "abess_binary_gic_nonzero",
      learner = "logistic"
    )
  )

  expect_identical(binary_fit$status, "success")
  expect_identical(binary_selected_fit$status, "success")
  expect_gte(length(binary_selected_fit$required_features), 2L)

  survival_data <- make_sb_fixture(n = 90L, seed = 8632L)
  horizon <- stats::median(survival_data$time)
  survival_task <- sb_surv_task(
    survival_data,
    "patient_id",
    "time",
    "event",
    c("gene_a", "gene_b", "gene_c"),
    "month",
    prediction_horizon = horizon
  )
  survival_fit <- sb_fit_final(
    survival_task,
    sb_pipeline(learner = "abess_cox_gic")
  )
  survival_selected_fit <- sb_fit_final(
    survival_task,
    sb_pipeline(
      selector = "abess_cox_gic_nonzero",
      learner = "coxph"
    )
  )

  expect_identical(survival_fit$status, "success")
  expect_identical(survival_selected_fit$status, "success")
  expect_gte(length(survival_selected_fit$required_features), 2L)
})

test_that("ABESS selection is refit inside every glmnet inner analysis fold", {
  skip_if_not_installed("abess", minimum_version = "0.4.11")
  skip_if_not_installed("glmnet", minimum_version = "4.1.2")

  binary_data <- make_sb_binary_fixture(n = 90L, seed = 8641L)
  binary_task <- sb_binary_task(
    binary_data,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b", "gene_c")
  )
  binary_result <- sb_run_cv(
    binary_task,
    sb_pipeline(
      selector = "abess_binary_gic_nonzero",
      learner = "glmnet_logistic",
      inner_folds = 3L,
      inner_seed = 8642L
    ),
    folds = 3L,
    seed = 8643L
  )

  survival_data <- make_sb_fixture(n = 90L, seed = 8644L)
  survival_task <- sb_surv_task(
    survival_data,
    "patient_id",
    "time",
    "event",
    c("gene_a", "gene_b", "gene_c"),
    "month"
  )
  survival_result <- sb_run_cv(
    survival_task,
    sb_pipeline(
      selector = "abess_cox_gic_nonzero",
      learner = "glmnet_cox",
      inner_folds = 3L,
      inner_seed = 8645L
    ),
    folds = 3L,
    seed = 8646L
  )

  for (result in list(binary_result, survival_result)) {
    expect_identical(result$status, "success")
    expect_true(result$complete_oof)
    expect_true(all(vapply(result$tuning, function(audit) {
      identical(audit$status, "success") &&
        identical(
          audit$inner_selection_method,
          "glmnet_exclude_function"
        ) &&
        length(audit$inner_selections) == 3L &&
        all(vapply(audit$inner_selections, function(selection) {
          length(selection$selected_features) >= 2L &&
            length(selection$analysis_patient_ids) > 0L
        }, logical(1)))
    }, logical(1))))
    expect_identical(
      result$provenance$abess_version,
      as.character(utils::packageVersion("abess"))
    )
  }
})
