test_that("classical binary boosting settings are fixed and explicit", {
  adaboost <- sb_pipeline(learner = "adaboost_binary")
  logitboost <- sb_pipeline(learner = "logitboost_binary")

  expect_identical(
    adaboost$id,
    "none__adaboost_binary__iterations50"
  )
  expect_identical(
    logitboost$id,
    "none__logitboost_binary__iterations51"
  )
  expect_identical(adaboost$learner_params$method, "AdaBoost.M1")
  expect_identical(adaboost$learner_params$rpart_cp, 0)
  expect_identical(logitboost$learner_params$method, "LogitBoost")
  expect_identical(logitboost$learner_params$base_learner, "decision_stump")
  expect_identical(
    adaboost$learner_params$iteration_semantics,
    "maximum_until_native_convergence"
  )
  expect_identical(logitboost$learner_params$iteration_semantics, "exact")
  expect_identical(adaboost$learner_params$tuning, "none")
  expect_identical(logitboost$learner_params$tuning, "none")

  expect_error(
    sb_pipeline(
      learner = "adaboost_binary",
      adaboost_iterations = 0L
    ),
    class = "stormbreaker_bad_integer"
  )
  expect_error(
    sb_pipeline(learner = "logistic", logitboost_iterations = 11L),
    class = "stormbreaker_irrelevant_classical_binary_boosting_parameter"
  )
  expect_error(
    sb_pipeline(
      learner = "adaboost_binary",
      logitboost_iterations = 11L
    ),
    class = "stormbreaker_irrelevant_classical_binary_boosting_parameter"
  )
})

test_that("classical binary boosting preserves native class-one probability", {
  skip_if_not_installed("fastAdaboost", minimum_version = "1.0.0")
  skip_if_not_installed("caTools", minimum_version = "1.18.3")

  d <- make_sb_binary_fixture(n = 54L, seed = 7911L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  new_x <- x[c(4L, 19L, 37L), , drop = FALSE]
  adapters <- list(
    adaboost_binary =
      sb_classical_binary_boosting_learner_adapter(
        "adaboost_binary",
        9L
      ),
    logitboost_binary =
      sb_classical_binary_boosting_learner_adapter(
        "logitboost_binary",
        11L
      )
  )

  for (learner in names(adapters)) {
    adapter <- adapters[[learner]]
    model <- adapter$fit(x, d$outcome)
    observed <- adapter$predict(model, new_x)
    expected <- if (identical(learner, "adaboost_binary")) {
      native <- stats::predict(model$engine, newdata = new_x)$prob
      unname(native[, match("1", unname(model$engine$classnames))])
    } else {
      native <- stats::predict(
        model$engine,
        as.matrix(new_x),
        type = "raw",
        nIter = model$iterations
      )
      unname(native[, match("1", as.character(model$engine$lablist))])
    }
    audit <- adapter$extract_tree_model(model)
    state <- unserialize(serialize(adapter$extract_state(model), NULL))
    restored <- adapter$predict_state(state, new_x)$probability

    expect_equal(observed, expected, tolerance = 0)
    expect_equal(restored, expected, tolerance = 0)
    expect_identical(state$learner, learner)
    expect_false("engine" %in% names(state))
    expect_type(state$model_raw, "raw")
    expect_gt(length(state$model_raw), 0L)
    expect_identical(audit$iterations, adapter$params$iterations)
    expect_identical(audit$fitted_iterations, model$fitted_iterations)
    expect_lte(audit$fitted_iterations, audit$iterations)
    expect_identical(audit$importance_type, "split_count")
    expect_identical(audit$feature_importance$feature, names(x))
    expect_true(all(audit$feature_importance$importance >= 0))
  }
})

test_that("AdaBoost records native convergence before its iteration maximum", {
  skip_if_not_installed("fastAdaboost", minimum_version = "1.0.0")

  outcome <- rep(c(0L, 1L), each = 10L)
  x <- data.frame(
    separating_feature = ifelse(outcome == 1L, 1, -1),
    noise = seq_along(outcome)
  )
  adapter <- sb_classical_binary_boosting_learner_adapter(
    "adaboost_binary",
    5L
  )
  model <- adapter$fit(x, outcome)
  audit <- adapter$extract_tree_model(model)

  expect_identical(model$iterations, 5L)
  expect_identical(model$fitted_iterations, 1L)
  expect_identical(audit$fitted_iterations, 1L)
  expect_true(all(adapter$predict(model, x) %in% c(0, 1)))
})

test_that("classical binary boosting runs CV and portable final prediction", {
  skip_if_not_installed("fastAdaboost", minimum_version = "1.0.0")
  skip_if_not_installed("caTools", minimum_version = "1.18.3")

  d <- make_sb_binary_fixture(n = 72L, seed = 7921L)
  predictors <- c("gene_a", "gene_b", "gene_c")
  task <- sb_binary_task(d, "patient_id", "outcome", predictors)
  pipelines <- list(
    sb_pipeline(
      learner = "adaboost_binary",
      adaboost_iterations = 9L
    ),
    sb_pipeline(
      selector = "univ_binary_auc_top_k",
      learner = "logitboost_binary",
      selector_n_features = 2L,
      logitboost_iterations = 11L
    )
  )
  version_fields <- c("fastadaboost_version", "catools_version")

  for (index in seq_along(pipelines)) {
    pipeline <- pipelines[[index]]
    result <- sb_run_cv(task, pipeline, folds = 3L, seed = 7922L)
    fit <- sb_fit_final(task, pipeline)
    new_data <- d[
      c(5L, 29L, 61L),
      c("patient_id", rev(predictors)),
      drop = FALSE
    ]
    prediction <- sb_predict_final(fit, new_data, "patient_id")
    direct <- sb_predict_classical_binary_boosting_learner_state(
      fit$learner_state,
      new_data[fit$required_features]
    )

    expect_identical(result$status, "success")
    expect_identical(result$complete_oof, TRUE)
    expect_length(result$tree_models, 3L)
    expect_identical(
      vapply(result$tree_models, `[[`, "", "status"),
      rep("success", 3L)
    )
    expect_identical(
      vapply(result$coefficients, `[[`, "", "status"),
      rep("not_applicable", 3L)
    )
    expect_true(version_fields[[index]] %in% names(result$provenance))
    expect_identical(fit$status, "success")
    expect_true(version_fields[[index]] %in% names(fit$provenance))
    expect_identical(prediction$status, "success")
    expect_equal(
      prediction$predictions$probability,
      direct$probability,
      tolerance = 0
    )
  }
})

test_that("classical binary boosting has no dependency fallback", {
  d <- make_sb_binary_fixture(n = 30L, seed = 7931L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  adaboost <- sb_classical_binary_boosting_learner_adapter(
    "adaboost_binary",
    3L
  )
  logitboost <- sb_classical_binary_boosting_learner_adapter(
    "logitboost_binary",
    3L
  )
  local_mocked_bindings(
    sb_fastadaboost_available = function() FALSE,
    sb_catools_available = function() FALSE
  )

  expect_error(
    adaboost$fit(x, d$outcome),
    class = "stormbreaker_missing_fastadaboost"
  )
  expect_error(
    logitboost$fit(x, d$outcome),
    class = "stormbreaker_missing_catools"
  )
})

test_that("classical binary boosting rejects invalid fits despite warnings", {
  skip_if_not_installed("fastAdaboost", minimum_version = "1.0.0")
  skip_if_not_installed("caTools", minimum_version = "1.18.3")

  d <- make_sb_binary_fixture(n = 30L, seed = 7941L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  adaboost <- sb_classical_binary_boosting_learner_adapter(
    "adaboost_binary",
    3L
  )
  logitboost <- sb_classical_binary_boosting_learner_adapter(
    "logitboost_binary",
    3L
  )
  local_mocked_bindings(
    sb_fastadaboost_fit = function(...) {
      warning("injected AdaBoost warning")
    },
    sb_catools_logitboost_fit = function(...) {
      warning("injected LogitBoost warning")
    }
  )

  expect_error(
    suppressWarnings(adaboost$fit(x, d$outcome)),
    class = "stormbreaker_invalid_classical_binary_boosting_fit"
  )
  expect_error(
    suppressWarnings(logitboost$fit(x, d$outcome)),
    class = "stormbreaker_invalid_classical_binary_boosting_fit"
  )
})
