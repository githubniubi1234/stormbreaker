test_that("survival CV uses each bagged-kNN assessment curve once", {
  skip_if_not_installed("bnnSurvival")
  skip_if_not_installed("survivalmodels")
  skip_if_not_installed("distr6")
  d <- sb_example_survival(90L, 82L)
  task <- sb_surv_task(d, "patient_id", "time", "event",
    c("gene_a", "gene_b", "gene_c"), "month",
    prediction_horizon = as.numeric(median(d$time)))
  bagged <- sb_predict_bagged_knn_survival_curves
  akritas <- sb_akritas_curves
  calls <- c(bagged = 0L, akritas = 0L)
  local_mocked_bindings(
    sb_predict_bagged_knn_survival_curves = function(...) {
      calls[["bagged"]] <<- calls[["bagged"]] + 1L
      bagged(...)
    },
    sb_akritas_curves = function(...) {
      calls[["akritas"]] <<- calls[["akritas"]] + 1L
      akritas(...)
    })
  for (id in c("bagged_knn_survival", "akritas_survival")) {
    result <- sb_run_cv(task, sb_pipeline("none", id), folds = 3L, seed = 82L)
    expect_identical(result$complete_oof, TRUE)
    expect_identical(result$complete_survival_probability_oof, TRUE)
  }
  expect_identical(calls, c(bagged = 3L, akritas = 6L))
})

test_that("shared assessment predictions preserve complete CV and horizon-only failures", {
  skip_if_not_installed("bnnSurvival")
  skip_if_not_installed("survivalmodels")
  skip_if_not_installed("distr6")
  withr::local_seed(17L)
  d <- sb_example_survival(90L, 82L)
  # One unique maximum makes exactly its assessment fold unable to support h.
  d$time[1L] <- max(d$time) + 1
  features <- c("gene_a", "gene_b", "gene_c")
  bagged <- sb_bagged_knn_survival_learner_adapter
  akritas <- sb_akritas_learner_adapter
  plain <- function(task, pipeline) {
    with_mocked_bindings(sb_run_cv(task, pipeline, folds = 3L, seed = 82L),
      sb_bagged_knn_survival_learner_adapter = function(...) {
        a <- bagged(...); a$predict_details <- NULL; a
      },
      sb_akritas_learner_adapter = function(...) {
        a <- akritas(...); a$predict_details <- NULL; a
      })
  }
  before <- .Random.seed
  for (h in c(median(d$time), max(d$time))) {
    task <- sb_surv_task(d, "patient_id", "time", "event", features, "month",
      prediction_horizon = as.numeric(h))
    for (id in c("bagged_knn_survival", "akritas_survival")) {
      pipelines <- list(sb_pipeline("none", id),
        sb_pipeline("univ_cox_top_k", id, selector_n_features = 2L),
        sb_pipeline("fixed_features", id, selector_features = features[1:2], feature_set_id = "two"))
      for (p in pipelines) {
        expected <- plain(task, p)
        actual <- sb_run_cv(task, p, folds = 3L, seed = 82L)
        expect_identical(actual, expected)
        expect_identical(actual$complete_oof, TRUE)
        expect_identical(actual$complete_survival_probability_oof, h < max(d$time))
        if (h == max(d$time)) {
          unsupported <- vapply(actual$baseline_hazards, function(item) {
            identical(item$error_class, "stormbreaker_prediction_horizon_unsupported")
          }, logical(1))
          expect_identical(sum(unsupported), 1L)
        }
        expect_identical(.Random.seed, before)
      }
    }
  }
  d$gene_b <- d$gene_a
  singular <- sb_surv_task(d, "patient_id", "time", "event", features[1:2], "month",
    prediction_horizon = as.numeric(median(d$time)))
  p <- sb_pipeline("none", "bagged_knn_survival")
  expected <- plain(singular, p)
  actual <- sb_run_cv(singular, p, folds = 3L, seed = 82L)
  expect_identical(actual, expected)
  expect_identical(actual$complete_oof, FALSE)
})

test_that("bagged-kNN shared curves preserve direct horizon queries", {
  skip_if_not_installed("bnnSurvival")
  d <- sb_example_survival(72L, 42L)
  d$time <- rep(seq_len(12L), length.out = nrow(d))
  for (p in c(1L, 3L)) {
    x <- d[c("gene_a", "gene_b", "gene_c")][seq_len(p)]
    adapter <- sb_bagged_knn_survival_learner_adapter("bagged_knn_survival")
    model <- adapter$fit(x, d$time, d$event)
    query <- x[c(1L, 9L, 30L), , drop = FALSE]
    details <- adapter$predict_details(model, query)
    expect_identical(details$risk, adapter$predict(model, query))
    for (h in c(0.5, 1, 1.5, 6, 6.75, 12)) {
      expected <- adapter$predict_survival(model, x, d$time, d$event, query, h)
      actual <- adapter$predict_survival(model, x, d$time, d$event, query, h,
        prediction_details = details)
      expect_identical(actual, expected)
    }
  }
})

test_that("nonfatal curve warnings are preserved and prevent prediction reuse", {
  skip_if_not_installed("bnnSurvival")
  d <- sb_example_survival(72L, 42L)
  task <- sb_surv_task(d, "patient_id", "time", "event",
    c("gene_a", "gene_b", "gene_c"), "month",
    prediction_horizon = as.numeric(median(d$time)))
  original <- sb_predict_bagged_knn_survival_curves
  factory <- sb_bagged_knn_survival_learner_adapter
  local_mocked_bindings(sb_predict_bagged_knn_survival_curves = function(...) {
    warning("prediction warning", call. = FALSE)
    original(...)
  })
  run <- function(shared) {
    warnings <- character()
    value <- withCallingHandlers(with_mocked_bindings(
      sb_run_cv(task, sb_pipeline("none", "bagged_knn_survival"), folds = 3L, seed = 82L),
      sb_bagged_knn_survival_learner_adapter = function(...) {
        a <- factory(...)
        if (!shared) a$predict_details <- NULL
        a
      }), warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      })
    list(value = value, warnings = warnings)
  }
  actual <- run(TRUE)
  expected <- run(FALSE)
  expect_identical(actual, expected)
  expect_identical(actual$warnings, rep("prediction warning", 6L))
})

test_that("risk-only CV does not request shared prediction details", {
  skip_if_not_installed("bnnSurvival")
  d <- sb_example_survival(72L, 42L)
  task <- sb_surv_task(d, "patient_id", "time", "event", c("gene_a", "gene_b"), "month")
  factory <- sb_bagged_knn_survival_learner_adapter
  local_mocked_bindings(sb_bagged_knn_survival_learner_adapter = function(...) {
    a <- factory(...)
    a$predict_details <- function(...) stop("Risk-only CV must not request details")
    a
  })
  result <- sb_run_cv(task, sb_pipeline("none", "bagged_knn_survival"), folds = 3L, seed = 82L)
  expect_identical(result$complete_oof, TRUE)
})

test_that("native bagged-kNN prediction errors are not converted to successful risks", {
  skip_if_not_installed("bnnSurvival")
  d <- sb_example_survival(72L, 42L)
  task <- sb_surv_task(d, "patient_id", "time", "event", c("gene_a", "gene_b"), "month",
    prediction_horizon = as.numeric(median(d$time)))
  p <- sb_pipeline("none", "bagged_knn_survival")
  factory <- sb_bagged_knn_survival_learner_adapter
  local_mocked_bindings(sb_predict_bagged_knn_engine = function(...) {
    warning("native engine warning", call. = FALSE)
    stop("native engine prediction failed")
  })
  actual <- sb_run_cv(task, p, folds = 3L, seed = 82L)
  expected <- with_mocked_bindings(sb_run_cv(task, p, folds = 3L, seed = 82L),
    sb_bagged_knn_survival_learner_adapter = function(...) {
      a <- factory(...); a$predict_details <- NULL; a
    })
  expect_identical(actual, expected)
  expect_identical(actual$complete_oof, FALSE)
  expect_identical(actual$runs$error_class, rep("simpleError", 3L))
  expect_identical(actual$runs$status, rep("failed", 3L))
  expect_true(all(grepl("native engine prediction failed", actual$runs$error_message, fixed = TRUE)))
})

test_that("ordinary bagged-kNN warnings preserve valid shared and plain predictions", {
  skip_if_not_installed("bnnSurvival")
  d <- sb_example_survival(72L, 42L)
  task <- sb_surv_task(d, "patient_id", "time", "event", c("gene_a", "gene_b"), "month",
    prediction_horizon = as.numeric(median(d$time)))
  p <- sb_pipeline("none", "bagged_knn_survival")
  native <- sb_predict_bagged_knn_engine
  factory <- sb_bagged_knn_survival_learner_adapter
  baseline <- sb_run_cv(task, p, folds = 3L, seed = 82L)
  local_mocked_bindings(sb_predict_bagged_knn_engine = function(...) {
    warning("native engine warning", call. = FALSE)
    native(...)
  })
  run <- function(shared) {
    messages <- character()
    value <- withCallingHandlers({
      if (shared) sb_run_cv(task, p, folds = 3L, seed = 82L) else
        with_mocked_bindings(sb_run_cv(task, p, folds = 3L, seed = 82L),
          sb_bagged_knn_survival_learner_adapter = function(...) {
            a <- factory(...); a$predict_details <- NULL; a
          })
    }, warning = function(w) {
      messages <<- c(messages, conditionMessage(w)); invokeRestart("muffleWarning")
    })
    list(value = value, messages = messages)
  }
  shared <- run(TRUE)
  plain <- run(FALSE)
  expect_identical(shared$value, plain$value)
  expect_identical(shared$value, baseline)
  expect_true(length(shared$messages) > 0L)
  expect_true(length(plain$messages) > 0L)
  expect_true(all(shared$messages == "native engine warning"))
  expect_true(all(plain$messages == "native engine warning"))
  expect_identical(shared$value$complete_oof, TRUE)
})
