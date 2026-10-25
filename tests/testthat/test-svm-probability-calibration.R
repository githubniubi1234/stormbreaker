test_that("SVM calibration keeps paired patients and fold-local Top-K", {
  skip_if_not_installed("e1071")
  d <- sb_example_binary(100L, 3800L)
  d$outcome <- rep(c(0L, 1L), 50L)
  d$group <- rep(paste0("pair", seq_len(50L)), each = 2L)
  predictors <- c("gene_a", "gene_b", "gene_c")
  task <- sb_binary_task(d, "patient_id", "outcome", predictors,
                         group_col = "group")
  pipeline <- sb_pipeline(selector = "univ_binary_auc_top_k",
    selector_n_features = 2L, learner = "svm_linear_binary")
  fit <- sb_fit_final(task, pipeline)
  expect_identical(fit$status, "success")
  expect_identical(fit$tuning$status, "success")
  expect_identical(fit$tuning$inner_folds, 5L)
  expect_length(fit$tuning$inner_selections, 5L)
  for (audit in fit$tuning$inner_selections) {
    train <- match(audit$analysis_patient_ids, d$patient_id)
    assessment <- match(audit$assessment_patient_ids, d$patient_id)
    expect_length(intersect(d$group[train], d$group[assessment]), 0L)
    direct <- sb_fit_univ_binary_auc_top_k(d[train, predictors],
                                           d$outcome[train], 2L)
    expect_identical(audit$selected_features, direct$selected_features)
  }
  cv <- sb_run_cv(task, pipeline, folds = 3L, seed = 3801L)
  expect_identical(cv$status, "success")
  expect_length(cv$tuning, 3L)
  for (audit in c(cv$tuning, list(fit$tuning))) {
    rows <- match(names(audit$fold_id), d$patient_id)
    expect_identical(anyNA(rows), FALSE)
    group_folds <- split(unname(audit$fold_id), d$group[rows])
    expect_identical(unname(vapply(group_folds, function(x) {
      length(unique(x))
    }, integer(1))), rep(1L, length(group_folds)))
    expect_identical(audit$calibration$method, "libsvm_platt_oof")
  }
})

test_that("unscaled unscreened calibration reproduces native LIBSVM", {
  skip_if_not_installed("e1071")
  d <- sb_example_binary(100L, 3802L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  for (first in c(0L, 1L)) for (learner in c(
    "svm_linear_binary", "svm_radial_binary",
    "svm_polynomial_binary", "svm_sigmoid_binary"
  )) {
    rows <- c(which(d$outcome == first), which(d$outcome != first))
    params <- sb_classical_binary_params(learner)
    params$standardize <- FALSE
    native <- sb_fit_classical_binary(x[rows, ], d$outcome[rows],
                                             learner, params)
    tuned <- sb_calibrate_svm_probability(x[rows, ], d$outcome[rows],
      learner, params, d$patient_id[rows])
    expect_equal(tuned$calibration$A, unname(native$engine$probA),
                 tolerance = 1e-8)
    expect_equal(tuned$calibration$B, unname(native$engine$probB),
                 tolerance = 1e-8)
    expect_identical(tuned$calibration$positive_class, as.character(first))
  }
})

test_that("all five SVM variants use native fold-local margins and portable probabilities", {
  skip_if_not_installed("e1071")
  d <- sb_example_binary(100L, 3803L)
  d$outcome <- rep(c(0L, 1L, 0L, 0L), 25L)
  d$group <- rep(paste0("pair", seq_len(50L)), each = 2L)
  predictors <- c("gene_a", "gene_b", "gene_c")
  task <- sb_binary_task(d, "patient_id", "outcome", predictors,
                         group_col = "group")
  native_fit <- function(rows, features, params) {
    args <- list(x = as.matrix(d[rows, features, drop = FALSE]),
      y = factor(d$outcome[rows], levels = c(0L, 1L)), type = "C-classification",
      kernel = params$kernel, cost = params$cost, scale = TRUE,
      probability = FALSE, fitted = FALSE)
    if (!is.null(params$gamma)) args$gamma <- 1 / length(features)
    if (!is.null(params$degree)) args$degree <- params$degree
    if (!is.null(params$coef0)) args$coef0 <- params$coef0
    if (!is.null(params$class_weights)) {
      counts <- table(args$y)
      args$class.weights <- setNames(length(rows) / (2 * as.numeric(counts)),
                                      names(counts))
    }
    do.call(e1071::svm, args)
  }
  set.seed(3804L)
  before <- .Random.seed
  for (learner in sb_svm_probability_learner_ids()) {
    pipeline <- sb_pipeline(selector = "univ_binary_auc_top_k",
      selector_n_features = 2L, learner = learner)
    fit <- sb_fit_final(task, pipeline)
    expect_identical(fit$status, "success")
    expect_identical(.Random.seed, before)
    expected_decision <- numeric(nrow(d))
    for (inner in fit$tuning$inner_selections) {
      train <- match(inner$analysis_patient_ids, d$patient_id)
      assess <- match(inner$assessment_patient_ids, d$patient_id)
      features <- sb_fit_univ_binary_auc_top_k(d[train, predictors],
                                                d$outcome[train], 2L)$selected_features
      native <- native_fit(train, features, pipeline$learner_params)
      decision <- attr(predict(native, as.matrix(d[assess, features]),
                                decision.values = TRUE), "decision.values")
      sign <- if (native$levels[native$labels[1L]] ==
                  fit$tuning$calibration$positive_class) 1 else -1
      expected_decision[assess] <- sign * as.numeric(decision)
      expect_identical(inner$selected_features, features)
      expect_equal(inner$scale_center, native$x.scale[["scaled:center"]],
                   tolerance = 1e-14)
      expect_equal(inner$scale_scale, native$x.scale[["scaled:scale"]],
                   tolerance = 1e-14)
      if (!is.null(pipeline$learner_params$class_weights)) {
        counts <- table(factor(d$outcome[train], levels = c(0L, 1L)))
        expect_equal(inner$class_weights,
          setNames(length(train) / (2 * as.numeric(counts)), names(counts)),
          tolerance = 0)
      }
    }
    expect_equal(fit$tuning$calibration$out_of_fold$decision, expected_decision,
                 tolerance = 1e-12)
    positive <- d$outcome == as.integer(fit$tuning$calibration$positive_class)
    target <- ifelse(positive, (sum(positive) + 1) / (sum(positive) + 2),
                      1 / (sum(!positive) + 2))
    # Some native SVMs have almost constant margins. Condition the independent
    # optimizer's coordinates; this does not change production calibration.
    centered <- (expected_decision - mean(expected_decision)) / sd(expected_decision)
    objective <- function(beta) {
      z <- beta[1L] * centered + beta[2L]
      sum(pmax(z, 0) - (1 - target) * z + log1p(exp(-abs(z))))
    }
    gradient <- function(beta) {
      z <- beta[1L] * centered + beta[2L]
      residual <- target - plogis(-z)
      c(sum(centered * residual), sum(residual))
    }
    optimum <- optim(c(0, 0), objective, gradient, method = "BFGS",
                      control = list(reltol = 1e-13))
    expect_identical(optimum$convergence, 0L)
    transformed <- c(fit$tuning$calibration$A * sd(expected_decision),
      fit$tuning$calibration$B + fit$tuning$calibration$A * mean(expected_decision))
    expect_equal(transformed, optimum$par, tolerance = 1e-5)
    expect_equal(objective(transformed), optimum$value, tolerance = 1e-9)
    native <- native_fit(seq_len(nrow(d)), fit$required_features,
                           pipeline$learner_params)
    native$probA <- fit$tuning$calibration$A
    native$probB <- fit$tuning$calibration$B
    native$compprob <- TRUE
    expected <- unname(attr(predict(native, as.matrix(d[fit$required_features]),
                                     probability = TRUE), "probabilities")[, "1"])
    portable <- unserialize(serialize(fit, NULL))
    actual <- sb_predict_final(portable, d[c("patient_id", predictors)],
                                "patient_id")
    expect_identical(actual$status, "success")
    expect_equal(actual$predictions$probability, expected, tolerance = 1e-14)
    expect_identical(portable$learner_state$params, pipeline$learner_params)
  }
})

test_that("calibration training is unchanged by its held-out fold", {
  skip_if_not_installed("e1071")
  d <- sb_example_binary(100L, 3805L)
  predictors <- c("gene_a", "gene_b", "gene_c")
  learner <- "svm_radial_binary"
  params <- sb_classical_binary_params(learner)
  selector <- sb_univ_binary_auc_top_k_selector_adapter(2L)
  declared <- sb_make_svm_calibration_folds(d$outcome, params$seed)
  fold <- if (declared$fold_id[1L] == 1L) 2L else 1L
  local_mocked_bindings(sb_make_svm_calibration_folds = function(...) declared)
  models <- list()
  original <- sb_fit_classical_binary
  local_mocked_bindings(sb_fit_classical_binary = function(...) {
    model <- original(...)
    models[[length(models) + 1L]] <<- model
    model
  })
  baseline <- sb_calibrate_svm_probability(d[predictors], d$outcome, learner,
    params, d$patient_id, selector = selector)
  changed <- d
  assess <- declared$fold_id == fold
  changed[assess, predictors] <- changed[assess, predictors] + 500
  changed$outcome[assess] <- 1L - changed$outcome[assess]
  perturbed <- sb_calibrate_svm_probability(changed[predictors], changed$outcome,
    learner, params, d$patient_id, selector = selector)
  expect_identical(perturbed$inner_selections[[fold]],
                    baseline$inner_selections[[fold]])
  first <- models[[fold]]
  second <- models[[fold + 5L]]
  probe <- as.matrix(d[first$features])
  expect_equal(attr(predict(first$engine, probe, decision.values = TRUE),
                     "decision.values"),
    attr(predict(second$engine, probe, decision.values = TRUE), "decision.values"),
    tolerance = 0)
})

test_that("SVM calibration preserves fixed blocks and singleton native folds", {
  skip_if_not_installed("e1071")
  d <- sb_example_binary(60L, 3806L)
  predictors <- c("gene_a", "gene_b", "gene_c")
  task <- sb_binary_task(d, "patient_id", "outcome", predictors)
  for (selector in c("none", "fixed_features")) {
    pipeline <- if (selector == "none") {
      sb_pipeline(learner = "svm_radial_binary")
    } else sb_pipeline(selector = "fixed_features",
      selector_features = c("gene_c", "gene_a"), feature_set_id = "predeclared",
      learner = "svm_radial_binary")
    fit <- sb_fit_final(task, pipeline)
    expect_identical(fit$status, "success")
    for (inner in fit$tuning$inner_selections) {
      expect_identical(inner$selected_features, fit$required_features)
    }
    expect_identical(fit$tuning$inner_selection_method,
      if (selector == "none") "none" else "predeclared_feature_subset")
  }
  plain <- sb_make_svm_calibration_folds(d$outcome, 1L)
  singleton <- sb_make_svm_calibration_folds(d$outcome, 1L, d$patient_id)
  expect_identical(singleton, plain)
})

test_that("single-class calibration subproblems retain native constant margins", {
  skip_if_not_installed("e1071")
  d <- sb_example_binary(10L, 3807L)
  d$outcome <- c(1L, rep(0L, 9L))
  learner <- "svm_linear_binary"
  params <- sb_classical_binary_params(learner)
  params$standardize <- FALSE
  x <- d[c("gene_a", "gene_b")]
  native <- sb_fit_classical_binary(x, d$outcome, learner, params)
  tuned <- sb_calibrate_svm_probability(x, d$outcome, learner, params,
                                        d$patient_id)
  expect_equal(tuned$calibration$A, unname(native$engine$probA), tolerance = 1e-8)
  expect_equal(tuned$calibration$B, unname(native$engine$probB), tolerance = 1e-8)
  expect_equal(sum(vapply(tuned$inner_selections, function(x) {
    x$status == "single_class_constant"
  }, logical(1))), 1L)
})

test_that("infeasible grouped calibration is not silently replaced by row folds", {
  skip_if_not_installed("e1071")
  d <- sb_example_binary(20L, 3808L)
  d$outcome <- rep(c(0L, 1L), 10L)
  d$group <- rep(paste0("group", seq_len(4L)), each = 5L)
  task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b"),
                         group_col = "group")
  fit <- sb_fit_final(task, sb_pipeline(learner = "svm_linear_binary"))
  expect_identical(fit$status, "failed")
  expect_identical(fit$error_stage, "tuning")
  expect_identical(fit$error_class, "stormbreaker_infeasible_inner_folds")
  expect_identical(fit$tuning$status, "failed")
  expect_identical(fit$tuning$inner_folds, 5L)
  expect_null(fit$learner_state)
})

test_that("SVM calibration failures remain visible in CV and final fitting", {
  skip_if_not_installed("e1071")
  d <- sb_example_binary(60L, 3809L)
  task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b"))
  pipeline <- sb_pipeline(learner = "svm_linear_binary")
  local_mocked_bindings(sb_calibrate_svm_probability = function(...) {
    sb_abort("Synthetic calibration failure", "stormbreaker_synthetic_calibration")
  })
  final <- sb_fit_final(task, pipeline)
  cv <- sb_run_cv(task, pipeline, folds = 3L, seed = 3810L)
  expect_identical(final$status, "failed")
  expect_identical(final$error_stage, "tuning")
  expect_identical(cv$status, "failed")
  expect_length(cv$tuning, 3L)
  for (audit in c(cv$tuning, list(final$tuning))) {
    expect_identical(audit$status, "failed")
    expect_identical(audit$error_class, "stormbreaker_synthetic_calibration")
    expect_identical(audit$inner_folds, 5L)
  }
})

test_that("rare classes do not block otherwise feasible grouped calibration", {
  skip_if_not_installed("e1071")
  d <- sb_example_binary(40L, 3811L)
  d$outcome <- c(1L, 1L, rep(0L, 38L))
  d$group <- rep(paste0("group", seq_len(20L)), each = 2L)
  task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b"),
                         group_col = "group")
  fit <- sb_fit_final(task, sb_pipeline(learner = "svm_linear_binary"))
  expect_identical(fit$status, "success")
  expect_identical(sort(unique(unname(fit$tuning$fold_id))), seq_len(5L))
  for (inner in fit$tuning$inner_selections) {
    train <- match(inner$analysis_patient_ids, d$patient_id)
    assess <- match(inner$assessment_patient_ids, d$patient_id)
    expect_length(intersect(d$group[train], d$group[assess]), 0L)
    if (length(unique(d$outcome[train])) == 1L) {
      expect_identical(inner$status, "single_class_constant")
      oof <- fit$tuning$calibration$out_of_fold
      expect_identical(oof$decision[match(inner$assessment_patient_ids, oof$patient_id)],
                         rep(-1, length(assess)))
    }
  }
  # Exactly K groups must also yield K nonempty folds when only one has class 1.
  for (seed in 1:10) {
    y <- c(1L, 1L, rep(0L, 8L))
    folds <- sb_make_svm_calibration_folds(y, seed,
      rep(paste0("g", seq_len(5L)), each = 2L))
    expect_identical(sort(unique(folds$fold_id)), seq_len(5L))
    expect_identical(unname(tabulate(folds$fold_id)), rep(2L, 5L))
  }
  error <- capture_sb_error(sb_make_group_stratified_fold_id(
    d$outcome, d$group, 5L, 1L, "stormbreaker_infeasible_inner_folds"))
  expect_s3_class(error, "stormbreaker_infeasible_inner_folds")
})
