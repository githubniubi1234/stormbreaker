test_that("grouped glmnet keeps row identities and patients at every fit level", {
  skip_if_not_installed("glmnet", minimum_version = "4.1.2")
  predictors <- c("gene_a", "gene_b", "gene_c")
  for (type in c("binary", "survival")) {
    d <- if (type == "binary") sb_example_binary(120L, 3750L) else {
      sb_example_survival(120L, 3751L)
    }
    d$group <- rep(paste0("patient", 1:60), each = 2)
    if (type == "binary") d$outcome <- rep(c(0L, 1L), 60)
    task <- if (type == "binary") {
      sb_binary_task(d, "patient_id", "outcome", predictors, group_col = "group")
    } else {
      sb_surv_task(d, "patient_id", "time", "event", predictors, "month",
                   group_col = "group")
    }
    expect_identical(sb_subset_task(task, 1:100)$group_col, "group")
    for (selected in c(FALSE, TRUE)) {
      selector <- if (!selected) "none" else if (type == "binary") {
        "univ_binary_auc_top_k"
      } else "univ_cox_top_k"
      pipeline <- if (selected) {
        sb_pipeline(selector = selector, selector_n_features = 2L,
          learner = if (type == "binary") "glmnet_logistic" else "glmnet_cox",
          inner_folds = 3L, inner_seed = 3752L)
      } else sb_pipeline(
        learner = if (type == "binary") "glmnet_logistic" else "glmnet_cox",
        inner_folds = 3L, inner_seed = 3752L)
      cv <- sb_run_cv(task, pipeline, folds = 3L, seed = 3753L)
      final <- sb_fit_final(task, pipeline)
      expect_identical(cv$status, "success")
      expect_identical(cv$complete_oof, TRUE)
      expect_identical(final$status, "success")
      expect_identical(final$development$group_col, "group")
      expect_identical(unname(final$development$group_ids), d$group)
      expect_identical(cv$predictions$patient_id, d$patient_id)
      outer <- cv$provenance$fold_assignment
      for (i in seq_along(cv$tuning)) {
        expect_identical(sort(names(cv$tuning[[i]]$fold_id)),
          sort(d$patient_id[outer$fold != i]))
      }
      audits <- c(cv$tuning, list(final$tuning))
      for (audit in audits) {
        rows <- match(names(audit$fold_id), d$patient_id)
        expect_identical(anyNA(rows), FALSE)
        expect_equal(anyDuplicated(names(audit$fold_id)), 0L)
        by_group <- split(unname(audit$fold_id), d$group[rows])
        expect_identical(unname(vapply(by_group, function(x) {
          length(unique(x))
        }, integer(1))), rep(1L, length(by_group)))
        expect_length(audit$inner_selections, if (selected) 3L else 0L)
        if (selected) for (inner in audit$inner_selections) {
          train <- rows[audit$fold_id != inner$fold]
          assess <- rows[audit$fold_id == inner$fold]
          expect_identical(inner$analysis_patient_ids, d$patient_id[train])
          expect_length(intersect(d$group[train], d$group[assess]), 0L)
          expected <- if (type == "binary") {
            sb_fit_univ_binary_auc_top_k(d[train, predictors], d$outcome[train], 2L)
          } else sb_fit_univ_cox_top_k(
            d[train, predictors], d$time[train], d$event[train], 2L)
          expect_identical(inner$selected_features, expected$selected_features)
        }
      }
      if (!selected) {
        args <- list(x = as.matrix(d[predictors]),
          y = if (type == "binary") d$outcome else survival::Surv(d$time, d$event),
          family = if (type == "binary") "binomial" else "cox",
          alpha = pipeline$learner_params$alpha,
          standardize = pipeline$learner_params$standardize,
          foldid = unname(final$tuning$fold_id), type.measure = "deviance")
        if (type == "survival") args$cox.ties <- "breslow"
        native <- do.call(glmnet::cv.glmnet, args)
        expect_identical(final$tuning$lambda, unname(native$lambda.min))
        expect_identical(final$tuning$curve$mean_deviance, native$cvm)
      }
    }
  }
})

test_that("Top-K glmnet combinations are admitted with inner CV", {
  survival_pipeline <- sb_pipeline(
    selector = "univ_cox_top_k",
    learner = "glmnet_cox",
    selector_n_features = 2L,
    inner_folds = 3L,
    inner_seed = 821L
  )
  binary_pipeline <- sb_pipeline(
    selector = "univ_binary_auc_top_k",
    learner = "glmnet_logistic",
    selector_n_features = 2L,
    inner_folds = 3L,
    inner_seed = 822L
  )

  expect_identical(survival_pipeline$task_type, "survival")
  expect_identical(binary_pipeline$task_type, "binary")
  expect_identical(survival_pipeline$learner_params$inner_folds, 3L)
  expect_identical(binary_pipeline$learner_params$inner_folds, 3L)
  expect_identical(
    sb_lint_pipeline(survival_pipeline)$status,
    "admitted"
  )
  expect_identical(
    sb_lint_pipeline(binary_pipeline)$status,
    "admitted"
  )
})

test_that("Cox Top-K is refitted inside every glmnet inner analysis fold", {
  skip_if_not_installed("glmnet", minimum_version = "4.1.2")

  d <- sb_example_survival(n = 90L, seed = 20260822L)
  predictors <- c("gene_a", "gene_b", "gene_c")
  x <- d[predictors]
  patient_ids <- d$patient_id
  fold_id <- sb_make_stratified_fold_id(d$event, 3L, 823L)
  selector <- sb_univ_cox_top_k_selector_adapter(2L)
  fit_outer <- function(data) {
    sb_fit_runtime_selector(
      selector,
      data[predictors],
      list(
        task_type = "survival",
        time = data$time,
        event = data$event
      )
    )
  }
  tune <- function(data) {
    sb_select_glmnet_cox_lambda(
      data[predictors],
      data$time,
      data$event,
      inner_folds = 3L,
      seed = 823L,
      selector = selector,
      fitted_selector = fit_outer(data),
      patient_ids = patient_ids
    )
  }
  local_mocked_bindings(
    sb_make_stratified_fold_id = function(...) fold_id
  )
  original_cv_glmnet <- glmnet::cv.glmnet
  observed_cox_ties <- character()
  local_mocked_bindings(
    cv.glmnet = function(..., cox.ties = NULL) {
      observed_cox_ties <<- c(observed_cox_ties, cox.ties)
      original_cv_glmnet(..., cox.ties = cox.ties)
    },
    .package = "glmnet"
  )

  original <- tune(d)
  expect_identical(original$inner_selection_method, "glmnet_exclude_function")
  expect_identical(original$alignment, "lambda")
  expect_length(original$inner_selections, 3L)
  for (fold in seq_len(3L)) {
    direct <- sb_fit_runtime_selector(
      selector,
      x[fold_id != fold, , drop = FALSE],
      list(
        task_type = "survival",
        time = d$time[fold_id != fold],
        event = d$event[fold_id != fold]
      )
    )
    audit <- original$inner_selections[[fold]]
    expect_identical(audit$fold, fold)
    expect_identical(
      audit$analysis_patient_ids,
      patient_ids[fold_id != fold]
    )
    expect_identical(audit$selected_features, direct$selected_features)
  }

  changed <- d
  assessment_rows <- fold_id == 1L
  changed$time[assessment_rows] <- changed$time[assessment_rows] + 1
  changed$event[assessment_rows] <- 1L - changed$event[assessment_rows]
  changed[assessment_rows, predictors] <-
    changed[assessment_rows, predictors] + 100
  perturbed <- tune(changed)
  expect_identical(
    perturbed$inner_selections[[1L]],
    original$inner_selections[[1L]]
  )
  expect_identical(observed_cox_ties, rep("breslow", 2L))
})

test_that("Top-K glmnet pipelines run and retain compact inner audits", {
  skip_if_not_installed("glmnet", minimum_version = "4.1.2")

  survival_data <- sb_example_survival(n = 90L, seed = 824L)
  survival_task <- sb_surv_task(
    survival_data,
    "patient_id",
    "time",
    "event",
    c("gene_a", "gene_b", "gene_c"),
    "month"
  )
  survival_pipeline <- sb_pipeline(
    selector = "univ_cox_top_k",
    learner = "glmnet_cox",
    selector_n_features = 2L,
    inner_folds = 3L,
    inner_seed = 825L
  )
  survival_cv <- sb_run_cv(
    survival_task,
    survival_pipeline,
    folds = 3L,
    seed = 826L
  )
  survival_fit <- sb_fit_final(survival_task, survival_pipeline)

  expect_identical(survival_cv$status, "success")
  expect_true(survival_cv$complete_oof)
  expect_true(all(vapply(survival_cv$tuning, function(audit) {
    identical(audit$status, "success") &&
      identical(audit$inner_selection_method, "glmnet_exclude_function") &&
      length(audit$inner_selections) == 3L
  }, logical(1))))
  expect_identical(survival_fit$status, "success")
  expect_length(survival_fit$tuning$inner_selections, 3L)
  expect_length(survival_fit$required_features, 2L)

  binary_data <- sb_example_binary(n = 90L, seed = 827L)
  binary_task <- sb_binary_task(
    binary_data,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b", "gene_c")
  )
  binary_pipeline <- sb_pipeline(
    selector = "univ_binary_auc_top_k",
    learner = "glmnet_logistic",
    selector_n_features = 2L,
    inner_folds = 3L,
    inner_seed = 828L
  )
  binary_cv <- sb_run_cv(
    binary_task,
    binary_pipeline,
    folds = 3L,
    seed = 829L
  )
  binary_fit <- sb_fit_final(binary_task, binary_pipeline)
  prediction <- sb_predict_final(
    binary_fit,
    binary_data[1:4, c("patient_id", "gene_c", "gene_a", "gene_b")],
    "patient_id"
  )

  expect_identical(binary_cv$status, "success")
  expect_true(binary_cv$complete_oof)
  expect_true(all(vapply(binary_cv$tuning, function(audit) {
    identical(audit$status, "success") &&
      identical(audit$inner_selection_method, "glmnet_exclude_function") &&
      length(audit$inner_selections) == 3L
  }, logical(1))))
  expect_identical(binary_fit$status, "success")
  expect_length(binary_fit$tuning$inner_selections, 3L)
  expect_length(binary_fit$required_features, 2L)
  expect_identical(prediction$status, "success")
  expect_true(all(is.finite(prediction$predictions$probability)))
})

test_that("an inner Top-K failure remains an explicit tuning failure", {
  skip_if_not_installed("glmnet", minimum_version = "4.1.2")

  n <- 30L
  outcome <- rep(c(0L, 1L), length.out = n)
  fold_id <- rep(seq_len(3L), length.out = n)
  x <- data.frame(
    signal = seq_len(n),
    fold_one_only = ifelse(fold_id == 1L, seq_len(n), 0)
  )
  patient_ids <- sprintf("P%03d", seq_len(n))
  selector <- sb_univ_binary_auc_top_k_selector_adapter(2L)
  outer <- sb_fit_runtime_selector(
    selector,
    x,
    list(task_type = "binary", outcome = outcome)
  )
  local_mocked_bindings(
    sb_make_stratified_fold_id = function(...) fold_id
  )

  error <- capture_sb_error(sb_select_glmnet_logistic_lambda(
    x,
    outcome,
    inner_folds = 3L,
    seed = 830L,
    selector = selector,
    fitted_selector = outer,
    patient_ids = patient_ids
  ))

  expect_s3_class(error, "stormbreaker_insufficient_valid_features")
  expect_match(conditionMessage(error), "only 1 features")
})

test_that("outer survival selection stays visible when inner tuning fails", {
  d <- sb_example_survival(n = 60L, seed = 831L)
  task <- sb_surv_task(
    d,
    "patient_id",
    "time",
    "event",
    c("gene_a", "gene_b", "gene_c"),
    "month"
  )
  pipeline <- sb_pipeline(
    selector = "univ_cox_top_k",
    learner = "glmnet_cox",
    selector_n_features = 2L,
    inner_folds = 3L,
    inner_seed = 832L
  )
  local_mocked_bindings(
    sb_select_glmnet_cox_lambda = function(...) {
      sb_abort(
        "synthetic inner-selection failure",
        "stormbreaker_synthetic_inner_failure"
      )
    }
  )

  result <- sb_run_cv(task, pipeline, folds = 3L, seed = 833L)

  expect_identical(result$status, "failed")
  expect_true(all(vapply(result$selections, function(audit) {
    identical(audit$status, "success") &&
      length(audit$selected_features) == 2L
  }, logical(1))))
  expect_true(all(vapply(result$tuning, function(audit) {
    identical(audit$status, "failed") &&
      identical(
        audit$error_class,
        "stormbreaker_synthetic_inner_failure"
      )
  }, logical(1))))
})
