test_that("PAMR selectors use native nonzero centroids for binary and survival tasks", {
  skip_if_not_installed("pamr")
  d <- sb_example_survival(180, 82)
  d$outcome <- d$event
  x <- d[c("gene_a", "gene_b", "gene_c")]
  names(x) <- c("gene A", "gene-B", "gene C")
  for (surv in c(FALSE, TRUE)) {
    id <- if (surv) "pamr_survival_nonzero" else "pamr_binary_nonzero"
    params <- sb_pamr_selector_params(id)
    a <- sb_pamr_selector_adapter(id, params)
    set.seed(19); before <- .Random.seed
    selected <- if (surv) a$fit(x, d$time, d$event) else a$fit(x, d$outcome)
    expect_identical(.Random.seed, before)
    data <- list(x = t(as.matrix(x)), geneid = names(x))
    if (surv) {
      data$survival.time <- d$time; data$censoring.status <- d$event
    } else data$y <- factor(d$outcome, levels = c(0, 1))
    invisible(capture.output(native <- pamr::pamr.train(data, threshold = 1,
      scale.sd = TRUE, offset.percent = 50, remove.zeros = FALSE,
      prior = if (surv) NULL else as.numeric(table(data$y) / nrow(x)),
      ngroup.survival = 2, sign.contrast = "both")))
    expected <- pamr::pamr.predict(native, data$x, threshold = 1, type = "nonzero")
    expect_identical(selected$selected_features, names(x)[expected])
    expect_identical(a$apply(unserialize(serialize(selected, NULL)), x), x[expected])
    expect_null(selected$engine)
    expect_null(selected$time)
    expect_null(selected$outcome)
    p <- sb_pipeline(id, if (surv) "coxph" else "logistic")
    expect_identical(sb_validate_pipeline_contract(p, p$task_type), p)
    features <- c("gene_a", "gene_b", "gene_c")
    task <- if (surv) sb_surv_task(d, "patient_id", "time", "event", features, "month") else
      sb_binary_task(d, "patient_id", "outcome", features)
    cv <- sb_run_cv(task, p, folds = 3, seed = 82)
    expect_identical(cv$status, "success")
    expect_identical(cv$provenance$pamr_version, as.character(packageVersion("pamr")))
    fit <- sb_fit_final(task, p)
    expect_identical(fit$status, "success")
    expect_identical(fit$provenance$pamr_version, as.character(packageVersion("pamr")))
    pred <- sb_predict_final(unserialize(serialize(fit, NULL)), d[c("patient_id", features)], "patient_id")
    expect_identical(pred$status, "success")
  }
})

test_that("PAMR selection is refitted within glmnet inner analysis folds", {
  skip_if_not_installed("pamr")
  skip_if_not_installed("glmnet")
  d <- withr::with_seed(82, {
    x <- matrix(rnorm(540), 180, 3)
    eta <- rowSums(x)
    t <- rexp(180, exp(eta) / 10)
    censor <- rexp(180, 0.02)
    data.frame(patient_id = paste0("p", seq_len(180)),
      time = pmin(t, censor), event = as.integer(t <= censor),
      outcome = rbinom(180, 1, plogis(eta)),
      gene_a = x[, 1], gene_b = x[, 2], gene_c = x[, 3])
  })
  features <- c("gene_a", "gene_b", "gene_c")
  original <- sb_pamr_selector_state
  seen <- list()
  local_mocked_bindings(sb_pamr_selector_state = function(engine, x, id, params) {
    seen[[length(seen) + 1L]] <<- d$patient_id[match(x[[1]], d$gene_a)]
    original(engine, x, id, params)
  })
  for (surv in c(FALSE, TRUE)) {
    task <- if (surv) sb_surv_task(d, "patient_id", "time", "event", features, "month") else
      sb_binary_task(d, "patient_id", "outcome", features)
    p <- sb_pipeline(if (surv) "pamr_survival_nonzero" else "pamr_binary_nonzero",
      if (surv) "glmnet_cox" else "glmnet_logistic", inner_folds = 3, inner_seed = 7)
    seen <- list()
    cv <- sb_run_cv(task, p, folds = 3, seed = 82)
    expect_gt(length(seen), 3L)
    expect_identical(any(lengths(seen) < 120L), TRUE)
    for (ids in seen) expect_identical(any(vapply(cv$folds,
      function(fold) all(ids %in% fold$analysis_ids), logical(1))), TRUE)
    expect_identical(cv$status, "success")
  }
})

test_that("PAMR does not pad a one-feature inner selection for glmnet", {
  skip_if_not_installed("pamr")
  skip_if_not_installed("glmnet")
  d <- sb_example_survival(180, 82); d$outcome <- d$event
  task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b", "gene_c"))
  cv <- sb_run_cv(task, sb_pipeline("pamr_binary_nonzero", "glmnet_logistic",
    inner_folds = 3, inner_seed = 7), folds = 3, seed = 82)
  expect_identical(cv$status, "failed")
  expect_identical(unique(cv$runs$error_class), "stormbreaker_glmnet_logistic_tuning_failure")
  expect_match(unique(cv$runs$error_message), "cannot retain 1 or less variables", fixed = TRUE)
})
