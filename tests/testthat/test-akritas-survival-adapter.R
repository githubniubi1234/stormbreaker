test_that("Akritas uses the native empirical CDF and native tied-time survival engine", {
  skip_if_not_installed("distr6")
  skip_if_not_installed("survivalmodels")
  withr::local_seed(296)
  x <- as.data.frame(matrix(rnorm(80 * 3), 80, 3))
  names(x) <- c("gene A", "gene-B", "gene C")
  rownames(x) <- paste0("PRIVATE_", 1:80)
  time <- ceiling(rexp(80, exp(x[[1]]) / 20))
  event <- rbinom(80, 1, .7)
  for (p in c(1, 3)) {
    a <- x[seq_len(p)]
    before <- .Random.seed
    fit <- sb_fit_akritas_survival(a, time, event, "akritas_survival", sb_akritas_params())
    expect_identical(.Random.seed, before)
    query <- a[c(2, 3, 7), , drop = FALSE]
    dist <- if (p == 1) distr6::Empirical$new(as.matrix(a)) else distr6::EmpiricalMV$new(as.matrix(a))
    ord <- order(time, decreasing = TRUE)
    expected <- getFromNamespace("C_Akritas", "survivalmodels")(
      cbind(time, event)[ord, ], fit$grid, fit$grid,
      as.numeric(dist$cdf(data = as.matrix(a)))[ord],
      as.numeric(dist$cdf(data = as.matrix(query))), .5)
    expect_identical(sb_akritas_curves(fit, query), expected)
    horizon <- as.numeric(median(time))
    state <- sb_extract_akritas_state(fit, list(prediction_horizon = horizon, time = time, time_unit = "month"))
    expect_null(rownames(state$basis))
    expect_null(rownames(state$truth))
    expect_null(state$engine)
    replay <- sb_predict_akritas_state(unserialize(serialize(state, NULL)), query)
    expect_identical(replay$risk, sb_predict_akritas_risk(fit, query))
    expect_identical(replay$survival_probability,
      as.numeric(sb_akritas_curves(fit, query, horizon)))
    expected_rmst <- as.numeric(cbind(1, expected)[,seq_along(fit$grid),drop=FALSE] %*% diff(c(0,fit$grid)))
    expect_identical(replay$risk, -expected_rmst)
  }
})

test_that("Akritas computes each distinct query CDF once and restores row order", {
  skip_if_not_installed("distr6")
  skip_if_not_installed("survivalmodels")
  withr::local_seed(961)
  x <- data.frame("gene A" = seq_len(60), "gene-B" = seq_len(60)^2,
    check.names = FALSE)
  time <- ceiling(rexp(60, 0.1))
  event <- rep(c(1L, 0L, 1L), 20L)
  native <- getFromNamespace("C_Akritas", "survivalmodels")
  original_lookup <- utils::getFromNamespace
  observed <- list()
  local_mocked_bindings(getFromNamespace = function(x, ns, ...) {
    f <- original_lookup(x, ns, ...)
    if (identical(x, "C_Akritas") && identical(ns, "survivalmodels")) {
      return(function(...) {
        arguments <- list(...)
        observed[[length(observed) + 1L]] <<- arguments[[5L]]
        f(...)
      })
    }
    f
  }, .package = "utils")
  for (p in c(1L, 2L)) {
    fit <- sb_fit_akritas_survival(x[seq_len(p)], time, event,
      "akritas_survival", sb_akritas_params())
    for (indices in list(c(7, 2, 7, 55, 2, 14), c(55, 7, 14), 7)) {
      query <- x[indices, seq_len(p), drop = FALSE]
      cdf <- sb_akritas_empirical_cdf(fit$basis, as.matrix(query))
      for (times in list(c(0, 2.5, max(time) + 1), 2.5)) {
        expected <- native(fit$truth, as.numeric(times), fit$grid,
          fit$cdf_train, cdf, fit$params$lambda)
        before <- .Random.seed
        actual <- sb_akritas_curves(fit, query, times)
        expect_identical(actual, expected)
        expect_identical(tail(observed, 1L)[[1L]], unique(cdf))
        expect_identical(.Random.seed, before)
      }
    }
  }
})

test_that("Akritas runs public CV and final state-only survival predictions", {
  skip_if_not_installed("distr6"); skip_if_not_installed("survivalmodels")
  d <- sb_example_survival(90, 82)
  features <- c("gene_a", "gene_b", "gene_c")
  task <- sb_surv_task(d, "patient_id", "time", "event", features, "month",
    prediction_horizon = as.numeric(median(d$time)))
  for (p in list(sb_pipeline("none", "akritas_survival"),
    sb_pipeline("univ_cox_top_k", "akritas_survival", selector_n_features = 2),
    sb_pipeline("fixed_features", "akritas_survival", selector_features = features[1:2], feature_set_id = "declared_two"))) {
    cv <- sb_run_cv(task, p, folds = 3, seed = 82)
    expect_identical(cv$status, "success")
    expect_identical(cv$complete_oof, TRUE)
    fit <- sb_fit_final(task, p)
    expect_identical(fit$status, "success")
    expect_identical(fit$provenance$distr6_version, as.character(packageVersion("distr6")))
    pred <- sb_predict_final(unserialize(serialize(fit, NULL)), d[c("patient_id", features)], "patient_id")
    expect_identical(pred$status, "success")
  }
})
