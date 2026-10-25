for (id in c("adaboost_binary", "logitboost_binary",
  "abess_binary_gic", "abess_cox_gic")) {
  test_that(paste(id, "preserves native results with ordinary warnings"), {
    package <- switch(id, adaboost_binary = "fastAdaboost",
      logitboost_binary = "caTools", "abess")
    skip_if_not_installed(package)
    survival_task <- identical(id, "abess_cox_gic")
    d <- if (survival_task) make_sb_fixture(n = 48L, seed = 8612L) else {
      make_sb_binary_fixture(n = 48L, seed = 8611L)
    }
    x <- d[c("gene_a", "gene_b", "gene_c")]
    adapter <- if (startsWith(id, "abess")) sb_abess_learner_adapter(id) else {
      sb_classical_binary_boosting_learner_adapter(id, 5L)
    }
    fit <- function() {
      if (survival_task) adapter$fit(x, d$time, d$event) else adapter$fit(x, d$outcome)
    }
    reference <- fit()
    expected <- adapter$predict(reference, x)
    if (startsWith(id, "abess")) {
      native_fit <- abess::abess
      local_mocked_bindings(abess = function(...) {
        warning("informational fit warning")
        native_fit(...)
      }, .package = "abess")
    } else if (id == "adaboost_binary") {
      native_fit <- sb_fastadaboost_fit
      local_mocked_bindings(sb_fastadaboost_fit = function(...) {
        warning("informational fit warning")
        native_fit(...)
      })
    } else {
      native_fit <- sb_catools_logitboost_fit
      local_mocked_bindings(sb_catools_logitboost_fit = function(...) {
        warning("informational fit warning")
        native_fit(...)
      })
    }
    native_predict <- stats::predict
    local_mocked_bindings(predict = function(...) {
      warning("informational prediction warning")
      native_predict(...)
    }, .package = "stats")
    seen <- character()
    observed <- withCallingHandlers({
      model <- fit()
      adapter$predict(model, x)
    }, warning = function(w) {
      seen <<- c(seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
    expect_identical(observed, expected)
    expect_identical(all(c("informational fit warning", "informational prediction warning") %in% seen), TRUE)
  })
}

test_that("ABESS coercion and changed-model warnings remain failures", {
  messages <- c(
    "c.max should be an integer. It is coerced to as.integer(c.max).",
    "x should be a *numeric* matrix/data.frame! The factor value are coerced to as.numeric(x).",
    "y is a binary variable and is not match to family = 'multinomial' or 'ordinal'. We change to family = 'binomial'",
    "The number of the category of y is relative large compare to nvars. The numerical result might be unstable."
  )
  for (message in messages) {
    error <- tryCatch(sb_abess_stop_on_warning({warning(message); 42},
      "ABESS fit", "stormbreaker_abess_fit_warning"), error = identity)
    expect_s3_class(error, "stormbreaker_abess_fit_warning")
  }
})
