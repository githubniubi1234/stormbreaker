test_that("saved score plots preserve native values and missing evidence", {
  ids <- c("rsf_minimal_depth", "sis_binary_nlogn", "mrmr_binary_nlogn",
    "jmi_binary_nlogn", "cmim_binary_nlogn", "disr_binary_nlogn", "jmim_binary_nlogn",
    "lpc_binary_nlogn", "lpc_cox_nlogn", "stabsel_binary_lasso", "stabsel_cox_lasso")
  for (id in ids) for (metric in names(sb_selector_score_views(id))) {
    scores <- data.frame(feature = c("a", "b", "c"), selected = c(TRUE, FALSE, FALSE))
    scores[[metric]] <- c(.9, .2, NA_real_)
    if (metric == "rank") scores[[metric]] <- c(1L, 2L, NA_integer_)
    s <- structure(list(id = id, scores = scores, params = list(cutoff = .75)), class = "stormbreaker_selector_fit")
    p <- sb_plot_selector_scores(s, metric)
    expect_identical(attr(p, "selector_scores"), scores)
    expect_setequal(p$data$value, scores[[metric]][1:2])
    expect_equal(nrow(p$data), 2)
    expect_s3_class(ggplot2::ggplotGrob(p), "gtable")
    if (metric == "selection_probability") expect_equal(p$layers[[2]]$data$xintercept, .75)
  }
})

test_that("doctor plots expose screening views without refitting", {
  d <- sb_example_binary(100, 11)
  t <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b", "gene_c"))
  f <- sb_fit_final(t, bi_pipeline("RF", bi_selector("mRMR")), capture_diagnostics = TRUE)
  expect_identical(f$status, "success")
  before <- f$selection$state
  p <- sb_plots(f)
  expect_s3_class(p$plots$selector_1_mrmr_score, "ggplot")
  expect_identical(f$selection$state, before)
  expect_null(p$unavailable$selector_1)
})

test_that("univariable selector plots expose the saved ranking statistics", {
  fits <- list()
  for (type in c("binary", "survival")) {
    d <- if (type == "binary") sb_example_binary(100, 21) else sb_example_survival(100, 21)
    genes <- c("gene_a", "gene_b", "gene_c")
    task <- if (type == "binary") sb_binary_task(d, "patient_id", "outcome", genes) else
      sb_surv_task(d, "patient_id", "time", "event", genes, "month")
    p <- sb_pipeline(selector = if (type == "binary") "univ_binary_auc_top_k" else "univ_cox_top_k",
      learner = if (type == "binary") "logistic" else "coxph", selector_n_features = 2L)
    fits[[type]] <- sb_fit_final(task, p, capture_diagnostics = TRUE)
  }
  local_mocked_bindings(sb_fit_final = function(...) stop("Plotting must not refit"))
  for (type in names(fits)) {
    fit <- fits[[type]]
    expect_identical(fit$status, "success")
    metric <- if (type == "binary") "discrimination" else "wald_z"
    key <- paste0("selector_1_", metric)
    out <- sb_plots(fit)
    expect_s3_class(out$plots[[key]], "ggplot")
    expect_identical(attr(out$plots[[key]], "selector_scores"), fit$selection$state$scores)
    expect_setequal(out$plots[[key]]$data$value, fit$selection$state$scores[[metric]])
    expect_null(out$unavailable$selector_1)
    expect_s3_class(ggplot2::ggplotGrob(out$plots[[key]]), "gtable")
  }
})
