test_that("captured selectors retain classic paths without changing selections", {
  withr::local_seed(81)
  x <- as.data.frame(matrix(rnorm(300 * 6), 300))
  names(x) <- paste0("gene", 1:6)
  y <- rbinom(300, 1, plogis(2 * x$gene1 - 2 * x$gene2))
  for (id in c("lasso_binary_nonzero", "enet_binary_nonzero", "rf_binary_importance")) {
    p <- sb_binary_model_selector_options(id)
    compact <- sb_fit_binary_model_selector(x, y, id, p)
    captured <- sb_captured_selector_adapter(id, p)$fit(x, y)
    expect_identical(compact$selected_features, captured$selected_features)
    expect_equal(compact$scores, captured$scores)
    if (id == "rf_binary_importance") {
      expect_equal(nrow(captured$evidence$oob_error), p$trees)
      expect_null(compact$evidence$oob_error)
      plot <- sb_selector_oob_plot(captured$evidence$oob_error)
    } else {
      path <- captured$evidence$coefficient_path
      expect_identical(rownames(path$beta), names(x))
      expect_gt(length(path$lambda), 10)
      expect_null(compact$evidence$coefficient_path)
      plot <- sb_selector_path_plot(path, captured$evidence$lambda)
    }
    expect_s3_class(ggplot2::ggplotGrob(plot), "gtable")
  }
})

test_that("Cox Lasso and Enet capture the fitted native path", {
  d <- sb_example_survival(180, 41)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  for (id in c("lasso_nonzero", "enet_cox_nonzero")) {
    s <- sb_captured_selector_adapter(id, list(inner_folds = 5L, seed = 1L))$fit(x, d$time, d$event)
    expect_gt(length(s$coefficient_path$lambda), 10)
    expect_identical(rownames(s$coefficient_path$beta), names(x))
    expect_s3_class(ggplot2::ggplotGrob(sb_selector_path_plot(s$coefficient_path, s$lambda)), "gtable")
  }
})

test_that("additional classic selector plots reach the doctor entrypoint", {
  withr::local_seed(81)
  d <- as.data.frame(matrix(rnorm(300 * 6), 300)); names(d) <- paste0("gene", 1:6)
  d$id <- paste0("p", 1:300); d$y <- rbinom(300, 1, plogis(2*d$gene1 - 2*d$gene2))
  task <- sb_binary_task(d, "id", "y", paste0("gene", 1:6))
  for (selector in c("ABESS", "StepGLM", "glmBoost")) {
    f <- sb_fit_final(task, bi_pipeline("RF", bi_selector(selector)), capture_diagnostics = TRUE)
    expect_identical(f$status, "success")
    z <- sb_plots(f)
    names <- grep("^selector_", names(z$plots), value = TRUE)
    expect_gte(length(names), 2)
    for (nm in names) expect_s3_class(ggplot2::ggplotGrob(z$plots[[nm]]), "gtable")
  }
  d <- sb_example_survival(180, 41)
  task <- sb_surv_task(d, "patient_id", "time", "event", c("gene_a", "gene_b", "gene_c"), "month")
  f <- sb_fit_final(task, sb_pipeline("coxboost_then_stepcox", "coxph", stepcox_direction = "both",
    selector_seed = 1L, selector_inner_folds = 5L), capture_diagnostics = TRUE)
  expect_identical(f$status, "success")
  z <- sb_plots(f)
  for (nm in c("selector_1_cv", "selector_1_coefficient_path")) {
    expect_s3_class(z$plots[[nm]], "ggplot")
    expect_s3_class(ggplot2::ggplotGrob(z$plots[[nm]]), "gtable")
  }
})
