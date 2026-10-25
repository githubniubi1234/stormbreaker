test_that("multiclass LPC retains native scores and selected genes", {
  skip_if_not_installed("lpc")
  for (k in c(3L, 4L)) {
    withr::local_seed(190)
    labels <- rep(seq_len(k), each = 20)
    d <- as.data.frame(matrix(rnorm(length(labels) * 60), ncol = 60))
    d[1:4] <- d[1:4] + labels
    genes <- names(d)
    d$id <- seq_along(labels)
    d$outcome <- paste0("Type ", labels)
    task <- sb_multiclass_task(d, "id", "outcome", genes,
      classes = rev(paste0("Type ", seq_len(k))))
    a <- sb_analyze_methods(task, "lpc")
    expect_identical(a$status$status, "success", info = a$status$reason)
    native <- NULL
    invisible(capture.output(native <- withr::with_seed(1,
      lpc::LPC(t(as.matrix(d[genes])), as.integer(factor(d$outcome, levels = task$classes)),
        type = "multiclass"))))
    tab <- sb_method_features(a, "lpc")
    expect_equal(tab$lpc_score, as.numeric(native$lpcscores))
    expect_equal(tab$base_score, as.numeric(native$tscores))
    ranked <- order(-abs(native$lpcscores), seq_along(genes))
    expected <- genes[head(ranked, sb_high_dimensional_nlogn_size(nrow(d), length(genes)))]
    expect_setequal(sb_selected_features(a, "lpc"), expected)
    expect_identical(sb_selected_features(unserialize(serialize(a, NULL)), "lpc"),
      sb_selected_features(a, "lpc"))
    expect_identical(a$results$lpc$plot_status$status, "available")
  }
})
