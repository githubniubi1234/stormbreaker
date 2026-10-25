test_that("multiclass information selectors reproduce native selected genes", {
  skip_if_not_installed("praznik")
  methods <- c("mrmr", "jmi", "cmim", "disr", "jmim")
  for (k in c(3L, 4L)) {
    withr::local_seed(122)
    y <- rep(seq_len(k), each = 15)
    d <- as.data.frame(matrix(rnorm(length(y) * 24), ncol = 24))
    d[[1]] <- y + rnorm(length(y), sd = .1)
    genes <- names(d)
    d$id <- seq_along(y)
    d$outcome <- paste0("Type ", y)
    task <- sb_multiclass_task(d, "id", "outcome", genes,
      classes = rev(paste0("Type ", seq_len(k))))
    a <- sb_analyze_methods(task, methods)
    expect_identical(a$status$status, rep("success", 5))
    for (method in methods) {
      f <- getExportedValue("praznik", toupper(method))
      args <- list(X = d[genes], Y = factor(d$outcome, levels = task$classes),
        k = sb_high_dimensional_nlogn_size(nrow(d), length(genes)), threads = 1L)
      if (method == "mrmr") args$positive <- FALSE
      native <- do.call(f, args)
      tab <- sb_method_features(a, method)
      expect_identical(tab$feature[order(tab$rank, na.last = NA)], genes[native$selection])
      expect_setequal(sb_selected_features(a, method), genes[native$selection])
      expect_identical(sb_selected_features(unserialize(serialize(a, NULL)), method),
        sb_selected_features(a, method))
      expect_identical(length(a$results[[method]]$plots), 1L)
    }
  }
})
