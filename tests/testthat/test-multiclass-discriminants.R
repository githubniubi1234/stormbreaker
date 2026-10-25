test_that("remaining classical engines preserve multiclass labels and saved predictions", {
  for (p in c("mda", "ipred", "sda", "klaR", "pamr", "kernlab", "kknn", "LiblineaR")) skip_if_not_installed(p)
  methods <- c("fda", "mda", "pda", "slda", "sda", "loclda", "pamr", "gaussian_process",
    "kernel_weighted_knn", "liblinear_l2_logistic", "liblinear_l1_logistic")
  for (k in c(3L, 4L)) {
    withr::local_seed(937)
    label <- rep(seq_len(k), each = 40)
    d <- data.frame(id = seq_along(label), outcome = paste0("Type ", label))
    for (j in seq_len(5)) d[[paste0("gene ", j)]] <- 2 * (label == j) + rnorm(length(label))
    task <- sb_multiclass_task(d, "id", "outcome", names(d)[-(1:2)],
      classes = rev(paste0("Type ", seq_len(k))))
    a <- sb_analyze_methods(task, methods)
    expect_identical(a$status$status, rep("success", length(methods)))
    for (method in methods) {
      fit <- a$results[[method]]$fit
      p <- sb_predict_multiclass(fit, d)
      expect_identical(colnames(p), task$classes)
      expect_equal(unname(rowSums(p)), rep(1, nrow(d)), tolerance = 1e-6)
      expect_gt(mean(task$classes[max.col(p)] == d$outcome), .7)
      rows <- c(nrow(d), 1, 41)
      expect_equal(unname(sb_predict_multiclass(fit, d[rows, ])), unname(p[rows, , drop = FALSE]))
      expect_equal(unname(sb_predict_multiclass(fit, d[rows[1], , drop = FALSE])),
        unname(p[rows[1], , drop = FALSE]))
      expect_equal(sb_predict_multiclass(unserialize(serialize(fit, NULL)), d[rows, ]),
        sb_predict_multiclass(fit, d[rows, ]))
    }
  }
  comparison <- sb_compare_multiclass(task, methods, folds = 3)
  expect_identical(comparison$comparison$status, rep("success", length(methods)))
})
