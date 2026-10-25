test_that("multinomial model families preserve class order and native saved predictions", {
  for (p in c("nnet", "brglm2", "mgcv", "pre")) skip_if_not_installed(p)
  methods <- c("logistic", "stepglm", "brglm2", "gam", "pre")
  for (k in c(3L, 4L)) {
    withr::local_seed(543)
    labels <- rep(seq_len(k), each = 35)
    d <- data.frame(id = seq_along(labels), label = paste0("Type ", labels))
    for (j in seq_len(3)) d[[paste0("gene ", j)]] <- 2 * (labels == j) + rnorm(length(labels))
    task <- sb_multiclass_task(d, "id", "label", names(d)[-(1:2)], classes = rev(paste0("Type ", seq_len(k))))
    a <- sb_analyze_methods(task, methods, inner_folds = 3)
    expect_identical(a$status$status, rep("success", length(methods)),
      info = paste(a$status$method, a$status$reason, collapse = "; "))
    for (method in methods) {
      fit <- a$results[[method]]$fit
      p <- sb_predict_multiclass(fit, d)
      expect_identical(colnames(p), task$classes)
      expect_equal(unname(rowSums(p)), rep(1, nrow(d)), tolerance = 1e-6)
      expect_gt(mean(task$classes[max.col(p)] == d$label), .65)
      rows <- c(nrow(d), 1, 36)
      expect_equal(unname(sb_predict_multiclass(fit, d[rows, ])), unname(p[rows, , drop = FALSE]))
      expect_equal(unname(sb_predict_multiclass(fit, d[1, , drop = FALSE])), unname(p[1, , drop = FALSE]))
      expect_equal(unname(sb_predict_multiclass(unserialize(serialize(fit, NULL)), d[rows, ])),
        unname(p[rows, , drop = FALSE]))
    }
  }
  cv <- sb_compare_multiclass(task, methods, folds = 3, inner_folds = 3)
  expect_identical(cv$comparison$status, rep("success", length(methods)),
    info = paste(cv$comparison$method, cv$comparison$reason, collapse = "; "))
})

test_that("unverified SLOPE is not advertised or fitted", {
  catalog <- sb_method_catalog("multiclass")
  expect_identical(catalog$supported[catalog$method == "slope"], FALSE)
  d <- datasets::iris
  d$id <- seq_len(nrow(d))
  task <- sb_multiclass_task(d, "id", "Species", names(d)[1:4])
  a <- sb_analyze_methods(task, "slope")
  expect_identical(a$status$status, "failed")
  expect_identical(a$results$slope$error_class, "stormbreaker_unsupported_method_task")
})
