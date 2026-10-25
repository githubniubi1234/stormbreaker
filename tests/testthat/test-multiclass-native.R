test_that("native probability engines support three and four classes", {
  for (p in c("e1071", "klaR", "nnet", "glmnet")) skip_if_not_installed(p)
  methods <- c("svm_linear", "svm_radial", "svm_radial_weighted", "svm_polynomial",
    "svm_sigmoid", "kernel_naive_bayes", "rda", "nnet", "knn", "ridge")
  for (k in c(3L, 4L)) {
    withr::local_seed(991)
    label <- rep(seq_len(k), each = 35L)
    d <- data.frame(id = seq_along(label), a = label + rnorm(length(label)),
      b = sin(label) + rnorm(length(label)), c = rnorm(length(label)),
      outcome = paste0("Type ", label))
    task <- sb_multiclass_task(d, "id", "outcome", c("a", "b", "c"),
      classes = rev(paste0("Type ", seq_len(k))))
    a <- sb_analyze_methods(task, methods, inner_folds = 3)
    expect_identical(a$status$status, rep("success", length(methods)))
    for (method in methods) {
      fit <- a$results[[method]]$fit
      rows <- c(nrow(d), 1, 36)
      p <- sb_predict_multiclass(fit, d[rows, ])
      expect_identical(colnames(p), task$classes)
      expect_equal(unname(rowSums(p)), rep(1, 3), tolerance = 1e-6)
      expect_equal(unname(sb_predict_multiclass(fit, d[rows[1], , drop = FALSE])),
        unname(p[1, , drop = FALSE]))
      expect_equal(sb_predict_multiclass(unserialize(serialize(fit, NULL)), d[rows, ]), p)
      expect_s3_class(tryCatch(sb_selected_features(a, method), error = identity),
        "stormbreaker_selection_not_defined")
      if (startsWith(method, "svm_")) {
        native <- attr(predict(fit$engine, d[rows, task$predictors], probability = TRUE), "probabilities")
        expect_equal(unname(p), unname(native[, fit$codes, drop = FALSE]))
        expect_equal(sum(a$results[[method]]$diagnostics$support_vectors$count), fit$engine$tot.nSV)
        counts <- table(factor(d$outcome[fit$engine$index], levels = task$classes))
        tab <- a$results[[method]]$diagnostics$support_vectors
        expect_equal(unname(tab$count[match(task$classes, as.character(tab$class))]), unname(as.integer(counts)))
      } else if (method == "rda") {
        native <- get("predict.rda", asNamespace("klaR"))(fit$engine, d[rows, task$predictors])$posterior
        expect_equal(unname(p), unname(native[, fit$codes, drop = FALSE]))
      } else if (method == "kernel_naive_bayes") {
        native <- predict(fit$engine, d[rows, task$predictors])$posterior
        expect_equal(unname(p), unname(native[, fit$codes, drop = FALSE]))
      } else if (method == "knn") {
        z <- sweep(sweep(as.matrix(d[rows, task$predictors]), 2, fit$center, "-"), 2, fit$scale, "/")
        for (i in seq_along(rows)) {
          distance <- rowSums((fit$engine$training - matrix(z[i, ], nrow(d), 3, byrow = TRUE))^2)
          neighbors <- order(distance, seq_along(distance))[1:5]
          expected <- as.numeric(table(factor(d$outcome[neighbors], levels = task$classes))) / 5
          expect_equal(unname(p[i, ]), expected)
        }
      } else if (method == "nnet") {
        z <- sweep(sweep(as.matrix(d[rows, task$predictors]), 2, fit$center, "-"), 2, fit$scale, "/")
        expect_equal(unname(p), unname(predict(fit$engine, z, type = "raw")))
        expect_equal(fit$center, vapply(d[task$predictors], mean, numeric(1)))
      }
    }
  }
  cv <- sb_compare_multiclass(task, c("svm_radial", "nnet", "knn", "ridge"),
    folds = 3, inner_folds = 3)
  expect_identical(cv$comparison$status, rep("success", 4))
})
