test_that("multiclass RDA uses klaR even when another package registers predict.rda", {
  skip_if_not_installed("klaR")
  d <- iris
  d$id <- seq_len(nrow(d))
  f <- sb_fit_multiclass(sb_multiclass_task(d, "id", "Species", names(iris)[1:4]), "rda", 3, 1)
  expected <- get("predict.rda", asNamespace("klaR"))(f$engine, d[1:5, 1:4])$posterior
  previous <- utils::getS3method("predict", "rda")
  withr::defer(registerS3method("predict", "rda", previous, envir = asNamespace("stats")))
  registerS3method("predict", "rda", function(...) stop("Unrelated RDA method"), envir = asNamespace("stats"))
  expect_error(stats::predict(f$engine, d[1:5, 1:4]), "Unrelated RDA")
  expect_equal(unname(sb_predict_multiclass(f, d[1:5, ])), unname(expected[, f$codes]))
})
