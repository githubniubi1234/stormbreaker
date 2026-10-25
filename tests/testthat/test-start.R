test_that("start page contains installed recipes and never runs them", {
  local_mocked_bindings(sb_fit_final = function(...) stop("Must not fit"),
    sb_compare_models = function(...) stop("Must not compare"),
    sb_analyze_methods = function(...) stop("Must not analyze"))
  file <- sb_start(view = FALSE)
  html <- paste(readLines(file, encoding = "UTF-8"), collapse = "\n")
  expect_identical(file.exists(file), TRUE)
  expect_identical(grepl("@@", html, fixed = TRUE), FALSE)
  for (action in c("compare", "validate")) {
    expect_match(html, paste0('data-action="', action, '"'), fixed = TRUE)
    expect_match(html, paste0('id="', action, '-code"'), fixed = TRUE)
  }
  expect_match(html, "eval_external(", fixed = TRUE)
  expect_identical(grepl('src="https://', html, fixed = TRUE), FALSE)
})

test_that("both doctor examples execute through real interfaces", {
  for (p in c("glmnet", "randomForest", "ggplot2")) skip_if_not_installed(p)
  withr::local_pdf(tempfile(fileext = ".pdf"))
  results <- lapply(c("compare", "validate"), function(action) {
    e <- new.env(parent = globalenv())
    sys.source(system.file("examples", paste0("doctor-", action, ".R"),
      package = "stormbreaker", mustWork = TRUE), envir = e)
    e
  })
  comparison <- results[[1L]]$comparison
  expect_length(results[[1L]]$run$pipelines, 3L)
  expect_identical(comparison$comparison$status, rep("success", 3))
  expect_identical(results[[1L]]$views$status$view, c("ROC", "Calibration", "Decision curve"))
  expect_identical(results[[2L]]$fit$status, "success")
  expect_identical(results[[2L]]$prediction$status, "success")
  expect_identical(results[[2L]]$evaluation$status, "success")
})
