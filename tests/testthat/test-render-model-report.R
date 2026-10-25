test_that("offline report reuses SHAP with accessible model and sample controls", {
  skip_if_not_installed("ggplot2")
  skip_if_not(capabilities("cairo"))
  d <- make_sb_binary_fixture(n = 50L, seed = 108L)
  task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b"))
  fit <- sb_fit_final(task, sb_pipeline(learner = "logistic"))
  new <- d[1:2, c("patient_id", "gene_a", "gene_b")]
  new$patient_id <- c("private_identifier_one", "private_identifier_two")
  file <- tempfile(fileext = ".html")
  path <- sb_render_model_report(fit, new, "patient_id", d, file,
    title = "Model <script>untrusted</script>"
  )
  html <- paste(readLines(path, warn = FALSE), collapse = "\n")
  expect_match(html, "Model &lt;script&gt;untrusted&lt;/script&gt;", fixed = TRUE)
  expect_match(html, "id='coefficient-view'", fixed = TRUE)
  expect_match(html, "id='shap-view'", fixed = TRUE)
  expect_match(html, "id='sample-2'", fixed = TRUE)
  expect_match(html, "replaceChildren", fixed = TRUE)
  expect_match(html, "Training reference: 50", fixed = TRUE)
  expect_match(html, "<table>", fixed = TRUE)
  expect_identical(grepl("private_identifier", html, fixed = TRUE), FALSE)
  expect_identical(grepl("src='https", html, fixed = TRUE), FALSE)
  expect_error(sb_render_model_report(fit, new, "patient_id", d, file),
    class = "stormbreaker_report_exists"
  )
})

test_that("survival report binds the same sample to waterfall and survival", {
  skip_if_not_installed("ggplot2")
  skip_if_not(capabilities("cairo"))
  d <- make_sb_fixture(n = 60L, seed = 107L)
  task <- sb_surv_task(d, "patient_id", "time", "event", c("gene_a", "gene_b"),
    "month",
    prediction_horizon = 6
  )
  fit <- sb_fit_final(task, sb_pipeline())
  file <- tempfile(fileext = ".html")
  path <- sb_render_model_report(fit, d[1:2, ], "patient_id", d, file, development = task)
  html <- paste(readLines(path, warn = FALSE), collapse = "\n")
  expect_match(html, "Predicted survival for observation 2", fixed = TRUE)
  expect_match(html, "SHAP waterfall for observation 2", fixed = TRUE)
  expect_match(html, "Survival at 6 month", fixed = TRUE)
  expect_match(html, "Baseline survival alone", fixed = TRUE)
  expect_error(sb_render_model_report(fit, d, "patient_id", d, tempfile(fileext = ".html")),
    class = "stormbreaker_missing_survival_curve"
  )
})
