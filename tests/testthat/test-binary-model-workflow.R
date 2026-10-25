test_that("one binary menu contains every existing prediction and feature method", {
  menu <- sb_models("binary")
  expect_equal(sum(menu$can_predict), nrow(sb_model_catalog("binary")))
  expect_equal(sum(menu$can_extract_features), length(sb_doctor_selector_menu("binary")))
  expect_identical(any(menu$model_id %in% sb_retired_learner_ids()), FALSE)
  expect_setequal(menu$model_id[menu$can_predict], sb_model_catalog("binary")$model_id)
  expect_identical(anyDuplicated(tolower(menu$name)), 0L)
  expect_identical(anyDuplicated(menu$model_id), 0L)
  expect_identical(menu$can_predict[menu$name == "Boruta"], FALSE)
  expect_identical(menu$can_extract_features[menu$name == "RF"], TRUE)
  expect_identical(menu$can_extract_features[menu$name == "Ridge"], FALSE)
  expect_identical(menu$r_packages[menu$name == "Boruta"], "Boruta, ranger")
})

test_that("ordered model names preserve native stage definitions and their order", {
  forward <- bi_pipeline(c("Lasso", "RF", "GBM"))
  old <- bi_pipeline("GBM", bi_selector(c("Lasso", "RF")))
  backward <- bi_pipeline(c("RF", "Lasso", "GBM"))
  expect_identical(forward, old)
  expect_identical(forward$id == backward$id, FALSE)
  expect_identical(bi_pipeline(c("PLSR", "Logistic"))$selector, "plsrglm_binary_weights")
  expect_identical(bi_pipeline(c("Enet", "Logistic"))$selector, "enet_binary_nonzero")
})

test_that("methods cannot silently acquire unsupported stage capabilities", {
  terminal <- tryCatch(bi_pipeline(c("Logistic", "Boruta")), error = identity)
  intermediate <- tryCatch(bi_pipeline(c("NaiveBayes", "Logistic")), error = identity)
  unknown <- tryCatch(bi_pipeline("not_a_method"), error = identity)
  expect_s3_class(terminal, "stormbreaker_method_cannot_predict")
  expect_s3_class(intermediate, "stormbreaker_method_cannot_extract")
  expect_s3_class(unknown, "stormbreaker_bad_model_choice")
})

test_that("three-stage enumeration covers the chosen method pool rather than presets", {
  plan <- sb_doctor_pipeline_set("binary", c("Lasso", "RF", "GBM"), c(1, 2, 3))
  expect_length(plan$pipelines, 39L)
  expect_length(plan$unavailable, 0L)
  expect_contains(names(plan$pipelines), "RF + Lasso + GBM")
  expect_contains(names(plan$pipelines), "Lasso + RF + GBM")
  expect_contains(names(plan$pipelines), "GBM + RF + Lasso")
  expect_identical(plan$pipelines[["RF + Lasso + GBM"]], bi_pipeline(c("RF", "Lasso", "GBM")))
})

test_that("feature-only methods are available in the same pool but not as predictors", {
  plan <- sb_doctor_pipeline_set("binary", c("Boruta", "Logistic"), c(1, 2))
  expect_named(plan$pipelines, c("Logistic", "Boruta + Logistic"))
  expect_length(plan$unavailable, 0L)
  expect_identical(plan$pipelines[[2L]]$selector, "boruta_binary_roughfix")
  expect_identical(plan$pipelines[[2L]]$learner, "logistic")
  defaults <- eval(formals(bi_run)$stages, list(models = NULL))
  expect_equal(defaults, c(1L, 2L))
})

test_that("unified ordered routes execute CV and report their saved genes", {
  data <- sb_example_binary(120, seed = 37)
  task <- sb_binary_task(data, "patient_id", "outcome", c("gene_a", "gene_b", "gene_c"))
  route <- bi_pipeline(c("Stepglm", "Logistic"))
  run <- bi_run(task, pipelines = list(route = route), folds = 3, seed = 7, plots = FALSE)
  completed <- sb_complete_run(run, "route", plots = FALSE)
  expect_identical(run$comparison$comparison$status, "success")
  expect_identical(completed$fit$status, "success")
  expect_setequal(completed$genes$gene, completed$fit$selection$state$selected_features)
  expect_identical(task$data, data)
})

test_that("the home page has one modeling entry and no separate selection task", {
  file <- sb_start(view = FALSE)
  html <- paste(readLines(file, encoding = "UTF-8"), collapse = "\n")
  expect_identical(grepl('data-action="select"', html, fixed = TRUE), FALSE)
  expect_identical(grepl('id="select"', html, fixed = TRUE), FALSE)
  expect_identical(grepl("@@", html, fixed = TRUE), FALSE)
  expect_match(html, "bi_run", fixed = TRUE)
  expect_match(html, "bi_pipeline", fixed = TRUE)
})

test_that("invalid pools retain explicit failures without inventing a model", {
  alone <- sb_doctor_pipeline_set("binary", "Boruta", 1)
  expect_length(alone$pipelines, 0L)
  expect_match(alone$unavailable$Boruta, "final prediction model", fixed = TRUE)
  repeated <- tryCatch(sb_doctor_pipeline_set("binary", c("RF", "random_forest_binary_fixed"), 1), error = identity)
  expect_s3_class(repeated, "stormbreaker_bad_model_choice")
  stages <- tryCatch(sb_doctor_pipeline_set("binary", "RF", 4), error = identity)
  expect_s3_class(stages, "error")
})

test_that("feature aliases and explicit native IDs resolve to the same route", {
  expect_identical(bi_pipeline(c("ElasticNet", "Logistic")), bi_pipeline(c("Enet", "Logistic")))
  expect_identical(bi_pipeline(c("boruta_binary_roughfix", "logistic")), bi_pipeline(c("Boruta", "Logistic")))
  expect_identical(bi_pipeline("random_forest_binary_fixed"), bi_pipeline("RF"))
})

test_that("previous native admission failures remain visible in the unified pool", {
  expected <- tryCatch(bi_pipeline(c("Boruta", "cancerclass")), error = identity)
  plan <- sb_doctor_pipeline_set("binary", c("Boruta", "cancerclass"), 2)
  expect_s3_class(expected, "stormbreaker_unadmitted_pipeline")
  expect_identical(plan$unavailable[["Boruta + cancerclass"]], conditionMessage(expected))
  expect_length(plan$pipelines, 0L)
})
