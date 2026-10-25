test_that("the menu includes all direct learners and three glmnet settings", {
  registry <- sb_pipeline_registry()
  for (type in c("binary", "survival")) {
    menu <- sb_model_catalog(type)
    expected <- unique(registry$learner[registry$task_type == type &
      registry$selector == "none" & registry$admission_status == "admitted"])
    expect_setequal(unique(menu$learner), expected)
    expect_setequal(menu$glmnet_alpha[!is.na(menu$glmnet_alpha)], c(0, .5, 1))
    pipelines <- sb_model_pipelines(type)
    expect_identical(names(pipelines), menu$model_id)
    expect_true(all(vapply(pipelines, function(p) p$selector == "none", logical(1))))
    expect_length(unique(vapply(pipelines, function(p) p$id, character(1))), nrow(menu))
  }
  expect_error(sb_model_catalog(NULL), "explicitly")
  expect_error(sb_model_pipelines("binary", "coxph"), "model IDs")
  expect_error(sb_model_pipelines("survival", character()), "model IDs")
  expect_error(sb_model_pipelines("binary", c("logistic", "logistic")), "distinct")
})

test_that("users choose task and genes and can request specific models", {
  genes <- c("gene_c", "gene_a")
  for (type in c("binary", "survival")) {
    d <- if (type == "binary") sb_example_binary(60, seed = 231) else
      sb_example_survival(60, seed = 231)
    task <- if (type == "binary") sb_binary_task(d, "patient_id", "outcome", genes) else
      sb_surv_task(d, "patient_id", "time", "event", genes, "month")
    models <- if (type == "binary") c("gaussian_naive_bayes", "logistic") else
      c("survreg_weibull", "coxph")
    result <- sb_compare_models(task, models = models, folds = 3L, seed = 91L, view = FALSE)
    expect_identical(result$input_genes, genes)
    expect_identical(result$requested_models, models)
    expect_identical(result$comparison$model_id, models)
    expect_true(all(result$comparison$complete_oof))
    expect_true(all(vapply(result$results, function(r)
      identical(r$provenance$fold_assignment, result$fold_assignment), logical(1))))
    expect_gte(result$elapsed_seconds, 0)
  }
})

test_that("all-model dispatch keeps missing dependencies and is independent of display", {
  task <- sb_binary_task(sb_example_binary(60, seed = 231), "patient_id", "outcome", "gene_a")
  menu <- sb_model_catalog("binary")
  menu$dependency_status <- "missing"
  dispatched <- list()
  local_mocked_bindings(
    sb_model_catalog = function(...) menu,
    sb_compare_cv = function(task, pipelines, ...) {
      dispatched[[length(dispatched) + 1L]] <<- pipelines
      list(comparison = data.frame(pipeline_id = vapply(pipelines, function(p) p$id, character(1))))
    },
    sb_attach_result_page = function(x, ...) x,
    sb_native_diagnostics = function(...) stop("Comparison must not extract diagnostics")
  )
  quiet <- sb_compare_models(task, models = NULL, view = FALSE)
  visible <- sb_compare_models(task, models = NULL, view = TRUE)
  expect_identical(quiet$requested_models, menu$model_id)
  expect_identical(visible$requested_models, menu$model_id)
  expect_identical(dispatched[[1L]], dispatched[[2L]])
  selected <- rev(menu$model_id[c(1L, 3L)])
  explicit <- sb_compare_models(task, models = selected, view = FALSE)
  expect_identical(explicit$requested_models, selected)
  expect_identical(names(dispatched[[3L]]), selected)
})
