test_that("explicit mixed-stage lists use one shared comparison", {
  for (type in c("binary", "survival")) {
    binary <- type == "binary"
    d <- if (binary) sb_example_binary(120, 37) else sb_example_survival(120, 37)
    build <- if (binary) bi_pipeline else surv_pipeline
    select <- if (binary) bi_selector else surv_selector
    method <- if (binary) "Stepglm" else "StepCox"
    model <- if (binary) "Logistic" else "Cox"
    pipelines <- list(single = build(model), pair = build(model, select(method)),
      chain = build(model, select(c(method, method))))
    task <- if (binary) sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b", "gene_c")) else
      sb_surv_task(d, "patient_id", "time", "event", c("gene_a", "gene_b", "gene_c"), "month")
    run <- if (binary) bi_run else surv_run
    result <- run(task, pipelines = pipelines, folds = 3, seed = 19, plots = FALSE)
    expected <- sb_compare_cv(task, unname(pipelines), folds = 3, seed = 19, progress = FALSE)
    expect_identical(result$pipelines, pipelines)
    expect_equal(result$comparison, expected)
    expect_length(result$comparison$results, 3L)
    expect_equal(result$comparison$comparison$status, rep("success", 3))
    completed <- sb_complete_run(result, "pair", plots = FALSE)
    expect_identical(completed$fit$status, "success")
  }
})

test_that("custom lists reject ambiguous names and incompatible task types", {
  d <- sb_example_binary(80, 11)
  task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b", "gene_c"))
  p <- bi_pipeline("Logistic")
  for (pipelines in list(list(p), setNames(list(p, p), c("same", "same")),
    list(first = p, duplicate = p), list(wrong = surv_pipeline("Cox")))) {
    problem <- tryCatch(bi_run(task, pipelines = pipelines, plots = FALSE), error = identity)
    expect_s3_class(problem, "error")
  }
  problem <- tryCatch(bi_run(task, models = "Logistic", pipelines = list(one = p)), error = identity)
  expect_match(conditionMessage(problem), "not both")
})

test_that("single external risk groups are reported and exported without plotting", {
  d <- sb_example_survival(120, 37)
  run <- surv_run(d, "time", "event", c("gene_a", "gene_b", "gene_c"),
    "patient_id", models = "Cox", folds = 3, plots = FALSE)
  run$details$Cox <- list(fit = sb_fit_final(run$development_task, run$pipelines$Cox))
  e <- d
  e$patient_id <- paste0("external", seq_len(nrow(e)))
  one <- e
  for (gene in c("gene_a", "gene_b", "gene_c")) one[[gene]] <- e[[gene]][1]
  local_mocked_bindings(sb_fit_final = function(...) stop("Unexpected refit"))
  for (plots in c(FALSE, TRUE)) {
    bundle <- sb_complete_run(run, "Cox", list(one = one, two = e), plots = plots)
    expect_identical(bundle$evaluations$one$status, "success")
    expect_identical(bundle$risk_groups$one$status, "single_group")
    expect_equal(bundle$risk_groups$one$n_low + bundle$risk_groups$one$n_high, nrow(e))
    expect_identical(bundle$risk_groups$two$status, "two_groups")
    expect_null(bundle$issues$risk$two)
    expect_match(bundle$issues$risk$one, "single_group")
    expect_equal(bundle$risk_groups$one$cutoff, bundle$risk_groups$two$cutoff)
    path <- tempfile("risk-export-")
    sb_export_run(bundle, path)
    expect_match(paste(read.csv(file.path(path, "issues.csv"))$message, collapse = " "), "single_group")
    expect_equal(read.csv(file.path(path, "risk-groups.csv"))$status, c("single_group", "two_groups"))
  }
})
