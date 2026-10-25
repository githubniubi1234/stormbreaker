test_that("explicit Boruta budgets retain defaults and survive contracts", {
  for (id in sb_boruta_selector_ids()) {
    learner <- if (grepl("survival", id)) "coxph" else "logistic"
    old <- sb_pipeline(id, learner)
    expect_identical(old$id, paste(id, learner, sep = "__"))
    p <- sb_pipeline(id, learner, selector_options = list(num_trees = 50, maximum_runs = 20))
    expect_identical(p$selector_params$num_trees, 50L)
    expect_identical(p$selector_params$maximum_runs, 20L)
    expect_identical(sb_validate_pipeline_contract(p, p$task_type), p)
    expect_false(identical(p$id, old$id))
  }
  expect_error(sb_pipeline("boruta_survival_roughfix", "coxph",
    selector_options = list(num_threads = 8)), class = "stormbreaker_bad_selector_options")
})
test_that("Boruta budget reaches the native call", {
  skip_if_not_installed("Boruta")
  seen <- NULL
  local_mocked_bindings(sb_call_boruta = function(x, response, params) {
    seen <<- params
    stop("observed native request")
  })
  d <- sb_example_survival(60, 11)
  p <- sb_pipeline("boruta_survival_roughfix", "coxph",
    selector_options = list(num_trees = 50L, maximum_runs = 20L))
  a <- sb_resolve_selector_adapter(p$selector, p$selector_params)
  expect_error(a$fit(d[c("gene_a", "gene_b")], d$time, d$event), "observed native request")
  expect_identical(seen$num_trees, 50L)
  expect_identical(seen$maximum_runs, 20L)
  expect_identical(seen$num_threads, 1L)
  captured <- sb_captured_selector_adapter(p$selector, p$selector_params)
  expect_error(captured$fit(d[c("gene_a", "gene_b")], d$time, d$event), "observed native request")
  expect_identical(seen$num_trees, 50L)
  default_captured <- sb_captured_selector_adapter(p$selector, NULL)
  expect_error(default_captured$fit(d[c("gene_a", "gene_b")], d$time, d$event), "observed native request")
  expect_identical(seen$num_trees, 500L)
})

test_that("retired PST selectors are rejected rather than substituted", {
  for (task in c("binary", "survival")) {
    id <- if (task == "binary") "vsolassobag_binary_pst" else "vsolassobag_cox_pst"
    learner <- if (task == "binary") "logistic" else "coxph"
    expect_false(id %in% sb_pipeline_registry(task_type = task)$selector)
    expect_false(id %in% unname(sb_doctor_selector_menu(task)))
    expect_error(sb_pipeline(id, learner), class = "stormbreaker_unsupported_selector")
    expect_error(sb_doctor_selector("VSOLassoBag_PST", task), "Unknown selector")
  }
})

