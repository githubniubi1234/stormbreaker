test_that("native information selectors share matching outer and inner analysis fits", {
  skip_if_not_installed("praznik")
  skip_if_not_installed("glmnet")
  native <- sb_fit_runtime_selector_uncached
  calls <- list()
  local_mocked_bindings(sb_fit_runtime_selector_uncached = function(selector, x, endpoint) {
    calls[[length(calls) + 1L]] <<- rownames(x)
    native(selector, x, endpoint)
  })
  withr::local_seed(1208L)
  d <- sb_example_binary(120L, 1208L)
  for (i in seq_len(21L)) d[[paste0("noise", i)]] <- rnorm(nrow(d))
  features <- c("gene_a", "gene_b", "gene_c", paste0("noise", seq_len(21L)))
  task <- sb_binary_task(d, "patient_id", "outcome", features)
  before <- .Random.seed
  for (selector in c("mrmr_binary_nlogn", sb_information_binary_selector_ids())) {
    pipelines <- lapply(c(1, .5), function(alpha) sb_pipeline(selector,
      "glmnet_logistic", glmnet_alpha = alpha, inner_folds = 3L, inner_seed = 101L))
    calls <- list()
    plain <- sb_compare_cv(task, pipelines, folds = 3L, seed = 1208L,
      reuse_selectors = FALSE)
    expect_length(calls, 24L)
    calls <- list()
    shared <- sb_compare_cv(task, pipelines, folds = 3L, seed = 1208L)
    expect_identical(shared, plain)
    expect_identical(unique(shared$comparison$status), "success")
    expect_length(calls, 12L)
    expected <- c(lapply(shared$results[[1]]$folds,
      function(f) rownames(task$data)[f$analysis_rows]),
      unlist(lapply(shared$results[[1]]$tuning, function(t)
        lapply(t$inner_selections, function(s) s$analysis_patient_ids)), recursive = FALSE))
    expect_identical(vapply(expected, function(ids) sum(vapply(calls,
      function(rows) identical(rows, ids), logical(1))), integer(1)), rep(1L, 12L))
    expect_identical(.Random.seed, before)
  }
})
