test_that("IECV doctor and stack final loops activate bounded model sharing", {
  skip_if_not_installed("survivalsvm")
  skip_if_not_installed("glmnet")
  d <- sb_example_survival(90L, 11032L)
  features <- c("gene_a", "gene_b", "gene_c")
  task <- sb_surv_task(d, "patient_id", "time", "event", features, "month")
  pp <- list(first = sb_pipeline("none", "survivalsvm"),
    second = sb_pipeline("fixed_features", "survivalsvm",
      selector_features = features, feature_set_id = "all"))
  native <- sb_fit_survivalsvm
  calls <- 0L
  final_counts <- integer()
  final <- sb_fit_final
  scope <- sb_with_selector_reuse
  reuse_enabled <- TRUE
  local_mocked_bindings(sb_fit_survivalsvm = function(...) {
    calls <<- calls + 1L; native(...)
  }, sb_fit_final = function(...) {
    before <- calls
    value <- final(...)
    final_counts <<- c(final_counts, calls - before)
    value
  }, sb_plots = function(...) list(),
  sb_with_selector_reuse = function(expr, enabled = TRUE, diagnostics = FALSE) {
    scope(expr, enabled = enabled && reuse_enabled, diagnostics = diagnostics)
  })
  for (entry in c("iecv", "doctor", "stack")) {
    run <- function() switch(entry,
      iecv = sb_fit_iecv_strategy(task, sb_validate_iecv_strategy(pp, "survival")),
      doctor = surv_run(task, pipelines = pp, folds = 3L, seed = 21L, plots = TRUE),
      stack = sb_fit_stack_final(task, sb_stack(pp, inner_folds = 3L, inner_seed = 21L)))
    reuse_enabled <- FALSE
    final_counts <- integer()
    plain <- run()
    expect_identical(final_counts, c(1L, 1L), info = entry)
    reuse_enabled <- TRUE
    final_counts <- integer()
    value <- run()
    expect_equal(value, plain, info = entry)
    fits <- switch(entry, iecv = value, doctor = lapply(value$details, `[[`, "fit"), stack = value$member_fits)
    expect_identical(vapply(fits, `[[`, character(1), "status"),
      stats::setNames(rep("success", 2L), names(fits)), info = entry)
    expect_identical(final_counts, c(1L, 0L), info = entry)
    expect_null(.sb_selector_reuse_context$active)
  }
})
