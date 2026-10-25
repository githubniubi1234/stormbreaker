test_that("sharing respects outer and inner glmnet analysis folds", {
  skip_if_not_installed("glmnet")
  skip_if_not_installed("Boruta")
  skip_if_not_installed("ranger")
  calls <- list()
  local_mocked_bindings(sb_fit_runtime_selector_uncached=function(selector,x,endpoint) {
    calls[[length(calls)+1L]] <<- list(rows=rownames(x),endpoint=endpoint)
    sb_new_forest_selector_fit(selector$id,names(x),names(x)[1:2],
      sb_forest_selector_params(selector$id),
      data.frame(feature=names(x),selected=names(x)%in%names(x)[1:2]))
  })
  d <- sb_example_binary(n=120,seed=1208)
  task <- sb_binary_task(d,"patient_id","outcome",c("gene_a","gene_b","gene_c"))
  pipelines <- lapply(c(1,.5),function(alpha)sb_pipeline(
    "boruta_binary_roughfix","glmnet_logistic",glmnet_alpha=alpha,
    inner_folds=3L,inner_seed=101L))
  plain <- sb_compare_cv(task,pipelines,folds=3L,seed=1208L,reuse_selectors=FALSE)
  plain_calls <- calls; calls <- list()
  shared <- sb_compare_cv(task,pipelines,folds=3L,seed=1208L)
  expect_true(all(shared$comparison$status=="success"))
  expect_equal(shared,plain,tolerance=1e-12)
  expect_length(plain_calls,24L)
  expect_length(calls,12L)
  # No full-development selector fit is allowed during CV, and inner fits
  # must have smaller analysis sets than the outer fits they belong to.
  sizes <- vapply(calls,function(x)length(x$rows),integer(1))
  expect_equal(sum(sizes==80L),3L)
  expected <- c(
    lapply(shared$results[[1]]$folds,function(f)rownames(task$data)[f$analysis_rows]),
    unlist(lapply(shared$results[[1]]$tuning,function(t)
      lapply(t$inner_selections,function(s)s$analysis_patient_ids)),recursive=FALSE))
  expect_length(expected,12L)
  # Stratifying each class can yield unequal inner-fold sizes. Compare the
  # exact declared analysis identities instead of assuming balanced sizes.
  expect_true(all(vapply(expected,function(ids)
    sum(vapply(calls,function(call)identical(call$rows,ids),logical(1)))==1L,
    logical(1))))
  expect_true(all(sizes<nrow(task$data)))
  expect_null(.sb_selector_reuse_context$active)
})

test_that("seeded selector reuse leaves caller random state unchanged", {
  withr::local_seed(891L)
  local_mocked_bindings(sb_fit_runtime_selector_uncached=function(selector,x,endpoint) {
    withr::with_seed(1L,list(selected_features=names(x),value=runif(3)))
  })
  selector <- sb_resolve_selector_adapter("boruta_binary_roughfix",list())
  x <- data.frame(a=1:4); endpoint <- list(task_type="binary",outcome=c(0,1,0,1))
  before <- .Random.seed
  out <- sb_with_selector_reuse({
    first <- sb_fit_runtime_selector(selector,x,endpoint)
    expect_identical(.Random.seed,before)
    second <- sb_fit_runtime_selector(selector,x,endpoint)
    expect_identical(first,second)
    expect_identical(.Random.seed,before)
  },diagnostics=TRUE)
  expect_identical(out$reuse$hits,1L)
  expect_identical(.Random.seed,before)
})
