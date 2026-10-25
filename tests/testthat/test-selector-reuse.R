test_that("reuse matches complete inputs, parameters and RNG kind", {
  calls <- 0L
  local_mocked_bindings(sb_fit_runtime_selector_uncached = function(selector,x,endpoint) {
    calls <<- calls + 1L
    list(selected_features=names(x), value=calls)
  })
  selector <- sb_resolve_selector_adapter("boruta_binary_roughfix",list())
  x <- data.frame(a=1:6,b=6:1,row.names=paste0("p",1:6))
  endpoint <- list(task_type="binary",outcome=rep(0:1,3))
  original_kind <- RNGkind()
  withr::defer(do.call(RNGkind,as.list(original_kind)))
  result <- sb_with_selector_reuse({
    first <- sb_fit_runtime_selector(selector,x,endpoint)
    expect_identical(sb_fit_runtime_selector(selector,x,endpoint),first)
    changed <- x; changed$a[1] <- 99L
    sb_fit_runtime_selector(selector,changed,endpoint)
    sb_fit_runtime_selector(selector,x[6:1,],endpoint)
    sb_fit_runtime_selector(selector,x[c("b","a")],endpoint)
    changed <- x; rownames(changed) <- paste0("q",1:6)
    sb_fit_runtime_selector(selector,changed,endpoint)
    changed_endpoint <- endpoint; changed_endpoint$outcome[1] <- 1L
    sb_fit_runtime_selector(selector,x,changed_endpoint)
    changed_selector <- selector; changed_selector$params <- list(seed=2L)
    sb_fit_runtime_selector(changed_selector,x,endpoint)
    RNGkind(normal.kind=if(original_kind[2]=="Inversion")"Box-Muller" else "Inversion")
    sb_fit_runtime_selector(selector,x,endpoint)
    do.call(RNGkind,as.list(original_kind))
    first
  },diagnostics=TRUE)
  expect_identical(calls,8L)
  expect_identical(result$reuse$hits,1L)
  expect_null(.sb_selector_reuse_context$active)
  sb_fit_runtime_selector(selector,x,endpoint)
  expect_identical(calls,9L)
})

test_that("warnings, errors and unoptimized selectors are executed again", {
  calls <- 0L
  local_mocked_bindings(sb_fit_runtime_selector_uncached=function(selector,x,endpoint) {
    calls <<- calls + 1L
    if(selector$params$mode=="warning")warning("same warning")
    if(selector$params$mode=="error")stop("same error")
    list(selected_features=names(x))
  })
  selector <- sb_resolve_selector_adapter("boruta_binary_roughfix",list())
  x <- data.frame(a=1:4); endpoint <- list(task_type="binary",outcome=c(0,1,0,1))
  sb_with_selector_reuse({
    selector$params <- list(mode="warning")
    for(i in 1:2)expect_warning(sb_fit_runtime_selector(selector,x,endpoint),"same warning")
    selector$params <- list(mode="error")
    for(i in 1:2)expect_error(sb_fit_runtime_selector(selector,x,endpoint),"same error")
    selector$id <- "none"; selector$params <- list(mode="none")
    for(i in 1:2)sb_fit_runtime_selector(selector,x,endpoint)
  })
  expect_identical(calls,6L)
  expect_error(sb_with_selector_reuse(stop("cleanup")),"cleanup")
  expect_null(.sb_selector_reuse_context$active)
})

test_that("disabled and nested contexts restore the outer context", {
  calls <- 0L
  local_mocked_bindings(sb_fit_runtime_selector_uncached=function(selector,x,endpoint) {
    calls <<- calls + 1L; list(selected_features=names(x))
  })
  selector <- sb_resolve_selector_adapter("boruta_binary_roughfix",list())
  x <- data.frame(a=1:4); endpoint <- list(task_type="binary",outcome=c(0,1,0,1))
  out <- sb_with_selector_reuse({
    sb_fit_runtime_selector(selector,x,endpoint)
    sb_with_selector_reuse({
      sb_fit_runtime_selector(selector,x,endpoint)
      sb_fit_runtime_selector(selector,x,endpoint)
    },enabled=FALSE)
    sb_with_selector_reuse({
      sb_fit_runtime_selector(selector,x,endpoint)
      sb_fit_runtime_selector(selector,x,endpoint)
    })
    sb_fit_runtime_selector(selector,x,endpoint)
  },diagnostics=TRUE)
  expect_identical(calls,4L)
  expect_identical(out$reuse$hits,1L)
  expect_null(.sb_selector_reuse_context$active)
  expect_error(sb_with_selector_reuse(1,enabled=NA),class="stormbreaker_bad_selector_reuse")
})
