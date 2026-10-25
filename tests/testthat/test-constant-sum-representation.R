constant_sum_fixture <- function() {
  withr::with_seed(41, {
    x <- as.data.frame(matrix(stats::rnorm(320), 80, 4))
    names(x) <- paste0("g", 1:4)
    x[] <- sweep(as.matrix(x), 1L, rowMeans(x), "-")
    list(x = x, outcome = rep(0:1, 40))
  })
}

test_that("formal-run constant-sum coordinates preserve native predictions", {
  for (package in c("brglm2", "MASS", "mda", "hda", "klaR")) skip_if_not_installed(package)
  d <- constant_sum_fixture()
  for (id in c("brglm2_logistic", "lda_binary", "fda_binary", "mda_binary", "hda_binary", "loclda_binary")) {
    pipeline <- sb_model_pipelines("binary", id)[[1L]]
    adapter <- sb_resolve_learner_adapter(id, pipeline$learner_params)
    model <- adapter$fit(d$x, d$outcome)
    state <- adapter$extract_state(model)
    expect_s3_class(state, "stormbreaker_learner_state")
    expect_identical(state$features, names(d$x), info = id)
    expect_identical(state$projection$features, names(d$x), info = id)
    prediction <- adapter$predict(model, d$x)
    replay <- adapter$predict_state(state, d$x)$probability
    expect_equal(replay, prediction, tolerance = 1e-12, info = id)
    family <- Filter(function(f) id %in% f[[1L]], sb_uniform_learner_adapter_families())[[1L]]
    native <- family[[2L]](id, state$inner_params)
    native_x <- as.data.frame(sweep(as.matrix(d$x), 2L, state$projection$center, "-") %*% state$projection$basis)
    names(native_x) <- paste0("contrast_", seq_len(ncol(native_x)))
    reference <- native$fit(native_x, d$outcome)
    expect_equal(prediction, native$predict(reference, native_x), tolerance = 1e-12, info = id)
    old <- if (id == "brglm2_logistic") {
      list(variant = "brglm2_constant_sum_v1", projection = state$projection, inner = state$inner)
    } else {
      list(round2_representation = "constant_sum_orthonormal_v1", projection = state$projection,
        inner_params = state$inner_params, inner = state$inner)
    }
    expect_equal(adapter$predict_state(old, d$x)$probability, prediction, tolerance = 1e-12, info = id)
    shifted <- d$x
    shifted$g1 <- shifted$g1 + 0.1
    error <- tryCatch(adapter$predict_state(state, shifted), error = identity)
    expect_s3_class(error, "stormbreaker_constant_sum_mismatch")
  }
})

test_that("ordinary inputs retain native learner states", {
  skip_if_not_installed("brglm2")
  d <- constant_sum_fixture()
  d$x$g1 <- seq_len(nrow(d$x)) / 100
  p <- sb_model_pipelines("binary", "brglm2_logistic")[[1L]]
  adapter <- sb_resolve_learner_adapter(p$learner, p$learner_params)
  model <- adapter$fit(d$x, d$outcome)
  expect_s3_class(adapter$extract_state(model), "stormbreaker_omics_specialist_binary_state")
  expect_null(model$projection)
})

test_that("constant-sum linear parameters stay in original gene units", {
  skip_if_not_installed("brglm2")
  d <- constant_sum_fixture()
  for (id in c("brglm2_logistic", "lda_binary")) {
    pipeline <- sb_model_pipelines("binary", id)[[1L]]
    adapter <- sb_resolve_learner_adapter(id, pipeline$learner_params)
    model <- adapter$fit(d$x, d$outcome)
    state <- adapter$extract_state(model)
    parameters <- if (id == "lda_binary") sb_discriminant_parameters(state) else sb_extended_logistic_parameters(state)
    expect_identical(names(parameters$coefficient), names(d$x))
    expect_equal(stats::plogis(parameters$intercept + as.vector(as.matrix(d$x) %*% parameters$coefficient)),
      adapter$predict(model, d$x), tolerance = 1e-10)
    wrong_view <- tryCatch({
      if (id == "lda_binary") sb_extended_logistic_parameters(state) else sb_discriminant_parameters(state)
    }, error = identity)
    expect_s3_class(wrong_view, if (id == "lda_binary") "stormbreaker_unsupported_coefficients" else
      "stormbreaker_unsupported_diagnostic")
  }
})
