test_that("specialized binary linear specifications are explicit", {
  expect_identical(
    sb_specialized_binary_linear_learner_ids(),
    c(
      "stepglm_binary",
      "glmboost_binary",
      "plsrglm_binary"
    )
  )
  expect_identical(
    sb_specialized_binary_linear_params(
      "stepglm_binary",
      stepglm_direction = "forward"
    ),
    list(
      direction = "forward",
      criterion = "AIC",
      k = 2,
      max_steps = 1000L
    )
  )
  expect_identical(
    sb_specialized_binary_linear_params("glmboost_binary"),
    list(
      iterations = 100L,
      shrinkage = 0.1,
      center = FALSE,
      family = "binomial_logit"
    )
  )
  expect_identical(
    sb_specialized_binary_linear_params("plsrglm_binary"),
    list(
      components = 2L,
      scale_x = TRUE,
      model = "pls-glm-logistic"
    )
  )
})

test_that("Stepglm direction is part of the pipeline specification", {
  default <- sb_pipeline(learner = "stepglm_binary")
  backward <- sb_pipeline(
    learner = "stepglm_binary",
    stepglm_direction = "backward"
  )
  bad <- capture_sb_error(sb_pipeline(
    learner = "stepglm_binary",
    stepglm_direction = "sideways"
  ))
  irrelevant <- capture_sb_error(sb_pipeline(
    learner = "logistic",
    stepglm_direction = "forward"
  ))

  expect_identical(default$learner_params$direction, "forward")
  expect_identical(
    default$id,
    "none__stepglm_binary__directionforward"
  )
  expect_identical(
    backward$id,
    "none__stepglm_binary__directionbackward"
  )
  expect_s3_class(bad, "stormbreaker_bad_stepglm_direction")
  expect_s3_class(
    irrelevant,
    "stormbreaker_irrelevant_stepglm_parameter"
  )
})

test_that("Stepglm normalizes equivalent character direction objects", {
  data <- sb_example_binary(n = 72L, seed = 8001L)
  x <- data[c("gene_a", "gene_b", "gene_c")]
  directions <- list(
    "forward",
    c(requested = "forward"),
    structure("forward", class = "declared_direction")
  )
  fits <- lapply(directions, function(direction) {
    params <- sb_specialized_binary_linear_params(
      "stepglm_binary",
      stepglm_direction = direction
    )
    adapter <- sb_specialized_binary_linear_learner_adapter(
      "stepglm_binary",
      params
    )
    model <- adapter$fit(x, data$outcome)
    list(
      params = params,
      model = model,
      probability = adapter$predict(model, x)
    )
  })

  expect_identical(
    lapply(fits, function(fit) fit$params$direction),
    rep(list("forward"), 3L)
  )
  expect_identical(
    lapply(fits, function(fit) fit$model$effective$selected_features),
    rep(list(fits[[1L]]$model$effective$selected_features), 3L)
  )
  expect_identical(
    lapply(fits, function(fit) fit$probability),
    rep(list(fits[[1L]]$probability), 3L)
  )
})

test_that("specialized binary linear adapters preserve fitted predictions", {
  skip_if_not_installed("mboost", minimum_version = "2.9.0")
  skip_if_not_installed("plsRglm", minimum_version = "1.6.0")

  data <- sb_example_binary(n = 72L, seed = 8001L)
  x <- data[c("gene_a", "gene_b", "gene_c")]
  new_x <- x[c(4L, 29L, 61L), , drop = FALSE]

  for (direction in c("forward", "backward", "both")) {
    params <- sb_specialized_binary_linear_params(
      "stepglm_binary",
      stepglm_direction = direction
    )
    adapter <- sb_specialized_binary_linear_learner_adapter(
      "stepglm_binary",
      params
    )
    model <- adapter$fit(x, data$outcome)
    probability <- adapter$predict(model, new_x)
    replay <- adapter$predict_state(adapter$extract_state(model), new_x)

    expect_identical(model$params$direction, direction)
    expect_equal(replay$probability, probability, tolerance = 0)
  }

  glmboost_params <- sb_specialized_binary_linear_params(
    "glmboost_binary"
  )
  glmboost_adapter <- sb_specialized_binary_linear_learner_adapter(
    "glmboost_binary",
    glmboost_params
  )
  glmboost_model <- glmboost_adapter$fit(x, data$outcome)
  glmboost_native <- mboost::glmboost(
    x = as.matrix(x),
    y = factor(data$outcome, levels = c(0L, 1L)),
    family = mboost::Binomial(type = "glm", link = "logit"),
    control = mboost::boost_control(
      mstop = glmboost_params$iterations,
      nu = glmboost_params$shrinkage,
      risk = "inbag",
      stopintern = FALSE,
      trace = FALSE
    ),
    center = glmboost_params$center
  )
  expect_equal(
    glmboost_adapter$predict(glmboost_model, new_x),
    as.numeric(stats::predict(
      glmboost_native,
      newdata = as.matrix(new_x),
      type = "response"
    )),
    tolerance = 1e-14
  )

  plsrglm_params <- sb_specialized_binary_linear_params(
    "plsrglm_binary"
  )
  plsrglm_adapter <- sb_specialized_binary_linear_learner_adapter(
    "plsrglm_binary",
    plsrglm_params
  )
  plsrglm_model <- plsrglm_adapter$fit(x, data$outcome)
  plsrglm_native <- plsRglm::PLS_glm(
    dataY = data$outcome,
    dataX = x,
    nt = plsrglm_params$components,
    modele = plsrglm_params$model,
    typeVC = "none",
    scaleX = plsrglm_params$scale_x,
    pvals.expli = FALSE,
    MClassed = FALSE,
    sparse = FALSE,
    naive = FALSE,
    verbose = FALSE
  )
  native_coefficients <- drop(plsrglm_native$Coeffs)
  native_probability <- stats::plogis(
    unname(native_coefficients[["Intercept"]]) +
      as.vector(
        as.matrix(new_x) %*% native_coefficients[names(new_x)]
      )
  )
  expect_equal(
    plsrglm_adapter$predict(plsrglm_model, new_x),
    native_probability,
    tolerance = 1e-14
  )

  for (model_adapter in list(
    list(model = glmboost_model, adapter = glmboost_adapter),
    list(model = plsrglm_model, adapter = plsrglm_adapter)
  )) {
    state <- model_adapter$adapter$extract_state(model_adapter$model)
    replay <- model_adapter$adapter$predict_state(
      unserialize(serialize(state, NULL)),
      new_x
    )
    direct <- model_adapter$adapter$predict(model_adapter$model, new_x)

    expect_identical("engine_raw" %in% names(state), FALSE)
    expect_equal(replay$probability, direct, tolerance = 0)
  }
})

test_that("specialized binary linear learners run CV and final prediction", {
  skip_if_not_installed("mboost", minimum_version = "2.9.0")
  skip_if_not_installed("plsRglm", minimum_version = "1.6.0")

  data <- sb_example_binary(n = 72L, seed = 8002L)
  task <- sb_binary_task(
    data,
    "patient_id",
    "outcome",
    predictors = c("gene_a", "gene_b", "gene_c")
  )
  learners <- sb_specialized_binary_linear_learner_ids()

  for (learner in learners) {
    pipeline <- if (identical(learner, "stepglm_binary")) {
      sb_pipeline(learner = learner, stepglm_direction = "both")
    } else {
      sb_pipeline(learner = learner)
    }
    result <- sb_run_cv(task, pipeline, folds = 3L, seed = 8003L)
    fit <- sb_fit_final(task, pipeline)
    prediction <- sb_predict_final(
      fit,
      data[1:6, c("patient_id", "gene_c", "gene_a", "gene_b")],
      "patient_id"
    )

    expect_identical(result$status, "success")
    expect_identical(result$complete_oof, TRUE)
    expect_identical(fit$status, "success")
    expect_identical(prediction$status, "success")
  }
})

test_that("specialized binary linear learners use legal feature selectors", {
  skip_if_not_installed("mboost", minimum_version = "2.9.0")

  data <- sb_example_binary(n = 72L, seed = 8004L)
  task <- sb_binary_task(
    data,
    "patient_id",
    "outcome",
    predictors = c("gene_a", "gene_b", "gene_c")
  )
  pipelines <- list(
    sb_pipeline(
      selector = "univ_binary_auc_top_k",
      learner = "stepglm_binary",
      selector_n_features = 2L,
      stepglm_direction = "forward"
    ),
    sb_pipeline(
      selector = "fixed_features",
      learner = "glmboost_binary",
      selector_features = c("gene_b", "gene_a"),
      feature_set_id = "two_genes"
    )
  )

  for (pipeline in pipelines) {
    result <- sb_run_cv(task, pipeline, folds = 3L, seed = 8005L)
    fit <- sb_fit_final(task, pipeline)

    expect_identical(result$status, "success")
    expect_identical(fit$status, "success")
    expect_length(fit$required_features, 2L)
  }
})
