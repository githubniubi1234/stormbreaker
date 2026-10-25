test_that("additional fixed binary learners keep distinct native methods", {
  learners <- c(
    "kernel_naive_bayes_binary",
    "rda_binary",
    "nnet_binary",
    "svm_polynomial_binary",
    "svm_sigmoid_binary",
    "bart_binary",
    "c50_binary",
    "kernel_weighted_knn_binary",
    "gaussian_process_binary",
    "pamr_binary",
    "neuralnet_multilayer_binary"
  )
  expect_setequal(
    intersect(sb_classical_binary_learner_ids(), learners),
    learners
  )
  expect_identical(
    sb_classical_binary_params(
      "kernel_naive_bayes_binary"
    ),
    list(
      prior = "empirical",
      usekernel = TRUE,
      fL = 0,
      adjust = 1
    )
  )
  expect_identical(
    sb_classical_binary_params("rda_binary"),
    list(
      prior = "empirical",
      gamma = 0.5,
      lambda = 0.5,
      estimate_error = FALSE
    )
  )
  expect_identical(
    sb_classical_binary_params("nnet_binary"),
    list(
      size = 5L,
      decay = 0.01,
      maxit = 2000L,
      rang = 0.1,
      entropy = TRUE,
      standardize = TRUE,
      max_weights = 1000000L,
      seed = 1L
    )
  )
  expect_identical(
    sb_classical_binary_params("svm_polynomial_binary"),
    list(
      kernel = "polynomial",
      cost = 1,
      gamma = "one_over_p",
      degree = 3L,
      coef0 = 0,
      standardize = TRUE,
      probability = TRUE,
      seed = 1L
    )
  )
  expect_identical(
    sb_classical_binary_params("svm_sigmoid_binary"),
    list(
      kernel = "sigmoid",
      cost = 1,
      gamma = "one_over_p",
      coef0 = 0,
      standardize = TRUE,
      probability = TRUE,
      seed = 1L
    )
  )
  expect_identical(
    sb_classical_binary_params("bart_binary"),
    list(
      num_trees = 50L,
      num_cutpoints = 100L,
      posterior_draws = 250L,
      burn_in = 100L,
      keep_every = 1L,
      prior_k = 2,
      tree_power = 2,
      tree_base = 0.95,
      sparse = FALSE,
      seed = 1L
    )
  )
  expect_identical(
    sb_classical_binary_params("c50_binary"),
    list(
      trials = 10L,
      rules = FALSE,
      winnow = FALSE,
      early_stopping = FALSE,
      seed = 1L
    )
  )
  expect_identical(
    sb_classical_binary_params("kernel_weighted_knn_binary"),
    list(k = 7L, distance = 2, kernel = "optimal", scale = TRUE)
  )
  expect_identical(
    sb_classical_binary_params("gaussian_process_binary"),
    list(
      kernel = "rbfdot",
      kernel_parameters = "automatic",
      scaled = TRUE,
      tolerance = 0.0005,
      cross_validation = 0L,
      keep_fitted = FALSE,
      seed = 1L
    )
  )
  expect_identical(
    sb_classical_binary_params("pamr_binary"),
    list(
      threshold = 1,
      scale_sd = TRUE,
      offset_percent = 50,
      remove_zero_thresholds = FALSE,
      prior = "empirical"
    )
  )
  expect_identical(
    sb_classical_binary_params("neuralnet_multilayer_binary"),
    list(
      hidden = c(3L, 2L),
      threshold = 0.01,
      stepmax = 500000L,
      algorithm = "rprop+",
      error_function = "sse",
      activation = "logistic",
      linear_output = FALSE,
      standardize = TRUE,
      seed = 1L
    )
  )
})

test_that("additional fixed binary probabilities match native engines", {
  data <- sb_example_binary(n = 72L, seed = 6810L)
  x <- data[c("gene_a", "gene_b", "gene_c")]
  assessment <- x[1:8, , drop = FALSE]
  learners <- c(
    "kernel_naive_bayes_binary",
    "rda_binary",
    "nnet_binary",
    "svm_polynomial_binary",
    "svm_sigmoid_binary"
  )

  for (learner in learners) {
    dependency <- sb_classical_binary_dependency(learner)
    if (!dependency$available) {
      next
    }
    params <- sb_classical_binary_params(learner)
    if (identical(learner, "cforest_binary")) {
      params$num_trees <- 20L
    }
    if (identical(learner, "bagging_binary")) {
      params$num_trees <- 20L
    }
    if (learner %in% c("gamboost_binary", "blackboost_binary")) {
      params$mstop <- 20L
    }
    if (identical(learner, "aorsf_binary")) {
      params$num_trees <- 20L
    }
    adapter <- sb_classical_binary_learner_adapter(learner, params)
    model <- adapter$fit(x, data$outcome)
    observed <- adapter$predict(model, assessment)
    native_engine <- if (identical(learner, "nnet_binary")) {
      model$engine$model
    } else {
      model$engine
    }
    expected_class <- switch(learner,
      kernel_naive_bayes_binary = "NaiveBayes",
      rda_binary = "rda",
      nnet_binary = "nnet",
      svm_polynomial_binary = "svm",
      svm_sigmoid_binary = "svm"
    )
    expect_s3_class(native_engine, expected_class)
    native <- switch(learner,
      kernel_naive_bayes_binary = stats::predict(
        model$engine,
        newdata = assessment,
        type = "raw"
      )$posterior[, "1"],
      rda_binary = get("predict.rda", envir = asNamespace("klaR"), inherits = FALSE)(
        model$engine,
        newdata = assessment
      )$posterior[, "1"],
      nnet_binary = {
        matrix_x <- as.matrix(assessment)
        scaled <- sweep(
          sweep(matrix_x, 2L, model$engine$center, "-"),
          2L,
          model$engine$scale,
          "/"
        )
        as.numeric(stats::predict(
          model$engine$model,
          newdata = scaled,
          type = "raw"
        ))
      },
      svm_polynomial_binary = {
        prediction <- stats::predict(
          model$engine,
          newdata = as.matrix(assessment),
          probability = TRUE
        )
        attr(prediction, "probabilities", exact = TRUE)[, "1"]
      },
      svm_sigmoid_binary = {
        prediction <- stats::predict(
          model$engine,
          newdata = as.matrix(assessment),
          probability = TRUE
        )
        attr(prediction, "probabilities", exact = TRUE)[, "1"]
      }
    )
    expect_equal(observed, unname(native), tolerance = 1e-14)
    if (identical(learner, "kernel_naive_bayes_binary")) {
      expect_identical(model$engine$usekernel, TRUE)
    }
    if (identical(learner, "rda_binary")) {
      expect_equal(
        unname(model$engine$regularization),
        c(0.5, 0.5),
        tolerance = 0
      )
    }
    if (identical(learner, "nnet_binary")) {
      expect_identical(model$engine$model$convergence, 0L)
    }
    if (identical(learner, "svm_polynomial_binary")) {
      expect_identical(model$engine$kernel, 1)
    }
    if (identical(learner, "svm_sigmoid_binary")) {
      expect_identical(model$engine$kernel, 3)
    }
  }
})

test_that("sixth-batch fixed binary probabilities match native engines", {
  data <- sb_example_binary(n = 72L, seed = 6815L)
  x <- data[c("gene_a", "gene_b", "gene_c")]
  assessment <- x[1:8, , drop = FALSE]
  learners <- c(
    "bart_binary",
    "c50_binary",
    "kernel_weighted_knn_binary",
    "gaussian_process_binary",
    "pamr_binary",
    "neuralnet_multilayer_binary"
  )

  for (learner in learners) {
    dependency <- sb_classical_binary_dependency(learner)
    if (!dependency$available) {
      next
    }
    params <- sb_classical_binary_params(learner)
    if (identical(learner, "bart_binary")) {
      params$posterior_draws <- 40L
      params$burn_in <- 20L
    }
    adapter <- sb_classical_binary_learner_adapter(learner, params)
    model <- adapter$fit(x, data$outcome)
    observed <- adapter$predict(model, assessment)
    grouping <- factor(data$outcome, levels = c(0L, 1L))

    native <- switch(learner,
      bart_binary = {
        prediction <- NULL
        method <- get(
          "predict.pbart",
          envir = asNamespace("BART"),
          inherits = FALSE
        )
        invisible(utils::capture.output(
          prediction <- method(
            model$engine,
            sb_bart_prediction_matrix(model$engine, as.matrix(assessment)),
            mc.cores = 1L,
            openmp = FALSE
          ),
          type = "output"
        ))
        expect_s3_class(model$engine, "pbart")
        prediction$prob.test.mean
      },
      c50_binary = {
        expect_s3_class(model$engine, "C5.0")
        expect_identical(
          as.integer(unname(model$engine$trials["Actual"])),
          params$trials
        )
        get(
          "predict.C5.0",
          envir = asNamespace("C50"),
          inherits = FALSE
        )(model$engine, assessment, type = "prob")[, "1"]
      },
      kernel_weighted_knn_binary = {
        training <- x
        training$.stormbreaker_outcome <- grouping
        fitted <- kknn::kknn(
          .stormbreaker_outcome ~ .,
          train = training,
          test = assessment,
          k = params$k,
          distance = params$distance,
          kernel = params$kernel,
          scale = params$scale
        )
        expect_s3_class(fitted, "kknn")
        fitted$prob[, "1"]
      },
      gaussian_process_binary = {
        expect_s4_class(model$engine, "gausspr")
        kernlab::predict(
          model$engine,
          as.matrix(assessment),
          type = "probabilities"
        )[, "1"]
      },
      pamr_binary = {
        expect_s3_class(model$engine, "pamrtrained")
        expect_equal(model$engine$prior, as.numeric(table(grouping) / 72))
        pamr::pamr.predict(
          model$engine,
          newx = t(as.matrix(assessment)),
          threshold = params$threshold,
          type = "posterior"
        )[, "1"]
      },
      neuralnet_multilayer_binary = {
        expect_s3_class(model$engine$model, "nn")
        expect_length(model$engine$model$weights[[1L]], 3L)
        matrix_x <- as.matrix(assessment)
        standardized <- sweep(
          sweep(matrix_x, 2L, model$engine$center, "-"),
          2L,
          model$engine$scale,
          "/"
        )
        standardized <- as.data.frame(standardized)
        names(standardized) <- names(assessment)
        as.numeric(get(
          "predict.nn",
          envir = asNamespace("neuralnet"),
          inherits = FALSE
        )(model$engine$model, standardized))
      }
    )
    expect_equal(observed, unname(native), tolerance = 1e-14)
  }
})

test_that("binary BART replays after native constant-column removal", {
  skip_if_not_installed("BART", minimum_version = "2.9.0")

  data <- sb_example_binary(n = 60L, seed = 6816L)
  data$constant <- 1
  x <- data[c("gene_a", "constant", "gene_b")]
  params <- sb_classical_binary_params("bart_binary")
  params$posterior_draws <- 30L
  params$burn_in <- 15L
  adapter <- sb_classical_binary_learner_adapter(
    "bart_binary",
    params
  )
  model <- adapter$fit(x, data$outcome)
  assessment <- x[1:6, , drop = FALSE]

  observed <- adapter$predict(model, assessment)
  replay <- adapter$predict_state(
    adapter$extract_state(model),
    assessment
  )$probability

  expect_true(length(model$engine$rm.const) > 0L)
  expect_equal(replay, observed, tolerance = 0)
})

test_that("fixed classical binary learners fit and preserve predictions", {
  data <- sb_example_binary(n = 72L, seed = 6801L)
  predictors <- c("gene_a", "gene_b", "gene_c")
  analysis <- data[predictors]
  assessment <- data[1:8, predictors, drop = FALSE]

  for (learner in sb_classical_binary_learner_ids()) {
    dependency <- sb_classical_binary_dependency(learner)
    if (!dependency$available) {
      next
    }
    adapter <- sb_classical_binary_learner_adapter(learner)
    model <- adapter$fit(analysis, data$outcome)
    probability <- adapter$predict(model, assessment)
    state_probability <- adapter$predict_state(
      adapter$extract_state(model),
      assessment
    )$probability

    expect_identical(adapter$id, learner)
    expect_length(probability, nrow(assessment))
    expect_identical(all(is.finite(probability)), TRUE)
    expect_identical(all(probability >= 0 & probability <= 1), TRUE)
    expect_equal(state_probability, probability, tolerance = 0)
  }
})

test_that("fixed classical binary learners execute through CV and final fit", {
  data <- sb_example_binary(n = 72L, seed = 6802L)
  predictors <- c("gene_a", "gene_b", "gene_c")
  task <- sb_binary_task(
    data,
    "patient_id",
    "outcome",
    predictors = predictors
  )
  original_params <- sb_classical_binary_params
  local_mocked_bindings(
    sb_classical_binary_params = function(learner) {
      params <- original_params(learner)
      if (identical(learner, "cforest_binary")) {
        params$num_trees <- 20L
      }
      if (identical(learner, "bagging_binary")) {
        params$num_trees <- 20L
      }
      if (learner %in% c(
        "gamboost_binary",
        "blackboost_binary"
      )) {
        params$mstop <- 20L
      }
      if (identical(learner, "aorsf_binary")) {
        params$num_trees <- 20L
      }
      params
    }
  )

  for (learner in sb_classical_binary_learner_ids()) {
    dependency <- sb_classical_binary_dependency(learner)
    if (!dependency$available) {
      next
    }
    pipeline <- sb_pipeline(learner = learner)
    result <- sb_run_cv(task, pipeline, folds = 3L, seed = 6803L)
    fit <- sb_fit_final(task, pipeline)
    prediction <- sb_predict_final(
      fit,
      data[1:6, c("patient_id", predictors), drop = FALSE],
      id_col = "patient_id"
    )

    expect_identical(result$status, "success")
    expect_identical(result$complete_oof, TRUE)
    expect_identical(result$runs$learner_used, rep(learner, 3L))
    if (!is.na(dependency$package)) {
      field <- paste0(tolower(dependency$package), "_version")
      expect_identical(
        result$provenance[[field]],
        as.character(utils::packageVersion(dependency$package))
      )
      expect_identical(
        fit$provenance[[field]],
        as.character(utils::packageVersion(dependency$package))
      )
    }
    expect_identical(fit$status, "success")
    expect_identical(fit$learner_state$learner, learner)
    expect_identical(prediction$status, "success")
    expect_identical(
      prediction$predictions$learner,
      rep(learner, 6L)
    )
  }
})

test_that("fixed classical binary learners enforce prediction schema", {
  data <- sb_example_binary(n = 60L, seed = 6804L)
  adapter <- sb_classical_binary_learner_adapter("lda_binary")
  model <- adapter$fit(data[c("gene_a", "gene_b")], data$outcome)

  error <- capture_sb_error(
    adapter$predict(model, data[c("gene_b", "gene_a")])
  )

  expect_s3_class(error, "stormbreaker_learner_schema_mismatch")
})

test_that("weighted radial SVM derives fold-local inverse-frequency weights", {
  skip_if_not_installed("e1071", minimum_version = "1.7.0")

  data <- sb_example_binary(n = 72L, seed = 6805L)
  data$outcome <- c(rep(0L, 48L), rep(1L, 24L))
  x <- data[c("gene_a", "gene_b", "gene_c")]
  adapter <- sb_classical_binary_learner_adapter(
    "svm_radial_weighted_binary"
  )
  original_svm <- e1071::svm
  observed_weights <- NULL

  testthat::with_mocked_bindings(
    svm = function(..., class.weights = NULL) {
      observed_weights <<- class.weights
      original_svm(..., class.weights = class.weights)
    },
    adapter$fit(x, data$outcome),
    .package = "e1071"
  )

  counts <- table(factor(data$outcome, levels = c(0L, 1L)))
  expected <- setNames(
    length(data$outcome) / (2 * as.numeric(counts)),
    names(counts)
  )
  expect_equal(observed_weights, expected, tolerance = 0)
})

test_that("LDA collinearity is a typed learner failure", {
  data <- sb_example_binary(n = 60L, seed = 6805L)
  data$duplicate_gene_a <- data$gene_a
  predictors <- c("gene_a", "gene_b", "duplicate_gene_a")
  adapter <- sb_classical_binary_learner_adapter("lda_binary")

  error <- capture_sb_error(
    adapter$fit(data[predictors], data$outcome)
  )
  task <- sb_binary_task(
    data,
    "patient_id",
    "outcome",
    predictors = predictors
  )
  pipeline <- sb_pipeline(learner = "lda_binary")
  result <- sb_run_cv(task, pipeline, folds = 3L, seed = 6806L)
  fit <- sb_fit_final(task, pipeline)

  expect_s3_class(error, "stormbreaker_lda_collinear_predictors")
  expect_identical(result$status, "failed")
  expect_identical(result$complete_oof, FALSE)
  expect_identical(
    unique(result$runs$error_class),
    "stormbreaker_lda_collinear_predictors"
  )
  expect_identical(fit$status, "failed")
  expect_identical(fit$error_stage, "learner_fit")
  expect_identical(
    fit$error_class,
    "stormbreaker_lda_collinear_predictors"
  )
})

test_that("fixed neural-network nonconvergence is a typed learner failure", {
  skip_if_not_installed("nnet", minimum_version = "7.3")

  data <- sb_example_binary(n = 60L, seed = 6811L)
  params <- sb_classical_binary_params("nnet_binary")
  params$maxit <- 1L
  error <- capture_sb_error(sb_fit_classical_binary(
    data[c("gene_a", "gene_b", "gene_c")],
    data$outcome,
    "nnet_binary",
    params
  ))

  expect_s3_class(error, "stormbreaker_nnet_nonconvergence")
})

test_that("multilayer neural-network nonconvergence is typed", {
  skip_if_not_installed("neuralnet", minimum_version = "1.44.0")

  data <- sb_example_binary(n = 60L, seed = 6817L)
  params <- sb_classical_binary_params(
    "neuralnet_multilayer_binary"
  )
  params$stepmax <- 1L
  error <- capture_sb_error(sb_fit_classical_binary(
    data[c("gene_a", "gene_b", "gene_c")],
    data$outcome,
    "neuralnet_multilayer_binary",
    params
  ))

  expect_s3_class(error, "stormbreaker_neuralnet_nonconvergence")
})

test_that("RDA prediction is isolated from conflicting S3 registrations", {
  skip_if_not_installed("klaR", minimum_version = "1.7.0")

  data <- sb_example_binary(n = 60L, seed = 6812L)
  x <- data[c("gene_a", "gene_b", "gene_c")]
  adapter <- sb_classical_binary_learner_adapter("rda_binary")
  model <- adapter$fit(x, data$outcome)
  expected <- get(
    "predict.rda",
    envir = asNamespace("klaR"),
    inherits = FALSE
  )(model$engine, newdata = x[1:6, , drop = FALSE])$posterior[, "1"]
  original <- getS3method("predict", "rda")
  withr::defer(registerS3method(
    "predict",
    "rda",
    original,
    envir = globalenv()
  ))
  registerS3method(
    "predict",
    "rda",
    function(object, ...) stop("conflicting predict.rda", call. = FALSE),
    envir = globalenv()
  )
  generic_error <- tryCatch(
    stats::predict(model$engine, newdata = x[1:6, , drop = FALSE]),
    error = identity
  )
  observed <- adapter$predict(model, x[1:6, , drop = FALSE])

  expect_s3_class(generic_error, "error")
  expect_equal(observed, unname(expected), tolerance = 1e-14)
})

test_that("saved fixed learner engines omit non-prediction training payloads", {
  data <- sb_example_binary(n = 60L, seed = 6807L)
  predictors <- c("gene_a", "gene_b", "gene_c")
  analysis <- data[predictors]
  engine_for <- function(learner) {
    adapter <- sb_classical_binary_learner_adapter(learner)
    model <- adapter$fit(analysis, data$outcome)
    unserialize(adapter$extract_state(model)$engine_raw)
  }

  knn <- engine_for("knn_binary")
  expect_identical(nrow(knn$training), nrow(analysis))
  expect_identical(knn$outcome, data$outcome)

  if (sb_classical_binary_dependency(
    "decision_tree_binary"
  )$available) {
    tree <- engine_for("decision_tree_binary")
    expect_identical(
      intersect(c("model", "x", "y", "where"), names(tree)),
      character()
    )
    expect_identical(environment(tree$terms), baseenv())
  }

  if (sb_classical_binary_dependency(
    "random_forest_binary"
  )$available) {
    forest <- engine_for("random_forest_binary")
    expect_identical(
      intersect(
        c(
          "predicted", "err.rate", "confusion", "votes", "oob.times",
          "y", "test", "inbag", "localImportance", "proximity"
        ),
        names(forest)
      ),
      character()
    )
  }

  if (sb_classical_binary_dependency(
    "svm_linear_binary"
  )$available) {
    svm <- engine_for("svm_linear_binary")
    expect_identical(
      intersect(
        c("call", "index", "na.action", "fitted", "decision.values"),
        names(svm)
      ),
      character()
    )
  }

  if (sb_classical_binary_dependency(
    "kernel_naive_bayes_binary"
  )$available) {
    kernel_nb <- engine_for("kernel_naive_bayes_binary")
    expect_identical(nrow(kernel_nb$x), 0L)
    expect_identical("call" %in% names(kernel_nb), FALSE)
  }

  if (sb_classical_binary_dependency("rda_binary")$available) {
    rda <- engine_for("rda_binary")
    expect_identical(
      intersect(c("call", "error.rate", "converged", "iter"), names(rda)),
      character()
    )
  }

  if (sb_classical_binary_dependency("nnet_binary")$available) {
    neural_net <- engine_for("nnet_binary")
    expect_identical(
      intersect(
        c("call", "fitted.values", "residuals"),
        names(neural_net$model)
      ),
      character()
    )
  }

  if (sb_classical_binary_dependency("bart_binary")$available) {
    bart <- engine_for("bart_binary")
    expect_true(all(names(bart) %in% c(
      "treedraws", "rm.const", "binaryOffset", "offset"
    )))
    expect_false(any(c("yhat.train", "prob.train", "varcount") %in%
      names(bart)))
  }

  if (sb_classical_binary_dependency("c50_binary")$available) {
    c50 <- engine_for("c50_binary")
    expect_identical(
      intersect(c("call", "output", "boostResults", "size"), names(c50)),
      character()
    )
  }

  if (sb_classical_binary_dependency(
    "kernel_weighted_knn_binary"
  )$available) {
    weighted_knn <- engine_for("kernel_weighted_knn_binary")
    expect_identical(nrow(weighted_knn$training), nrow(analysis))
    expect_identical(as.integer(weighted_knn$outcome), data$outcome + 1L)
  }

  if (sb_classical_binary_dependency(
    "gaussian_process_binary"
  )$available) {
    gaussian_process <- engine_for("gaussian_process_binary")
    expect_s4_class(gaussian_process, "gausspr")
    expect_identical(methods::slot(gaussian_process, "fitted"), NA)
  }

  if (sb_classical_binary_dependency("pamr_binary")$available) {
    pamr <- engine_for("pamr_binary")
    expect_identical(
      intersect(c("y", "yhat", "prob", "errors", "call"), names(pamr)),
      character()
    )
  }

  if (sb_classical_binary_dependency(
    "neuralnet_multilayer_binary"
  )$available) {
    neuralnet <- engine_for("neuralnet_multilayer_binary")
    expect_setequal(
      names(neuralnet$model),
      c("model.list", "weights", "act.fct", "linear.output")
    )
  }

  if (sb_classical_binary_dependency("gbm_binary")$available) {
    gbm <- engine_for("gbm_binary")
    expect_identical(
      intersect(
        c("fit", "train.error", "valid.error", "oobag.improve"),
        names(gbm)
      ),
      character()
    )
  }
})

test_that("CForest safely predicts more rows than its training data", {
  skip_if_not_installed("partykit", minimum_version = "1.2.0")

  data <- sb_example_binary(n = 72L, seed = 6821L)
  x <- data[c("gene_a", "gene_b", "gene_c")]
  params <- sb_classical_binary_params("cforest_binary")
  params$num_trees <- 20L
  adapter <- sb_classical_binary_learner_adapter(
    "cforest_binary", params
  )
  model <- adapter$fit(x, data$outcome)
  training_data <- model$engine$data
  rows <- rep(rev(seq_len(nrow(x))), length.out = 2L * nrow(x) + 5L)
  assessment <- x[rows, , drop = FALSE]
  chunks <- split(seq_len(nrow(assessment)), ceiling(
    seq_len(nrow(assessment)) / nrow(x)
  ))
  native_predict <- stats::predict
  expected <- unlist(lapply(chunks, function(index) {
    native_predict(
      model$engine, newdata = assessment[index, , drop = FALSE], type = "prob"
    )[, "1"]
  }), use.names = FALSE)
  native_batch_sizes <- integer()
  local_mocked_bindings(
    predict = function(object, newdata, ...) {
      if (inherits(object, "cforest")) {
        native_batch_sizes <<- c(native_batch_sizes, nrow(newdata))
        if (nrow(newdata) > nrow(object$data)) {
          stop("The synthetic test intercepted an unsafe native CForest batch.")
        }
      }
      native_predict(object, newdata = newdata, ...)
    },
    .package = "stats"
  )

  observed <- adapter$predict(model, assessment)
  state <- sb_extract_classical_binary_learner_state(model)
  replayed <- sb_predict_classical_binary_learner_state(
    state, assessment
  )$probability
  expect_equal(observed, expected, tolerance = 0)
  expect_equal(replayed, expected, tolerance = 0)
  expect_identical(native_batch_sizes, rep(unname(lengths(chunks)), 2L))
  expect_identical(model$params, params)
  expect_identical(model$engine$data, training_data)
  expect_equal(
    adapter$predict(model, assessment[1:8, , drop = FALSE]),
    expected[1:8],
    tolerance = 0
  )
  expect_identical(tail(native_batch_sizes, 1L), 8L)
})

test_that("new fixed learners run with learned and fixed feature selectors", {
  data <- sb_example_binary(n = 72L, seed = 6808L)
  task <- sb_binary_task(
    data,
    "patient_id",
    "outcome",
    predictors = c("gene_a", "gene_b", "gene_c")
  )
  pipelines <- list(
    sb_pipeline(
      selector = "univ_binary_auc_top_k",
      learner = "lda_binary",
      selector_n_features = 2L
    ),
    sb_pipeline(
      selector = "fixed_features",
      learner = "decision_tree_binary",
      selector_features = c("gene_b", "gene_a"),
      feature_set_id = "two_genes"
    ),
    sb_pipeline(
      selector = "univ_binary_auc_top_k",
      learner = "kernel_naive_bayes_binary",
      selector_n_features = 2L
    ),
    sb_pipeline(
      selector = "fixed_features",
      learner = "nnet_binary",
      selector_features = c("gene_b", "gene_a"),
      feature_set_id = "neural_two_genes"
    ),
    sb_pipeline(
      selector = "univ_binary_auc_top_k",
      learner = "bart_binary",
      selector_n_features = 2L
    ),
    sb_pipeline(
      selector = "fixed_features",
      learner = "pamr_binary",
      selector_features = c("gene_b", "gene_a"),
      feature_set_id = "pamr_two_genes"
    )
  )

  for (pipeline in pipelines) {
    dependency <- sb_classical_binary_dependency(pipeline$learner)
    if (!dependency$available) {
      next
    }
    result <- sb_run_cv(task, pipeline, folds = 3L, seed = 6809L)
    fit <- sb_fit_final(task, pipeline)

    expect_identical(result$status, "success")
    expect_identical(fit$status, "success")
    expect_length(fit$required_features, 2L)
  }
})
