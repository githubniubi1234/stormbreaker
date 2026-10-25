test_that("random ferns, robust logistic, and LS-SVM are distinct learners", {
  learners <- c(
    "random_ferns_binary",
    "robust_logistic_mqle",
    "lssvm_radial_binary"
  )

  expect_true(all(learners %in% sb_classical_binary_learner_ids()))
  expect_identical(
    sb_classical_binary_params("random_ferns_binary"),
    list(
      depth = 5L,
      num_ferns = 500L,
      importance = "none",
      save_forest = TRUE,
      threads = 1L,
      seed = 1L,
      score_mapping = "binary_softmax"
    )
  )
  expect_identical(
    sb_classical_binary_params("robust_logistic_mqle"),
    list(
      method = "Mqle",
      weights_on_x = "none",
      tuning_constant = 3.5,
      convergence_tolerance = 1e-4,
      convergence_test = "coef",
      max_iterations = 50L,
      family = "binomial_logit"
    )
  )
  expect_identical(
    sb_classical_binary_params("lssvm_radial_binary"),
    list(
      kernel = "rbfdot",
      kernel_parameters = "automatic",
      sigma_seed = 1L,
      tau = 0.01,
      scaled = TRUE,
      reduced = FALSE,
      cross_validation = 0L,
      keep_fitted = FALSE,
      score_mapping = "binary_softmax"
    )
  )

  dependencies <- vapply(
    learners,
    function(learner) {
      sb_classical_binary_dependency(learner)$package
    },
    character(1)
  )
  expect_identical(
    unname(dependencies),
    c("rFerns", "robustbase", "kernlab")
  )

  registry <- sb_pipeline_registry()
  expected_selectors <- expected_fixed_learner_selectors("binary")
  for (learner in learners) {
    rows <- registry[registry$learner == learner, , drop = FALSE]
    expect_identical(unique(rows$task_type), "binary")
    expect_setequal(rows$selector, expected_selectors)
  }
  expect_identical(
    unique(registry$tuning_semantics[
      registry$learner == "lssvm_radial_binary"
    ]),
    "analysis_fold_predictor_only_sigest"
  )
})

test_that("new binary learners replay their declared native predictions", {
  skip_if_not_installed("rFerns", minimum_version = "2.0.0")
  skip_if_not_installed("robustbase", minimum_version = "0.99.0")
  skip_if_not_installed("kernlab", minimum_version = "0.9.0")

  data <- make_sb_binary_fixture(n = 84L, seed = 8811L)
  x <- data[c("gene_a", "gene_b", "gene_c")]
  names(x) <- c("gene one", "gene-two", "gene_c")
  assessment <- x[c(3L, 19L, 47L, 72L), , drop = FALSE]

  learners <- c(
    "random_ferns_binary",
    "robust_logistic_mqle",
    "lssvm_radial_binary"
  )
  for (learner in learners) {
    params <- sb_classical_binary_params(learner)
    if (identical(learner, "random_ferns_binary")) {
      params$num_ferns <- 120L
    }
    adapter <- sb_classical_binary_learner_adapter(learner, params)
    model <- adapter$fit(x, data$outcome)
    observed <- adapter$predict(model, assessment)
    state <- unserialize(serialize(adapter$extract_state(model), NULL))
    replay <- adapter$predict_state(state, assessment)$probability

    native <- switch(learner,
      random_ferns_binary = {
        fit <- withr::with_seed(
          params$seed,
          rFerns::rFerns(
            x = x,
            y = factor(data$outcome, levels = c(0L, 1L)),
            depth = params$depth,
            ferns = params$num_ferns,
            importance = params$importance,
            saveForest = params$save_forest,
            threads = params$threads
          )
        )
        score <- as.matrix(stats::predict(
          fit,
          x = assessment,
          scores = TRUE
        ))
        stats::plogis(score[, "1"] - score[, "0"])
      },
      robust_logistic_mqle = {
        internal_names <- paste0(
          ".stormbreaker_robust_x",
          seq_len(ncol(x))
        )
        analysis <- x
        names(analysis) <- internal_names
        analysis$.stormbreaker_robust_outcome <- factor(
          data$outcome,
          levels = c(0L, 1L)
        )
        fit <- robustbase::glmrob(
          stats::reformulate(
            internal_names,
            response = ".stormbreaker_robust_outcome",
            env = baseenv()
          ),
          data = analysis,
          family = stats::binomial(link = "logit"),
          method = params$method,
          weights.on.x = params$weights_on_x,
          control = robustbase::glmrobMqle.control(
            acc = params$convergence_tolerance,
            test.acc = params$convergence_test,
            maxit = params$max_iterations,
            tcc = params$tuning_constant
          ),
          model = FALSE,
          x = FALSE
        )
        native_new <- assessment
        names(native_new) <- internal_names
        unname(stats::predict(fit, native_new, type = "response"))
      },
      lssvm_radial_binary = {
        fit <- withr::with_seed(params$sigma_seed, {
          invisible(utils::capture.output(
            native_fit <- kernlab::lssvm(
              x = as.matrix(x),
              y = factor(data$outcome, levels = c(0L, 1L)),
              scaled = params$scaled,
              kernel = params$kernel,
              kpar = params$kernel_parameters,
              type = "classification",
              tau = params$tau,
              reduced = params$reduced,
              cross = params$cross_validation,
              fit = params$keep_fitted
            )
          ))
          native_fit
        })
        score <- kernlab::predict(
          fit,
          as.matrix(assessment),
          type = "decision"
        )
        stats::plogis(score[, 2L] - score[, 1L])
      }
    )

    expect_equal(observed, unname(as.numeric(native)), tolerance = 1e-12)
    expect_equal(replay, observed, tolerance = 0)
    expect_true(all(observed > 0 & observed < 1))
  }

  fern_params <- sb_classical_binary_params(
    "random_ferns_binary"
  )
  fern_params$num_ferns <- 120L
  fern_adapter <- sb_classical_binary_learner_adapter(
    "random_ferns_binary",
    fern_params
  )
  fern_state <- fern_adapter$extract_state(
    fern_adapter$fit(x, data$outcome)
  )
  fern_engine <- unserialize(fern_state$engine_raw)
  expect_setequal(
    names(fern_engine),
    c("model", "classLabels", "isStruct", "parameters", "type")
  )
  expect_false("training_data_retention" %in% names(fern_state))

  robust_adapter <- sb_classical_binary_learner_adapter(
    "robust_logistic_mqle"
  )
  robust_state <- robust_adapter$extract_state(
    robust_adapter$fit(x, data$outcome)
  )
  robust_engine <- unserialize(robust_state$engine_raw)
  expect_setequal(
    names(robust_engine),
    c("intercept", "coefficient", "iterations")
  )
  expect_false("training_data_retention" %in% names(robust_state))

  lssvm_adapter <- sb_classical_binary_learner_adapter(
    "lssvm_radial_binary"
  )
  lssvm_state <- lssvm_adapter$extract_state(
    lssvm_adapter$fit(x, data$outcome)
  )
  lssvm_engine <- unserialize(lssvm_state$engine_raw)
  expect_identical(
    lssvm_state$training_data_retention,
    "required_scaled_analysis_predictors_for_radial_kernel_prediction"
  )
  expect_identical(nrow(lssvm_engine$scaled_training), nrow(x))
  expect_false("outcome" %in% names(lssvm_engine))
})

test_that("fixed Random Ferns isolates its deterministic seed", {
  skip_if_not_installed("rFerns", minimum_version = "2.0.0")

  data <- make_sb_binary_fixture(n = 72L, seed = 8812L)
  x <- data[c("gene_a", "gene_b", "gene_c")]
  params <- sb_classical_binary_params("random_ferns_binary")
  params$num_ferns <- 100L
  adapter <- sb_classical_binary_learner_adapter(
    "random_ferns_binary",
    params
  )
  withr::local_seed(8813L)
  rng_before <- .Random.seed
  first <- adapter$fit(x, data$outcome)
  expect_identical(.Random.seed, rng_before)
  second <- adapter$fit(x, data$outcome)
  expect_identical(.Random.seed, rng_before)
  expect_identical(first$engine$model, second$engine$model)
  expect_equal(
    adapter$predict(first, x[1:8, , drop = FALSE]),
    adapter$predict(second, x[1:8, , drop = FALSE]),
    tolerance = 0
  )
})

test_that("fixed radial LS-SVM isolates its predictor-only sigest seed", {
  skip_if_not_installed("kernlab", minimum_version = "0.9.0")

  data <- make_sb_binary_fixture(n = 72L, seed = 8814L)
  x <- data[c("gene_a", "gene_b", "gene_c")]
  adapter <- sb_classical_binary_learner_adapter(
    "lssvm_radial_binary"
  )
  withr::local_seed(8815L)
  rng_before <- .Random.seed
  first <- adapter$fit(x, data$outcome)
  expect_identical(.Random.seed, rng_before)
  second <- adapter$fit(x, data$outcome)
  expect_identical(.Random.seed, rng_before)
  expect_equal(first$engine$sigma, second$engine$sigma, tolerance = 0)
  expect_gt(first$engine$sigma, 0)
  expect_equal(
    adapter$predict(first, x[1:8, , drop = FALSE]),
    adapter$predict(second, x[1:8, , drop = FALSE]),
    tolerance = 0
  )
})

test_that("new binary learners execute through fold-local CV and final fit", {
  skip_if_not_installed("rFerns", minimum_version = "2.0.0")
  skip_if_not_installed("robustbase", minimum_version = "0.99.0")
  skip_if_not_installed("kernlab", minimum_version = "0.9.0")

  data <- make_sb_binary_fixture(n = 72L, seed = 8821L)
  task <- sb_binary_task(
    data,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b", "gene_c")
  )
  dependencies <- c(
    random_ferns_binary = "rferns_version",
    robust_logistic_mqle = "robustbase_version",
    lssvm_radial_binary = "kernlab_version"
  )
  for (learner in names(dependencies)) {
    pipeline <- sb_pipeline(
      selector = "univ_binary_auc_top_k",
      learner = learner,
      selector_n_features = 2L
    )
    comparison <- sb_run_cv(task, pipeline, folds = 3L, seed = 8822L)
    fit <- sb_fit_final(task, pipeline)
    prediction <- sb_predict_final(
      fit,
      data[1:6, c("patient_id", "gene_c", "gene_a", "gene_b")],
      "patient_id"
    )

    expect_identical(comparison$status, "success")
    expect_identical(comparison$complete_oof, TRUE)
    expect_true(all(comparison$runs$status == "success"))
    expect_true(unname(dependencies[[learner]]) %in% names(
      comparison$provenance
    ))
    expect_identical(fit$status, "success")
    expect_identical(prediction$status, "success")
    expect_length(prediction$predictions$probability, 6L)
  }
})
