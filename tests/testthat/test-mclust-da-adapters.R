make_mclust_da_fixture <- function(n = 160L, seed = 9201L) {
  withr::with_seed(seed, {
    gene_a <- stats::rnorm(n)
    gene_b <- stats::rnorm(n)
    gene_c <- stats::rnorm(n)
    gene_d <- stats::rnorm(n)
    probability <- stats::plogis(
      1.1 * gene_a - 0.7 * gene_b + 0.35 * gene_c
    )
    data.frame(
      patient_id = sprintf("MCLUST%03d", seq_len(n)),
      outcome = stats::rbinom(n, 1L, probability),
      check.names = FALSE,
      "gene A" = gene_a,
      "gene-B" = gene_b,
      "gene.C" = gene_c,
      "gene D" = gene_d
    )
  })
}

native_mclust_da_fit <- function(x, outcome, params) {
  environment <- new.env(parent = asNamespace("mclust"))
  environment$data <- data.matrix(x)
  environment$class <- factor(outcome, levels = c(0L, 1L))
  environment$G <- params$components
  environment$modelNames <- if (ncol(x) == 1L) {
    params$univariate_covariance_models
  } else {
    params$multivariate_covariance_models
  }
  environment$modelType <- params$model_type
  environment$control <- do.call(mclust::emControl, params$em_control)
  eval(
    quote(MclustDA(
      data = data,
      class = class,
      G = G,
      modelNames = modelNames,
      modelType = modelType,
      prior = NULL,
      control = control,
      warn = FALSE,
      verbose = FALSE
    )),
    envir = environment
  )
}

test_that("mclust DA learners declare one analysis-fold BIC specification", {
  skip_if_not_installed("mclust", minimum_version = "6.1.0")

  expect_identical(
    sb_mclust_da_learner_ids(),
    c("edda_binary_bic", "mclustda_binary_bic")
  )
  edda <- sb_mclust_da_params("edda_binary_bic")
  mixture <- sb_mclust_da_params("mclustda_binary_bic")

  expect_identical(edda$model_type, "EDDA")
  expect_identical(edda$components, 1L)
  expect_identical(mixture$model_type, "MclustDA")
  expect_identical(mixture$components, 1:3)
  expect_identical(edda$univariate_covariance_models, c("E", "V"))
  expect_identical(
    edda$multivariate_covariance_models,
    mclust::mclust.options("emModelNames")
  )
  expect_identical(edda$class_prior, "analysis_fold_empirical")
  expect_identical(edda$parameter_prior, "none")
  expect_identical(edda$scale, "original")
  expect_identical(edda$criterion, "BIC")
  expect_identical(edda$native_candidate_warnings, FALSE)
  expect_identical(edda$all_model_failure, "typed")
  expect_identical(
    edda$initialization,
    list(method = "known_class_labels")
  )
  expect_identical(
    mixture$initialization,
    list(
      method = "native_hierarchical",
      subset_max = 2000L,
      hc_model_name = "VVV",
      hc_use = "SVD",
      seed = 1L
    )
  )

  expect_identical(
    sb_pipeline(learner = "edda_binary_bic")$learner_params,
    edda
  )
  expect_identical(
    sb_pipeline(learner = "mclustda_binary_bic")$learner_params,
    mixture
  )
})

test_that("mclust DA probabilities match the native posterior", {
  skip_if_not_installed("mclust", minimum_version = "6.1.0")
  data <- make_mclust_da_fixture()
  features <- c("gene A", "gene-B", "gene.C", "gene D")
  x <- data[features]
  new_x <- x[1:13, , drop = FALSE]

  for (learner in sb_mclust_da_learner_ids()) {
    params <- sb_mclust_da_params(learner)
    adapter <- sb_mclust_da_learner_adapter(learner)
    model <- adapter$fit(x, data$outcome)
    native <- native_mclust_da_fit(x, data$outcome, params)
    expected <- unname(stats::predict(
      native,
      newdata = data.matrix(new_x)
    )$z[, "1"])

    expect_equal(adapter$predict(model, new_x), expected, tolerance = 1e-12)
    expect_identical(model$audit$class_prior, native$prop)
    expect_identical(model$audit$selected$class, c("0", "1"))
  }
})

test_that("mclust DA uses the native univariate model grid", {
  skip_if_not_installed("mclust", minimum_version = "6.1.0")
  data <- make_mclust_da_fixture(seed = 9202L)
  x <- data["gene A"]

  for (learner in sb_mclust_da_learner_ids()) {
    adapter <- sb_mclust_da_learner_adapter(learner)
    model <- adapter$fit(x, data$outcome)
    probability <- adapter$predict(model, x[1:9, , drop = FALSE])

    expect_length(probability, 9L)
    expect_equal(range(probability), pmin(pmax(range(probability), 0), 1))
    expect_identical(model$params$active_covariance_models, c("E", "V"))
  }
})

test_that("mclust DA rejects a native silent reduction of the BIC grid", {
  skip_if_not_installed("mclust", minimum_version = "6.1.0")
  data <- make_mclust_da_fixture(seed = 9203L)
  x <- data[c("gene A", "gene-B", "gene.C", "gene D")]
  params <- sb_mclust_da_params("mclustda_binary_bic")
  engine <- native_mclust_da_fit(x, data$outcome, params)
  engine$models[[1L]]$BIC <- engine$models[[1L]]$BIC[-3L, , drop = FALSE]

  error <- tryCatch(
    {
      sb_validate_mclust_da_search_grid(engine, params)
      NULL
    },
    error = identity
  )
  expect_s3_class(error, "stormbreaker_mclust_da_search_reduction")
})

test_that("EDDA rejects its real dimension-based covariance-grid reduction", {
  skip_if_not_installed("mclust", minimum_version = "6.1.0")
  x <- withr::with_seed(92031L, {
    as.data.frame(matrix(
      stats::rnorm(20L * 10L),
      nrow = 20L,
      ncol = 10L
    ))
  })
  names(x) <- sprintf("gene_%02d", seq_len(ncol(x)))
  outcome <- rep(0:1, each = 10L)

  error <- tryCatch(
    {
      sb_mclust_da_learner_adapter("edda_binary_bic")$fit(x, outcome)
      NULL
    },
    error = identity
  )

  expect_s3_class(error, "stormbreaker_mclust_da_search_reduction")
})

test_that("EDDA keeps the full grid when one class exceeds the dimension", {
  skip_if_not_installed("mclust", minimum_version = "6.1.0")
  outcome <- c(rep(0L, 10L), rep(1L, 20L))
  x <- withr::with_seed(92032L, {
    matrix_x <- matrix(
      stats::rnorm(30L * 15L),
      nrow = 30L,
      ncol = 15L
    )
    matrix_x[, 1L] <- matrix_x[, 1L] + 0.6 * outcome
    as.data.frame(matrix_x)
  })
  names(x) <- sprintf("gene_%02d", seq_len(ncol(x)))

  model <- sb_mclust_da_learner_adapter("edda_binary_bic")$fit(x, outcome)

  expect_s3_class(model, "stormbreaker_mclust_da")
  expect_identical(
    model$audit$candidate_summary$requested,
    c(14L, 14L)
  )
})

test_that("mclust DA portable states omit analysis rows and replay", {
  skip_if_not_installed("mclust", minimum_version = "6.1.0")
  data <- make_mclust_da_fixture(seed = 9204L)
  features <- c("gene A", "gene-B", "gene.C", "gene D")
  x <- data[features]
  new_x <- x[1:11, , drop = FALSE]

  for (learner in sb_mclust_da_learner_ids()) {
    adapter <- sb_mclust_da_learner_adapter(learner)
    model <- adapter$fit(x, data$outcome)
    state <- unserialize(serialize(adapter$extract_state(model, list()), NULL))

    expect_equal(
      adapter$predict_state(state, new_x)$probability,
      adapter$predict(model, new_x),
      tolerance = 1e-12
    )
    expect_setequal(names(state$engine), c("type", "prop", "models"))
    for (class_model in state$engine$models) {
      expect_setequal(names(class_model), c("modelName", "parameters"))
    }
    expect_null(state$data)
    expect_null(state$class)
    expect_null(state$outcome)
    expect_null(state$engine$data)
    expect_null(state$engine$class)
    expect_null(state$engine$call)
  }
})

test_that("mclust DA does not alter the caller RNG state", {
  skip_if_not_installed("mclust", minimum_version = "6.1.0")
  data <- make_mclust_da_fixture(seed = 9205L)
  x <- data[c("gene A", "gene-B", "gene.C", "gene D")]

  set.seed(9206L)
  before <- .Random.seed
  for (learner in sb_mclust_da_learner_ids()) {
    sb_mclust_da_learner_adapter(learner)$fit(x, data$outcome)
  }
  expect_identical(.Random.seed, before)
})

test_that("MclustDA contains and audits sampled hierarchical initialization", {
  skip_if_not_installed("mclust", minimum_version = "6.1.0")
  if (!"package:mclust" %in% search()) {
    suppressWarnings(withr::local_package("mclust"))
  }
  original_options <- list(
    subset = mclust::mclust.options("subset"),
    hcModelName = mclust::mclust.options("hcModelName"),
    hcUse = mclust::mclust.options("hcUse")
  )
  withr::defer(do.call(mclust::mclust.options, original_options))

  n <- 240L
  outcome <- rep(0:1, each = n / 2L)
  x <- withr::with_seed(92061L, {
    data.frame(
      gene_a = stats::rnorm(n) + 0.8 * outcome,
      gene_b = stats::rnorm(n) - 0.5 * outcome,
      gene_c = stats::rnorm(n),
      gene_d = stats::rnorm(n)
    )
  })
  params <- sb_mclust_da_params("mclustda_binary_bic")
  params$initialization$subset_max <- 50L
  adapter <- sb_mclust_da_learner_adapter(
    "mclustda_binary_bic",
    params = params
  )
  new_x <- x[1:15, , drop = FALSE]

  mclust::mclust.options(
    subset = 40,
    hcModelName = "EII",
    hcUse = "VARS"
  )
  set.seed(92062L)
  rng_before_first <- .Random.seed
  first <- adapter$fit(x, outcome)
  first_probability <- adapter$predict(first, new_x)

  expect_identical(.Random.seed, rng_before_first)
  expect_identical(mclust::mclust.options("subset"), 40)
  expect_identical(mclust::mclust.options("hcModelName"), "EII")
  expect_identical(mclust::mclust.options("hcUse"), "VARS")
  expect_identical(
    first$audit$candidate_summary$initialization_subset_size,
    c(50L, 50L)
  )

  mclust::mclust.options(
    subset = 30,
    hcModelName = "VVV",
    hcUse = "SVD"
  )
  set.seed(92063L)
  rng_before_second <- .Random.seed
  second <- adapter$fit(x, outcome)
  second_probability <- adapter$predict(second, new_x)

  expect_identical(.Random.seed, rng_before_second)
  expect_identical(mclust::mclust.options("subset"), 30)
  expect_identical(mclust::mclust.options("hcModelName"), "VVV")
  expect_identical(mclust::mclust.options("hcUse"), "SVD")
  expect_identical(
    second$audit$candidate_summary$initialization_subset_size,
    c(50L, 50L)
  )
  expect_identical(first$audit$selected, second$audit$selected)
  expect_equal(first_probability, second_probability, tolerance = 1e-12)
})

test_that("mclust DA learners run through CV and outcome-blind prediction", {
  skip_if_not_installed("mclust", minimum_version = "6.1.0")
  data <- make_mclust_da_fixture(n = 180L, seed = 9207L)
  features <- c("gene A", "gene-B", "gene.C", "gene D")
  task <- sb_binary_task(data, "patient_id", "outcome", features)
  new_data <- data[1:9, c("patient_id", features), drop = FALSE]

  for (learner in sb_mclust_da_learner_ids()) {
    pipeline <- sb_pipeline(
      selector = "univ_binary_auc_top_k",
      learner = learner,
      selector_n_features = 3L
    )
    cv <- sb_run_cv(task, pipeline, folds = 3L, seed = 9208L)
    fit <- sb_fit_final(task, pipeline)
    prediction <- sb_predict_final(fit, new_data, "patient_id")

    expect_identical(cv$status, "success")
    expect_identical(cv$complete_oof, TRUE)
    expect_identical(length(cv$selections), 3L)
    expect_identical(
      vapply(cv$selections, function(selection) {
        length(selection$selected_features)
      }, integer(1)),
      rep.int(3L, 3L)
    )
    expect_identical(
      unique(vapply(cv$coefficients, `[[`, "", "status")),
      "not_applicable"
    )
    expect_identical(
      cv$provenance$mclust_version,
      as.character(utils::packageVersion("mclust"))
    )
    expect_identical(fit$status, "success")
    expect_identical(prediction$status, "success")
    expect_identical(
      fit$provenance$mclust_version,
      as.character(utils::packageVersion("mclust"))
    )
    expect_equal(nrow(prediction$predictions), nrow(new_data))
    expect_equal(
      range(prediction$predictions$probability),
      pmin(pmax(range(prediction$predictions$probability), 0), 1)
    )
  }
})

test_that("mclust DA learners retain the supported selector routes", {
  registry <- sb_pipeline_registry()
  rows <- registry[
    registry$learner %in% c("edda_binary_bic", "mclustda_binary_bic"), ,
    drop = FALSE
  ]

  expected_selectors <- expected_fixed_learner_selectors("binary")
  expect_identical(nrow(rows), 2L * length(expected_selectors))
  for (learner in c("edda_binary_bic", "mclustda_binary_bic")) {
    expect_setequal(rows$selector[rows$learner == learner], expected_selectors)
  }
  expect_identical(unique(rows$tuning_semantics), "analysis_fold_bic")
  expect_identical(
    all(grepl("mclust", rows$optional_dependency, fixed = TRUE)),
    TRUE
  )
})
