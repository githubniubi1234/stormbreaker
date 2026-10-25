make_forest_selector_binary_fixture <- function(n = 100L, seed = 8711L) {
  withr::with_seed(seed, {
    outcome <- sample(rep(c(0L, 1L), length.out = n))
    data.frame(
      patient_id = sprintf("FSB%03d", seq_len(n)),
      outcome = outcome,
      gene_a = stats::rnorm(n, 1.2 * outcome),
      gene_b = stats::rnorm(n, -0.9 * outcome),
      gene_c = stats::rnorm(n, 0.7 * outcome),
      gene_d = stats::rnorm(n),
      gene_e = stats::rnorm(n),
      gene_f = stats::rnorm(n),
      check.names = FALSE
    )
  })
}

make_forest_selector_survival_fixture <- function(n = 100L, seed = 8712L) {
  withr::with_seed(seed, {
    gene_a <- stats::rnorm(n)
    gene_b <- stats::rnorm(n)
    gene_c <- stats::rnorm(n)
    gene_d <- stats::rnorm(n)
    gene_e <- stats::rnorm(n)
    gene_f <- stats::rnorm(n)
    event_time <- stats::rexp(
      n,
      rate = 0.06 * exp(1.1 * gene_a - 0.9 * gene_b + 0.5 * gene_c)
    ) + 0.1
    cutoff <- stats::median(event_time)
    data.frame(
      patient_id = sprintf("FSS%03d", seq_len(n)),
      time = round(pmin(event_time, cutoff), 4L),
      event = as.integer(event_time <= cutoff),
      gene_a = gene_a,
      gene_b = gene_b,
      gene_c = gene_c,
      gene_d = gene_d,
      gene_e = gene_e,
      gene_f = gene_f,
      check.names = FALSE
    )
  })
}

test_that("forest selectors expose fixed native specifications", {
  expect_identical(
    sb_forest_selector_ids(),
    c(
      "boruta_binary_roughfix",
      "boruta_survival_roughfix",
      "varselrf_binary_oob_1se"
    )
  )

  boruta <- sb_forest_selector_params("boruta_binary_roughfix")
  expect_identical(boruta$method, "Boruta")
  expect_identical(boruta$p_value, 0.01)
  expect_true(boruta$multiple_comparison_adjustment)
  expect_identical(boruta$maximum_runs, 100L)
  expect_identical(boruta$num_trees, 500L)
  expect_identical(boruta$num_threads, 1L)
  expect_identical(boruta$tentative_rule, "TentativeRoughFix")
  expect_identical(boruta$seed, 1L)

  varselrf <- sb_forest_selector_params("varselrf_binary_oob_1se")
  expect_identical(varselrf$method, "varSelRF")
  expect_identical(varselrf$c_sd, 1)
  expect_identical(varselrf$initial_trees, 5000L)
  expect_identical(varselrf$iteration_trees, 2000L)
  expect_identical(varselrf$drop_fraction, 0.2)
  expect_true(varselrf$whole_range)
  expect_false(varselrf$recompute_importance)
  expect_identical(varselrf$seed, 1L)
  expect_identical(
    names(sb_forest_selector_provenance("boruta_binary_roughfix")),
    c("boruta_version", "ranger_version")
  )
  expect_identical(
    names(sb_forest_selector_provenance("varselrf_binary_oob_1se")),
    c("varselrf_version", "randomforest_version")
  )
})

test_that("Boruta calls pin the ranger normalized-permutation backend", {
  skip_if_not_installed("Boruta", minimum_version = "9.0.0")

  captured_importance <- NULL
  observed <- testthat::with_mocked_bindings(
    sb_call_boruta(
      data.frame(gene_a = 1:4, gene_b = 4:1),
      factor(c(0L, 0L, 1L, 1L)),
      sb_forest_selector_params("boruta_binary_roughfix")
    ),
    Boruta = function(..., getImp) {
      captured_importance <<- getImp
      "captured"
    },
    .package = "Boruta"
  )

  expect_identical(observed, "captured")
  expect_identical(captured_importance, Boruta::getImpRfZ)
})

test_that("Boruta selectors reproduce native rough-fixed decisions", {
  skip_if_not_installed("Boruta", minimum_version = "9.0.0")
  skip_if_not_installed("ranger", minimum_version = "0.18.0")

  binary <- make_forest_selector_binary_fixture()
  x <- binary[paste0("gene_", letters[1:6])]
  adapter <- sb_forest_selector_adapter("boruta_binary_roughfix")
  params <- sb_forest_selector_params("boruta_binary_roughfix")

  set.seed(8713L)
  caller_rng <- .Random.seed
  fitted <- adapter$fit(x, binary$outcome)
  expect_identical(.Random.seed, caller_rng)

  native <- withr::with_seed(params$seed, Boruta::Boruta(
    x,
    factor(binary$outcome, levels = c(0L, 1L)),
    getImp = Boruta::getImpRfZ,
    pValue = params$p_value,
    mcAdj = params$multiple_comparison_adjustment,
    maxRuns = params$maximum_runs,
    doTrace = 0L,
    holdHistory = TRUE,
    num.trees = params$num_trees,
    num.threads = params$num_threads
  ))
  if (any(native$finalDecision == "Tentative")) {
    native <- Boruta::TentativeRoughFix(native)
  }

  expect_identical(
    fitted$selected_features,
    Boruta::getSelectedAttributes(native, withTentative = FALSE)
  )
  expect_identical(
    fitted$selected_features,
    fitted$scores$feature[fitted$scores$selected]
  )
  expect_null(fitted$engine)
  expect_null(fitted$call)

  survival <- make_forest_selector_survival_fixture()
  survival_x <- survival[paste0("gene_", letters[1:6])]
  survival_fit <- sb_forest_selector_adapter(
    "boruta_survival_roughfix"
  )$fit(survival_x, survival$time, survival$event)
  expect_gte(length(survival_fit$selected_features), 1L)
  expect_identical(
    names(sb_apply_forest_selector(survival_fit, survival_x)),
    survival_fit$selected_features
  )
  expect_null(survival_fit$engine)
  expect_null(survival_fit$call)
})

test_that("varSelRF selector reproduces its native one-SE OOB subset", {
  skip_if_not_installed("varSelRF", minimum_version = "0.7-9")
  skip_if_not_installed("randomForest", minimum_version = "4.7")

  binary <- make_forest_selector_binary_fixture(seed = 8721L)
  x <- binary[paste0("gene_", letters[1:6])]
  params <- sb_forest_selector_params("varselrf_binary_oob_1se")
  adapter <- sb_forest_selector_adapter("varselrf_binary_oob_1se")

  set.seed(8722L)
  caller_rng <- .Random.seed
  fitted <- adapter$fit(x, binary$outcome)
  expect_identical(.Random.seed, caller_rng)

  native <- withr::with_seed(params$seed, varSelRF::varSelRF(
    xdata = x,
    Class = factor(binary$outcome, levels = c(0L, 1L)),
    c.sd = params$c_sd,
    mtryFactor = params$mtry_factor,
    ntree = params$initial_trees,
    ntreeIterat = params$iteration_trees,
    vars.drop.num = NULL,
    vars.drop.frac = params$drop_fraction,
    whole.range = params$whole_range,
    recompute.var.imp = params$recompute_importance,
    verbose = FALSE,
    returnFirstForest = FALSE,
    keep.forest = FALSE
  ))

  expect_identical(fitted$selected_features, native$selected.vars)
  expect_gte(length(fitted$selected_features), 2L)
  expect_identical(nrow(fitted$selection_history), nrow(native$selec.history))
  expect_false("Vars.in.Forest" %in% names(fitted$selection_history))
  expect_null(fitted$engine)
  expect_null(fitted$call)
})

test_that("forest selectors add only legal generated routes", {
  registry <- sb_pipeline_registry()
  restricted_internal_cv <- c(
    setdiff(sb_rule_aft_learner_ids(), "bart_aft_survival"),
    sb_slope_binary_learner_ids(),
    "cancerclass_binary"
  )
  binary_learners <- registry$learner[
    registry$task_type == "binary" & registry$selector == "none" &
      !registry$learner %in% restricted_internal_cv
  ]
  survival_learners <- registry$learner[
    registry$task_type == "survival" & registry$selector == "none" &
      !registry$learner %in% c(
        sb_analysis_tuned_survival_learner_ids(),
        restricted_internal_cv
      )
  ]

  expect_setequal(
    registry$learner[registry$selector == "boruta_binary_roughfix"],
    binary_learners
  )
  expect_setequal(
    registry$learner[registry$selector == "varselrf_binary_oob_1se"],
    binary_learners
  )
  expect_setequal(
    registry$learner[registry$selector == "boruta_survival_roughfix"],
    survival_learners
  )
  expect_identical(
    nrow(registry[
      registry$selector == "boruta_survival_roughfix" &
        registry$learner %in% sb_analysis_tuned_survival_learner_ids(),
      ,
      drop = FALSE
    ]),
    0L
  )
})

test_that("forest selectors are repeated inside glmnet inner folds", {
  skip_if_not_installed("Boruta", minimum_version = "9.0.0")
  skip_if_not_installed("ranger", minimum_version = "0.18.0")
  skip_if_not_installed("varSelRF", minimum_version = "0.7-9")
  skip_if_not_installed("randomForest", minimum_version = "4.7")
  skip_if_not_installed("glmnet", minimum_version = "4.1.2")

  binary <- make_forest_selector_binary_fixture(n = 110L, seed = 8731L)
  binary_task <- sb_binary_task(
    binary,
    "patient_id",
    "outcome",
    paste0("gene_", letters[1:6])
  )
  survival <- make_forest_selector_survival_fixture(n = 110L, seed = 8732L)
  survival_task <- sb_surv_task(
    survival,
    "patient_id",
    "time",
    "event",
    paste0("gene_", letters[1:6]),
    "month"
  )

  fits <- list(
    sb_fit_final(binary_task, sb_pipeline(
      selector = "boruta_binary_roughfix",
      learner = "glmnet_logistic",
      inner_folds = 3L,
      inner_seed = 8733L
    )),
    sb_fit_final(binary_task, sb_pipeline(
      selector = "varselrf_binary_oob_1se",
      learner = "glmnet_logistic",
      inner_folds = 3L,
      inner_seed = 8734L
    )),
    sb_fit_final(survival_task, sb_pipeline(
      selector = "boruta_survival_roughfix",
      learner = "glmnet_cox",
      inner_folds = 3L,
      inner_seed = 8735L
    ))
  )

  for (fit in fits) {
    expect_identical(fit$status, "success")
    expect_identical(
      fit$tuning$inner_selection_method,
      "glmnet_exclude_function"
    )
    expect_length(fit$tuning$inner_selections, 3L)
    expect_true(all(vapply(fit$tuning$inner_selections, function(selection) {
      length(selection$analysis_patient_ids) > 0L &&
        length(selection$selected_features) > 0L
    }, logical(1))))
  }
  expect_identical(
    fits[[1L]]$provenance$boruta_version,
    as.character(utils::packageVersion("Boruta"))
  )
  expect_identical(
    fits[[1L]]$provenance$ranger_version,
    as.character(utils::packageVersion("ranger"))
  )
  expect_identical(
    fits[[2L]]$provenance$varselrf_version,
    as.character(utils::packageVersion("varSelRF"))
  )
  expect_identical(
    fits[[2L]]$provenance$randomforest_version,
    as.character(utils::packageVersion("randomForest"))
  )
  expect_identical(
    fits[[3L]]$provenance$boruta_version,
    as.character(utils::packageVersion("Boruta"))
  )
})
