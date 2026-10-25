test_that("the admitted pipeline resolves explicit component adapters", {
  runtime <- sb_resolve_pipeline_runtime(sb_pipeline())

  expect_s3_class(runtime, "stormbreaker_pipeline_runtime")
  expect_identical(runtime$id, "none__coxph")
  expect_identical(runtime$selector$id, "none")
  expect_type(runtime$selector$fit, "closure")
  expect_type(runtime$selector$apply, "closure")
  expect_identical(runtime$learner$id, "coxph")
  expect_type(runtime$learner$fit, "closure")
  expect_type(runtime$learner$predict, "closure")
  expect_type(runtime$learner$predict_survival, "closure")
  expect_type(runtime$learner$extract_coefficients, "closure")
})

test_that("uniform learner adapter families are distinct and admitted", {
  families <- sb_uniform_learner_adapter_families()
  learner_ids <- unlist(
    lapply(families, function(family) family[[1L]]),
    use.names = FALSE
  )
  admitted <- sb_pipeline_registry(level = "learner")$learner

  expect_identical(anyDuplicated(learner_ids), 0L)
  expect_identical(setdiff(learner_ids, admitted), character())
})

test_that("glmnet logistic resolves an explicit tuning runtime", {
  pipeline <- sb_pipeline(
    learner = "glmnet_logistic",
    glmnet_alpha = 0.5,
    inner_folds = 3L,
    inner_seed = 101L
  )
  runtime <- sb_resolve_pipeline_runtime(pipeline)

  expect_s3_class(runtime, "stormbreaker_pipeline_runtime")
  expect_identical(runtime$id, pipeline$id)
  expect_identical(runtime$selector$id, "none")
  expect_identical(runtime$learner$id, "glmnet_logistic")
  expect_identical(runtime$learner$params, pipeline$learner_params)
  expect_identical(runtime$learner$requires_tuning, TRUE)
  expect_identical(
    sb_validate_pipeline_runtime(runtime, pipeline),
    runtime
  )
})

test_that("fixed CatBoost resolves an executable tree runtime", {
  pipeline <- sb_pipeline(
    learner = "catboost_binary",
    catboost_iterations = 7L,
    catboost_depth = 3L,
    catboost_seed = 101L
  )
  runtime <- sb_resolve_pipeline_runtime(pipeline)

  expect_s3_class(runtime, "stormbreaker_pipeline_runtime")
  expect_identical(runtime$id, pipeline$id)
  expect_identical(runtime$learner$id, "catboost_binary")
  expect_identical(runtime$learner$params, pipeline$learner_params)
  expect_type(runtime$learner$fit, "closure")
  expect_type(runtime$learner$predict, "closure")
  expect_type(runtime$learner$extract_tree_model, "closure")
  expect_identical(
    sb_validate_pipeline_runtime(runtime, pipeline),
    runtime
  )
})

test_that("the runtime resolver rejects unadmitted component identifiers", {
  pipeline <- sb_pipeline()
  pipeline$selector <- "unknown_selector"
  selector_error <- capture_sb_error(
    sb_resolve_pipeline_runtime(pipeline)
  )
  expect_s3_class(selector_error, "stormbreaker_unresolved_selector_adapter")

  pipeline <- sb_pipeline()
  pipeline$learner <- "unknown_learner"
  learner_error <- capture_sb_error(
    sb_resolve_pipeline_runtime(pipeline)
  )
  expect_s3_class(learner_error, "stormbreaker_unresolved_learner_adapter")
})

test_that("an inconsistent runtime is rejected as a visible fold failure", {
  d <- make_sb_fixture(n = 20)
  task <- sb_surv_task(
    d, "patient_id", "time", "event", "gene_a", "month"
  )
  split <- list(
    fold = 1L,
    analysis_rows = 1:10,
    assessment_rows = 11:20,
    analysis_ids = d$patient_id[1:10],
    assessment_ids = d$patient_id[11:20]
  )
  runtime <- sb_resolve_pipeline_runtime(sb_pipeline())
  runtime$learner$id <- "different_learner"

  result <- sb_run_one_fold(
    task, sb_pipeline(), runtime, split, seed = 3L
  )

  expect_identical(result$run$status, "failed")
  expect_identical(
    result$run$error_class,
    "stormbreaker_bad_pipeline_runtime"
  )
  expect_identical(all(is.na(result$predictions$risk)), TRUE)
})

test_that("the none selector locks names without retaining patient data", {
  x <- data.frame(
    gene_b = c(4, 5, 6),
    gene_a = c(1, 2, 3),
    check.names = FALSE
  )
  selector <- sb_none_selector_adapter()
  fitted <- selector$fit(x)

  expect_s3_class(fitted, "stormbreaker_selector_fit")
  expect_named(fitted, c("id", "selected_features"))
  expect_identical(fitted$id, "none")
  expect_identical(fitted$selected_features, names(x))
  expect_identical(any(vapply(fitted, is.data.frame, logical(1))), FALSE)
  expect_identical(any(vapply(fitted, is.matrix, logical(1))), FALSE)
  expect_identical(selector$apply(fitted, x), x)
})

test_that("the none selector enforces its locked predictor schema", {
  x <- data.frame(gene_a = 1:3, gene_b = 4:6)
  selector <- sb_none_selector_adapter()
  fitted <- selector$fit(x)

  reordered_error <- capture_sb_error(
    selector$apply(fitted, x[c("gene_b", "gene_a")])
  )
  expect_s3_class(
    reordered_error,
    "stormbreaker_selector_schema_mismatch"
  )

  missing_error <- capture_sb_error(
    selector$apply(fitted, x["gene_a"])
  )
  expect_s3_class(missing_error, "stormbreaker_selector_schema_mismatch")

  extra <- x
  extra$gene_c <- 7:9
  extra_error <- capture_sb_error(selector$apply(fitted, extra))
  expect_s3_class(extra_error, "stormbreaker_selector_schema_mismatch")
})

test_that("the Cox adapter audits coefficients and reconstructs risk", {
  d <- make_sb_fixture(n = 60)
  analysis_rows <- seq_len(40)
  assessment_rows <- 41:60
  x_analysis <- d[
    analysis_rows,
    c("gene_b", "gene_a"),
    drop = FALSE
  ]
  x_assessment <- d[
    assessment_rows,
    c("gene_b", "gene_a"),
    drop = FALSE
  ]
  learner <- sb_coxph_learner_adapter()
  model <- learner$fit(
    x_analysis,
    d$time[analysis_rows],
    d$event[analysis_rows]
  )
  observed <- learner$predict(model, x_assessment)
  coefficients <- learner$extract_coefficients(model)

  manual <- survival::coxph(
    survival::Surv(time, event) ~ gene_b + gene_a,
    data = d[analysis_rows, , drop = FALSE],
    ties = "efron",
    singular.ok = FALSE,
    x = TRUE,
    y = TRUE,
    model = TRUE
  )
  expected <- unname(stats::predict(
    manual,
    newdata = d[assessment_rows, , drop = FALSE],
    type = "lp",
    reference = "zero"
  ))
  reconstructed <- unname(as.vector(
    as.matrix(x_assessment) %*% coefficients$coefficient
  ))

  expect_named(
    coefficients,
    c("feature", "coefficient", "nonzero", "direction"),
    ignore.order = FALSE
  )
  expect_identical(coefficients$feature, c("gene_b", "gene_a"))
  expect_equal(
    coefficients$coefficient,
    unname(stats::coef(model)[coefficients$feature]),
    tolerance = 0
  )
  expect_identical(coefficients$nonzero, coefficients$coefficient != 0)
  expect_identical(
    coefficients$direction,
    ifelse(
      coefficients$coefficient > 0,
      "higher_hazard",
      ifelse(coefficients$coefficient < 0, "lower_hazard", "zero")
    )
  )
  expect_equal(observed, expected, tolerance = 1e-12)
  expect_equal(observed, reconstructed, tolerance = 1e-12)

  schema_error <- capture_sb_error(
    learner$predict(model, x_assessment[c("gene_a", "gene_b")])
  )
  expect_s3_class(schema_error, "stormbreaker_learner_schema_mismatch")
})

test_that("ordinary Cox survival probabilities retain Efron baseline semantics", {
  d <- make_sb_fixture(n = 60L, seed = 20260819L)
  analysis_rows <- seq_len(40L)
  assessment_rows <- 41:60
  predictors <- c("gene_b", "gene_a")
  x_analysis <- d[analysis_rows, predictors, drop = FALSE]
  x_assessment <- d[assessment_rows, predictors, drop = FALSE]
  horizon <- 12
  learner <- sb_coxph_learner_adapter()
  model <- learner$fit(
    x_analysis,
    d$time[analysis_rows],
    d$event[analysis_rows]
  )

  observed <- learner$predict_survival(
    model,
    x_analysis,
    d$time[analysis_rows],
    d$event[analysis_rows],
    x_assessment,
    horizon
  )
  basehaz <- survival::basehaz(model, centered = FALSE)
  expected_hazard <- tail(basehaz$hazard[basehaz$time <= horizon], 1L)
  expected_survival <- drop(summary(
    survival::survfit(model, newdata = x_assessment, se.fit = FALSE),
    times = horizon,
    extend = TRUE
  )$surv)
  risk <- learner$predict(model, x_assessment)
  reconstructed <- exp(-expected_hazard * exp(risk))

  expect_named(
    observed,
    c(
      "prediction_horizon", "cumulative_baseline_hazard",
      "baseline_survival", "source", "ties", "probability"
    ),
    ignore.order = FALSE
  )
  expect_identical(observed$prediction_horizon, horizon)
  expect_identical(observed$source, "coxph_basehaz_zero_reference")
  expect_identical(observed$ties, "efron")
  expect_equal(
    observed$cumulative_baseline_hazard,
    expected_hazard,
    tolerance = 1e-14
  )
  expect_equal(
    observed$baseline_survival,
    exp(-expected_hazard),
    tolerance = 1e-14
  )
  expect_equal(
    observed$probability,
    unname(expected_survival),
    tolerance = 1e-12
  )
  expect_equal(observed$probability, reconstructed, tolerance = 1e-12)
})

test_that("Cox coefficient audits preserve non-syntactic feature names", {
  d <- make_sb_fixture(n = 40L, seed = 20260819L)
  x_analysis <- d[1:30, c("gene_a", "gene_b"), drop = FALSE]
  x_assessment <- d[31:40, c("gene_a", "gene_b"), drop = FALSE]
  names(x_analysis) <- c("gene a", "gene-b")
  names(x_assessment) <- names(x_analysis)
  learner <- sb_coxph_learner_adapter()

  model <- learner$fit(
    x_analysis,
    d$time[1:30],
    d$event[1:30]
  )
  coefficients <- learner$extract_coefficients(model)
  observed <- learner$predict(model, x_assessment)
  survival_prediction <- learner$predict_survival(
    model,
    x_analysis,
    d$time[1:30],
    d$event[1:30],
    x_assessment,
    prediction_horizon = 12
  )
  reconstructed <- unname(as.vector(
    as.matrix(x_assessment[coefficients$feature]) %*%
      coefficients$coefficient
  ))

  expect_identical(coefficients$feature, c("gene a", "gene-b"))
  expect_true(all(is.finite(coefficients$coefficient)))
  expect_equal(reconstructed, observed, tolerance = 1e-12)
  expect_length(survival_prediction$probability, nrow(x_assessment))
  expect_true(all(is.finite(survival_prediction$probability)))
})

test_that("the unpenalized Cox adapter rejects p greater than or equal to n", {
  x <- data.frame(
    gene_a = c(0.1, 0.4, 0.2, 0.8),
    gene_b = c(1.2, 0.9, 1.5, 1.1),
    gene_c = c(2.1, 2.4, 2.3, 2.8),
    gene_d = c(3.2, 3.4, 3.1, 3.9)
  )

  error <- capture_sb_error(
    sb_fit_coxph(
      x,
      time = c(1, 2, 3, 4),
      event = c(1, 0, 1, 0)
    )
  )

  expect_s3_class(error, "stormbreaker_high_dimensional_cox")
  expect_match(conditionMessage(error), "4 predictors with 4 patients")
  expect_match(conditionMessage(error), "fold-local selector or glmnet Cox")
})

test_that("the fixed glmnet Cox adapter matches a direct path fit", {
  skip_if_not_installed("glmnet")

  d <- make_sb_fixture(n = 60)
  analysis_rows <- seq_len(40)
  assessment_rows <- 41:60
  x_analysis <- d[
    analysis_rows,
    c("gene_a", "gene_b", "gene_c"),
    drop = FALSE
  ]
  x_assessment <- d[
    assessment_rows,
    c("gene_a", "gene_b", "gene_c"),
    drop = FALSE
  ]
  learner <- sb_glmnet_cox_learner_adapter(
    alpha = 1,
    lambda = 0.05,
    standardize = TRUE
  )
  model <- learner$fit(
    x_analysis,
    d$time[analysis_rows],
    d$event[analysis_rows]
  )
  observed <- learner$predict(model, x_assessment)

  manual <- glmnet::glmnet(
    x = as.matrix(x_analysis),
    y = survival::Surv(
      d$time[analysis_rows],
      d$event[analysis_rows]
    ),
    family = "cox",
    alpha = 1,
    standardize = TRUE
  )
  expected <- unname(drop(stats::predict(
    manual,
    newx = as.matrix(x_assessment),
    s = 0.05,
    type = "link"
  )))

  expect_identical(learner$id, "glmnet_cox")
  expect_equal(observed, expected, tolerance = 1e-12)
  expect_length(observed, nrow(x_assessment))
  expect_identical(all(is.finite(observed)), TRUE)

  schema_error <- capture_sb_error(
    learner$predict(model, x_assessment[c("gene_b", "gene_a", "gene_c")])
  )
  expect_s3_class(
    schema_error,
    "stormbreaker_learner_schema_mismatch"
  )
})

test_that("glmnet Cox survival probabilities use fixed-beta Breslow ties", {
  skip_if_not_installed("glmnet")

  d <- make_sb_fixture(n = 60L, seed = 20260819L)
  d$time <- rep(seq_len(10L), each = 6L)
  d$event <- rep(c(1L, 0L, 1L, 0L, 1L, 0L), 10L)
  analysis_rows <- seq_len(45L)
  assessment_rows <- 46:60
  predictors <- c("gene_a", "gene_b", "gene_c")
  x_analysis <- d[analysis_rows, predictors, drop = FALSE]
  x_assessment <- d[assessment_rows, predictors, drop = FALSE]
  time <- d$time[analysis_rows]
  event <- d$event[analysis_rows]
  horizon <- 6
  learner <- sb_glmnet_cox_learner_adapter(
    alpha = 1,
    lambda = 0.05,
    standardize = TRUE
  )
  model <- learner$fit(x_analysis, time, event)

  observed <- learner$predict_survival(
    model,
    x_analysis,
    time,
    event,
    x_assessment,
    horizon
  )
  analysis_risk <- learner$predict(model, x_analysis)
  oracle_data <- data.frame(time = time, event = event, eta = analysis_risk)
  oracle <- survival::coxph(
    survival::Surv(time, event) ~ eta,
    data = oracle_data,
    ties = "breslow",
    init = 1,
    iter = 0,
    singular.ok = FALSE
  )
  oracle_basehaz <- survival::basehaz(oracle, centered = FALSE)
  expected_hazard <- tail(
    oracle_basehaz$hazard[oracle_basehaz$time <= horizon],
    1L
  )
  assessment_risk <- learner$predict(model, x_assessment)
  expected_probability <- exp(-expected_hazard * exp(assessment_risk))
  efron_oracle <- survival::coxph(
    survival::Surv(time, event) ~ eta,
    data = oracle_data,
    ties = "efron",
    init = 1,
    iter = 0,
    singular.ok = FALSE
  )
  efron_basehaz <- survival::basehaz(efron_oracle, centered = FALSE)
  efron_hazard <- tail(
    efron_basehaz$hazard[efron_basehaz$time <= horizon],
    1L
  )

  expect_identical(observed$source, "glmnet_fixed_beta_breslow")
  expect_identical(observed$ties, "breslow")
  expect_equal(
    observed$cumulative_baseline_hazard,
    expected_hazard,
    tolerance = 1e-12
  )
  expect_equal(
    observed$probability,
    expected_probability,
    tolerance = 1e-12
  )
  expect_gt(abs(expected_hazard - efron_hazard), 1e-8)
  expect_equal(
    learner$extract_coefficients(model)$coefficient,
    unname(as.matrix(stats::coef(model, s = 0.05))[, 1L]),
    tolerance = 0
  )
})

test_that("glmnet Cox reports every fitted coefficient at its lambda", {
  skip_if_not_installed("glmnet")

  d <- make_sb_fixture(n = 60)
  predictors <- c("gene_a", "gene_b", "gene_c")
  learner <- sb_glmnet_cox_learner_adapter(
    alpha = 1,
    lambda = 0.05,
    standardize = TRUE
  )
  model <- learner$fit(
    d[predictors],
    d$time,
    d$event
  )
  observed <- learner$extract_coefficients(model)
  expected <- drop(as.matrix(stats::coef(model, s = 0.05)))

  expect_named(
    observed,
    c("feature", "coefficient", "nonzero", "direction"),
    ignore.order = FALSE
  )
  expect_identical(observed$feature, predictors)
  expect_equal(observed$coefficient, unname(expected), tolerance = 0)
  expect_identical(observed$nonzero, unname(expected != 0))
  expect_identical(
    observed$direction,
    unname(ifelse(expected > 0, "higher_hazard", ifelse(
      expected < 0,
      "lower_hazard",
      "zero"
    )))
  )
  expect_identical(any(observed$nonzero), TRUE)
  expect_identical(any(!observed$nonzero), TRUE)
})

test_that("the fixed glmnet Cox adapter rejects invalid settings", {
  alpha_error <- capture_sb_error(
    sb_glmnet_cox_learner_adapter(alpha = 1.2, lambda = 0.05)
  )
  expect_s3_class(alpha_error, "stormbreaker_bad_glmnet_alpha")

  lambda_error <- capture_sb_error(
    sb_glmnet_cox_learner_adapter(alpha = 1, lambda = 0)
  )
  expect_s3_class(lambda_error, "stormbreaker_bad_glmnet_lambda")

  standardize_error <- capture_sb_error(
    sb_glmnet_cox_learner_adapter(
      alpha = 1,
      lambda = 0.05,
      standardize = NA
    )
  )
  expect_s3_class(
    standardize_error,
    "stormbreaker_bad_glmnet_standardize"
  )
})

test_that("the fixed glmnet Cox adapter reports a missing dependency", {
  local_mocked_bindings(
    sb_glmnet_available = function() FALSE
  )
  learner <- sb_glmnet_cox_learner_adapter(
    alpha = 1,
    lambda = 0.05,
    standardize = TRUE
  )
  d <- make_sb_fixture(n = 20)

  error <- capture_sb_error(learner$fit(
    d[c("gene_a", "gene_b")],
    d$time,
    d$event
  ))

  expect_s3_class(error, "stormbreaker_missing_glmnet")
})

test_that("the fixed glmnet Cox adapter rejects non-positive time", {
  skip_if_not_installed("glmnet")

  d <- make_sb_fixture(n = 60)
  d$time[1] <- 0
  learner <- sb_glmnet_cox_learner_adapter(
    alpha = 1,
    lambda = 0.05,
    standardize = TRUE
  )

  error <- capture_sb_error(learner$fit(
    d[c("gene_a", "gene_b", "gene_c")],
    d$time,
    d$event
  ))

  expect_s3_class(error, "stormbreaker_glmnet_nonpositive_time")
})

test_that("glmnet Cox refuses lambda values outside its fitted path", {
  skip_if_not_installed("glmnet")

  x <- data.frame(
    g1 = c(
      -1.941518376040617833, 1.181180891438521785,
      1.859910861543160987, 1.074012260281989928,
      -0.027346966413741798
    ),
    g2 = c(
      -0.033330342275919272, -1.516067623151099220,
      0.790385343888228009, -0.210734180691775658,
      -0.656742925361075369
    )
  )
  time <- c(
    2.34298038283798737, 0.62733672885224223,
    1.22186382110686331, 0.60436541587114334,
    1.96968345521929300
  )
  event <- c(1, 1, 1, 0, 0)
  learner <- sb_glmnet_cox_learner_adapter(
    alpha = 1,
    lambda = 4.0337729269749e-12,
    standardize = TRUE
  )

  error <- capture_sb_error(learner$fit(x, time, event))

  expect_s3_class(error, "stormbreaker_glmnet_lambda_outside_path")
})

test_that("assessment endpoints cannot change fixed-split adapter risks", {
  d <- make_sb_fixture(n = 60)
  task <- sb_surv_task(
    d, "patient_id", "time", "event", c("gene_a", "gene_b"),
    "month"
  )
  split <- list(
    fold = 1L,
    analysis_rows = seq_len(48),
    assessment_rows = 49:60,
    analysis_ids = d$patient_id[seq_len(48)],
    assessment_ids = d$patient_id[49:60]
  )
  runtime <- sb_resolve_pipeline_runtime(sb_pipeline())
  original <- sb_run_one_fold(
    task, sb_pipeline(), runtime, split, seed = 17L
  )

  changed <- d
  changed$time[split$assessment_rows] <-
    changed$time[split$assessment_rows] + 1000
  changed$event[split$assessment_rows] <-
    1L - changed$event[split$assessment_rows]
  changed_task <- sb_surv_task(
    changed, "patient_id", "time", "event", c("gene_a", "gene_b"),
    "month"
  )
  perturbed <- sb_run_one_fold(
    changed_task, sb_pipeline(), runtime, split, seed = 17L
  )

  expect_equal(
    original$predictions$risk,
    perturbed$predictions$risk,
    tolerance = 0
  )
  expect_identical(original$coefficients$status, "success")
  expect_identical(original$coefficients, perturbed$coefficients)
})

test_that("each adapter operation receives only its permitted fold data", {
  d <- make_sb_fixture(n = 20)
  task <- sb_surv_task(
    d, "patient_id", "time", "event", c("gene_a", "gene_b"),
    "month"
  )
  split <- list(
    fold = 1L,
    analysis_rows = 1:12,
    assessment_rows = 13:20,
    analysis_ids = d$patient_id[1:12],
    assessment_ids = d$patient_id[13:20]
  )
  received <- new.env(parent = emptyenv())
  runtime <- structure(list(
    id = "none__coxph",
    selector = list(
      id = "none",
      params = list(),
      needs_outcome = FALSE,
      fit = function(x) {
        received$selector_fit <- x
        list(features = names(x))
      },
      apply = function(fitted, x) {
        if (is.null(received$selector_apply)) {
          received$selector_apply <- list()
        }
        received$selector_apply[[length(received$selector_apply) + 1L]] <- x
        x[fitted$features]
      }
    ),
    learner = list(
      id = "coxph",
      fit = function(x, time, event) {
        received$learner_fit <- list(x = x, time = time, event = event)
        list(features = names(x))
      },
      predict = function(model, x) {
        received$learner_predict <- x
        rep.int(0, nrow(x))
      },
      predict_survival = function(model, x_analysis, time, event,
                                  x_assessment, prediction_horizon) {
        stop("not requested by this fixture")
      },
      extract_coefficients = function(model) {
        data.frame(
          feature = model$features,
          coefficient = rep(0, length(model$features)),
          nonzero = rep(FALSE, length(model$features)),
          direction = rep("zero", length(model$features)),
          stringsAsFactors = FALSE
        )
      },
      extract_state = function(model, context = list()) {
        list(features = model$features)
      },
      predict_state = function(state, x) rep.int(0, nrow(x))
    )
  ), class = "stormbreaker_pipeline_runtime")

  result <- sb_run_one_fold(
    task, sb_pipeline(), runtime, split, seed = 23L
  )

  expected_analysis <- d[1:12, c("gene_a", "gene_b"), drop = FALSE]
  expected_assessment <- d[13:20, c("gene_a", "gene_b"), drop = FALSE]
  expect_identical(received$selector_fit, expected_analysis)
  expect_identical(received$selector_apply[[1]], expected_analysis)
  expect_identical(received$selector_apply[[2]], expected_assessment)
  expect_identical(received$learner_fit$x, expected_analysis)
  expect_identical(received$learner_fit$time, d$time[1:12])
  expect_identical(received$learner_fit$event, d$event[1:12])
  expect_identical(received$learner_predict, expected_assessment)
  expect_identical(result$run$status, "success")
})

test_that("an adapter failure is recorded without invoking a fallback", {
  d <- make_sb_fixture(n = 20)
  task <- sb_surv_task(
    d, "patient_id", "time", "event", "gene_a", "month"
  )
  split <- list(
    fold = 1L,
    analysis_rows = seq_len(10),
    assessment_rows = 11:20,
    analysis_ids = d$patient_id[seq_len(10)],
    assessment_ids = d$patient_id[11:20]
  )
  learner_fit_calls <- 0L
  runtime <- structure(list(
    id = "none__coxph",
    selector = list(
      id = "none",
      params = list(),
      needs_outcome = FALSE,
      fit = function(x) {
        sb_abort(
          "Synthetic selector failure.",
          "stormbreaker_test_selector_failure"
        )
      },
      apply = function(fitted, x) x
    ),
    learner = list(
      id = "coxph",
      fit = function(x, time, event) {
        learner_fit_calls <<- learner_fit_calls + 1L
        NULL
      },
      predict = function(model, x) rep.int(0, nrow(x)),
      predict_survival = function(model, x_analysis, time, event,
                                  x_assessment, prediction_horizon) {
        stop("not requested by this fixture")
      },
      extract_coefficients = function(model) NULL,
      extract_state = function(model, context = list()) list(),
      predict_state = function(state, x) rep.int(0, nrow(x))
    )
  ), class = "stormbreaker_pipeline_runtime")

  result <- sb_run_one_fold(
    task, sb_pipeline(), runtime, split, seed = 19L
  )

  expect_identical(result$run$status, "failed")
  expect_identical(
    result$run$error_class,
    "stormbreaker_test_selector_failure"
  )
  expect_identical(result$run$learner_used, NA_character_)
  expect_identical(learner_fit_calls, 0L)
  expect_identical(all(is.na(result$predictions$risk)), TRUE)
})

test_that("a prediction failure retains an already fitted Cox audit", {
  d <- make_sb_fixture(n = 30L, seed = 20260819L)
  task <- sb_surv_task(
    d,
    "patient_id",
    "time",
    "event",
    c("gene_a", "gene_b"),
    "month"
  )
  split <- list(
    fold = 1L,
    analysis_rows = 1:20,
    assessment_rows = 21:30,
    analysis_ids = d$patient_id[1:20],
    assessment_ids = d$patient_id[21:30]
  )
  pipeline <- sb_pipeline()
  runtime <- sb_resolve_pipeline_runtime(pipeline)
  runtime$learner$predict <- function(model, x) {
    sb_abort(
      "Synthetic Cox prediction failure.",
      "stormbreaker_test_cox_prediction_failure"
    )
  }

  result <- sb_run_one_fold(task, pipeline, runtime, split, seed = 1L)

  expect_identical(result$run$status, "failed")
  expect_identical(
    result$run$error_class,
    "stormbreaker_test_cox_prediction_failure"
  )
  expect_true(all(is.na(result$predictions$risk)))
  expect_identical(result$coefficients$status, "success")
  expect_identical(result$coefficients$values$feature, task$predictors)
  expect_true(all(is.finite(result$coefficients$values$coefficient)))
})

test_that("univariable Cox Top-K matches independent manual scores", {
  d <- make_sb_fixture(n = 60)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  selector <- sb_univ_cox_top_k_selector_adapter(n_features = 2L)
  fitted <- selector$fit(x, d$time, d$event)

  manual_p <- vapply(names(x), function(feature) {
    fit <- survival::coxph(
      stats::reformulate(
        feature,
        response = "survival::Surv(time, event)"
      ),
      data = d,
      ties = "efron",
      singular.ok = FALSE
    )
    unname(summary(fit)$coefficients[1, "Pr(>|z|)"])
  }, numeric(1))
  expected <- names(x)[order(manual_p, seq_along(manual_p))][1:2]

  expect_s3_class(fitted, "stormbreaker_selector_fit")
  expect_identical(fitted$id, "univ_cox_top_k")
  expect_identical(fitted$selected_features, expected)
  expect_named(
    fitted$scores,
    c(
      "feature", "original_order", "coefficient", "standard_error",
      "wald_z", "p_value", "status", "error_class", "error_message",
      "rank", "selected"
    ),
    ignore.order = FALSE
  )
  expect_equal(
    fitted$scores$p_value,
    unname(manual_p[match(fitted$scores$feature, names(manual_p))]),
    tolerance = 1e-12
  )
  expect_identical(sum(fitted$scores$selected), 2L)
  expect_identical(selector$apply(fitted, x), x[expected])
})

test_that("Top-K breaks exact score ties by original predictor order", {
  d <- make_sb_fixture(n = 60)
  x <- data.frame(
    first = d$gene_a,
    second = d$gene_a,
    third = d$gene_b
  )
  selector <- sb_univ_cox_top_k_selector_adapter(n_features = 1L)
  fitted <- selector$fit(x, d$time, d$event)

  expect_identical(fitted$selected_features, "first")
})

test_that("Top-K records invalid features and fails when K is infeasible", {
  d <- make_sb_fixture(n = 60)
  x <- data.frame(valid = d$gene_a, constant = 1)
  selector <- sb_univ_cox_top_k_selector_adapter(n_features = 1L)
  fitted <- selector$fit(x, d$time, d$event)

  expect_identical(fitted$selected_features, "valid")
  expect_identical(
    fitted$scores$status[fitted$scores$feature == "constant"],
    "failed"
  )
  expect_identical(
    is.na(fitted$scores$p_value[fitted$scores$feature == "constant"]),
    TRUE
  )

  infeasible <- sb_univ_cox_top_k_selector_adapter(n_features = 2L)
  error <- capture_sb_error(infeasible$fit(x, d$time, d$event))
  expect_s3_class(error, "stormbreaker_insufficient_valid_features")
})

test_that("Top-K does not let a single patient define a feature score", {
  d <- make_sb_fixture(n = 60)
  x <- data.frame(
    stable = d$gene_a,
    singleton = c(1, rep(0, 59))
  )
  selector <- sb_univ_cox_top_k_selector_adapter(n_features = 1L)

  fitted <- selector$fit(x, d$time, d$event)
  singleton_score <- fitted$scores[
    fitted$scores$feature == "singleton",
    ,
    drop = FALSE
  ]

  expect_identical(fitted$selected_features, "stable")
  expect_identical(singleton_score$status, "failed")
  expect_identical(
    singleton_score$error_class,
    "stormbreaker_selector_insufficient_feature_support"
  )
  expect_match(singleton_score$error_message, "single patient")
})

test_that("Top-K rejects K larger than the candidate feature count", {
  d <- make_sb_fixture(n = 20)
  x <- d[c("gene_a", "gene_b")]
  selector <- sb_univ_cox_top_k_selector_adapter(n_features = 3L)

  error <- capture_sb_error(selector$fit(x, d$time, d$event))

  expect_s3_class(error, "stormbreaker_insufficient_candidate_features")
})

test_that("a Top-K selector failure never reaches learner fitting", {
  d <- make_sb_fixture(n = 20)
  d$constant_a <- 1
  d$constant_b <- 1
  task <- sb_surv_task(
    d, "patient_id", "time", "event",
    c("constant_a", "constant_b"), "month"
  )
  pipeline <- sb_pipeline(
    selector = "univ_cox_top_k",
    selector_n_features = 1L
  )
  runtime <- sb_resolve_pipeline_runtime(pipeline)
  learner_fit_calls <- 0L
  runtime$learner$fit <- function(x, time, event) {
    learner_fit_calls <<- learner_fit_calls + 1L
    NULL
  }
  split <- list(
    fold = 1L,
    analysis_rows = 1:12,
    assessment_rows = 13:20,
    analysis_ids = d$patient_id[1:12],
    assessment_ids = d$patient_id[13:20]
  )

  result <- sb_run_one_fold(task, pipeline, runtime, split, seed = 7L)

  expect_identical(result$run$status, "failed")
  expect_identical(
    result$run$error_class,
    "stormbreaker_insufficient_valid_features"
  )
  expect_identical(learner_fit_calls, 0L)
  expect_identical(all(is.na(result$predictions$risk)), TRUE)
})

test_that("Top-K apply enforces the fitted source schema", {
  d <- make_sb_fixture(n = 60)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  selector <- sb_univ_cox_top_k_selector_adapter(n_features = 2L)
  fitted <- selector$fit(x, d$time, d$event)

  reordered_error <- capture_sb_error(
    selector$apply(fitted, x[c("gene_b", "gene_a", "gene_c")])
  )
  expect_s3_class(
    reordered_error,
    "stormbreaker_selector_schema_mismatch"
  )
})

test_that("Top-K receives only analysis endpoints in fold execution", {
  d <- make_sb_fixture(n = 60)
  task <- sb_surv_task(
    d, "patient_id", "time", "event",
    c("gene_a", "gene_b", "gene_c"), "month"
  )
  split <- list(
    fold = 1L,
    analysis_rows = 1:48,
    assessment_rows = 49:60,
    analysis_ids = d$patient_id[1:48],
    assessment_ids = d$patient_id[49:60]
  )
  pipeline <- sb_pipeline(
    selector = "univ_cox_top_k",
    selector_n_features = 2L
  )
  runtime <- sb_resolve_pipeline_runtime(pipeline)
  captured <- new.env(parent = emptyenv())
  original_fit <- runtime$selector$fit
  runtime$selector$fit <- function(x, time, event) {
    captured$x <- x
    captured$time <- time
    captured$event <- event
    original_fit(x, time, event)
  }

  result <- sb_run_one_fold(task, pipeline, runtime, split, seed = 17L)

  expect_identical(captured$x, d[1:48, task$predictors, drop = FALSE])
  expect_identical(captured$time, d$time[1:48])
  expect_identical(captured$event, d$event[1:48])
  expect_identical(result$run$status, "success")
  expect_length(result$selection$selected_features, 2L)
})

test_that("assessment endpoints cannot change fixed-split Top-K selection", {
  d <- make_sb_fixture(n = 60)
  pipeline <- sb_pipeline(
    selector = "univ_cox_top_k",
    selector_n_features = 2L
  )
  split <- list(
    fold = 1L,
    analysis_rows = 1:48,
    assessment_rows = 49:60,
    analysis_ids = d$patient_id[1:48],
    assessment_ids = d$patient_id[49:60]
  )
  make_task <- function(data) {
    sb_surv_task(
      data, "patient_id", "time", "event",
      c("gene_a", "gene_b", "gene_c"), "month"
    )
  }
  original <- sb_run_one_fold(
    make_task(d), pipeline, sb_resolve_pipeline_runtime(pipeline), split,
    seed = 17L
  )
  changed <- d
  changed$time[49:60] <- changed$time[49:60] + 1000
  changed$event[49:60] <- 1L - changed$event[49:60]
  perturbed <- sb_run_one_fold(
    make_task(changed), pipeline, sb_resolve_pipeline_runtime(pipeline),
    split, seed = 17L
  )

  expect_identical(
    original$selection$selected_features,
    perturbed$selection$selected_features
  )
  expect_identical(original$selection$scores, perturbed$selection$scores)
  expect_equal(
    original$predictions$risk,
    perturbed$predictions$risk,
    tolerance = 0
  )
})
