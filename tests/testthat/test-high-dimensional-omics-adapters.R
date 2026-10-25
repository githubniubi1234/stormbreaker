make_high_dimensional_omics_fixture <- function(n = 72L, p = 48L,
                                                seed = 9201L) {
  withr::with_seed(seed, {
    outcome <- sample(rep(c(0L, 1L), length.out = n))
    x <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
    x[, 1L] <- x[, 1L] + 1.4 * outcome
    x[, 2L] <- x[, 2L] - outcome
    x[, 3L] <- x[, 3L] + 0.7 * outcome
    x <- as.data.frame(x, optional = TRUE)
    names(x) <- paste0("gene_", seq_len(p))
    list(x = x, outcome = outcome)
  })
}

make_high_dimensional_survival_fixture <- function(n = 84L, p = 48L,
                                                   seed = 9203L) {
  withr::with_seed(seed, {
    x <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
    event_time <- stats::rexp(
      n,
      rate = exp(1.1 * x[, 1L] - 0.9 * x[, 2L] + 0.6 * x[, 3L])
    )
    censor_time <- stats::rexp(n, rate = 0.3)
    x <- as.data.frame(x, optional = TRUE)
    names(x) <- paste0("gene_", seq_len(p))
    list(
      x = x,
      time = pmin(event_time, censor_time) + 0.01,
      event = as.integer(event_time <= censor_time)
    )
  })
}

test_that("SIS screening uses its fold-local native marginal ranking", {
  skip_if_not_installed("SIS", minimum_version = "1.5")

  fixture <- make_high_dimensional_omics_fixture(p = 80L)
  params <- sb_sis_binary_selector_params()
  expected_size <- sb_sis_binary_screen_size(
    nrow(fixture$x),
    ncol(fixture$x)
  )
  native <- withr::with_seed(params$seed, SIS::SIS(
    x = data.matrix(fixture$x),
    y = fixture$outcome,
    family = params$family,
    penalty = params$penalty,
    tune = params$tuning,
    type.measure = "deviance",
    nsis = expected_size,
    iter = FALSE,
    varISIS = "vanilla",
    perm = FALSE,
    greedy = FALSE,
    seed = params$seed,
    standardize = params$standardize,
    boot_ci = FALSE,
    parallel = FALSE
  ))

  set.seed(9202L)
  caller_rng <- .Random.seed
  fitted <- sb_fit_sis_binary_selector(fixture$x, fixture$outcome)

  expect_identical(.Random.seed, caller_rng)
  expect_identical(
    fitted$selected_features,
    names(fixture$x)[native$sis.ix0]
  )
  expect_identical(length(fitted$selected_features), expected_size)
  expect_identical(
    names(sb_apply_high_dimensional_selector(fitted, fixture$x)),
    fitted$selected_features
  )
})

test_that("VSURF exposes its native prediction-stage subset", {
  skip_if_not_installed("VSURF", minimum_version = "1.2.1")

  fixture <- make_high_dimensional_omics_fixture(p = 30L, seed = 9211L)
  params <- sb_vsurf_binary_selector_params()
  mtry <- max(floor(ncol(fixture$x) / 3), 1L)
  native <- withr::with_seed(params$seed, VSURF::VSURF(
    x = data.matrix(fixture$x),
    y = factor(fixture$outcome, levels = c(0L, 1L)),
    mtry = mtry,
    ntree.thres = params$threshold_trees,
    nfor.thres = params$threshold_forests,
    nmin = params$nmin,
    ntree.interp = params$interpretation_trees,
    nfor.interp = params$interpretation_forests,
    nsd = params$nsd,
    ntree.pred = params$prediction_trees,
    nfor.pred = params$prediction_forests,
    nmj = params$nmj,
    RFimplem = params$random_forest_engine,
    parallel = FALSE,
    verbose = FALSE
  ))

  set.seed(9212L)
  caller_rng <- .Random.seed
  fitted <- sb_fit_vsurf_binary_selector(fixture$x, fixture$outcome)
  diagnostic <- sb_selector_diagnostics(fitted, "vsurf_importance")
  expect_equal(diagnostic$table$importance[native$imp.mean.dec.ind], as.numeric(native$imp.mean.dec), tolerance = 0)
  expect_equal(diagnostic$table$importance_sd[native$imp.mean.dec.ind], as.numeric(native$imp.sd.dec), tolerance = 0)

  expect_identical(.Random.seed, caller_rng)
  expect_identical(
    fitted$selected_features,
    names(fixture$x)[native$varselect.pred]
  )
  expect_identical(
    fitted$selection_history$n_features,
    as.integer(c(
      length(native$varselect.thres),
      length(native$varselect.interp),
      length(native$varselect.pred)
    ))
  )
})

test_that("mRMR preserves its native ordered high-dimensional subset", {
  skip_if_not_installed("praznik", minimum_version = "13.0.0")

  fixture <- make_high_dimensional_omics_fixture(p = 90L, seed = 9213L)
  params <- sb_mrmr_binary_selector_params()
  expected_size <- sb_high_dimensional_nlogn_size(
    nrow(fixture$x),
    ncol(fixture$x)
  )
  native <- praznik::MRMR(
    X = fixture$x,
    Y = factor(fixture$outcome, levels = c(0L, 1L)),
    k = expected_size,
    positive = params$positive_scores_only,
    threads = params$threads
  )

  set.seed(9214L)
  caller_rng <- .Random.seed
  fitted <- sb_fit_mrmr_binary_selector(fixture$x, fixture$outcome)

  expect_identical(.Random.seed, caller_rng)
  expect_identical(
    fitted$selected_features,
    names(fixture$x)[as.integer(native$selection)]
  )
  expect_equal(
    fitted$scores$mrmr_score[
      order(fitted$scores$rank, na.last = NA)
    ],
    as.numeric(native$score),
    tolerance = 0
  )
})

test_that("information filters preserve their native ordered subsets", {
  skip_if_not_installed("praznik", minimum_version = "13.0.0")

  fixture <- make_high_dimensional_omics_fixture(p = 90L, seed = 9215L)
  expected_size <- sb_high_dimensional_nlogn_size(
    nrow(fixture$x),
    ncol(fixture$x)
  )
  native_methods <- c(
    jmi_binary_nlogn = "JMI",
    cmim_binary_nlogn = "CMIM",
    disr_binary_nlogn = "DISR",
    jmim_binary_nlogn = "JMIM"
  )

  for (selector in names(native_methods)) {
    params <- sb_information_binary_selector_params(selector)
    native <- do.call(
      getExportedValue("praznik", native_methods[[selector]]),
      list(
        X = fixture$x,
        Y = factor(fixture$outcome, levels = c(0L, 1L)),
        k = expected_size,
        threads = params$threads
      )
    )

    set.seed(9216L)
    caller_rng <- .Random.seed
    fitted <- sb_fit_information_binary_selector(
      fixture$x,
      fixture$outcome,
      selector
    )

    expect_identical(.Random.seed, caller_rng)
    expect_identical(
      fitted$selected_features,
      names(fixture$x)[as.integer(native$selection)]
    )
    expect_equal(
      fitted$scores$information_score[
        order(fitted$scores$rank, na.last = NA)
      ],
      as.numeric(native$score),
      tolerance = 0
    )
    expect_identical(fitted$params$native_method, native_methods[[selector]])
  }
})

test_that("LPC selectors preserve native gene scores and absolute ranking", {
  skip_if_not_installed("lpc", minimum_version = "1.0.2.1")

  binary <- make_high_dimensional_omics_fixture(
    n = 80L,
    p = 60L,
    seed = 9233L
  )
  survival <- make_high_dimensional_survival_fixture(
    n = 84L,
    p = 60L,
    seed = 9234L
  )
  specifications <- list(
    list(
      selector = "lpc_binary_nlogn",
      x = binary$x,
      response = binary$outcome + 1L,
      censoring_status = NULL,
      fit = function() {
        sb_fit_lpc_binary_selector(binary$x, binary$outcome)
      }
    ),
    list(
      selector = "lpc_cox_nlogn",
      x = survival$x,
      response = survival$time,
      censoring_status = survival$event,
      fit = function() {
        sb_fit_lpc_cox_selector(
          survival$x,
          survival$time,
          survival$event
        )
      }
    )
  )

  for (specification in specifications) {
    params <- sb_lpc_selector_params(specification$selector)
    native <- NULL
    invisible(utils::capture.output(
      native <- withr::with_seed(
        params$seed,
        lpc::LPC(
          x = t(data.matrix(specification$x)),
          y = specification$response,
          type = params$native_type,
          soft.thresh = NULL,
          censoring.status = specification$censoring_status
        )
      )
    ))
    native_score <- as.numeric(native$lpcscores)
    expected_size <- sb_high_dimensional_nlogn_size(
      nrow(specification$x),
      ncol(specification$x)
    )
    expected_order <- order(
      -abs(native_score),
      seq_along(native_score),
      method = "radix"
    )

    set.seed(9235L)
    caller_rng <- .Random.seed
    fitted <- specification$fit()

    expect_identical(.Random.seed, caller_rng)
    expect_identical(
      fitted$selected_features,
      names(specification$x)[expected_order[seq_len(expected_size)]]
    )
    expect_equal(fitted$scores$lpc_score, native_score, tolerance = 0)
    expect_equal(
      fitted$scores$base_score,
      as.numeric(native$tscores),
      tolerance = 0
    )
    expect_identical(
      fitted$params$realized_soft_threshold,
      as.numeric(native$soft.thresh)
    )
  }
})

test_that("LPC fails transparently when its native threshold is unavailable", {
  skip_if_not_installed("lpc", minimum_version = "1.0.2.1")

  fixture <- make_high_dimensional_omics_fixture(
    n = 60L,
    p = 20L,
    seed = 9236L
  )

  expect_error(
    sb_fit_lpc_binary_selector(fixture$x, fixture$outcome),
    class = "stormbreaker_invalid_lpc_selector"
  )
})

test_that("stability selectors reproduce native selected variables", {
  skip_if_not_installed("stabs", minimum_version = "0.7-1")
  skip_if_not_installed("glmnet", minimum_version = "4.1-2")

  binary <- make_high_dimensional_omics_fixture(p = 60L, seed = 9216L)
  survival <- make_high_dimensional_survival_fixture(
    p = 60L,
    seed = 9217L
  )
  specifications <- list(
    list(
      selector = "stabsel_binary_lasso",
      x = binary$x,
      response = binary$outcome,
      fit = function() {
        sb_fit_stabsel_binary_selector(binary$x, binary$outcome)
      }
    ),
    list(
      selector = "stabsel_cox_lasso",
      x = survival$x,
      response = survival::Surv(survival$time, survival$event),
      fit = function() {
        sb_fit_stabsel_cox_selector(
          survival$x,
          survival$time,
          survival$event
        )
      }
    )
  )

  for (specification in specifications) {
    params <- sb_stabsel_selector_params(specification$selector)
    subsample_q <- sb_stabsel_subsample_q(
      nrow(specification$x),
      ncol(specification$x)
    )
    fitfun_args <- list(
      type = params$base_selection_type,
      family = params$family,
      standardize = params$standardize
    )
    if (identical(params$family, "cox")) {
      fitfun_args$cox.ties <- params$cox_ties
    }
    native <- withr::with_seed(params$seed, stabs::stabsel(
      x = data.matrix(specification$x),
      y = specification$response,
      fitfun = stabs::glmnet.lasso,
      args.fitfun = fitfun_args,
      cutoff = params$cutoff,
      q = subsample_q,
      B = params$complementary_pairs,
      sampling.type = params$sampling_type,
      assumption = params$assumption,
      papply = parallel::mclapply,
      mc.cores = 1L,
      verbose = FALSE
    ))

    set.seed(9218L)
    caller_rng <- .Random.seed
    fitted <- specification$fit()

    expect_identical(.Random.seed, caller_rng)
    expect_identical(
      fitted$selected_features,
      names(specification$x)[unname(as.integer(stabs::selected(native)))]
    )
    expect_equal(
      fitted$scores$selection_probability,
      unname(native$max[names(specification$x)]),
      tolerance = 0
    )
  }
})

test_that("Cox stability selection locks glmnet ties to Breslow", {
  skip_if_not_installed("stabs", minimum_version = "0.7-1")
  skip_if_not_installed("glmnet", minimum_version = "4.1-2")

  fixture <- make_high_dimensional_survival_fixture(
    n = 30L,
    p = 6L,
    seed = 9219L
  )
  captured <- new.env(parent = emptyenv())
  local_mocked_bindings(
    stabsel = function(x, y, ..., args.fitfun) {
      captured$args.fitfun <- args.fitfun
      structure(
        list(
          selected = 1L,
          max = stats::setNames(
            c(1, rep(0, ncol(x) - 1L)),
            colnames(x)
          ),
          PFER = 0
        ),
        class = "stabsel"
      )
    },
    .package = "stabs"
  )

  sb_fit_stabsel_cox_selector(
    fixture$x,
    fixture$time,
    fixture$event
  )

  expect_identical(
    captured$args.fitfun,
    list(
      type = "conservative",
      family = "cox",
      standardize = TRUE,
      cox.ties = "breslow"
    )
  )
})

test_that("HDDA predicts native class-one posteriors from p greater than n", {
  skip_if_not_installed("HDclassif", minimum_version = "2.2.2")

  fixture <- make_high_dimensional_omics_fixture(n = 60L, p = 90L)
  adapter <- sb_hdda_binary_learner_adapter("hdda_binary_cattell")
  model <- adapter$fit(fixture$x, fixture$outcome)
  assessment <- fixture$x[1:7, , drop = FALSE]
  native <- stats::predict(model$engine, data = data.matrix(assessment))
  observed <- adapter$predict(model, assessment)

  expect_equal(observed, unname(native$posterior[, "1"]), tolerance = 1e-12)
  state <- unserialize(serialize(adapter$extract_state(model), NULL))
  expect_equal(
    adapter$predict_state(state, assessment)$probability,
    observed,
    tolerance = 1e-12
  )
})

test_that("VSURF and HDDA complete final fit and outcome-blind prediction", {
  skip_if_not_installed("VSURF", minimum_version = "1.2.1")
  skip_if_not_installed("HDclassif", minimum_version = "2.2.2")

  fixture <- make_high_dimensional_omics_fixture(p = 24L, seed = 9215L)
  data <- data.frame(
    patient_id = sprintf("HF%03d", seq_len(nrow(fixture$x))),
    outcome = fixture$outcome,
    fixture$x,
    check.names = FALSE
  )
  task <- sb_binary_task(
    data,
    "patient_id",
    "outcome",
    names(fixture$x)
  )
  new_data <- data[1:6, c("patient_id", names(fixture$x))]
  new_data$patient_id <- paste0("NEW_", new_data$patient_id)
  fits <- list(
    sb_fit_final(task, sb_pipeline(
      selector = "vsurf_binary_prediction",
      learner = "logistic"
    )),
    sb_fit_final(task, sb_pipeline(learner = "hdda_binary_cattell"))
  )

  for (fit in fits) {
    expect_identical(fit$status, "success")
    prediction <- sb_predict_final(fit, new_data, "patient_id")
    expect_identical(prediction$status, "success")
    expect_identical(nrow(prediction$predictions), nrow(new_data))
    expect_gte(min(prediction$predictions$probability), 0)
    expect_lte(max(prediction$predictions$probability), 1)
  }
  expect_identical(
    fits[[1L]]$provenance$vsurf_version,
    as.character(utils::packageVersion("VSURF"))
  )
  expect_identical(
    fits[[2L]]$provenance$hdclassif_version,
    as.character(utils::packageVersion("HDclassif"))
  )
})

test_that("information and stability selectors complete final prediction", {
  skip_if_not_installed("praznik", minimum_version = "13.0.0")
  skip_if_not_installed("stabs", minimum_version = "0.7-1")
  skip_if_not_installed("glmnet", minimum_version = "4.1-2")

  binary <- make_high_dimensional_omics_fixture(
    n = 84L,
    p = 40L,
    seed = 9223L
  )
  binary_data <- data.frame(
    patient_id = sprintf("MB%03d", seq_len(nrow(binary$x))),
    outcome = binary$outcome,
    binary$x,
    check.names = FALSE
  )
  binary_task <- sb_binary_task(
    binary_data,
    "patient_id",
    "outcome",
    names(binary$x)
  )
  binary_new <- binary_data[1:6, c("patient_id", names(binary$x))]
  binary_new$patient_id <- paste0("NEW_", binary_new$patient_id)

  survival <- make_high_dimensional_survival_fixture(
    n = 90L,
    p = 40L,
    seed = 9224L
  )
  survival_data <- data.frame(
    patient_id = sprintf("MS%03d", seq_len(nrow(survival$x))),
    time = survival$time,
    event = survival$event,
    survival$x,
    check.names = FALSE
  )
  survival_task <- sb_surv_task(
    survival_data,
    "patient_id",
    "time",
    "event",
    names(survival$x),
    "month",
    prediction_horizon = stats::median(survival$time)
  )
  survival_new <- survival_data[1:6, c("patient_id", names(survival$x))]
  survival_new$patient_id <- paste0("NEW_", survival_new$patient_id)

  binary_selectors <- c(
    "mrmr_binary_nlogn",
    sb_information_binary_selector_ids(),
    "stabsel_binary_lasso"
  )
  binary_fits <- lapply(binary_selectors, function(selector) {
    sb_fit_final(binary_task, sb_pipeline(
      selector = selector,
      learner = "gaussian_naive_bayes"
    ))
  })
  fits <- c(
    binary_fits,
    list(
      sb_fit_final(survival_task, sb_pipeline(
        selector = "stabsel_cox_lasso",
        learner = "coxph"
      ))
    )
  )
  new_data <- c(
    rep(list(binary_new), length(binary_fits)),
    list(survival_new)
  )

  for (index in seq_along(fits)) {
    expect_identical(fits[[index]]$status, "success")
    prediction <- sb_predict_final(
      fits[[index]],
      new_data[[index]],
      "patient_id"
    )
    expect_identical(prediction$status, "success")
    expect_identical(nrow(prediction$predictions), 6L)
  }
  praznik_indexes <- which(binary_selectors != "stabsel_binary_lasso")
  expect_identical(
    vapply(
      fits[praznik_indexes],
      function(fit) fit$provenance$praznik_version,
      character(1)
    ),
    rep(as.character(utils::packageVersion("praznik")), 5L)
  )
  expect_identical(
    fits[[length(binary_fits)]]$provenance$stabs_version,
    as.character(utils::packageVersion("stabs"))
  )
  expect_identical(
    fits[[length(fits)]]$provenance$glmnet_version,
    as.character(utils::packageVersion("glmnet"))
  )

  cv <- sb_run_cv(
    survival_task,
    pipeline = sb_pipeline(
      selector = "stabsel_cox_lasso",
      learner = "coxph"
    ),
    folds = 3L,
    seed = 9225L
  )
  expect_identical(cv$status, "success")
  expect_identical(
    cv$provenance$stabs_version,
    as.character(utils::packageVersion("stabs"))
  )
})

test_that("LPC selectors complete serializable outcome-blind prediction", {
  skip_if_not_installed("lpc", minimum_version = "1.0.2.1")
  skip_if_not_installed("e1071", minimum_version = "1.7-13")

  binary <- make_high_dimensional_omics_fixture(
    n = 84L,
    p = 60L,
    seed = 9237L
  )
  binary_data <- data.frame(
    patient_id = sprintf("LB%03d", seq_len(nrow(binary$x))),
    outcome = binary$outcome,
    binary$x,
    check.names = FALSE
  )
  binary_task <- sb_binary_task(
    binary_data,
    "patient_id",
    "outcome",
    names(binary$x)
  )

  survival <- make_high_dimensional_survival_fixture(
    n = 90L,
    p = 60L,
    seed = 9238L
  )
  survival_data <- data.frame(
    patient_id = sprintf("LS%03d", seq_len(nrow(survival$x))),
    time = survival$time,
    event = survival$event,
    survival$x,
    check.names = FALSE
  )
  survival_task <- sb_surv_task(
    survival_data,
    "patient_id",
    "time",
    "event",
    names(survival$x),
    "month",
    prediction_horizon = stats::median(survival$time)
  )

  fits <- list(
    sb_fit_final(binary_task, sb_pipeline(
      selector = "lpc_binary_nlogn",
      learner = "gaussian_naive_bayes"
    )),
    sb_fit_final(survival_task, sb_pipeline(
      selector = "lpc_cox_nlogn",
      learner = "coxph"
    ))
  )
  new_data <- list(
    binary_data[1:6, c("patient_id", names(binary$x))],
    survival_data[1:6, c("patient_id", names(survival$x))]
  )

  for (index in seq_along(fits)) {
    fit <- unserialize(serialize(fits[[index]], NULL))
    prediction <- sb_predict_final(fit, new_data[[index]], "patient_id")

    expect_identical(fit$status, "success")
    expect_identical(prediction$status, "success")
    expect_identical(nrow(prediction$predictions), 6L)
    prediction_values <- if (index == 1L) {
      prediction$predictions$probability
    } else {
      prediction$predictions$risk
    }
    expect_true(all(is.finite(prediction_values)))
    if (index == 2L) {
      expect_true(all(is.finite(
        prediction$predictions$survival_probability
      )))
    }
    expect_identical(
      fit$provenance$lpc_version,
      as.character(utils::packageVersion("lpc"))
    )
  }
  expect_true(all(
    fits[[1L]]$selection$selected_features %in% names(binary$x)
  ))
  expect_true(all(
    fits[[2L]]$selection$selected_features %in% names(survival$x)
  ))
})

test_that("high-dimensional methods resolve through admitted binary routes", {
  registry <- sb_pipeline_registry()
  unrestricted <- registry$learner[
    registry$task_type == "binary" & registry$selector == "none" &
      !registry$learner %in% c(
        "pre_binary_rule_ensemble",
        sb_slope_binary_learner_ids(),
        "cancerclass_binary"
      )
  ]

  for (selector in setdiff(
    sb_high_dimensional_binary_selector_ids(),
    "stabsel_binary_lasso"
  )) {
    expect_setequal(
      registry$learner[registry$selector == selector],
      unrestricted
    )
  }
  expect_setequal(
    registry$learner[registry$selector == "stabsel_binary_lasso"],
    setdiff(unrestricted, "glmnet_logistic")
  )
  expect_setequal(
    registry$learner[registry$selector == "lpc_cox_nlogn"],
    registry$learner[registry$selector == "boruta_survival_roughfix"]
  )
  expect_setequal(
    registry$learner[registry$selector == "stabsel_cox_lasso"],
    setdiff(
      registry$learner[registry$selector == "boruta_survival_roughfix"],
      "glmnet_cox"
    )
  )
  expect_setequal(
    registry$selector[registry$learner == "hdda_binary_cattell"],
    expected_fixed_learner_selectors("binary")
  )

  binary_specifications <- lapply(
    sb_high_dimensional_binary_selector_ids(),
    function(selector) {
      sb_pipeline(selector = selector, learner = "logistic")
    }
  )
  specifications <- c(
    binary_specifications,
    list(
      sb_pipeline(selector = "lpc_cox_nlogn", learner = "coxph"),
      sb_pipeline(selector = "stabsel_cox_lasso", learner = "coxph"),
      sb_pipeline(learner = "hdda_binary_cattell")
    )
  )
  for (pipeline in specifications) {
    runtime <- sb_resolve_pipeline_runtime(pipeline)
    expect_identical(sb_validate_pipeline_runtime(runtime, pipeline), runtime)
  }
})

test_that("SIS is refitted inside glmnet inner analysis folds", {
  skip_if_not_installed("SIS", minimum_version = "1.5")
  skip_if_not_installed("glmnet", minimum_version = "4.1.2")

  fixture <- make_high_dimensional_omics_fixture(n = 84L, p = 36L)
  data <- data.frame(
    patient_id = sprintf("HD%03d", seq_len(nrow(fixture$x))),
    outcome = fixture$outcome,
    fixture$x,
    check.names = FALSE
  )
  task <- sb_binary_task(
    data,
    "patient_id",
    "outcome",
    names(fixture$x)
  )
  fit <- sb_fit_final(task, sb_pipeline(
    selector = "sis_binary_nlogn",
    learner = "glmnet_logistic",
    inner_folds = 3L,
    inner_seed = 9221L
  ))

  expect_identical(fit$status, "success")
  expect_identical(
    fit$tuning$inner_selection_method,
    "glmnet_exclude_function"
  )
  expect_length(fit$tuning$inner_selections, 3L)
  expect_identical(
    fit$provenance$sis_version,
    as.character(utils::packageVersion("SIS"))
  )
})

test_that("praznik selectors are refitted inside glmnet inner folds", {
  skip_if_not_installed("praznik", minimum_version = "13.0.0")
  skip_if_not_installed("glmnet", minimum_version = "4.1-2")

  binary <- make_high_dimensional_omics_fixture(
    n = 90L,
    p = 36L,
    seed = 9231L
  )
  binary_data <- data.frame(
    patient_id = sprintf("IB%03d", seq_len(nrow(binary$x))),
    outcome = binary$outcome,
    binary$x,
    check.names = FALSE
  )
  binary_task <- sb_binary_task(
    binary_data,
    "patient_id",
    "outcome",
    names(binary$x)
  )
  selectors <- c("mrmr_binary_nlogn", sb_information_binary_selector_ids())
  for (selector in selectors) {
    fit <- sb_fit_final(binary_task, sb_pipeline(
      selector = selector,
      learner = "glmnet_logistic",
      inner_folds = 3L,
      inner_seed = 9232L
    ))

    expect_identical(fit$status, "success")
    expect_identical(
      fit$tuning$inner_selection_method,
      "glmnet_exclude_function"
    )
    expect_length(fit$tuning$inner_selections, 3L)
    expect_gt(
      min(vapply(
        fit$tuning$inner_selections,
        function(selection) length(selection$selected_features),
        integer(1)
      )),
      0L
    )
  }
})

test_that("LPC selectors are refitted inside glmnet inner folds", {
  skip_if_not_installed("lpc", minimum_version = "1.0.2.1")
  skip_if_not_installed("glmnet", minimum_version = "4.1-2")

  binary <- make_high_dimensional_omics_fixture(
    n = 90L,
    p = 60L,
    seed = 9239L
  )
  binary_data <- data.frame(
    patient_id = sprintf("LI%03d", seq_len(nrow(binary$x))),
    outcome = binary$outcome,
    binary$x,
    check.names = FALSE
  )
  binary_task <- sb_binary_task(
    binary_data,
    "patient_id",
    "outcome",
    names(binary$x)
  )

  survival <- make_high_dimensional_survival_fixture(
    n = 90L,
    p = 60L,
    seed = 9240L
  )
  survival_data <- data.frame(
    patient_id = sprintf("LC%03d", seq_len(nrow(survival$x))),
    time = survival$time,
    event = survival$event,
    survival$x,
    check.names = FALSE
  )
  survival_task <- sb_surv_task(
    survival_data,
    "patient_id",
    "time",
    "event",
    names(survival$x),
    "month",
    prediction_horizon = stats::median(survival$time)
  )

  fits <- list(
    sb_fit_final(binary_task, sb_pipeline(
      selector = "lpc_binary_nlogn",
      learner = "glmnet_logistic",
      inner_folds = 3L,
      inner_seed = 9241L
    )),
    sb_fit_final(survival_task, sb_pipeline(
      selector = "lpc_cox_nlogn",
      learner = "glmnet_cox",
      inner_folds = 3L,
      inner_seed = 9242L
    ))
  )

  for (fit in fits) {
    expect_identical(fit$status, "success")
    expect_identical(
      fit$tuning$inner_selection_method,
      "glmnet_exclude_function"
    )
    expect_length(fit$tuning$inner_selections, 3L)
    expect_gt(
      min(vapply(
        fit$tuning$inner_selections,
        function(selection) length(selection$selected_features),
        integer(1)
      )),
      0L
    )
    expect_identical(
      fit$provenance$lpc_version,
      as.character(utils::packageVersion("lpc"))
    )
  }
})
