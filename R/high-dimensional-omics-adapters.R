sb_high_dimensional_binary_selector_ids <- function() {
  c(
    "sis_binary_nlogn",
    "vsurf_binary_prediction",
    "mrmr_binary_nlogn",
    sb_information_binary_selector_ids(),
    "lpc_binary_nlogn",
    "stabsel_binary_lasso"
  )
}

sb_information_binary_selector_ids <- function() {
  c(
    "jmi_binary_nlogn",
    "cmim_binary_nlogn",
    "disr_binary_nlogn",
    "jmim_binary_nlogn"
  )
}

sb_high_dimensional_survival_selector_ids <- function() {
  c("lpc_cox_nlogn", "stabsel_cox_lasso")
}

sb_lpc_selector_ids <- function() {
  c("lpc_binary_nlogn", "lpc_cox_nlogn")
}

sb_high_dimensional_selector_ids <- function() {
  c(
    sb_high_dimensional_binary_selector_ids(),
    sb_high_dimensional_survival_selector_ids()
  )
}

sb_hdda_binary_learner_ids <- function() {
  "hdda_binary_cattell"
}

sb_high_dimensional_method_ids <- function() {
  c(sb_high_dimensional_selector_ids(), sb_hdda_binary_learner_ids())
}

sb_high_dimensional_method_packages <- function(method) {
  switch(method,
    sis_binary_nlogn = "SIS",
    vsurf_binary_prediction = "VSURF",
    mrmr_binary_nlogn = "praznik",
    jmi_binary_nlogn = "praznik",
    cmim_binary_nlogn = "praznik",
    disr_binary_nlogn = "praznik",
    jmim_binary_nlogn = "praznik",
    lpc_binary_nlogn = "lpc",
    lpc_cox_nlogn = "lpc",
    stabsel_binary_lasso = c("stabs", "glmnet"),
    stabsel_cox_lasso = c("stabs", "glmnet"),
    hdda_binary_cattell = "HDclassif",
    sb_abort(
      paste0("Unknown high-dimensional method `", method, "`."),
      "stormbreaker_unsupported_high_dimensional_method"
    )
  )
}

sb_require_high_dimensional_method <- function(method) {
  packages <- sb_high_dimensional_method_packages(method)
  missing <- packages[
    !vapply(packages, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing) > 0L) {
    sb_abort(
      paste0(
        "Package", if (length(missing) == 1L) " `" else "s `",
        paste(missing, collapse = "`, `"), "` ",
        if (length(missing) == 1L) "is" else "are",
        " required for method `", method, "`."
      ),
      "stormbreaker_missing_high_dimensional_dependency"
    )
  }
  invisible(packages)
}

sb_high_dimensional_method_provenance <- function(method) {
  packages <- sb_high_dimensional_method_packages(method)
  fields <- c(
    SIS = "sis_version",
    VSURF = "vsurf_version",
    praznik = "praznik_version",
    lpc = "lpc_version",
    stabs = "stabs_version",
    glmnet = "glmnet_version",
    HDclassif = "hdclassif_version"
  )
  versions <- lapply(packages, function(package) {
    if (requireNamespace(package, quietly = TRUE)) {
      as.character(utils::packageVersion(package))
    } else {
      NA_character_
    }
  })
  stats::setNames(versions, unname(fields[packages]))
}

sb_high_dimensional_native_call <- function(code, method, stage) {
  tryCatch(
    withCallingHandlers(
      force(code),
      warning = function(warning) {
        sb_abort(
          paste0(
            "Method `", method, "` raised a ", stage, " warning: ",
            conditionMessage(warning)
          ),
          "stormbreaker_high_dimensional_native_warning",
          parent = warning
        )
      }
    ),
    error = function(error) {
      if (inherits(error, "stormbreaker_error")) {
        stop(error)
      }
      sb_abort(
        paste0(
          "Method `", method, "` failed during ", stage, ": ",
          conditionMessage(error)
        ),
        "stormbreaker_high_dimensional_native_error",
        parent = error
      )
    }
  )
}

sb_validate_high_dimensional_binary_inputs <- function(x, outcome, method) {
  x <- sb_validate_predictor_frame(x, "x")
  outcome <- sb_validate_binary_outcome_vector(outcome, nrow(x))
  if (!all(c(0L, 1L) %in% outcome)) {
    sb_abort(
      paste0("A `", method, "` analysis fold must contain both classes."),
      "stormbreaker_high_dimensional_single_class"
    )
  }
  list(x = x, outcome = outcome)
}

sb_sis_binary_selector_params <- function() {
  list(
    family = "binomial",
    penalty = "lasso",
    tuning = "bic",
    selection_output = "sis.ix0_before_penalized_fit",
    screen_size_rule = "min(p, max(2, floor(n/log(n))))",
    iterative = FALSE,
    standardize = TRUE,
    parallel = FALSE,
    seed = 1L
  )
}

sb_high_dimensional_nlogn_size <- function(n, p) {
  as.integer(min(p, max(2L, floor(n / log(n)))))
}

sb_sis_binary_screen_size <- function(n, p) {
  sb_high_dimensional_nlogn_size(n, p)
}

sb_fit_sis_binary_selector <- function(x, outcome) {
  method <- "sis_binary_nlogn"
  sb_require_high_dimensional_method(method)
  inputs <- sb_validate_high_dimensional_binary_inputs(x, outcome, method)
  if (ncol(inputs$x) < 2L) {
    sb_abort(
      "SIS screening requires at least two candidate predictors.",
      "stormbreaker_sis_too_few_predictors"
    )
  }

  params <- sb_sis_binary_selector_params()
  screen_size <- sb_sis_binary_screen_size(nrow(inputs$x), ncol(inputs$x))
  native <- sb_high_dimensional_native_call(
    withr::with_seed(params$seed, SIS::SIS(
      x = data.matrix(inputs$x),
      y = inputs$outcome,
      family = params$family,
      penalty = params$penalty,
      tune = params$tuning,
      type.measure = "deviance",
      nsis = screen_size,
      iter = params$iterative,
      varISIS = "vanilla",
      perm = FALSE,
      greedy = FALSE,
      seed = params$seed,
      standardize = params$standardize,
      boot_ci = FALSE,
      parallel = params$parallel
    )),
    method,
    "fit"
  )

  selected_index <- as.integer(native$sis.ix0)
  ranking <- if (is.list(native$ix_list) && length(native$ix_list) >= 1L) {
    as.integer(native$ix_list[[1L]])
  } else {
    integer()
  }
  p <- ncol(inputs$x)
  valid <- length(selected_index) == screen_size &&
    !anyNA(selected_index) && !anyDuplicated(selected_index) &&
    all(selected_index >= 1L & selected_index <= p) &&
    length(ranking) == p && !anyNA(ranking) &&
    !anyDuplicated(ranking) && setequal(ranking, seq_len(p))
  if (!valid) {
    sb_abort(
      "SIS did not return a complete marginal screening ranking.",
      "stormbreaker_invalid_sis_selector"
    )
  }

  source_features <- names(inputs$x)
  selected_features <- source_features[selected_index]
  active_params <- params
  active_params$screen_size <- screen_size
  structure(
    list(
      id = method,
      source_features = source_features,
      selected_features = selected_features,
      params = active_params,
      scores = data.frame(
        feature = source_features,
        rank = as.integer(match(seq_len(p), ranking)),
        selected = source_features %in% selected_features,
        stringsAsFactors = FALSE
      )
    ),
    class = "stormbreaker_selector_fit"
  )
}

sb_vsurf_binary_selector_params <- function() {
  list(
    method = "VSURF",
    selected_subset = "prediction",
    random_forest_engine = "randomForest",
    mtry_rule = "max(floor(p/3), 1)",
    threshold_trees = 500L,
    threshold_forests = 20L,
    interpretation_trees = 100L,
    interpretation_forests = 10L,
    prediction_trees = 100L,
    prediction_forests = 10L,
    nmin = 1L,
    nsd = 1,
    nmj = 1,
    parallel = FALSE,
    seed = 1L
  )
}

sb_fit_vsurf_binary_selector <- function(x, outcome) {
  method <- "vsurf_binary_prediction"
  inputs <- sb_validate_high_dimensional_binary_inputs(x, outcome, method)
  sb_fit_vsurf_classification_selector(inputs$x, factor(inputs$outcome, levels = c(0L, 1L)),
    "vsurf_binary_prediction")
}

sb_fit_vsurf_classification_selector <- function(x, response, method) {
  sb_require_high_dimensional_method("vsurf_binary_prediction")
  inputs <- list(x = x)
  if (ncol(inputs$x) < 2L) {
    sb_abort(
      "VSURF selection requires at least two candidate predictors.",
      "stormbreaker_vsurf_too_few_predictors"
    )
  }

  params <- sb_vsurf_binary_selector_params()
  mtry <- max(floor(ncol(inputs$x) / 3), 1L)
  native <- sb_high_dimensional_native_call(
    withr::with_seed(params$seed, VSURF::VSURF(
      x = data.matrix(inputs$x),
      y = response,
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
      parallel = params$parallel,
      verbose = FALSE
    )),
    method,
    "fit"
  )

  selected_index <- as.integer(native$varselect.pred)
  ranking <- as.integer(native$imp.mean.dec.ind)
  importance <- as.numeric(native$imp.mean.dec)
  importance_sd <- as.numeric(native$imp.sd.dec)
  p <- ncol(inputs$x)
  valid <- length(selected_index) > 0L && !anyNA(selected_index) &&
    !anyDuplicated(selected_index) &&
    all(selected_index >= 1L & selected_index <= p) &&
    length(ranking) == p && !anyNA(ranking) &&
    !anyDuplicated(ranking) && setequal(ranking, seq_len(p)) &&
    length(importance) == p && length(importance_sd) == p &&
    !anyNA(importance) && !anyNA(importance_sd) &&
    all(is.finite(importance)) && all(is.finite(importance_sd))
  if (!valid) {
    sb_abort(
      "VSURF did not return a valid prediction-stage feature subset.",
      "stormbreaker_invalid_vsurf_selector"
    )
  }

  source_features <- names(inputs$x)
  selected_features <- source_features[selected_index]
  rank_by_feature <- integer(p)
  mean_by_feature <- numeric(p)
  sd_by_feature <- numeric(p)
  rank_by_feature[ranking] <- seq_len(p)
  mean_by_feature[ranking] <- importance
  sd_by_feature[ranking] <- importance_sd
  active_params <- params
  active_params$mtry <- as.integer(mtry)
  structure(
    list(
      id = method,
      source_features = source_features,
      selected_features = selected_features,
      params = active_params,
      scores = data.frame(
        feature = source_features,
        importance = mean_by_feature,
        importance_sd = sd_by_feature,
        rank = as.integer(rank_by_feature),
        selected = source_features %in% selected_features,
        stringsAsFactors = FALSE
      ),
      selection_history = data.frame(
        stage = c("threshold", "interpretation", "prediction"),
        n_features = as.integer(c(
          length(native$varselect.thres),
          length(native$varselect.interp),
          length(native$varselect.pred)
        )),
        stringsAsFactors = FALSE
      )
    ),
    class = "stormbreaker_selector_fit"
  )
}

sb_mrmr_binary_selector_params <- function() {
  list(
    method = "minimum_redundancy_maximum_relevance",
    screen_size_rule = "min(p, max(2, floor(n/log(n))))",
    positive_scores_only = FALSE,
    threads = 1L
  )
}

sb_fit_mrmr_binary_selector <- function(x, outcome) {
  method <- "mrmr_binary_nlogn"
  sb_require_high_dimensional_method(method)
  inputs <- sb_validate_high_dimensional_binary_inputs(x, outcome, method)
  if (ncol(inputs$x) < 2L) {
    sb_abort(
      "mRMR screening requires at least two candidate predictors.",
      "stormbreaker_mrmr_too_few_predictors"
    )
  }

  params <- sb_mrmr_binary_selector_params()
  screen_size <- sb_high_dimensional_nlogn_size(
    nrow(inputs$x),
    ncol(inputs$x)
  )
  native <- sb_high_dimensional_native_call(
    praznik::MRMR(
      X = inputs$x,
      Y = factor(inputs$outcome, levels = c(0L, 1L)),
      k = screen_size,
      positive = params$positive_scores_only,
      threads = params$threads
    ),
    method,
    "fit"
  )

  selected_index <- as.integer(native$selection)
  selected_score <- as.numeric(native$score)
  p <- ncol(inputs$x)
  valid <- length(selected_index) == screen_size &&
    !anyNA(selected_index) && !anyDuplicated(selected_index) &&
    all(selected_index >= 1L & selected_index <= p) &&
    length(selected_score) == screen_size && !anyNA(selected_score) &&
    all(is.finite(selected_score))
  if (!valid) {
    sb_abort(
      "mRMR did not return a valid ordered feature subset.",
      "stormbreaker_invalid_mrmr_selector"
    )
  }

  source_features <- names(inputs$x)
  selected_features <- source_features[selected_index]
  score_by_feature <- rep(NA_real_, p)
  score_by_feature[selected_index] <- selected_score
  active_params <- params
  active_params$screen_size <- screen_size
  structure(
    list(
      id = method,
      source_features = source_features,
      selected_features = selected_features,
      params = active_params,
      scores = data.frame(
        feature = source_features,
        mrmr_score = score_by_feature,
        rank = as.integer(match(seq_len(p), selected_index)),
        selected = source_features %in% selected_features,
        stringsAsFactors = FALSE
      )
    ),
    class = "stormbreaker_selector_fit"
  )
}

sb_information_binary_selector_params <- function(selector) {
  native_methods <- c(
    jmi_binary_nlogn = "JMI",
    cmim_binary_nlogn = "CMIM",
    disr_binary_nlogn = "DISR",
    jmim_binary_nlogn = "JMIM"
  )
  criteria <- c(
    jmi_binary_nlogn = "joint_mutual_information",
    cmim_binary_nlogn = "conditional_mutual_information_maximization",
    disr_binary_nlogn = "double_input_symmetrical_relevance",
    jmim_binary_nlogn = "minimal_joint_mutual_information_maximization"
  )
  if (!selector %in% names(native_methods)) {
    sb_abort(
      paste0("Unknown information selector `", selector, "`."),
      "stormbreaker_unsupported_selector"
    )
  }
  list(
    method = unname(criteria[[selector]]),
    native_method = unname(native_methods[[selector]]),
    screen_size_rule = "min(p, max(2, floor(n/log(n))))",
    threads = 1L
  )
}

sb_run_information_binary_selector <- function(selector, x, outcome,
                                               screen_size, threads) {
  switch(selector,
    jmi_binary_nlogn = praznik::JMI(
      X = x,
      Y = outcome,
      k = screen_size,
      threads = threads
    ),
    cmim_binary_nlogn = praznik::CMIM(
      X = x,
      Y = outcome,
      k = screen_size,
      threads = threads
    ),
    disr_binary_nlogn = praznik::DISR(
      X = x,
      Y = outcome,
      k = screen_size,
      threads = threads
    ),
    jmim_binary_nlogn = praznik::JMIM(
      X = x,
      Y = outcome,
      k = screen_size,
      threads = threads
    )
  )
}

sb_fit_information_binary_selector <- function(x, outcome, selector) {
  params <- sb_information_binary_selector_params(selector)
  sb_require_high_dimensional_method(selector)
  inputs <- sb_validate_high_dimensional_binary_inputs(x, outcome, selector)
  if (ncol(inputs$x) < 2L) {
    sb_abort(
      paste0(params$native_method, " screening requires at least two candidate predictors."),
      "stormbreaker_information_selector_too_few_predictors"
    )
  }

  screen_size <- sb_high_dimensional_nlogn_size(
    nrow(inputs$x),
    ncol(inputs$x)
  )
  native <- sb_high_dimensional_native_call(
    sb_run_information_binary_selector(
      selector,
      inputs$x,
      factor(inputs$outcome, levels = c(0L, 1L)),
      screen_size,
      params$threads
    ),
    selector,
    "fit"
  )

  selected_index <- as.integer(native$selection)
  selected_score <- as.numeric(native$score)
  p <- ncol(inputs$x)
  valid <- length(selected_index) == screen_size &&
    !anyNA(selected_index) && !anyDuplicated(selected_index) &&
    all(selected_index >= 1L & selected_index <= p) &&
    length(selected_score) == screen_size && !anyNA(selected_score) &&
    all(is.finite(selected_score))
  if (!valid) {
    sb_abort(
      paste0(
        params$native_method,
        " did not return a valid ordered feature subset."
      ),
      "stormbreaker_invalid_information_selector"
    )
  }

  source_features <- names(inputs$x)
  selected_features <- source_features[selected_index]
  score_by_feature <- rep(NA_real_, p)
  score_by_feature[selected_index] <- selected_score
  active_params <- params
  active_params$screen_size <- screen_size
  structure(
    list(
      id = selector,
      source_features = source_features,
      selected_features = selected_features,
      params = active_params,
      scores = data.frame(
        feature = source_features,
        information_score = score_by_feature,
        rank = as.integer(match(seq_len(p), selected_index)),
        selected = source_features %in% selected_features,
        stringsAsFactors = FALSE
      )
    ),
    class = "stormbreaker_selector_fit"
  )
}

sb_lpc_selector_params <- function(selector) {
  if (!selector %in% sb_lpc_selector_ids()) {
    sb_abort(
      paste0("Unknown LPC selector `", selector, "`."),
      "stormbreaker_unsupported_selector"
    )
  }
  list(
    method = "lassoed_principal_components",
    native_type = if (identical(selector, "lpc_binary_nlogn")) {
      "two class"
    } else {
      "survival"
    },
    ranking = "decreasing_absolute_lpc_score",
    screen_size_rule = "min(p, max(2, floor(n/log(n))))",
    soft_threshold = "native_10_half_split_predictive_advantage",
    seed = 1L
  )
}

sb_fit_lpc_selector <- function(x, endpoint, selector) {
  sb_require_high_dimensional_method(selector)
  params <- sb_lpc_selector_params(selector)
  inputs <- if (identical(selector, "lpc_binary_nlogn")) {
    binary <- sb_validate_high_dimensional_binary_inputs(
      x,
      endpoint$outcome,
      selector
    )
    list(
      x = binary$x,
      response = binary$outcome + 1L,
      censoring_status = NULL
    )
  } else {
    survival <- sb_validate_classical_survival_inputs(
      x,
      endpoint$time,
      endpoint$event,
      selector
    )
    list(
      x = survival$x,
      response = survival$time,
      censoring_status = survival$event
    )
  }
  screen_size <- sb_high_dimensional_nlogn_size(
    nrow(inputs$x),
    ncol(inputs$x)
  )
  native <- sb_high_dimensional_native_call(
    {
      fitted <- NULL
      invisible(utils::capture.output(
        fitted <- withr::with_seed(
          params$seed,
          lpc::LPC(
            x = t(data.matrix(inputs$x)),
            y = inputs$response,
            type = params$native_type,
            soft.thresh = NULL,
            censoring.status = inputs$censoring_status
          )
        ),
        type = "output"
      ))
      fitted
    },
    selector,
    "fit"
  )

  lpc_score <- as.numeric(native$lpcscores)
  base_score <- as.numeric(native$tscores)
  soft_threshold <- as.numeric(native$soft.thresh)
  p <- ncol(inputs$x)
  valid <- length(lpc_score) == p && !anyNA(lpc_score) &&
    all(is.finite(lpc_score)) &&
    length(base_score) == p && !anyNA(base_score) &&
    all(is.finite(base_score)) &&
    length(soft_threshold) == 1L && !is.na(soft_threshold) &&
    is.finite(soft_threshold) && soft_threshold >= 0
  if (!valid) {
    sb_abort(
      if (p < 50L && length(soft_threshold) == 0L) {
        paste0("Native LPC automatic threshold search requires at least 50 predictors; ",
               "the supplied ", p, " predictors produced an empty threshold. ",
               "No padding or substitute threshold was used.")
      } else {
        "LPC did not return finite gene scores and one soft threshold."
      },
      if (p < 50L && length(soft_threshold) == 0L) {
        c("stormbreaker_lpc_dimension_unsupported", "stormbreaker_invalid_lpc_selector")
      } else "stormbreaker_invalid_lpc_selector"
    )
  }

  source_features <- names(inputs$x)
  ranking <- order(-abs(lpc_score), seq_len(p), method = "radix")
  selected_index <- ranking[seq_len(screen_size)]
  selected_features <- source_features[selected_index]
  active_params <- params
  active_params$screen_size <- screen_size
  active_params$realized_soft_threshold <- soft_threshold
  structure(
    list(
      id = selector,
      source_features = source_features,
      selected_features = selected_features,
      params = active_params,
      scores = data.frame(
        feature = source_features,
        lpc_score = lpc_score,
        base_score = base_score,
        rank = as.integer(match(seq_len(p), ranking)),
        selected = source_features %in% selected_features,
        stringsAsFactors = FALSE
      )
    ),
    class = "stormbreaker_selector_fit"
  )
}

sb_fit_lpc_binary_selector <- function(x, outcome) {
  sb_fit_lpc_selector(
    x,
    list(outcome = outcome),
    "lpc_binary_nlogn"
  )
}

sb_fit_lpc_cox_selector <- function(x, time, event) {
  sb_fit_lpc_selector(
    x,
    list(time = time, event = event),
    "lpc_cox_nlogn"
  )
}

sb_stabsel_selector_params <- function(selector) {
  if (!selector %in% c("stabsel_binary_lasso", "stabsel_cox_lasso")) {
    sb_abort(
      paste0("Unknown stability selector `", selector, "`."),
      "stormbreaker_unsupported_selector"
    )
  }
  params <- list(
    base_selector = "glmnet_lasso",
    family = if (identical(selector, "stabsel_binary_lasso")) {
      "binomial"
    } else {
      "cox"
    },
    subsample_q_rule = paste0(
      "min(max(2, floor(n/log(n))), floor(p/2)) with p >= 4"
    ),
    base_selection_type = "conservative",
    cutoff = 0.75,
    complementary_pairs = 50L,
    sampling_type = "SS",
    assumption = "unimodal",
    standardize = TRUE,
    parallel = FALSE,
    seed = 1L
  )
  if (identical(selector, "stabsel_cox_lasso")) {
    params$cox_ties <- "breslow"
  }
  params
}

sb_stabsel_subsample_q <- function(n, p) {
  as.integer(min(max(2L, floor(n / log(n))), floor(p / 2)))
}

sb_fit_stabsel_selector <- function(x, endpoint, selector) {
  sb_require_high_dimensional_method(selector)
  inputs <- if (identical(selector, "stabsel_binary_lasso")) {
    binary <- sb_validate_high_dimensional_binary_inputs(
      x,
      endpoint$outcome,
      selector
    )
    list(
      x = binary$x,
      response = binary$outcome
    )
  } else {
    survival <- sb_validate_classical_survival_inputs(
      x,
      endpoint$time,
      endpoint$event,
      selector
    )
    list(
      x = survival$x,
      response = survival::Surv(survival$time, survival$event)
    )
  }
  if (ncol(inputs$x) < 4L) {
    sb_abort(
      "Stability selection requires at least four candidate predictors.",
      "stormbreaker_stabsel_too_few_predictors"
    )
  }

  params <- sb_stabsel_selector_params(selector)
  subsample_q <- sb_stabsel_subsample_q(
    nrow(inputs$x),
    ncol(inputs$x)
  )
  fitfun_args <- list(
    type = params$base_selection_type,
    family = params$family,
    standardize = params$standardize
  )
  if (identical(params$family, "cox")) {
    fitfun_args$cox.ties <- params$cox_ties
  }
  native <- sb_high_dimensional_native_call(
    withr::with_seed(params$seed, stabs::stabsel(
      x = data.matrix(inputs$x),
      y = inputs$response,
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
    )),
    selector,
    "fit"
  )

  selected_index <- unname(as.integer(stabs::selected(native)))
  source_features <- names(inputs$x)
  probabilities <- native$max[source_features]
  p <- ncol(inputs$x)
  valid <- length(probabilities) == p && !anyNA(probabilities) &&
    all(is.finite(probabilities)) &&
    all(probabilities >= 0 & probabilities <= 1) &&
    !anyNA(selected_index) && !anyDuplicated(selected_index) &&
    all(selected_index >= 1L & selected_index <= p) &&
    is.numeric(native$PFER) && length(native$PFER) == 1L &&
    is.finite(native$PFER) && native$PFER >= 0
  if (!valid) {
    sb_abort(
      "Stability selection did not return valid selection probabilities.",
      "stormbreaker_invalid_stabsel_selector"
    )
  }
  if (length(selected_index) == 0L) {
    sb_abort(
      "Stability selection did not retain any predictors at cutoff 0.75.",
      "stormbreaker_empty_stabsel_selector"
    )
  }

  selected_features <- source_features[selected_index]
  active_params <- params
  active_params$subsample_q <- subsample_q
  active_params$realized_pfer <- unname(as.numeric(native$PFER))
  structure(
    list(
      id = selector,
      source_features = source_features,
      selected_features = selected_features,
      params = active_params,
      scores = data.frame(
        feature = source_features,
        selection_probability = unname(as.numeric(probabilities)),
        rank = as.integer(rank(-probabilities, ties.method = "min")),
        selected = source_features %in% selected_features,
        stringsAsFactors = FALSE
      )
    ),
    class = "stormbreaker_selector_fit"
  )
}

sb_fit_stabsel_binary_selector <- function(x, outcome) {
  sb_fit_stabsel_selector(
    x,
    list(outcome = outcome),
    "stabsel_binary_lasso"
  )
}

sb_fit_stabsel_cox_selector <- function(x, time, event) {
  sb_fit_stabsel_selector(
    x,
    list(time = time, event = event),
    "stabsel_cox_lasso"
  )
}

sb_apply_high_dimensional_selector <- function(fitted, x) {
  valid <- inherits(fitted, "stormbreaker_selector_fit") &&
    fitted$id %in% sb_high_dimensional_selector_ids() &&
    is.character(fitted$source_features) &&
    is.character(fitted$selected_features) &&
    length(fitted$selected_features) > 0L &&
    !anyNA(fitted$selected_features) &&
    !anyDuplicated(fitted$selected_features) &&
    all(fitted$selected_features %in% fitted$source_features)
  if (!valid) {
    sb_abort(
      "The fitted high-dimensional selector state is invalid.",
      "stormbreaker_bad_selector_fit"
    )
  }
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), fitted$source_features)) {
    sb_abort(
      "Predictor columns do not match the fitted selector schema.",
      "stormbreaker_selector_schema_mismatch"
    )
  }
  x[fitted$selected_features]
}

sb_high_dimensional_selector_adapter <- function(selector) {
  if (!selector %in% sb_high_dimensional_selector_ids()) {
    sb_abort(
      paste0("Unknown high-dimensional selector `", selector, "`."),
      "stormbreaker_unsupported_selector"
    )
  }
  fit <- if (selector %in% sb_information_binary_selector_ids()) {
    function(x, outcome) {
      sb_fit_information_binary_selector(x, outcome, selector)
    }
  } else {
    switch(selector,
      sis_binary_nlogn = sb_fit_sis_binary_selector,
      vsurf_binary_prediction = sb_fit_vsurf_binary_selector,
      mrmr_binary_nlogn = sb_fit_mrmr_binary_selector,
      lpc_binary_nlogn = sb_fit_lpc_binary_selector,
      lpc_cox_nlogn = sb_fit_lpc_cox_selector,
      stabsel_binary_lasso = sb_fit_stabsel_binary_selector,
      stabsel_cox_lasso = sb_fit_stabsel_cox_selector
    )
  }
  list(
    id = selector,
    params = list(),
    needs_outcome = TRUE,
    fit = fit,
    apply = sb_apply_high_dimensional_selector
  )
}

sb_hdda_binary_params <- function(learner = "hdda_binary_cattell") {
  if (!learner %in% sb_hdda_binary_learner_ids()) {
    sb_abort(
      paste0("Unknown HDDA binary learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  }
  list(
    model = "AkjBkQkDk",
    dimension_selection = "Cattell",
    threshold = 0.2,
    scaling = TRUE,
    noise_control = 1e-8,
    probability_mapping = "native_class_one_posterior"
  )
}

sb_hdda_binary_parameter_semantics <- function(
  learner = "hdda_binary_cattell"
) {
  sb_hdda_binary_params(learner)
  paste0(
    "HDDA AkjBkQkDk class-specific subspaces, Cattell dimension rule, ",
    "threshold 0.2, analysis-fold scaling, and native posterior probabilities"
  )
}

sb_fit_hdda_binary <- function(x, outcome, learner, params) {
  sb_require_high_dimensional_method(learner)
  inputs <- sb_validate_high_dimensional_binary_inputs(x, outcome, learner)
  engine <- sb_high_dimensional_native_call(
    HDclassif::hdda(
      data = data.matrix(inputs$x),
      cls = factor(inputs$outcome, levels = c(0L, 1L)),
      model = params$model,
      graph = FALSE,
      d_select = params$dimension_selection,
      threshold = params$threshold,
      show = FALSE,
      scaling = params$scaling,
      noise.ctrl = params$noise_control
    ),
    learner,
    "fit"
  )
  valid <- inherits(engine, "hdc") && identical(engine$K, 2L) &&
    identical(engine$kname, c("0", "1"))
  if (!valid) {
    sb_abort(
      "HDDA did not return a fitted two-class discriminant model.",
      "stormbreaker_invalid_hdda_fit"
    )
  }
  structure(
    list(
      learner = learner,
      features = names(inputs$x),
      params = params,
      engine = engine
    ),
    class = "stormbreaker_hdda_binary"
  )
}

sb_predict_hdda_binary <- function(model, x) {
  if (!inherits(model, "stormbreaker_hdda_binary") ||
    !model$learner %in% sb_hdda_binary_learner_ids()) {
    sb_abort(
      "The fitted HDDA binary model is invalid.",
      "stormbreaker_bad_hdda_model"
    )
  }
  sb_require_high_dimensional_method(model$learner)
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), model$features)) {
    sb_abort(
      "Prediction columns do not match the fitted HDDA model.",
      "stormbreaker_hdda_schema_mismatch"
    )
  }
  prediction <- sb_high_dimensional_native_call(
    stats::predict(model$engine, data = data.matrix(x)),
    model$learner,
    "prediction"
  )
  posterior <- prediction$posterior
  probability <- if (is.matrix(posterior) && "1" %in% colnames(posterior)) {
    unname(as.numeric(posterior[, "1"]))
  } else {
    numeric()
  }
  valid <- length(probability) == nrow(x) && !anyNA(probability) &&
    all(is.finite(probability)) &&
    all(probability >= 0 & probability <= 1)
  if (!valid) {
    sb_abort(
      "HDDA produced invalid class-one posterior probabilities.",
      "stormbreaker_invalid_binary_prediction"
    )
  }
  probability
}

sb_extract_hdda_binary_state <- function(model, context = list()) {
  if (!inherits(model, "stormbreaker_hdda_binary")) {
    sb_bad_learner_state("The fitted HDDA binary state is invalid.")
  }
  structure(
    unclass(model),
    class = c(
      "stormbreaker_hdda_binary_learner_state",
      "stormbreaker_learner_state"
    )
  )
}

sb_predict_hdda_binary_state <- function(state, x) {
  learner <- if (is.list(state)) state[["learner"]] else NULL
  if (!is.character(learner) || length(learner) != 1L || is.na(learner) ||
    !learner %in% sb_hdda_binary_learner_ids()) {
    sb_bad_learner_state("The saved HDDA binary state is invalid.")
  }
  x <- sb_validate_learner_state_x(state, x, learner)
  model <- structure(unclass(state), class = "stormbreaker_hdda_binary")
  list(probability = sb_predict_hdda_binary(model, x))
}

sb_hdda_binary_learner_adapter <- function(
  learner,
  params = sb_hdda_binary_params(learner)
) {
  if (!learner %in% sb_hdda_binary_learner_ids()) {
    sb_abort(
      paste0("Unknown HDDA binary learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  }
  list(
    id = learner,
    params = params,
    fit = function(x, outcome) {
      sb_fit_hdda_binary(x, outcome, learner, params)
    },
    predict = sb_predict_hdda_binary,
    extract_state = sb_extract_hdda_binary_state,
    predict_state = sb_predict_hdda_binary_state
  )
}
