sb_rule_aft_learner_ids <- function() {
  c(
    "bart_aft_survival",
    "pre_binary_rule_ensemble",
    "pre_survival_rule_ensemble",
    "penaft_gehan_elastic_net"
  )
}

sb_pre_rule_ensemble_learner_ids <- function() {
  c("pre_binary_rule_ensemble", "pre_survival_rule_ensemble")
}

sb_rule_aft_survival_learner_ids <- function() {
  c(
    "bart_aft_survival",
    "pre_survival_rule_ensemble",
    "penaft_gehan_elastic_net"
  )
}

sb_rule_aft_params <- function(learner) {
  switch(learner,
    bart_aft_survival = list(
      type = "abart",
      error_distribution = "normal",
      num_trees = 50L,
      num_cutpoints = 100L,
      time_grid_size = 50L,
      posterior_draws = 250L,
      burn_in = 100L,
      keep_every = 1L,
      prior_k = 2,
      tree_power = 2,
      tree_base = 0.95,
      sparse = FALSE,
      seed = 1L,
      cores = 1L
    ),
    pre_binary_rule_ensemble = list(
      family = "binomial",
      ensemble_type = "both",
      sample_fraction = 0.5,
      maximum_depth = 3L,
      learning_rate = 0.01,
      mtry_rule = "ceiling_sqrt_p",
      num_trees = 500L,
      normalize = TRUE,
      standardize = FALSE,
      internal_folds = 5L,
      penalty_rule = "lambda.1se",
      seed = 1L
    ),
    pre_survival_rule_ensemble = list(
      family = "cox",
      ensemble_type = "both",
      sample_fraction = 0.5,
      maximum_depth = 3L,
      learning_rate = 0.01,
      mtry_rule = "ceiling_sqrt_p",
      num_trees = 500L,
      normalize = TRUE,
      standardize = FALSE,
      internal_folds = 5L,
      penalty_rule = "lambda.1se",
      seed = 1L
    ),
    penaft_gehan_elastic_net = list(
      penalty = "EN",
      alpha = 0.5,
      num_lambda = 50L,
      minimum_lambda_ratio = 0.1,
      standardize = TRUE,
      internal_folds = 5L,
      absolute_tolerance = 1e-8,
      relative_tolerance = 2.5e-4,
      maximum_iterations = 10000L,
      selection_rule = "minimum_cv_linear_predictor_score",
      seed = 1L
    ),
    sb_abort(
      paste0("Unknown rule/AFT learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  )
}

sb_rule_aft_parameter_semantics <- function(learner) {
  params <- sb_rule_aft_params(learner)
  switch(learner,
    bart_aft_survival = paste0(
      "normal-error AFT BART with ", params$num_trees,
      " trees, ", params$posterior_draws, " posterior draws, ",
      params$burn_in, " burn-in draws, seed ", params$seed,
      ", and one-core execution fixed before resampling"
    ),
    pre_binary_rule_ensemble = paste0(
      "binary prediction-rule ensemble with linear terms, ",
      params$num_trees, " trees, depth ", params$maximum_depth,
      ", and ceiling-square-root mtry; the rule library is learned from the ",
      "complete analysis fold before seeded ", params$internal_folds,
      "-fold coefficient-only lambda.1se CV on that fixed library"
    ),
    pre_survival_rule_ensemble = paste0(
      "Cox prediction-rule ensemble with linear terms, ",
      params$num_trees, " trees, depth ", params$maximum_depth,
      ", and ceiling-square-root mtry; the rule library is learned from the ",
      "complete analysis fold before seeded ", params$internal_folds,
      "-fold coefficient-only lambda.1se CV on that fixed library"
    ),
    penaft_gehan_elastic_net = paste0(
      "semiparametric Gehan AFT elastic net with alpha ", params$alpha,
      ", ", params$num_lambda, " lambdas, and seeded event-stratified ",
      params$internal_folds,
      "-fold linear-predictor-score selection inside each analysis fold"
    )
  )
}

sb_rule_aft_dependency <- function(learner) {
  packages <- switch(learner,
    bart_aft_survival = "BART",
    pre_binary_rule_ensemble = "pre",
    pre_survival_rule_ensemble = c("pre", "mboost"),
    penaft_gehan_elastic_net = "penAFT",
    sb_abort(
      paste0("Unknown rule/AFT learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  )
  installed <- stats::setNames(
    vapply(packages, requireNamespace, logical(1), quietly = TRUE),
    packages
  )
  list(
    package = packages[[1L]],
    packages = packages,
    installed = installed,
    available = all(installed)
  )
}

sb_require_rule_aft_dependency <- function(learner) {
  dependency <- sb_rule_aft_dependency(learner)
  if (!isTRUE(dependency$available)) {
    missing <- dependency$packages[!dependency$installed]
    sb_abort(
      paste0(
        "Package", if (length(missing) > 1L) "s" else "", " `",
        paste(missing, collapse = "`, `"), "` ",
        if (length(missing) > 1L) "are" else "is",
        " required for learner `", learner, "`."
      ),
      "stormbreaker_missing_rule_aft_dependency"
    )
  }
  invisible(dependency)
}

sb_rule_aft_stop_on_warning <- function(
  code, learner, stage, allow_partial_glmnet_path = FALSE
) {
  withCallingHandlers(
    force(code),
    warning = function(warning) {
      message <- conditionMessage(warning)
      partial_glmnet_path <- isTRUE(allow_partial_glmnet_path) &&
        grepl(
          paste0(
            "^from glmnet C\\+\\+ code \\(error code -[0-9]+\\); ",
            ".*solutions for larger (values of lambda|lambdas) returned$"
          ),
          message
        )
      if (partial_glmnet_path) {
        invokeRestart("muffleWarning")
      }
      sb_abort(
        paste0(
          "Learner `", learner, "` raised a ", stage, " warning: ",
          message
        ),
        "stormbreaker_rule_aft_warning",
        parent = warning
      )
    }
  )
}

sb_new_rule_aft_model <- function(learner, features, params, ...) {
  structure(
    c(
      list(learner = learner, features = features, params = params),
      list(...)
    ),
    class = "stormbreaker_rule_aft_model"
  )
}

sb_validate_rule_aft_model_x <- function(model, x) {
  valid <- inherits(model, "stormbreaker_rule_aft_model") &&
    is.list(model) && model$learner %in% sb_rule_aft_learner_ids() &&
    is.character(model$features) && length(model$features) > 0L &&
    !anyNA(model$features) && !anyDuplicated(model$features)
  if (!isTRUE(valid)) {
    sb_abort(
      "The supplied rule/AFT model is incomplete or invalid.",
      "stormbreaker_bad_rule_aft_fit"
    )
  }
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), model$features)) {
    sb_abort(
      "Prediction columns do not match the fitted rule/AFT learner schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  x
}

sb_compact_pre_engine <- function(engine) {
  engine[c("call", "weights", "data", "modmat")] <- NULL
  environment(engine$formula) <- baseenv()
  engine
}

sb_fit_pre_rule_ensemble <- function(x, endpoint, learner, params) {
  if (identical(learner, "pre_binary_rule_ensemble")) {
    inputs <- sb_validate_classical_binary_inputs(x, endpoint, learner)
    formula_data <- sb_classical_binary_formula_data(
      inputs$x,
      inputs$outcome
    )
  } else {
    inputs <- sb_validate_additional_survival_inputs(
      x,
      endpoint$time,
      endpoint$event,
      learner
    )
    formula_data <- sb_additional_survival_formula_data(
      inputs$x,
      inputs$time,
      inputs$event
    )
  }
  sb_require_rule_aft_dependency(learner)
  effective_mtry <- as.integer(ceiling(sqrt(ncol(inputs$x))))
  engine <- NULL
  invisible(utils::capture.output(
    engine <- withr::with_seed(
      params$seed,
      sb_rule_aft_stop_on_warning(
        pre::pre(
          formula = formula_data$formula,
          data = formula_data$data,
          family = params$family,
          type = params$ensemble_type,
          sampfrac = params$sample_fraction,
          maxdepth = params$maximum_depth,
          learnrate = params$learning_rate,
          mtry = effective_mtry,
          ntrees = params$num_trees,
          normalize = params$normalize,
          standardize = params$standardize,
          nfolds = params$internal_folds,
          tree.unbiased = TRUE,
          removecomplements = TRUE,
          removeduplicates = TRUE,
          verbose = FALSE
        ),
        learner,
        "fitting",
        allow_partial_glmnet_path = TRUE
      )
    ),
    type = "output"
  ))
  if (!inherits(engine, "pre")) {
    sb_abort(
      "The PRE fit did not return a prediction-rule ensemble.",
      "stormbreaker_incomplete_pre_fit"
    )
  }
  sb_new_rule_aft_model(
    learner,
    names(inputs$x),
    params,
    engine = sb_compact_pre_engine(engine),
    effective_mtry = effective_mtry
  )
}

sb_fit_bart_aft <- function(x, time, event, learner, params) {
  inputs <- sb_validate_additional_survival_inputs(x, time, event, learner)
  sb_require_rule_aft_dependency(learner)
  engine <- NULL
  invisible(utils::capture.output(
    engine <- withr::with_seed(
      params$seed,
      sb_rule_aft_stop_on_warning(
        BART::abart(
          x.train = as.matrix(inputs$x),
          times = inputs$time,
          delta = inputs$event,
          K = params$time_grid_size,
          type = params$type,
          ntype = 1L,
          sparse = params$sparse,
          k = params$prior_k,
          power = params$tree_power,
          base = params$tree_base,
          ntree = params$num_trees,
          numcut = params$num_cutpoints,
          ndpost = params$posterior_draws,
          nskip = params$burn_in,
          keepevery = params$keep_every,
          printevery = .Machine$integer.max,
          mc.cores = params$cores,
          seed = params$seed
        ),
        learner,
        "fitting"
      )
    ),
    type = "output"
  ))
  post_sigma <- utils::tail(
    as.numeric(engine$sigma),
    params$posterior_draws
  )
  kept_columns <- as.integer(engine$rm.const)
  if (length(kept_columns) == 0L) {
    kept_columns <- seq_len(ncol(inputs$x))
  }
  valid <- inherits(engine, "abart") && is.list(engine$treedraws) &&
    is.numeric(engine$offset) && length(engine$offset) == 1L &&
    is.finite(engine$offset) &&
    length(post_sigma) == params$posterior_draws &&
    !anyNA(post_sigma) && all(is.finite(post_sigma)) &&
    all(post_sigma > 0) && length(kept_columns) > 0L &&
    !anyNA(kept_columns) && !anyDuplicated(kept_columns) &&
    all(kept_columns >= 1L & kept_columns <= ncol(inputs$x))
  if (!isTRUE(valid)) {
    sb_abort(
      "AFT-BART did not return a portable posterior prediction state.",
      "stormbreaker_incomplete_bart_aft_fit"
    )
  }
  compact_engine <- structure(
    list(
      treedraws = engine$treedraws,
      offset = unname(as.numeric(engine$offset)),
      sigma = unname(post_sigma),
      kept_columns = kept_columns
    ),
    class = "stormbreaker_bart_aft_engine"
  )
  sb_new_rule_aft_model(
    learner,
    names(inputs$x),
    params,
    engine = compact_engine
  )
}

sb_fit_penaft <- function(x, time, event, learner, params) {
  inputs <- sb_validate_additional_survival_inputs(x, time, event, learner)
  sb_require_rule_aft_dependency(learner)
  fold_id <- sb_make_stratified_fold_id(
    inputs$event,
    params$internal_folds,
    params$seed,
    stratum_label = "event"
  )
  cv_index <- lapply(seq_len(params$internal_folds), function(fold) {
    which(fold_id == fold)
  })
  fit <- NULL
  invisible(utils::capture.output(
    fit <- withr::with_seed(
      params$seed,
      sb_rule_aft_stop_on_warning(
        penAFT::penAFT.cv(
          X = as.matrix(inputs$x),
          logY = log(inputs$time),
          delta = inputs$event,
          nlambda = params$num_lambda,
          lambda.ratio.min = params$minimum_lambda_ratio,
          penalty = params$penalty,
          alpha = params$alpha,
          tol.abs = params$absolute_tolerance,
          tol.rel = params$relative_tolerance,
          standardize = params$standardize,
          cv.index = cv_index,
          admm.max.iter = params$maximum_iterations,
          quiet = TRUE
        ),
        learner,
        "fitting"
      )
    ),
    type = "output"
  ))
  cv_error <- as.numeric(fit$cv.err.linPred)
  selected_index <- if (length(cv_error) > 0L &&
    !anyNA(cv_error) && all(is.finite(cv_error))) {
    min(which(cv_error == min(cv_error)))
  } else {
    integer()
  }
  full <- fit$full.fit
  valid <- inherits(fit, "penAFT.cv") && length(selected_index) == 1L &&
    is.numeric(full$lambda) && length(full$lambda) == length(cv_error) &&
    !anyNA(full$lambda) && all(is.finite(full$lambda)) &&
    nrow(full$beta) == ncol(inputs$x) &&
    ncol(full$beta) == length(full$lambda)
  if (!isTRUE(valid)) {
    sb_abort(
      "penAFT did not return a complete cross-validated coefficient path.",
      "stormbreaker_incomplete_penaft_fit"
    )
  }
  coefficient <- as.numeric(as.matrix(full$beta[, selected_index]))
  if (isTRUE(full$standardize)) {
    coefficient <- coefficient / as.numeric(full$X.sd)
    intercept <- -sum(as.numeric(full$X.mean) * coefficient)
  } else {
    intercept <- 0
  }
  if (length(coefficient) != ncol(inputs$x) || anyNA(coefficient) ||
    any(!is.finite(coefficient)) || !is.finite(intercept)) {
    sb_abort(
      "penAFT did not return portable selected coefficients.",
      "stormbreaker_incomplete_penaft_fit"
    )
  }
  sb_new_rule_aft_model(
    learner,
    names(inputs$x),
    params,
    intercept = unname(intercept),
    coefficient = unname(coefficient),
    selected_lambda = unname(as.numeric(full$lambda[[selected_index]])),
    selected_cv_error = unname(cv_error[[selected_index]]),
    diagnostic_path = list(full_fit = full, cv_error = cv_error,
      selected_index = selected_index)
  )
}

sb_predict_pre_rule_ensemble <- function(model, x) {
  x <- sb_validate_rule_aft_model_x(model, x)
  sb_require_rule_aft_dependency(model$learner)
  prediction <- NULL
  prediction_type <- if (identical(
    model$learner,
    "pre_binary_rule_ensemble"
  )) {
    "response"
  } else {
    "link"
  }
  invisible(utils::capture.output(
    prediction <- sb_rule_aft_stop_on_warning(
      stats::predict(
        model$engine,
        newdata = x,
        type = prediction_type,
        penalty.par.val = model$params$penalty_rule
      ),
      model$learner,
      "prediction"
    ),
    type = "output"
  ))
  prediction <- unname(as.numeric(prediction))
  valid <- length(prediction) == nrow(x) && !anyNA(prediction) &&
    all(is.finite(prediction))
  if (identical(model$learner, "pre_binary_rule_ensemble")) {
    valid <- valid && all(prediction >= 0 & prediction <= 1)
  }
  if (!isTRUE(valid)) {
    sb_abort(
      "PRE produced incomplete or invalid predictions.",
      "stormbreaker_invalid_pre_prediction"
    )
  }
  prediction
}

sb_bart_aft_log_time_draws <- function(model, x) {
  x <- sb_validate_rule_aft_model_x(model, x)
  sb_require_rule_aft_dependency(model$learner)
  engine <- model$engine
  matrix_x <- as.matrix(x)[, engine$kept_columns, drop = FALSE]
  draws <- NULL
  invisible(utils::capture.output(
    draws <- sb_rule_aft_stop_on_warning(
      BART::pwbart(
        matrix_x,
        engine$treedraws,
        mu = engine$offset,
        mc.cores = model$params$cores,
        dodraws = TRUE
      ),
      model$learner,
      "prediction"
    ),
    type = "output"
  ))
  draws <- as.matrix(draws)
  valid <- nrow(draws) == length(engine$sigma) &&
    ncol(draws) == nrow(x) && !anyNA(draws) && all(is.finite(draws))
  if (!isTRUE(valid)) {
    sb_abort(
      "AFT-BART produced incomplete posterior log-time draws.",
      "stormbreaker_invalid_bart_aft_prediction"
    )
  }
  draws
}

sb_penaft_log_time_link <- function(model, x) {
  x <- sb_validate_rule_aft_model_x(model, x)
  link <- unname(
    model$intercept + as.vector(as.matrix(x) %*% model$coefficient)
  )
  if (length(link) != nrow(x) || anyNA(link) || any(!is.finite(link))) {
    sb_abort(
      "penAFT produced incomplete predicted log times.",
      "stormbreaker_invalid_penaft_prediction"
    )
  }
  link
}

sb_predict_rule_aft <- function(model, x) {
  learner <- model$learner
  if (learner %in% sb_pre_rule_ensemble_learner_ids()) {
    return(sb_predict_pre_rule_ensemble(model, x))
  }
  if (identical(learner, "bart_aft_survival")) {
    return(-colMeans(sb_bart_aft_log_time_draws(model, x)))
  }
  if (identical(learner, "penaft_gehan_elastic_net")) {
    return(-sb_penaft_log_time_link(model, x))
  }
  sb_abort(
    "The supplied rule/AFT model is not supported.",
    "stormbreaker_bad_rule_aft_fit"
  )
}

sb_bart_aft_survival_probability <- function(model, x, prediction_horizon) {
  draws <- sb_bart_aft_log_time_draws(model, x)
  standardized <- (log(prediction_horizon) - draws) / model$engine$sigma
  probability <- unname(colMeans(stats::pnorm(
    standardized,
    lower.tail = FALSE
  )))
  if (length(probability) != ncol(draws) || anyNA(probability) ||
    any(!is.finite(probability)) || any(probability < 0 | probability > 1)) {
    sb_abort(
      "AFT-BART produced invalid horizon survival probabilities.",
      "stormbreaker_invalid_bart_aft_prediction"
    )
  }
  probability
}

sb_rule_aft_direct_survival_source <- function(learner) {
  switch(learner,
    bart_aft_survival = "BART_abart_normal_posterior_survival_direct",
    penaft_gehan_elastic_net =
      "penAFT_gehan_analysis_residual_kaplan_meier",
    NULL
  )
}

sb_predict_rule_aft_survival <- function(
  model, x_analysis, time, event, x_assessment, prediction_horizon
) {
  inputs <- sb_validate_additional_survival_inputs(
    x_analysis,
    time,
    event,
    model$learner
  )
  if (!identical(names(inputs$x), model$features)) {
    sb_abort(
      "Analysis columns do not match the fitted rule/AFT learner schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  prediction_horizon <- sb_validate_prediction_horizon(
    prediction_horizon,
    inputs$time
  )
  if (identical(model$learner, "bart_aft_survival")) {
    return(list(
      prediction_horizon = prediction_horizon,
      source = sb_rule_aft_direct_survival_source(model$learner),
      ties = "not_applicable",
      probability = sb_bart_aft_survival_probability(
        model,
        x_assessment,
        prediction_horizon
      )
    ))
  }
  if (identical(model$learner, "pre_survival_rule_ensemble")) {
    analysis_risk <- sb_predict_pre_rule_ensemble(model, inputs$x)
    assessment_risk <- sb_predict_pre_rule_ensemble(model, x_assessment)
    baseline <- sb_breslow_baseline_hazard(
      inputs$time,
      inputs$event,
      analysis_risk,
      prediction_horizon
    )
    return(sb_survival_prediction_result(
      prediction_horizon,
      baseline,
      source = "pre_rule_ensemble_analysis_breslow",
      ties = "breslow",
      risk = assessment_risk
    ))
  }
  analysis_link <- sb_penaft_log_time_link(model, inputs$x)
  assessment_link <- sb_penaft_log_time_link(model, x_assessment)
  curve <- sb_fit_mboost_gehan_residual_survival(
    inputs$time,
    inputs$event,
    analysis_link
  )
  list(
    prediction_horizon = prediction_horizon,
    source = sb_rule_aft_direct_survival_source(model$learner),
    ties = "not_applicable",
    probability = sb_predict_mboost_gehan_residual_survival(
      curve,
      assessment_link,
      prediction_horizon
    ),
    residual_survival = curve
  )
}

sb_extract_rule_aft_coefficients <- function(model) {
  if (!identical(model$learner, "penaft_gehan_elastic_net")) {
    sb_abort(
      "Only the penAFT learner exposes portable linear coefficients.",
      "stormbreaker_rule_aft_coefficients_not_available"
    )
  }
  values <- model$coefficient
  data.frame(
    feature = model$features,
    coefficient = values,
    nonzero = values != 0,
    direction = ifelse(
      values > 0,
      "longer_predicted_survival",
      ifelse(values < 0, "shorter_predicted_survival", "zero")
    ),
    stringsAsFactors = FALSE
  )
}

sb_extract_rule_aft_state <- function(model, context = list()) {
  learner <- model$learner
  is_pre <- learner %in% sb_pre_rule_ensemble_learner_ids()
  state <- list(
    learner = learner,
    features = model$features,
    params = model$params,
    baseline_hazard = NULL,
    survival_probability = NULL
  )
  if (is_pre || identical(learner, "bart_aft_survival")) {
    state$engine_raw <- serialize(
      model$engine,
      connection = NULL,
      ascii = FALSE,
      xdr = TRUE,
      version = 3L
    )
  }
  if (is_pre) {
    state$effective_mtry <- model$effective_mtry
  }
  if (identical(learner, "penaft_gehan_elastic_net")) {
    state$intercept <- model$intercept
    state$coefficient <- model$coefficient
    state$selected_lambda <- model$selected_lambda
    state$selected_cv_error <- model$selected_cv_error
  }
  if (identical(learner, "pre_binary_rule_ensemble")) {
    state$prediction_reference <- "pre_native_rule_ensemble"
    state$prediction_formula <- "positive_class_probability"
  } else {
    risk <- switch(learner,
      bart_aft_survival = list(
        reference = "BART_abart_normal_native",
        formula = "negative_posterior_mean_log_time"
      ),
      pre_survival_rule_ensemble = list(
        reference = "pre_native_rule_ensemble",
        formula = "pre_cox_link_prediction"
      ),
      penaft_gehan_elastic_net = list(
        reference = "penAFT_native",
        formula = "negative_centered_elastic_net_log_time_prediction"
      )
    )
    state$risk_reference <- risk$reference
    state$risk_formula <- risk$formula
  }
  if (learner %in% sb_rule_aft_survival_learner_ids() &&
    "prediction_horizon" %in% names(context)) {
    probability <- sb_predict_rule_aft_survival(
      model,
      context$x,
      context$time,
      context$event,
      context$x[1L, , drop = FALSE],
      context$prediction_horizon
    )
    time_unit <- if ("time_unit" %in% names(context)) {
      context$time_unit
    } else {
      NA_character_
    }
    if (identical(learner, "pre_survival_rule_ensemble")) {
      state$baseline_hazard <- list(
        prediction_horizon = probability$prediction_horizon,
        time_unit = time_unit,
        cumulative_baseline_hazard = probability$cumulative_baseline_hazard,
        baseline_survival = probability$baseline_survival,
        source = probability$source,
        ties = probability$ties
      )
    } else {
      state$survival_probability <- list(
        prediction_horizon = probability$prediction_horizon,
        time_unit = time_unit,
        source = probability$source,
        ties = probability$ties
      )
      if (identical(learner, "penaft_gehan_elastic_net")) {
        state$survival_probability$residual_survival <-
          probability$residual_survival
      }
    }
  }
  structure(
    state,
    class = c("stormbreaker_rule_aft_state", "stormbreaker_learner_state")
  )
}

sb_restore_rule_aft_model <- function(state) {
  learner <- state$learner
  if (learner %in% sb_pre_rule_ensemble_learner_ids() ||
    identical(learner, "bart_aft_survival")) {
    if (!is.raw(state$engine_raw) || length(state$engine_raw) == 0L) {
      sb_bad_learner_state("The saved rule/AFT engine is invalid.")
    }
    engine <- tryCatch(
      unserialize(state$engine_raw),
      error = function(error) {
        sb_abort(
          paste0(
            "The saved `", learner, "` engine could not be loaded: ",
            conditionMessage(error)
          ),
          "stormbreaker_invalid_rule_aft_state",
          parent = error
        )
      }
    )
    valid_engine <- if (learner %in% sb_pre_rule_ensemble_learner_ids()) {
      inherits(engine, "pre")
    } else {
      inherits(engine, "stormbreaker_bart_aft_engine")
    }
    if (!isTRUE(valid_engine)) {
      sb_bad_learner_state("The saved rule/AFT engine is invalid.")
    }
    return(sb_new_rule_aft_model(
      learner,
      state$features,
      state$params,
      engine = engine,
      effective_mtry = state$effective_mtry
    ))
  }
  valid <- is.numeric(state$intercept) && length(state$intercept) == 1L &&
    is.finite(state$intercept) && is.numeric(state$coefficient) &&
    length(state$coefficient) == length(state$features) &&
    !anyNA(state$coefficient) && all(is.finite(state$coefficient))
  if (!isTRUE(valid)) {
    sb_bad_learner_state("The saved penAFT coefficients are invalid.")
  }
  sb_new_rule_aft_model(
    learner,
    state$features,
    state$params,
    intercept = state$intercept,
    coefficient = state$coefficient,
    selected_lambda = state$selected_lambda,
    selected_cv_error = state$selected_cv_error
  )
}

sb_validate_rule_aft_horizon_state <- function(state, learner) {
  probability <- if (identical(learner, "pre_survival_rule_ensemble")) {
    state$baseline_hazard
  } else {
    state$survival_probability
  }
  valid <- is.list(probability) &&
    is.numeric(probability$prediction_horizon) &&
    length(probability$prediction_horizon) == 1L &&
    is.finite(probability$prediction_horizon) &&
    probability$prediction_horizon > 0 &&
    is.character(probability$time_unit) &&
    length(probability$time_unit) == 1L &&
    !is.na(probability$time_unit) && nzchar(probability$time_unit)
  if (!isTRUE(valid)) {
    sb_bad_learner_state("The saved rule/AFT horizon state is invalid.")
  }
  probability
}

sb_predict_rule_aft_state <- function(state, x) {
  learner <- if (is.list(state)) state$learner else NULL
  if (!inherits(state, "stormbreaker_rule_aft_state") ||
    !is.character(learner) || length(learner) != 1L || is.na(learner) ||
    !learner %in% sb_rule_aft_learner_ids()) {
    sb_bad_learner_state("The saved rule/AFT learner state is invalid.")
  }
  x <- sb_validate_learner_state_x(state, x, learner)
  model <- sb_restore_rule_aft_model(state)
  prediction <- sb_predict_rule_aft(model, x)
  if (identical(learner, "pre_binary_rule_ensemble")) {
    return(list(probability = prediction))
  }
  result <- list(risk = prediction)
  if (!is.null(state$baseline_hazard)) {
    baseline <- sb_validate_rule_aft_horizon_state(state, learner)
    result$prediction_horizon <- baseline$prediction_horizon
    result$time_unit <- baseline$time_unit
    result$survival_probability <- sb_survival_probability_from_hazard(
      prediction,
      baseline$cumulative_baseline_hazard
    )
  } else if (!is.null(state$survival_probability)) {
    probability <- sb_validate_rule_aft_horizon_state(state, learner)
    result$prediction_horizon <- probability$prediction_horizon
    result$time_unit <- probability$time_unit
    if (identical(learner, "bart_aft_survival")) {
      result$survival_probability <- sb_bart_aft_survival_probability(
        model,
        x,
        probability$prediction_horizon
      )
    } else {
      result$survival_probability <-
        sb_predict_mboost_gehan_residual_survival(
          probability$residual_survival,
          -prediction,
          probability$prediction_horizon
        )
    }
  }
  result
}

sb_rule_aft_learner_adapter <- function(
  learner,
  params = sb_rule_aft_params(learner)
) {
  if (!learner %in% sb_rule_aft_learner_ids()) {
    sb_abort(
      paste0("Unknown rule/AFT learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  }
  fit <- if (identical(learner, "pre_binary_rule_ensemble")) {
    function(x, outcome) {
      sb_fit_pre_rule_ensemble(x, outcome, learner, params)
    }
  } else if (identical(learner, "pre_survival_rule_ensemble")) {
    function(x, time, event) {
      sb_fit_pre_rule_ensemble(
        x,
        list(time = time, event = event),
        learner,
        params
      )
    }
  } else if (identical(learner, "bart_aft_survival")) {
    function(x, time, event) {
      sb_fit_bart_aft(x, time, event, learner, params)
    }
  } else {
    function(x, time, event) {
      sb_fit_penaft(x, time, event, learner, params)
    }
  }
  adapter <- list(
    id = learner,
    params = params,
    fit = fit,
    predict = sb_predict_rule_aft,
    extract_state = sb_extract_rule_aft_state,
    predict_state = sb_predict_rule_aft_state
  )
  if (learner %in% sb_rule_aft_survival_learner_ids()) {
    adapter$predict_survival <- sb_predict_rule_aft_survival
  }
  if (identical(learner, "penaft_gehan_elastic_net")) {
    adapter$extract_coefficients <- sb_extract_rule_aft_coefficients
  }
  adapter
}
