sb_survivalsvm_linear_learner_ids <- function() {
  c(
    "survivalsvm",
    "survivalsvm_vanbelle1",
    "survivalsvm_vanbelle2",
    "survivalsvm_hybrid"
  )
}

sb_survivalsvm_additive_learner_ids <- function() {
  c(
    "survivalsvm_additive",
    "survivalsvm_vanbelle1_additive",
    "survivalsvm_vanbelle2_additive",
    "survivalsvm_hybrid_additive"
  )
}

sb_survivalsvm_rbf_learner_ids <- function() {
  c(
    "survivalsvm_rbf",
    "survivalsvm_vanbelle1_rbf",
    "survivalsvm_vanbelle2_rbf",
    "survivalsvm_hybrid_rbf"
  )
}

sb_survivalsvm_kernel_learner_ids <- function() {
  c(
    sb_survivalsvm_additive_learner_ids(),
    sb_survivalsvm_rbf_learner_ids()
  )
}

sb_survivalsvm_learner_ids <- function() {
  c(
    sb_survivalsvm_linear_learner_ids(),
    sb_survivalsvm_kernel_learner_ids()
  )
}

sb_survivalsvm_type <- function(learner) {
  type <- switch(learner,
    survivalsvm = "regression",
    survivalsvm_vanbelle1 = "vanbelle1",
    survivalsvm_vanbelle2 = "vanbelle2",
    survivalsvm_hybrid = "hybrid",
    survivalsvm_additive = "regression",
    survivalsvm_vanbelle1_additive = "vanbelle1",
    survivalsvm_vanbelle2_additive = "vanbelle2",
    survivalsvm_hybrid_additive = "hybrid",
    survivalsvm_rbf = "regression",
    survivalsvm_vanbelle1_rbf = "vanbelle1",
    survivalsvm_vanbelle2_rbf = "vanbelle2",
    survivalsvm_hybrid_rbf = "hybrid",
    NULL
  )
  if (is.null(type)) {
    sb_abort(
      paste0("Unknown fixed survival-SVM learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  }
  type
}

sb_survivalsvm_kernel <- function(learner) {
  if (learner %in% sb_survivalsvm_linear_learner_ids()) {
    return("lin_kernel")
  }
  if (learner %in% sb_survivalsvm_additive_learner_ids()) {
    return("add_kernel")
  }
  if (learner %in% sb_survivalsvm_rbf_learner_ids()) {
    return("rbf_kernel")
  }
  sb_abort(
    paste0("Unknown fixed survival-SVM learner `", learner, "`."),
    "stormbreaker_unsupported_learner"
  )
}

sb_classical_survival_learner_ids <- function() {
  c(
    "coxboost",
    "plsrcox",
    "superpc",
    sb_survivalsvm_linear_learner_ids()
  )
}

sb_analysis_tuned_survival_learner_ids <- function() {
  c("coxboost_cv", "plsrcox_cv", "superpc_cv", "gbm_cox_cv")
}

sb_classical_survival_linear_learner_ids <- function() {
  c(
    sb_classical_survival_learner_ids(),
    setdiff(sb_analysis_tuned_survival_learner_ids(), "gbm_cox_cv")
  )
}

sb_linear_survival_learner_ids <- function() {
  c(
    "coxph",
    sb_coxphf_learner_ids(),
    "glmnet_cox",
    sb_classical_survival_linear_learner_ids(),
    sb_abess_survival_learner_ids(),
    sb_ncvreg_survival_learner_ids()
  )
}

sb_survival_learner_ids <- function() {
  c(
    sb_linear_survival_learner_ids(),
    sb_survivalsvm_kernel_learner_ids(),
    "catboost_cox",
    "catboost_aft_survival",
    "gbm_cox",
    "gbm_cox_cv",
    sb_mboost_cox_learner_ids(),
    sb_mboost_aft_learner_ids(),
    sb_mboost_ranking_learner_ids(),
    sb_ranger_survival_learner_ids(),
    sb_rlt_survival_learner_ids(),
    sb_bagged_knn_survival_learner_ids(), sb_akritas_learner_ids(),
    sb_ahaz_survival_learner_ids(),
    sb_bujar_linear_learner_ids(),
    sb_rule_aft_survival_learner_ids(),
    sb_xgboost_survival_learner_ids(),
    sb_flexsurv_parametric_learner_ids(),
    sb_flexsurvspline_learner_ids(),
    sb_mtlr_survival_learner_ids(),
    sb_pamr_survival_learner_ids(),
    sb_additional_survival_learner_ids(),
    sb_timereg_survival_learner_ids(),
    sb_mgcv_gam_survival_learner_ids()
  )
}

sb_coxboost_available <- function() {
  requireNamespace("CoxBoost", quietly = TRUE)
}

sb_plsrcox_available <- function() {
  requireNamespace("plsRcox", quietly = TRUE)
}

sb_superpc_available <- function() {
  requireNamespace("superpc", quietly = TRUE)
}

sb_survivalsvm_available <- function() {
  requireNamespace("survivalsvm", quietly = TRUE)
}

sb_validate_classical_survival_inputs <- function(x, time, event, learner) {
  x <- sb_validate_predictor_frame(x, "x")
  time <- sb_validate_endpoint_vector(time, "time", nrow(x))
  event <- sb_validate_endpoint_vector(event, "event", nrow(x))
  if (any(time < 0)) {
    sb_abort(
      "`time` must contain non-negative values.",
      "stormbreaker_bad_adapter_time"
    )
  }
  if (!all(event %in% c(0, 1))) {
    sb_abort(
      "`event` must contain only 0/1 values.",
      "stormbreaker_bad_adapter_event"
    )
  }
  if (!any(event == 1L)) {
    sb_abort(
      paste0("A fixed ", learner, " analysis fold must contain an event."),
      "stormbreaker_classical_survival_no_event"
    )
  }
  list(x = x, time = as.numeric(time), event = as.integer(event))
}

sb_validate_positive_number <- function(value, argument, subclass) {
  valid <- is.numeric(value) && length(value) == 1L &&
    !is.na(value) && is.finite(value) && value > 0
  if (!valid) {
    sb_abort(
      paste0("`", argument, "` must be one finite positive number."),
      subclass
    )
  }
  as.numeric(value)
}

sb_classical_survival_stop_on_warning <- function(code, label, subclass) {
  withCallingHandlers(
    force(code),
    warning = function(warning) {
      sb_abort(
        paste0(
          "The ", label, " raised a warning: ",
          conditionMessage(warning)
        ),
        subclass,
        parent = warning
      )
    }
  )
}

sb_coxboost_params <- function(stepno = 100L, penalty = NULL) {
  if (!is.null(penalty)) {
    penalty <- sb_validate_positive_number(
      penalty,
      "coxboost_penalty",
      "stormbreaker_bad_coxboost_penalty"
    )
  }
  list(
    stepno = sb_check_positive_integer(
      stepno,
      "coxboost_stepno",
      minimum = 1L
    ),
    penalty = penalty,
    penalty_rule = if (is.null(penalty)) {
      "9_times_analysis_event_count"
    } else {
      "fixed_numeric"
    },
    standardize = TRUE,
    criterion = "pscore"
  )
}

sb_coxboost_cv_params <- function(inner_folds = 10L, seed = 1L,
                                  maximum_steps = 500L,
                                  start_penalty = 500) {
  list(
    inner_folds = sb_check_positive_integer(
      inner_folds,
      "learner_tuning_folds",
      minimum = 3L
    ),
    seed = sb_check_positive_integer(seed, "learner_seed", minimum = 1L),
    maximum_steps = sb_check_positive_integer(
      maximum_steps,
      "maximum_steps",
      minimum = 1L
    ),
    start_penalty = sb_validate_positive_number(
      start_penalty,
      "start_penalty",
      "stormbreaker_bad_coxboost_penalty"
    ),
    standardize = TRUE,
    tuning_type = "verweij"
  )
}

sb_plsrcox_params <- function(components = 2L) {
  list(
    components = sb_check_positive_integer(
      components,
      "plsrcox_components",
      minimum = 1L
    ),
    scale_x = TRUE,
    scale_y = FALSE,
    sparse_stop = FALSE
  )
}

sb_plsrcox_cv_params <- function(inner_folds = 5L, seed = 1L,
                                 maximum_components = 10L) {
  list(
    inner_folds = sb_check_positive_integer(
      inner_folds,
      "learner_tuning_folds",
      minimum = 3L
    ),
    seed = sb_check_positive_integer(seed, "learner_seed", minimum = 1L),
    maximum_components = sb_check_positive_integer(
      maximum_components,
      "maximum_components",
      minimum = 1L
    ),
    selection_rule = "minimum_cv_error",
    scale_x = TRUE,
    scale_y = FALSE,
    sparse_stop = FALSE
  )
}

sb_superpc_params <- function(threshold = 1, components = 1L) {
  valid_threshold <- is.numeric(threshold) && length(threshold) == 1L &&
    !is.na(threshold) && is.finite(threshold) && threshold >= 0
  if (!valid_threshold) {
    sb_abort(
      "`superpc_threshold` must be one finite non-negative number.",
      "stormbreaker_bad_superpc_threshold"
    )
  }
  components <- sb_check_positive_integer(
    components,
    "superpc_components",
    minimum = 1L
  )
  if (components > 3L) {
    sb_abort(
      "`superpc_components` must be an integer from 1 through 3.",
      "stormbreaker_bad_superpc_components"
    )
  }
  list(
    threshold = as.numeric(threshold),
    components = components,
    s0_perc = 0.5
  )
}

sb_superpc_cv_params <- function(inner_folds = 10L, seed = 1L) {
  list(
    inner_folds = sb_check_positive_integer(
      inner_folds,
      "learner_tuning_folds",
      minimum = 3L
    ),
    seed = sb_check_positive_integer(seed, "learner_seed", minimum = 1L),
    threshold_count = 20L,
    tuning_components = 3L,
    prediction_components = 1L,
    minimum_features = 2L,
    s0_perc = 0.5
  )
}

sb_survivalsvm_params <- function(
  gamma_mu = 1,
  learner = "survivalsvm"
) {
  type <- sb_survivalsvm_type(learner)
  params <- list(
    gamma_mu = sb_validate_positive_number(
      gamma_mu,
      "survivalsvm_gamma_mu",
      "stormbreaker_bad_survivalsvm_gamma_mu"
    ),
    type = type,
    kernel = sb_survivalsvm_kernel(learner),
    opt_meth = "quadprog",
    standardize = TRUE,
    calibration_ties = "breslow"
  )
  if (identical(params$kernel, "rbf_kernel")) {
    params["kernel_pars"] <- list(NULL)
    params$kernel_parameter_rule <-
      "analysis_median_squared_distance"
  }
  if (!identical(type, "regression")) {
    params <- append(params, list(diff_meth = "makediff3"), after = 2L)
  }
  params
}

sb_new_classical_survival_model <- function(
  learner, features, coefficient, params, audit, engine,
  baseline_source
) {
  valid <- learner %in% sb_classical_survival_linear_learner_ids() &&
    is.character(features) && length(features) > 0L &&
    !anyNA(features) && !anyDuplicated(features) &&
    is.numeric(coefficient) && length(coefficient) == length(features) &&
    !anyNA(coefficient) && all(is.finite(coefficient))
  if (!valid) {
    sb_abort(
      paste0("Fixed ", learner, " did not yield portable coefficients."),
      "stormbreaker_invalid_classical_survival_fit"
    )
  }
  structure(
    list(
      learner = learner,
      features = features,
      coefficient = unname(coefficient),
      params = params,
      audit = audit,
      baseline_source = baseline_source,
      engine = engine
    ),
    class = c(
      paste0("stormbreaker_", learner),
      "stormbreaker_classical_survival_fixed"
    )
  )
}

sb_predict_classical_survival <- function(model, x) {
  if (!inherits(model, "stormbreaker_classical_survival_fixed")) {
    sb_abort(
      "The supplied model is not a fitted classical survival learner.",
      "stormbreaker_bad_classical_survival_fit"
    )
  }
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), model$features)) {
    sb_abort(
      "Prediction columns do not match the fitted survival learner schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  risk <- unname(as.vector(as.matrix(x) %*% model$coefficient))
  if (length(risk) != nrow(x) || anyNA(risk) || any(!is.finite(risk))) {
    sb_abort(
      "The fitted survival learner produced invalid risk predictions.",
      "stormbreaker_invalid_classical_survival_prediction"
    )
  }
  risk
}

sb_extract_classical_survival_coefficients <- function(model) {
  if (!inherits(model, "stormbreaker_classical_survival_fixed")) {
    sb_abort(
      "The supplied model is not a fitted classical survival learner.",
      "stormbreaker_bad_classical_survival_fit"
    )
  }
  data.frame(
    feature = model$features,
    coefficient = model$coefficient,
    nonzero = model$coefficient != 0,
    direction = ifelse(
      model$coefficient > 0,
      "higher_hazard",
      ifelse(model$coefficient < 0, "lower_hazard", "zero")
    ),
    stringsAsFactors = FALSE
  )
}

sb_predict_classical_survival_survival <- function(
  model, x_analysis, time, event, x_assessment, prediction_horizon
) {
  inputs <- sb_validate_classical_survival_inputs(
    x_analysis,
    time,
    event,
    model$learner
  )
  prediction_horizon <- sb_validate_prediction_horizon(
    prediction_horizon,
    inputs$time
  )
  analysis_risk <- sb_predict_classical_survival(model, inputs$x)
  assessment_risk <- sb_predict_classical_survival(
    model,
    x_assessment
  )
  cumulative_hazard <- sb_breslow_baseline_hazard(
    inputs$time,
    inputs$event,
    analysis_risk,
    prediction_horizon
  )
  sb_survival_prediction_result(
    prediction_horizon,
    cumulative_hazard,
    source = model$baseline_source,
    ties = "breslow",
    risk = assessment_risk
  )
}

sb_fit_coxboost_cv_native <- function(x, time, event, params,
                                      label = "CoxBoost CV") {
  inputs <- sb_validate_classical_survival_inputs(x, time, event, label)
  if (!sb_coxboost_available()) {
    sb_abort(
      "Package `CoxBoost` is required for the CV-tuned CoxBoost learner.",
      "stormbreaker_missing_coxboost"
    )
  }
  fold_id <- sb_make_stratified_fold_id(
    inputs$event,
    params$inner_folds,
    params$seed,
    stratum_label = "event"
  )
  folds <- lapply(seq_len(params$inner_folds), function(fold) {
    which(fold_id == fold)
  })
  x_matrix <- as.matrix(inputs$x)
  penalty_search <- sb_classical_survival_stop_on_warning(
    withr::with_seed(
      params$seed,
      CoxBoost::optimCoxBoostPenalty(
        time = inputs$time,
        status = inputs$event,
        x = x_matrix,
        start.penalty = params$start_penalty,
        parallel = FALSE,
        trace = FALSE
      )
    ),
    paste0(label, " penalty search"),
    "stormbreaker_coxboost_cv_warning"
  )
  penalty <- as.numeric(penalty_search$penalty)
  if (length(penalty) != 1L || !is.finite(penalty) || penalty <= 0) {
    sb_abort(
      "CV-tuned CoxBoost returned an invalid penalty.",
      "stormbreaker_invalid_coxboost_cv"
    )
  }
  cross_validation <- sb_classical_survival_stop_on_warning(
    withr::with_seed(
      params$seed,
      CoxBoost::cv.CoxBoost(
        time = inputs$time,
        status = inputs$event,
        x = x_matrix,
        maxstepno = params$maximum_steps,
        K = params$inner_folds,
        type = params$tuning_type,
        penalty = penalty,
        folds = folds,
        parallel = FALSE,
        multicore = FALSE,
        trace = FALSE
      )
    ),
    paste0(label, " step search"),
    "stormbreaker_coxboost_cv_warning"
  )
  optimal_step <- as.integer(cross_validation$optimal.step)
  if (length(optimal_step) != 1L || is.na(optimal_step) ||
    optimal_step < 1L || optimal_step > params$maximum_steps) {
    sb_abort(
      "CV-tuned CoxBoost returned an invalid optimal step.",
      "stormbreaker_invalid_coxboost_cv"
    )
  }
  engine <- sb_classical_survival_stop_on_warning(
    CoxBoost::CoxBoost(
      time = inputs$time,
      status = inputs$event,
      x = x_matrix,
      stepno = optimal_step,
      penalty = penalty,
      standardize = params$standardize,
      trace = FALSE
    ),
    paste0(label, " final fit"),
    "stormbreaker_coxboost_cv_warning"
  )
  native_coefficient <- as.numeric(stats::coef(
    engine,
    at.step = optimal_step,
    scaled = FALSE
  ))
  scale <- as.numeric(engine$sdx)
  valid <- length(native_coefficient) == ncol(inputs$x) &&
    length(scale) == ncol(inputs$x) && !anyNA(native_coefficient) &&
    all(is.finite(native_coefficient)) && !anyNA(scale) &&
    all(is.finite(scale)) && all(scale > 0)
  if (!valid) {
    sb_abort(
      "CV-tuned CoxBoost returned invalid coefficients.",
      "stormbreaker_invalid_coxboost_cv"
    )
  }
  names(native_coefficient) <- names(inputs$x)
  list(
    inputs = inputs,
    engine = engine,
    native_coefficient = native_coefficient,
    coefficient = native_coefficient / scale,
    penalty = penalty,
    optimal_step = optimal_step,
    fold_id = fold_id,
    curve = data.frame(
      step = seq_along(cross_validation$mean.logplik) - 1L,
      mean_log_likelihood = as.numeric(cross_validation$mean.logplik),
      standard_error = as.numeric(cross_validation$se.logplik),
      stringsAsFactors = FALSE
    )
  )
}

sb_fit_coxboost_cv <- function(x, time, event, params) {
  native <- sb_fit_coxboost_cv_native(x, time, event, params)
  sb_new_classical_survival_model(
    learner = "coxboost_cv",
    features = names(native$inputs$x),
    coefficient = native$coefficient,
    params = params,
    audit = list(
      effective_penalty = native$penalty,
      optimal_step = native$optimal_step,
      fold_id = native$fold_id,
      tuning_curve = native$curve,
      standardize = params$standardize,
      tuning_type = params$tuning_type
    ),
    engine = native$engine,
    baseline_source = "coxboost_cv_breslow_log_hazard"
  )
}

sb_fit_coxboost <- function(x, time, event, params) {
  inputs <- sb_validate_classical_survival_inputs(
    x,
    time,
    event,
    "CoxBoost"
  )
  if (!sb_coxboost_available()) {
    sb_abort(
      "Package `CoxBoost` is required for the fixed CoxBoost learner.",
      "stormbreaker_missing_coxboost"
    )
  }
  effective_penalty <- if (is.null(params$penalty)) {
    9 * sum(inputs$event == 1L)
  } else {
    params$penalty
  }
  engine <- sb_classical_survival_stop_on_warning(
    CoxBoost::CoxBoost(
      time = inputs$time,
      status = inputs$event,
      x = as.matrix(inputs$x),
      standardize = params$standardize,
      stepno = params$stepno,
      penalty = effective_penalty,
      criterion = params$criterion,
      trace = FALSE
    ),
    "fixed CoxBoost learner",
    "stormbreaker_coxboost_fit_warning"
  )
  internal <- as.numeric(stats::coef(
    engine,
    at.step = params$stepno,
    scaled = FALSE
  ))
  scale <- as.numeric(engine$sdx)
  coefficient <- internal / scale
  sb_new_classical_survival_model(
    learner = "coxboost",
    features = names(inputs$x),
    coefficient = coefficient,
    params = params,
    audit = list(
      effective_penalty = as.numeric(effective_penalty),
      analysis_events = as.integer(sum(inputs$event == 1L)),
      fitted_steps = as.integer(engine$stepno),
      standardize = isTRUE(engine$standardize),
      criterion = params$criterion
    ),
    engine = engine,
    baseline_source = "coxboost_breslow_log_hazard"
  )
}

sb_fit_plsrcox_cv <- function(x, time, event, params) {
  inputs <- sb_validate_classical_survival_inputs(x, time, event, "plsRcox CV")
  if (!sb_plsrcox_available()) {
    sb_abort(
      "Package `plsRcox` is required for the CV-tuned plsRcox learner.",
      "stormbreaker_missing_plsrcox"
    )
  }
  single_predictor <- ncol(inputs$x) == 1L
  # Small late-component effects can trigger the native infinite-coefficient
  # warning when likelihood convergence precedes score convergence. Tighten
  # likelihood accuracy; keep the native Cholesky and infinite-coefficient
  # warning thresholds. toler.inf otherwise changes implicitly with eps.
  cox_control <- survival::coxph.control(
    eps = 1e-11,
    iter.max = 50L,
    toler.inf = survival::coxph.control()$toler.inf
  )
  search_warnings <- list()
  warning_indices <- list()
  abort_search_warning <- function(warning) {
    sb_abort(
      paste0("The plsRcox component search raised a warning: ",
             conditionMessage(warning)),
      "stormbreaker_plsrcox_cv_warning",
      parent = warning
    )
  }
  cross_validation <- if (single_predictor) {
    # There is exactly one feasible component, so the native component-count
    # search has no choice to make. cv.plsRcox drops one-column folds to a
    # vector; retain its native final learner at the unique component count.
    list(
      lambda.min5 = 1L,
      folds = list(),
      cv.error5 = NA_real_,
      cv.se5 = NA_real_
    )
  } else {
    withCallingHandlers(
      withr::with_seed(
        params$seed,
        plsRcox::cv.plsRcox(
          list(x = inputs$x, time = inputs$time, status = inputs$event),
          nfold = params$inner_folds,
          nt = params$maximum_components,
          plot.it = FALSE,
          folddetails = TRUE,
          verbose = FALSE,
          control = cox_control
        )
      ),
      warning = function(warning) {
        pattern <- paste0(
          "^Loglik converged before variable[[:space:]]+",
          "([0-9,[:space:]]+); coefficient may be infinite\\.[[:space:]]*$"
        )
        message <- conditionMessage(warning)
        matched <- regmatches(message, regexec(pattern, message))[[1L]]
        if (length(matched) != 2L) abort_search_warning(warning)
        indices <- as.integer(strsplit(
          gsub("[[:space:]]", "", matched[[2L]]), ",", fixed = TRUE
        )[[1L]])
        if (anyNA(indices) || any(indices < 1L |
                                 indices > params$maximum_components)) {
          abort_search_warning(warning)
        }
        search_warnings[[length(search_warnings) + 1L]] <<- warning
        warning_indices[[length(warning_indices) + 1L]] <<- indices
        invokeRestart("muffleWarning")
      }
    )
  }
  selected_components <- as.integer(cross_validation$lambda.min5)
  if (length(selected_components) != 1L || is.na(selected_components) ||
    selected_components < 1L ||
    selected_components > params$maximum_components) {
    sb_abort(
      "CV-tuned plsRcox returned an invalid component count.",
      "stormbreaker_invalid_plsrcox_cv"
    )
  }
  if (length(search_warnings) > 0L) {
    # Native component k fits Cox models with at most k columns. A warning
    # naming only later columns cannot come from the selected prefix.
    affected <- which(vapply(warning_indices, function(indices) {
      any(indices <= selected_components)
    }, logical(1)))
    if (length(affected) > 0L) {
      abort_search_warning(search_warnings[[affected[[1L]]]])
    }
    completed <- cross_validation$completed.cv5
    valid <- is.matrix(completed) &&
      nrow(completed) >= selected_components + 1L &&
      ncol(completed) == params$inner_folds &&
      !anyNA(completed[selected_components + 1L, ]) &&
      all(completed[selected_components + 1L, ])
    if (!valid) {
      sb_abort(
        "CV-tuned plsRcox did not complete every fold at the selected component count.",
        "stormbreaker_invalid_plsrcox_cv"
      )
    }
  }
  search_messages <- vapply(search_warnings, conditionMessage, character(1))
  engine <- sb_classical_survival_stop_on_warning(
    plsRcox::plsRcox(
      inputs$x,
      time = inputs$time,
      event = inputs$event,
      nt = selected_components,
      control = cox_control,
      scaleX = params$scale_x,
      scaleY = params$scale_y,
      sparseStop = params$sparse_stop,
      allres = FALSE,
      verbose = FALSE
    ),
    "CV-tuned plsRcox final fit",
    "stormbreaker_plsrcox_cv_warning"
  )
  coefficient_matrix <- engine$Coeffs
  coefficient <- if (is.matrix(coefficient_matrix) &&
    ncol(coefficient_matrix) >= 1L) {
    as.numeric(coefficient_matrix[, 1L])
  } else {
    numeric()
  }
  if (!identical(rownames(coefficient_matrix), names(inputs$x)) ||
    !identical(as.integer(engine$computed_nt), selected_components)) {
    coefficient <- numeric()
  }
  model <- sb_new_classical_survival_model(
    learner = "plsrcox_cv",
    features = names(inputs$x),
    coefficient = coefficient,
    params = params,
    audit = list(
      selected_components = selected_components,
      maximum_components = params$maximum_components,
      component_selection = if (single_predictor) {
        "single_predictor_has_one_feasible_component"
      } else {
        params$selection_rule
      },
      folds = cross_validation$folds,
      cv_error = as.numeric(cross_validation$cv.error5),
      cv_standard_error = as.numeric(cross_validation$cv.se5),
      component_search_warnings = search_messages,
      cox_control = cox_control
    ),
    engine = engine,
    baseline_source = "plsrcox_cv_breslow_log_hazard"
  )
  if (length(search_messages) > 0L) {
    attr(model, "stormbreaker_fit_warnings") <- search_messages
  }
  model
}

sb_fit_plsrcox <- function(x, time, event, params) {
  inputs <- sb_validate_classical_survival_inputs(
    x,
    time,
    event,
    "plsRcox"
  )
  if (!sb_plsrcox_available()) {
    sb_abort(
      "Package `plsRcox` is required for the fixed plsRcox learner.",
      "stormbreaker_missing_plsrcox"
    )
  }
  if (params$components > ncol(inputs$x)) {
    sb_abort(
      "`plsrcox_components` cannot exceed the fitted predictor count.",
      "stormbreaker_plsrcox_too_many_components"
    )
  }
  engine <- sb_classical_survival_stop_on_warning(
    plsRcox::plsRcox(
      inputs$x,
      time = inputs$time,
      event = inputs$event,
      nt = params$components,
      scaleX = params$scale_x,
      scaleY = params$scale_y,
      sparseStop = params$sparse_stop,
      allres = FALSE,
      verbose = FALSE
    ),
    "fixed plsRcox learner",
    "stormbreaker_plsrcox_fit_warning"
  )
  coefficient_matrix <- engine$Coeffs
  coefficient <- if (is.matrix(coefficient_matrix) &&
    ncol(coefficient_matrix) >= 1L) {
    as.numeric(coefficient_matrix[, 1L])
  } else {
    numeric()
  }
  coefficient_names <- rownames(coefficient_matrix)
  if (!is.null(coefficient_names) &&
    !identical(coefficient_names, names(inputs$x))) {
    coefficient <- numeric()
  }
  if (!identical(as.integer(engine$computed_nt), params$components)) {
    coefficient <- numeric()
  }
  sb_new_classical_survival_model(
    learner = "plsrcox",
    features = names(inputs$x),
    coefficient = coefficient,
    params = params,
    audit = list(
      computed_components = as.integer(engine$computed_nt),
      scale_x = params$scale_x,
      scale_y = params$scale_y,
      sparse_stop = params$sparse_stop
    ),
    engine = engine,
    baseline_source = "plsrcox_breslow_log_hazard"
  )
}

sb_superpc_portable_coefficient <- function(predicted, features,
                                            components) {
  selected <- as.logical(predicted$which.features)
  u <- predicted$u
  if (is.null(dim(u))) {
    u <- matrix(u, ncol = 1L)
  }
  d <- as.numeric(predicted$d)
  outcome_coefficient <- as.numeric(predicted$coef)
  valid <- length(selected) == length(features) &&
    nrow(u) == sum(selected) && ncol(u) == components &&
    length(d) == components && length(outcome_coefficient) == components
  coefficient <- numeric()
  if (valid) {
    scaling <- colSums(sweep(abs(u), 2L, d, "/"))
    selected_coefficient <- as.vector(
      u %*% (outcome_coefficient / (scaling * d))
    )
    coefficient <- numeric(length(features))
    coefficient[selected] <- selected_coefficient
  }
  list(coefficient = coefficient, selected = selected)
}

sb_fit_superpc_cv <- function(x, time, event, params) {
  inputs <- sb_validate_classical_survival_inputs(x, time, event, "SuperPC CV")
  if (!sb_superpc_available()) {
    sb_abort(
      "Package `superpc` is required for the CV-tuned SuperPC learner.",
      "stormbreaker_missing_superpc"
    )
  }
  if (ncol(inputs$x) < params$minimum_features) {
    sb_abort(
      paste0(
        "SuperPC CV threshold tuning requires at least ",
        params$minimum_features, " predictors but received ", ncol(inputs$x),
        "."
      ),
      "stormbreaker_superpc_cv_too_few_features"
    )
  }
  superpc_data <- list(
    x = t(as.matrix(inputs$x)),
    y = inputs$time,
    censoring.status = inputs$event,
    featurenames = names(inputs$x),
    data.names = as.character(seq_len(nrow(inputs$x)))
  )
  trained <- sb_classical_survival_stop_on_warning(
    withr::with_seed(
      params$seed,
      superpc::superpc.train(
        data = superpc_data,
        type = "survival",
        s0.perc = params$s0_perc
      )
    ),
    "SuperPC training call",
    "stormbreaker_superpc_cv_warning"
  )
  cross_validation <- sb_classical_survival_stop_on_warning(
    withr::with_seed(
      params$seed,
      superpc::superpc.cv(
        fit = trained,
        data = superpc_data,
        n.threshold = params$threshold_count,
        n.fold = params$inner_folds,
        n.components = params$tuning_components,
        min.features = params$minimum_features,
        max.features = nrow(superpc_data$x),
        compute.fullcv = TRUE,
        compute.preval = TRUE
      )
    ),
    "SuperPC threshold search",
    "stormbreaker_superpc_cv_warning"
  )
  best_index <- which.max(cross_validation$scor[1L, ])
  threshold <- as.numeric(cross_validation$thresholds[[best_index]])
  if (length(threshold) != 1L || !is.finite(threshold) || threshold < 0) {
    sb_abort(
      "CV-tuned SuperPC returned an invalid threshold.",
      "stormbreaker_invalid_superpc_cv"
    )
  }
  predicted <- sb_classical_survival_stop_on_warning(
    superpc::superpc.predict(
      object = trained,
      data = superpc_data,
      newdata = superpc_data,
      threshold = threshold,
      n.components = params$prediction_components,
      prediction.type = "continuous"
    ),
    "CV-tuned SuperPC final prediction call",
    "stormbreaker_superpc_cv_warning"
  )
  portable <- sb_superpc_portable_coefficient(
    predicted,
    names(inputs$x),
    params$prediction_components
  )
  sb_new_classical_survival_model(
    learner = "superpc_cv",
    features = names(inputs$x),
    coefficient = portable$coefficient,
    params = params,
    audit = list(
      selected_features = names(inputs$x)[portable$selected],
      selected_threshold = threshold,
      selected_threshold_index = as.integer(best_index),
      threshold_count = params$threshold_count,
      tuning_components = params$tuning_components,
      prediction_components = params$prediction_components,
      requested_folds = params$inner_folds,
      effective_folds = as.integer(cross_validation$n.fold),
      cv_scores = cross_validation$scor
    ),
    engine = list(
      trained = trained,
      cross_validation = cross_validation,
      prediction = predicted
    ),
    baseline_source = "superpc_cv_breslow_log_hazard"
  )
}

sb_superpc_single_feature_prediction <- function(trained, data, selected) {
  # The sole principal component is the centered selected predictor, up to
  # sign. Preserve matrix shape and use SuperPC's own outcome fit unchanged.
  x <- data$x[selected, , drop = FALSE]
  centered <- sweep(x, 1L, rowMeans(x), "-")
  decomposition <- svd(centered, nu = 1L, nv = 0L)
  u <- decomposition$u[, 1L, drop = FALSE]
  d <- decomposition$d[1L]
  scaling <- sum(abs(u) / d)
  score <- (t(centered) %*% u) / (scaling * d)
  coefficient <- superpc::superpc.fit.to.outcome(
    trained, data, score, print = FALSE
  )$results$coef
  out <- scale(score, center = FALSE, scale = sign(coefficient))
  list(
    u = u, d = d, coef = coefficient, which.features = selected,
    v.pred = out,
    v.pred.1df = apply(
      scale(out, center = FALSE, scale = 1 / abs(coefficient)), 1L, sum
    ),
    n.components = 1L, prediction.type = "continuous"
  )
}

sb_fit_superpc <- function(x, time, event, params) {
  inputs <- sb_validate_classical_survival_inputs(
    x,
    time,
    event,
    "SuperPC"
  )
  if (!sb_superpc_available()) {
    sb_abort(
      "Package `superpc` is required for the fixed SuperPC learner.",
      "stormbreaker_missing_superpc"
    )
  }
  if (ncol(inputs$x) < 2L) {
    sb_abort(
      paste0(
        "Fixed SuperPC training requires at least 2 predictors because ",
        "the native survival-score path does not preserve a one-row ",
        "feature matrix."
      ),
      "stormbreaker_superpc_too_few_features"
    )
  }
  superpc_data <- list(
    x = t(as.matrix(inputs$x)),
    y = inputs$time,
    censoring.status = inputs$event,
    featurenames = names(inputs$x),
    data.names = as.character(seq_len(nrow(inputs$x)))
  )
  trained <- sb_classical_survival_stop_on_warning(
    superpc::superpc.train(
      data = superpc_data,
      type = "survival",
      s0.perc = params$s0_perc
    ),
    "fixed SuperPC training call",
    "stormbreaker_superpc_train_warning"
  )
  predicted <- tryCatch(
    sb_classical_survival_stop_on_warning(
      {
        selected <- abs(trained$feature.scores) >= params$threshold
        if (sum(selected) == 1L && params$components == 1L) {
          sb_superpc_single_feature_prediction(trained, superpc_data, selected)
        } else {
          superpc::superpc.predict(
            object = trained,
            data = superpc_data,
            newdata = superpc_data,
            threshold = params$threshold,
            n.components = params$components,
            prediction.type = "continuous"
          )
        }
      },
      "fixed SuperPC prediction call",
      "stormbreaker_superpc_predict_warning"
    ),
    error = function(error) {
      if (inherits(error, "stormbreaker_error")) stop(error)
      selected_count <- sum(abs(trained$feature.scores) >= params$threshold)
      sb_abort(
        paste0(
          "Native SuperPC prediction failed after selecting ", selected_count,
          " of ", ncol(inputs$x), " predictors (threshold=", params$threshold,
          ", components=", params$components, "): ", conditionMessage(error)
        ),
        "stormbreaker_superpc_predict_error",
        parent = error
      )
    }
  )
  portable <- sb_superpc_portable_coefficient(
    predicted,
    names(inputs$x),
    params$components
  )
  sb_new_classical_survival_model(
    learner = "superpc",
    features = names(inputs$x),
    coefficient = portable$coefficient,
    params = params,
    audit = list(
      selected_features = names(inputs$x)[portable$selected],
      threshold = params$threshold,
      components = params$components,
      s0_perc = params$s0_perc
    ),
    engine = list(trained = trained, prediction = predicted),
    baseline_source = "superpc_breslow_log_hazard"
  )
}

sb_fit_survivalsvm_calibration <- function(time, event, raw_risk, ties) {
  calibration_data <- data.frame(
    time = time,
    event = event,
    raw_risk = raw_risk
  )
  sb_classical_survival_stop_on_warning(
    survival::coxph(
      survival::Surv(time, event) ~ raw_risk,
      data = calibration_data,
      ties = ties,
      singular.ok = FALSE
    ),
    "survival-SVM Cox calibration",
    "stormbreaker_survivalsvm_calibration_warning"
  )
}

sb_survivalsvm_native_components <- function(engine, type, feature_count) {
  fit <- engine$model.fit
  basis <- matrix(numeric(), nrow = 0L, ncol = 0L)
  row_weight <- numeric()
  intercept <- numeric()
  valid <- FALSE
  if (identical(type, "regression")) {
    basis <- as.matrix(fit$SV)
    row_weight <- as.numeric(fit$Beta)
    intercept <- as.numeric(fit$b0)
    valid <- ncol(basis) == feature_count &&
      nrow(basis) == length(row_weight) && length(intercept) == 1L &&
      is.finite(intercept) && !anyNA(basis) && all(is.finite(basis)) &&
      !anyNA(row_weight) && all(is.finite(row_weight))
  } else {
    basis <- as.matrix(fit$Xtrain)
    differences <- as.matrix(fit$DifMat$Mat)
    alpha <- as.numeric(fit$Alpha)
    valid <- ncol(basis) == feature_count &&
      nrow(basis) == ncol(differences) &&
      nrow(differences) == length(alpha) &&
      !anyNA(basis) && all(is.finite(basis)) &&
      !anyNA(differences) && all(is.finite(differences)) &&
      !anyNA(alpha) && all(is.finite(alpha))
    if (valid) {
      row_weight <- as.vector(crossprod(differences, alpha))
      intercept <- 0
      if (identical(type, "hybrid")) {
        beta <- as.numeric(fit$Beta)
        beta_star <- as.numeric(fit$Betarstar)
        delta <- as.numeric(fit$Delta)
        intercept <- as.numeric(fit$b0)
        valid <- length(beta) == nrow(basis) &&
          length(beta_star) == nrow(basis) &&
          length(delta) == nrow(basis) && length(intercept) == 1L &&
          is.finite(intercept) && !anyNA(beta) && all(is.finite(beta)) &&
          !anyNA(beta_star) && all(is.finite(beta_star)) &&
          !anyNA(delta) && all(is.finite(delta))
        if (valid) {
          row_weight <- row_weight + beta - delta * beta_star
        }
      }
    }
  }
  valid <- valid && nrow(basis) > 0L && length(row_weight) == nrow(basis) &&
    !anyNA(row_weight) && all(is.finite(row_weight))
  if (!valid) {
    sb_abort(
      "Fixed survival-SVM returned invalid native prediction components.",
      "stormbreaker_invalid_survivalsvm_fit"
    )
  }
  list(basis = basis, row_weight = row_weight, intercept = intercept)
}

sb_survivalsvm_linear_components <- function(engine, type, feature_count) {
  native <- sb_survivalsvm_native_components(engine, type, feature_count)
  weight <- as.vector(crossprod(native$basis, native$row_weight))
  if (length(weight) != feature_count || anyNA(weight) ||
    any(!is.finite(weight))) {
    sb_abort(
      "Fixed survival-SVM returned an invalid linear model.",
      "stormbreaker_invalid_survivalsvm_fit"
    )
  }
  list(weight = weight, intercept = native$intercept)
}

sb_survivalsvm_additive_kernel_matrix <- function(basis, x) {
  basis <- as.matrix(basis)
  x <- as.matrix(x)
  ranges <- apply(basis, 2L, max) - apply(basis, 2L, min)
  if (any(ranges == 0)) {
    sb_abort(
      "The fitted additive survival-SVM basis contains a constant feature.",
      "stormbreaker_invalid_survivalsvm_fit"
    )
  }
  transposed_x <- t(x)
  rows <- lapply(seq_len(nrow(basis)), function(index) {
    contribution <- t(
      (ranges - abs(transposed_x - basis[index, ])) / ranges
    )
    rowSums(contribution)
  })
  do.call(rbind, rows)
}

sb_survivalsvm_rbf_kernel_matrix <- function(basis, x, kernel_pars) {
  basis <- as.matrix(basis)
  x <- as.matrix(x)
  valid <- is.numeric(kernel_pars) && length(kernel_pars) > 0L &&
    !anyNA(kernel_pars) && all(is.finite(kernel_pars)) &&
    all(kernel_pars > 0)
  if (!valid) {
    sb_abort(
      "The fitted RBF survival-SVM has invalid kernel parameters.",
      "stormbreaker_invalid_survivalsvm_fit"
    )
  }
  basis_norm <- crossprod(
    t(rowSums(basis^2)),
    t(rep(1, nrow(x)))
  )
  x_norm <- crossprod(
    t(rowSums(x^2)),
    t(rep(1, nrow(basis)))
  )
  squared_distance <- basis_norm + t(x_norm) -
    2 * tcrossprod(basis, x)
  exp(-squared_distance / (2 * kernel_pars[[1L]]))
}

sb_survivalsvm_kernel_matrix <- function(basis, x, kernel, kernel_pars) {
  switch(kernel,
    add_kernel = sb_survivalsvm_additive_kernel_matrix(basis, x),
    rbf_kernel = sb_survivalsvm_rbf_kernel_matrix(
      basis,
      x,
      kernel_pars
    ),
    sb_abort(
      paste0("Unsupported portable survival-SVM kernel `", kernel, "`."),
      "stormbreaker_invalid_survivalsvm_fit"
    )
  )
}

sb_survivalsvm_kernel_reference <- function(kernel) {
  switch(kernel,
    add_kernel = "survivalsvm_additive_kernel_native",
    rbf_kernel = "survivalsvm_rbf_kernel_native",
    sb_abort(
      paste0("Unsupported portable survival-SVM kernel `", kernel, "`."),
      "stormbreaker_invalid_survivalsvm_fit"
    )
  )
}

sb_survivalsvm_kernel_model_class <- function(kernel) {
  switch(kernel,
    add_kernel = "stormbreaker_survivalsvm_additive",
    rbf_kernel = "stormbreaker_survivalsvm_rbf",
    sb_abort(
      paste0("Unsupported portable survival-SVM kernel `", kernel, "`."),
      "stormbreaker_invalid_survivalsvm_fit"
    )
  )
}

sb_survivalsvm_kernel_state_class <- function(kernel) {
  switch(kernel,
    add_kernel = "stormbreaker_survivalsvm_additive_learner_state",
    rbf_kernel = "stormbreaker_survivalsvm_rbf_learner_state",
    sb_abort(
      paste0("Unsupported portable survival-SVM kernel `", kernel, "`."),
      "stormbreaker_invalid_survivalsvm_fit"
    )
  )
}

sb_new_survivalsvm_kernel_model <- function(
  learner, features, basis, row_weight, intercept,
  calibration_coefficient, kernel_pars, params, audit, engine,
  baseline_source
) {
  kernel <- sb_survivalsvm_kernel(learner)
  valid_kernel_pars <- if (identical(kernel, "rbf_kernel")) {
    is.numeric(kernel_pars) && length(kernel_pars) > 0L &&
      !anyNA(kernel_pars) && all(is.finite(kernel_pars)) &&
      all(kernel_pars > 0)
  } else {
    is.null(kernel_pars)
  }
  valid <- learner %in% sb_survivalsvm_kernel_learner_ids() &&
    identical(params$kernel, kernel) && valid_kernel_pars &&
    is.character(features) && length(features) > 0L &&
    !anyNA(features) && !anyDuplicated(features) &&
    is.matrix(basis) && is.numeric(basis) && ncol(basis) == length(features) &&
    nrow(basis) > 0L && !anyNA(basis) && all(is.finite(basis)) &&
    is.numeric(row_weight) && length(row_weight) == nrow(basis) &&
    !anyNA(row_weight) && all(is.finite(row_weight)) &&
    is.numeric(intercept) && length(intercept) == 1L &&
    !is.na(intercept) && is.finite(intercept) &&
    is.numeric(calibration_coefficient) &&
    length(calibration_coefficient) == 1L &&
    !is.na(calibration_coefficient) && is.finite(calibration_coefficient)
  if (!valid) {
    sb_abort(
      paste0("Fixed ", learner, " did not yield a portable kernel model."),
      "stormbreaker_invalid_survivalsvm_fit"
    )
  }
  structure(
    list(
      learner = learner,
      features = features,
      basis = basis,
      row_weight = unname(row_weight),
      intercept = unname(intercept),
      calibration_coefficient = unname(calibration_coefficient),
      kernel_pars = kernel_pars,
      params = params,
      audit = audit,
      baseline_source = baseline_source,
      engine = engine
    ),
    class = c(
      paste0("stormbreaker_", learner),
      sb_survivalsvm_kernel_model_class(kernel),
      "stormbreaker_survivalsvm_kernel_fixed"
    )
  )
}

sb_predict_survivalsvm_kernel <- function(model, x) {
  if (!inherits(model, "stormbreaker_survivalsvm_kernel_fixed")) {
    sb_abort(
      "The supplied model is not a fitted kernel survival-SVM learner.",
      "stormbreaker_bad_survivalsvm_fit"
    )
  }
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), model$features)) {
    sb_abort(
      "Prediction columns do not match the fitted survival-SVM schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  scaled_x <- sweep(
    sweep(as.matrix(x), 2L, model$audit$center, "-"),
    2L,
    model$audit$scale,
    "/"
  )
  kernel <- sb_survivalsvm_kernel_matrix(
    model$basis,
    scaled_x,
    model$params$kernel,
    model$kernel_pars
  )
  native_score <- as.vector(crossprod(model$row_weight, kernel)) +
    model$intercept
  risk <- -model$calibration_coefficient * native_score
  if (length(risk) != nrow(x) || anyNA(risk) || any(!is.finite(risk))) {
    sb_abort(
      "The fitted kernel survival-SVM produced invalid risk predictions.",
      "stormbreaker_invalid_survivalsvm_prediction"
    )
  }
  unname(risk)
}

sb_predict_survivalsvm_kernel_survival <- function(
  model, x_analysis, time, event, x_assessment, prediction_horizon
) {
  inputs <- sb_validate_classical_survival_inputs(
    x_analysis,
    time,
    event,
    model$learner
  )
  prediction_horizon <- sb_validate_prediction_horizon(
    prediction_horizon,
    inputs$time
  )
  analysis_risk <- sb_predict_survivalsvm_kernel(model, inputs$x)
  assessment_risk <- sb_predict_survivalsvm_kernel(
    model,
    x_assessment
  )
  # A common log-risk offset cancels in H0(t) * exp(risk). Store the
  # baseline at the training mean to avoid overflow from native intercepts.
  risk_reference <- mean(analysis_risk)
  cumulative_hazard <- sb_breslow_baseline_hazard(
    inputs$time,
    inputs$event,
    analysis_risk - risk_reference,
    prediction_horizon
  )
  result <- sb_survival_prediction_result(
    prediction_horizon,
    cumulative_hazard,
    source = model$baseline_source,
    ties = "breslow",
    risk = assessment_risk - risk_reference
  )
  result$risk_reference <- risk_reference
  result
}

sb_extract_survivalsvm_kernel_state <- function(model, context = list()) {
  if (!inherits(model, "stormbreaker_survivalsvm_kernel_fixed")) {
    sb_bad_learner_state(
      "The fitted kernel survival-SVM learner state is invalid."
    )
  }
  kernel <- sb_survivalsvm_kernel(model$learner)
  basis <- model$basis
  rownames(basis) <- NULL
  baseline_hazard <- NULL
  if ("prediction_horizon" %in% names(context)) {
    probability <- sb_predict_survivalsvm_kernel_survival(
      model,
      context$x,
      context$time,
      context$event,
      context$x[1L, , drop = FALSE],
      context$prediction_horizon
    )
    baseline_hazard <- list(
      prediction_horizon = probability$prediction_horizon,
      time_unit = if ("time_unit" %in% names(context)) {
        context$time_unit
      } else {
        NA_character_
      },
      cumulative_baseline_hazard = probability$cumulative_baseline_hazard,
      risk_reference = probability$risk_reference,
      baseline_survival = probability$baseline_survival,
      source = probability$source,
      ties = probability$ties
    )
  }
  structure(
    list(
      learner = model$learner,
      features = model$features,
      params = model$params,
      audit = model$audit,
      center = model$audit$center,
      scale = model$audit$scale,
      basis = basis,
      row_weight = model$row_weight,
      intercept = model$intercept,
      calibration_coefficient = model$calibration_coefficient,
      kernel_pars = model$kernel_pars,
      risk_reference = sb_survivalsvm_kernel_reference(kernel),
      risk_formula =
        "negative_native_kernel_score_times_analysis_cox_coefficient",
      training_data_retention =
        "required_scaled_kernel_basis_without_outcomes",
      baseline_hazard = baseline_hazard
    ),
    class = c(
      sb_survivalsvm_kernel_state_class(kernel),
      "stormbreaker_survivalsvm_kernel_learner_state",
      "stormbreaker_learner_state"
    )
  )
}

sb_predict_survivalsvm_kernel_state <- function(state, x) {
  learner <- if (is.list(state)) state[["learner"]] else NULL
  if (!is.character(learner) || length(learner) != 1L || is.na(learner) ||
    !learner %in% sb_survivalsvm_kernel_learner_ids()) {
    sb_bad_learner_state("The saved kernel survival-SVM state is invalid.")
  }
  kernel_type <- sb_survivalsvm_kernel(learner)
  x <- sb_validate_learner_state_x(state, x, learner)
  basis <- state[["basis"]]
  row_weight <- state[["row_weight"]]
  intercept <- state[["intercept"]]
  calibration_coefficient <- state[["calibration_coefficient"]]
  kernel_pars <- state[["kernel_pars"]]
  center <- state[["center"]]
  scale <- state[["scale"]]
  valid <- is.matrix(basis) && is.numeric(basis) &&
    identical(dim(basis), c(length(row_weight), ncol(x))) &&
    nrow(basis) > 0L && !anyNA(basis) && all(is.finite(basis)) &&
    is.numeric(row_weight) && !anyNA(row_weight) &&
    all(is.finite(row_weight)) &&
    is.numeric(intercept) && length(intercept) == 1L &&
    !is.na(intercept) && is.finite(intercept) &&
    is.numeric(calibration_coefficient) &&
    length(calibration_coefficient) == 1L &&
    !is.na(calibration_coefficient) && is.finite(calibration_coefficient) &&
    is.numeric(center) && length(center) == ncol(x) &&
    !anyNA(center) && all(is.finite(center)) &&
    is.numeric(scale) && length(scale) == ncol(x) &&
    !anyNA(scale) && all(is.finite(scale)) && all(scale > 0)
  valid_kernel_pars <- if (identical(kernel_type, "rbf_kernel")) {
    is.numeric(kernel_pars) && length(kernel_pars) > 0L &&
      !anyNA(kernel_pars) && all(is.finite(kernel_pars)) &&
      all(kernel_pars > 0)
  } else {
    is.null(kernel_pars)
  }
  valid <- valid && valid_kernel_pars
  if (!valid) {
    sb_bad_learner_state("The saved kernel survival-SVM state is invalid.")
  }
  scaled_x <- sweep(
    sweep(as.matrix(x), 2L, center, "-"),
    2L,
    scale,
    "/"
  )
  kernel <- sb_survivalsvm_kernel_matrix(
    basis,
    scaled_x,
    kernel_type,
    kernel_pars
  )
  native_score <- as.vector(crossprod(row_weight, kernel)) + intercept
  risk <- unname(-calibration_coefficient * native_score)
  if (anyNA(risk) || any(!is.finite(risk))) {
    sb_abort(
      "The saved kernel survival-SVM state produced invalid risk scores.",
      "stormbreaker_invalid_survivalsvm_prediction"
    )
  }
  result <- list(risk = risk)
  baseline <- state[["baseline_hazard"]]
  if (!is.null(baseline)) {
    valid_baseline <- is.list(baseline) &&
      is.numeric(baseline[["prediction_horizon"]]) &&
      length(baseline[["prediction_horizon"]]) == 1L &&
      !is.na(baseline[["prediction_horizon"]]) &&
      is.finite(baseline[["prediction_horizon"]]) &&
      baseline[["prediction_horizon"]] > 0
    if (!valid_baseline) {
      sb_bad_learner_state(
        "The saved kernel survival-SVM baseline-hazard state is invalid."
      )
    }
    result$prediction_horizon <- baseline[["prediction_horizon"]]
    result$time_unit <- baseline[["time_unit"]]
    # Legacy states used the uncentered baseline (reference zero).
    risk_reference <- baseline[["risk_reference"]]
    if (is.null(risk_reference)) {
      risk_reference <- 0
    }
    if (!is.numeric(risk_reference) || length(risk_reference) != 1L ||
      is.na(risk_reference) || !is.finite(risk_reference)) {
      sb_bad_learner_state(
        "The saved kernel survival-SVM baseline risk reference is invalid."
      )
    }
    result$survival_probability <- sb_survival_probability_from_hazard(
      risk - risk_reference,
      baseline[["cumulative_baseline_hazard"]]
    )
  }
  result
}

sb_fit_survivalsvm <- function(
  x,
  time,
  event,
  params,
  learner = "survivalsvm"
) {
  inputs <- sb_validate_classical_survival_inputs(
    x,
    time,
    event,
    "survival-SVM"
  )
  if (!sb_survivalsvm_available()) {
    sb_abort(
      "Package `survivalsvm` is required for the fixed survival-SVM learner.",
      "stormbreaker_missing_survivalsvm"
    )
  }
  minimum_features <- if (
    identical(params$kernel, "add_kernel") ||
      (identical(params$kernel, "rbf_kernel") &&
        !identical(params$type, "regression"))
  ) {
    2L
  } else {
    1L
  }
  if (ncol(inputs$x) < minimum_features) {
    sb_abort(
      paste0(
        "Native ", params$type, " survival-SVM with ", params$kernel,
        " requires at least ", minimum_features, " predictors but received ",
        ncol(inputs$x), "."
      ),
      "stormbreaker_survivalsvm_too_few_features"
    )
  }
  center <- vapply(inputs$x, mean, numeric(1))
  scale <- vapply(inputs$x, stats::sd, numeric(1))
  scale[!is.finite(scale) | scale == 0] <- 1
  scaled_x <- as.data.frame(
    sweep(sweep(as.matrix(inputs$x), 2L, center, "-"), 2L, scale, "/"),
    optional = TRUE,
    stringsAsFactors = FALSE
  )
  names(scaled_x) <- names(inputs$x)
  time_name <- ".stormbreaker_time"
  while (time_name %in% names(scaled_x)) {
    time_name <- paste0(time_name, "_")
  }
  event_name <- ".stormbreaker_event"
  while (event_name %in% c(names(scaled_x), time_name)) {
    event_name <- paste0(event_name, "_")
  }
  analysis <- scaled_x
  analysis[[time_name]] <- inputs$time
  analysis[[event_name]] <- inputs$event
  response <- paste0(
    "survival::Surv(",
    sb_quote_name(time_name), ", ",
    sb_quote_name(event_name), ")"
  )
  rhs <- paste(vapply(names(inputs$x), sb_quote_name, ""), collapse = " + ")
  formula <- stats::as.formula(
    paste(response, "~", rhs),
    env = parent.frame()
  )
  if ("diff_meth" %in% names(params)) {
    # survivalsvm 0.0.6 resolves this exported function in the caller frame.
    makediff3 <- survivalsvm::makediff3
  }
  kernel_parameters <- params$kernel_pars
  if (identical(params$kernel_parameter_rule, "analysis_median_squared_distance")) {
    kernel_parameters <- stats::median(as.numeric(stats::dist(scaled_x))^2)
  }
  engine <- sb_classical_survival_stop_on_warning(
    survivalsvm::survivalsvm(
      formula = formula,
      data = analysis,
      type = params$type,
      diff.meth = if ("diff_meth" %in% names(params)) {
        params$diff_meth
      } else {
        NULL
      },
      gamma.mu = if (identical(params$type, "hybrid")) {
        rep(params$gamma_mu, 2L)
      } else {
        params$gamma_mu
      },
      opt.meth = params$opt_meth,
      kernel = params$kernel,
      kernel.pars = kernel_parameters
    ),
    "fixed survival-SVM learner",
    "stormbreaker_survivalsvm_fit_warning"
  )
  kernel_pars <- NULL
  if (params$kernel %in% c("add_kernel", "rbf_kernel")) {
    components <- sb_survivalsvm_native_components(
      engine,
      params$type,
      ncol(inputs$x)
    )
    if (identical(params$kernel, "rbf_kernel")) {
      kernel_pars <- as.numeric(engine$model.fit$Kernel$Kernpar)
    }
    kernel <- sb_survivalsvm_kernel_matrix(
      components$basis,
      scaled_x,
      params$kernel,
      kernel_pars
    )
    native_score <- as.vector(crossprod(components$row_weight, kernel)) +
      components$intercept
    raw_risk <- -native_score
  } else {
    components <- sb_survivalsvm_linear_components(
      engine,
      params$type,
      ncol(inputs$x)
    )
    scaled_weight <- components$weight
    raw_risk <- -(
      as.vector(as.matrix(scaled_x) %*% scaled_weight) + components$intercept
    )
  }
  calibration <- sb_fit_survivalsvm_calibration(
    inputs$time,
    inputs$event,
    raw_risk,
    params$calibration_ties
  )
  calibration_coefficient <- as.numeric(stats::coef(calibration))
  audit <- list(
    calibration_coefficient = calibration_coefficient,
    raw_score = switch(params$type,
      regression = "negative_predicted_survival_time",
      vanbelle1 = "negative_native_vanbelle1_ranking_score",
      vanbelle2 = "negative_native_vanbelle2_ranking_score",
      hybrid = "negative_native_hybrid_score"
    ),
    center = center,
    scale = scale,
    type = params$type,
    diff_meth = if ("diff_meth" %in% names(params)) {
      params$diff_meth
    } else {
      NULL
    },
    kernel = params$kernel,
    kernel_pars = kernel_pars,
    opt_meth = params$opt_meth
  )
  baseline_source <- paste0(
    learner,
    "_breslow_calibrated_log_hazard"
  )
  if (params$kernel %in% c("add_kernel", "rbf_kernel")) {
    return(sb_new_survivalsvm_kernel_model(
      learner = learner,
      features = names(inputs$x),
      basis = components$basis,
      row_weight = components$row_weight,
      intercept = components$intercept,
      calibration_coefficient = calibration_coefficient,
      kernel_pars = kernel_pars,
      params = params,
      audit = audit,
      engine = engine,
      baseline_source = baseline_source
    ))
  }
  coefficient <- -calibration_coefficient * scaled_weight / scale
  sb_new_classical_survival_model(
    learner = learner,
    features = names(inputs$x),
    coefficient = coefficient,
    params = params,
    audit = audit,
    engine = engine,
    baseline_source = baseline_source
  )
}

sb_classical_survival_adapter <- function(id, params, fit_function) {
  list(
    id = id,
    params = params,
    fit = function(x, time, event) {
      fit_function(x, time, event, params)
    },
    predict = sb_predict_classical_survival,
    predict_survival = sb_predict_classical_survival_survival,
    extract_coefficients = sb_extract_classical_survival_coefficients,
    extract_state = function(model, context = list()) {
      sb_extract_cox_learner_state(
        model,
        learner = id,
        coefficients = sb_extract_classical_survival_coefficients(model),
        predict_survival = sb_predict_classical_survival_survival,
        context = context,
        params = model$params,
        audit = model$audit
      )
    },
    predict_state = sb_predict_cox_learner_state
  )
}

sb_coxboost_learner_adapter <- function(stepno = 100L,
                                              penalty = NULL) {
  sb_classical_survival_adapter(
    "coxboost",
    sb_coxboost_params(stepno, penalty),
    sb_fit_coxboost
  )
}

sb_coxboost_cv_learner_adapter <- function(inner_folds = 10L,
                                           seed = 1L,
                                           maximum_steps = 500L,
                                           start_penalty = 500) {
  sb_classical_survival_adapter(
    "coxboost_cv",
    sb_coxboost_cv_params(
      inner_folds = inner_folds,
      seed = seed,
      maximum_steps = maximum_steps,
      start_penalty = start_penalty
    ),
    sb_fit_coxboost_cv
  )
}

sb_plsrcox_learner_adapter <- function(components = 2L) {
  sb_classical_survival_adapter(
    "plsrcox",
    sb_plsrcox_params(components),
    sb_fit_plsrcox
  )
}

sb_plsrcox_cv_learner_adapter <- function(inner_folds = 5L,
                                          seed = 1L,
                                          maximum_components = 10L) {
  sb_classical_survival_adapter(
    "plsrcox_cv",
    sb_plsrcox_cv_params(
      inner_folds = inner_folds,
      seed = seed,
      maximum_components = maximum_components
    ),
    sb_fit_plsrcox_cv
  )
}

sb_superpc_learner_adapter <- function(threshold = 1,
                                             components = 1L) {
  sb_classical_survival_adapter(
    "superpc",
    sb_superpc_params(threshold, components),
    sb_fit_superpc
  )
}

sb_superpc_cv_learner_adapter <- function(inner_folds = 10L,
                                          seed = 1L) {
  sb_classical_survival_adapter(
    "superpc_cv",
    sb_superpc_cv_params(inner_folds = inner_folds, seed = seed),
    sb_fit_superpc_cv
  )
}

sb_survivalsvm_learner_adapter <- function(
  gamma_mu = 1,
  learner = "survivalsvm",
  params = sb_survivalsvm_params(gamma_mu, learner)
) {
  if (learner %in% sb_survivalsvm_kernel_learner_ids()) {
    return(list(
      id = learner,
      params = params,
      fit = function(x, time, event) {
        sb_fit_survivalsvm(x, time, event, params, learner)
      },
      predict = sb_predict_survivalsvm_kernel,
      predict_survival = sb_predict_survivalsvm_kernel_survival,
      extract_state = sb_extract_survivalsvm_kernel_state,
      predict_state = sb_predict_survivalsvm_kernel_state
    ))
  }
  sb_classical_survival_adapter(
    learner,
    params,
    function(x, time, event, params) {
      sb_fit_survivalsvm(x, time, event, params, learner)
    }
  )
}
