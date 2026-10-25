sb_rlt_binary_learner_ids <- function() {
  "rlt_reinforcement_binary"
}

sb_rlt_survival_learner_ids <- function() {
  "rlt_reinforcement_survival"
}

sb_rlt_learner_ids <- function() {
  c(sb_rlt_binary_learner_ids(), sb_rlt_survival_learner_ids())
}

sb_rlt_available <- function() {
  requireNamespace("RLT", quietly = TRUE)
}

sb_rlt_dependency <- function() {
  list(package = "RLT", available = sb_rlt_available())
}

sb_require_rlt_dependency <- function() {
  if (!sb_rlt_available()) {
    sb_abort(
      "Package `RLT` is required for reinforcement-learning forests.",
      "stormbreaker_missing_rlt"
    )
  }
  invisible(TRUE)
}

sb_rlt_learner_params <- function(learner) {
  if (!learner %in% sb_rlt_learner_ids()) {
    sb_abort(
      paste0("Unknown RLT learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  }
  list(
    model = if (learner %in% sb_rlt_binary_learner_ids()) {
      "classification"
    } else {
      "survival"
    },
    num_trees = 100L,
    reinforcement = TRUE,
    combination_split_size = 1L,
    native_nmin = 5L,
    embedded_stop_size = 1L,
    cores = 1L,
    importance = FALSE,
    track_observations = FALSE,
    seed = 1L,
    split_rule = "axis_aligned",
    remaining_controls = "RLT_native_defaults"
  )
}

sb_rlt_parameter_semantics <- function(learner) {
  params <- sb_rlt_learner_params(learner)
  paste0(
    "reinforcement-learning forest with ", params$num_trees,
    " trees, native internal-node split threshold nmin ", params$native_nmin,
    ", axis-aligned splits, embedded-model stop size ",
    params$embedded_stop_size,
    ", one-core execution, and seed ", params$seed,
    " fixed; remaining split, resampling, and muting controls use native ",
    "RLT defaults"
  )
}

sb_rlt_stop_on_warning <- function(code, learner, stage) {
  withCallingHandlers(
    force(code),
    warning = function(warning) {
      message <- conditionMessage(warning)
      resets_weights <- grepl("^(Subject|Variable) weights .*reset to equal weights$", message)
      mismatches_features <- identical(message,
        "test data variables names does not match training data...")
      if (!resets_weights && !mismatches_features) return(invisible(NULL))
      sb_abort(
        paste0(
          "Learner `", learner, "` raised a ", stage, " warning: ",
          message
        ),
        "stormbreaker_rlt_warning",
        parent = warning
      )
    }
  )
}

sb_validate_rlt_survival_native_domain <- function(inputs, learner) {
  if (any(inputs$time <= 0)) {
    sb_abort(
      "RLT survival fitting requires strictly positive follow-up times.",
      "stormbreaker_bad_rlt_time"
    )
  }
  distinct_event_times <- unique(inputs$time[inputs$event == 1L])
  if (length(distinct_event_times) < 2L) {
    # RLT 3.2.6 terminates the Windows R process in compiled code when its
    # survival time grid has one point. An R condition cannot catch that native
    # process exit, so this check belongs immediately before the unsafe call.
    sb_abort(
      paste0(
        "Learner `", learner, "` requires at least two distinct event times ",
        "because the native RLT survival engine cannot safely fit a one-point ",
        "event-time grid."
      ),
      "stormbreaker_rlt_insufficient_event_times"
    )
  }
  inputs
}

sb_validate_rlt_engine <- function(engine, learner, features, params) {
  expected_model <- params$model
  valid <- inherits(engine, "RLT") && is.list(engine) &&
    identical(engine$model, expected_model) &&
    is.numeric(engine$p) && length(engine$p) == 1L &&
    identical(as.integer(engine$p), as.integer(length(features))) &&
    identical(engine$variablenames, features) &&
    is.numeric(engine$ntrees) && length(engine$ntrees) == 1L &&
    identical(as.integer(engine$ntrees), params$num_trees) &&
    is.list(engine$FittedTrees) &&
    length(engine$FittedTrees) == params$num_trees &&
    is.numeric(engine$ncat) && length(engine$ncat) == length(features) &&
    is.numeric(engine$parameters.int) &&
    is.numeric(engine$parameters.double)
  if (valid && learner %in% sb_rlt_binary_learner_ids()) {
    valid <- identical(engine$ylevels, c("control", "case")) &&
      identical(as.integer(engine$nclass), 2L)
  }
  if (valid && learner %in% sb_rlt_survival_learner_ids()) {
    valid <- is.list(engine$FittedSurv) &&
      length(engine$FittedSurv) == params$num_trees &&
      is.numeric(engine$timepoints) && length(engine$timepoints) >= 2L &&
      !anyNA(engine$timepoints) && all(is.finite(engine$timepoints)) &&
      all(engine$timepoints > 0) && all(diff(engine$timepoints) > 0)
  }
  if (!isTRUE(valid)) {
    sb_abort(
      paste0("Learner `", learner, "` returned an invalid native RLT fit."),
      "stormbreaker_invalid_rlt_fit"
    )
  }
  engine
}

sb_rlt_prediction_engine <- function(engine, learner, features, params) {
  engine <- sb_validate_rlt_engine(engine, learner, features, params)
  keep <- c(
    "FittedTrees", "FittedSurv", "timepoints", "ylevels", "nclass",
    "model", "p", "variablenames", "ncat", "xlevels", "ntrees",
    "parameters.int", "parameters.double"
  )
  compact <- engine[intersect(keep, names(engine))]
  class(compact) <- "RLT"
  sb_validate_rlt_engine(compact, learner, features, params)
}

sb_fit_rlt_learner <- function(x, endpoint, learner, params) {
  sb_require_rlt_dependency()
  inputs <- if (learner %in% sb_rlt_binary_learner_ids()) {
    sb_validate_classical_binary_inputs(x, endpoint$outcome, learner)
  } else {
    sb_validate_classical_survival_inputs(
      x,
      endpoint$time,
      endpoint$event,
      learner
    ) |>
      sb_validate_rlt_survival_native_domain(learner)
  }
  matrix_x <- as.matrix(inputs$x)
  storage.mode(matrix_x) <- "double"
  arguments <- list(
    x = matrix_x,
    model = params$model,
    print.summary = 0L,
    use.cores = params$cores,
    ntrees = params$num_trees,
    nmin = params$native_nmin,
    reinforcement = params$reinforcement,
    combsplit = params$combination_split_size,
    embed.n.th = params$embedded_stop_size,
    importance = params$importance,
    track.obs = params$track_observations
  )
  if (learner %in% sb_rlt_binary_learner_ids()) {
    arguments$y <- factor(
      inputs$outcome,
      levels = c(0L, 1L),
      labels = c("control", "case")
    )
  } else {
    arguments$y <- inputs$time
    arguments$censor <- inputs$event
  }
  engine <- sb_rlt_stop_on_warning(
    withr::with_seed(params$seed, do.call(RLT::RLT, arguments)),
    learner,
    "fitting"
  )
  engine <- sb_rlt_prediction_engine(
    engine,
    learner,
    names(inputs$x),
    params
  )
  structure(
    list(
      learner = learner,
      features = names(inputs$x),
      params = params,
      engine = engine
    ),
    class = "stormbreaker_rlt_model"
  )
}

sb_validate_rlt_model_x <- function(model, x) {
  valid <- inherits(model, "stormbreaker_rlt_model") && is.list(model) &&
    model$learner %in% sb_rlt_learner_ids() &&
    is.character(model$features) && length(model$features) > 0L &&
    !anyNA(model$features) && !anyDuplicated(model$features) &&
    is.list(model$params)
  if (!isTRUE(valid)) {
    sb_abort(
      "The supplied model is not a valid fitted RLT learner.",
      "stormbreaker_bad_rlt_fit"
    )
  }
  sb_validate_rlt_engine(
    model$engine,
    model$learner,
    model$features,
    model$params
  )
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), model$features)) {
    sb_abort(
      "Prediction columns do not match the fitted RLT learner schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  x
}

sb_predict_rlt_native <- function(model, x) {
  x <- sb_validate_rlt_model_x(model, x)
  matrix_x <- as.matrix(x)
  storage.mode(matrix_x) <- "double"
  predict_method <- get(
    "predict.RLT",
    envir = asNamespace("RLT"),
    inherits = FALSE
  )
  prediction <- sb_rlt_stop_on_warning(
    predict_method(model$engine, testx = matrix_x),
    model$learner,
    "prediction"
  )
  list(x = x, prediction = prediction)
}

sb_predict_rlt_binary <- function(model, x) {
  native <- sb_predict_rlt_native(model, x)
  if (!model$learner %in% sb_rlt_binary_learner_ids()) {
    sb_abort(
      "A survival RLT model cannot produce binary probabilities.",
      "stormbreaker_bad_rlt_fit"
    )
  }
  probabilities <- as.matrix(native$prediction$ProbPrediction)
  valid <- nrow(probabilities) == nrow(native$x) &&
    ncol(probabilities) == 2L && !anyNA(probabilities) &&
    all(is.finite(probabilities)) &&
    all(probabilities >= 0 & probabilities <= 1) &&
    all(abs(rowSums(probabilities) - 1) <= 1e-10)
  if (!isTRUE(valid)) {
    sb_abort(
      "The RLT classification forest produced invalid class probabilities.",
      "stormbreaker_invalid_rlt_prediction"
    )
  }
  unname(probabilities[, 2L])
}

sb_predict_rlt_survival_curves <- function(model, x) {
  native <- sb_predict_rlt_native(model, x)
  if (!model$learner %in% sb_rlt_survival_learner_ids()) {
    sb_abort(
      "A binary RLT model cannot produce survival curves.",
      "stormbreaker_bad_rlt_fit"
    )
  }
  hazard <- as.matrix(native$prediction$SurvPred)
  timepoints <- as.numeric(native$prediction$timepoints)
  valid <- length(timepoints) >= 2L && !anyNA(timepoints) &&
    all(is.finite(timepoints)) && all(timepoints > 0) &&
    all(diff(timepoints) > 0) && nrow(hazard) == nrow(native$x) &&
    ncol(hazard) == length(timepoints) && !anyNA(hazard) &&
    all(is.finite(hazard)) && all(hazard >= 0 & hazard <= 1)
  if (!isTRUE(valid)) {
    sb_abort(
      "The RLT survival forest produced invalid discrete hazard increments.",
      "stormbreaker_invalid_rlt_prediction"
    )
  }
  survival <- matrix(
    NA_real_,
    nrow = nrow(hazard),
    ncol = ncol(hazard)
  )
  for (row in seq_len(nrow(hazard))) {
    survival[row, ] <- cumprod(1 - hazard[row, ])
  }
  valid_survival <- !anyNA(survival) && all(is.finite(survival)) &&
    all(survival >= 0 & survival <= 1) &&
    all(apply(survival, 1L, function(values) all(diff(values) <= 0)))
  if (!isTRUE(valid_survival)) {
    sb_abort(
      "The reconstructed RLT survival curves are invalid.",
      "stormbreaker_invalid_rlt_prediction"
    )
  }
  interval <- diff(c(0, timepoints))
  previous <- cbind(
    1,
    survival[, -ncol(survival), drop = FALSE]
  )
  risk <- -as.vector(previous %*% interval)
  if (length(risk) != nrow(native$x) || anyNA(risk) ||
    any(!is.finite(risk))) {
    sb_abort(
      "The RLT survival forest produced invalid restricted-mean risk scores.",
      "stormbreaker_invalid_rlt_prediction"
    )
  }
  list(
    risk = unname(risk),
    time = timepoints,
    hazard = hazard,
    survival = survival
  )
}

sb_predict_rlt_learner <- function(model, x) {
  if (model$learner %in% sb_rlt_binary_learner_ids()) {
    return(sb_predict_rlt_binary(model, x))
  }
  sb_predict_rlt_survival_curves(model, x)$risk
}

sb_rlt_direct_survival_source <- function() {
  "rlt_reconstructed_survival_from_native_discrete_hazard"
}

sb_predict_rlt_survival <- function(
  model, x_analysis, time, event, x_assessment, prediction_horizon
) {
  if (!model$learner %in% sb_rlt_survival_learner_ids()) {
    sb_abort(
      "A binary RLT learner cannot predict survival probability.",
      "stormbreaker_bad_rlt_fit"
    )
  }
  inputs <- sb_validate_classical_survival_inputs(
    x_analysis,
    time,
    event,
    model$learner
  )
  if (!identical(names(inputs$x), model$features)) {
    sb_abort(
      "Analysis columns do not match the fitted RLT learner schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  prediction_horizon <- sb_validate_prediction_horizon(
    prediction_horizon,
    inputs$time
  )
  curves <- sb_predict_rlt_survival_curves(model, x_assessment)
  eligible <- which(curves$time <= prediction_horizon)
  probability <- if (length(eligible) == 0L) {
    rep(1, length(curves$risk))
  } else {
    unname(curves$survival[, max(eligible)])
  }
  if (length(probability) != nrow(x_assessment) || anyNA(probability) ||
    any(!is.finite(probability)) || any(probability < 0 | probability > 1)) {
    sb_abort(
      "The RLT survival forest produced invalid horizon probabilities.",
      "stormbreaker_invalid_rlt_prediction"
    )
  }
  list(
    prediction_horizon = prediction_horizon,
    source = sb_rlt_direct_survival_source(),
    ties = "not_applicable",
    probability = probability
  )
}

sb_extract_rlt_state <- function(model, context = list()) {
  sb_validate_rlt_model_x(
    model,
    as.data.frame(
      matrix(
        0,
        nrow = 1L,
        ncol = length(model$features),
        dimnames = list(NULL, model$features)
      ),
      check.names = FALSE
    )
  )
  survival_probability <- NULL
  if (model$learner %in% sb_rlt_survival_learner_ids() &&
    "prediction_horizon" %in% names(context)) {
    probability <- sb_predict_rlt_survival(
      model,
      context$x,
      context$time,
      context$event,
      context$x[1L, , drop = FALSE],
      context$prediction_horizon
    )
    survival_probability <- list(
      prediction_horizon = probability$prediction_horizon,
      time_unit = if ("time_unit" %in% names(context)) {
        context$time_unit
      } else {
        NA_character_
      },
      source = probability$source,
      ties = probability$ties
    )
  }
  structure(
    list(
      learner = model$learner,
      features = model$features,
      params = model$params,
      engine_serialization = "R_xdr_v3_RLT_compact_prediction_engine",
      engine_raw = serialize(
        model$engine,
        connection = NULL,
        ascii = FALSE,
        xdr = TRUE,
        version = 3L
      ),
      risk_reference = if (
        model$learner %in% sb_rlt_survival_learner_ids()
      ) {
        "RLT_native_discrete_hazard"
      } else {
        NULL
      },
      risk_formula = if (
        model$learner %in% sb_rlt_survival_learner_ids()
      ) {
        "negative_step_RMST_over_native_event_grid"
      } else {
        NULL
      },
      baseline_hazard = NULL,
      survival_probability = survival_probability
    ),
    class = c("stormbreaker_rlt_learner_state", "stormbreaker_learner_state")
  )
}

sb_restore_rlt_model <- function(state, learner) {
  if (!is.raw(state$engine_raw) || length(state$engine_raw) == 0L) {
    sb_bad_learner_state("The saved RLT prediction engine is invalid.")
  }
  engine <- tryCatch(
    unserialize(state$engine_raw),
    error = function(error) {
      sb_abort(
        paste0(
          "The saved RLT prediction engine could not be loaded: ",
          conditionMessage(error)
        ),
        "stormbreaker_invalid_rlt_state",
        parent = error
      )
    }
  )
  engine <- sb_validate_rlt_engine(
    engine,
    learner,
    state$features,
    state$params
  )
  structure(
    list(
      learner = learner,
      features = state$features,
      params = state$params,
      engine = engine
    ),
    class = "stormbreaker_rlt_model"
  )
}

sb_validate_rlt_probability_state <- function(state) {
  probability <- state$survival_probability
  valid <- is.list(probability) &&
    is.numeric(probability$prediction_horizon) &&
    length(probability$prediction_horizon) == 1L &&
    !is.na(probability$prediction_horizon) &&
    is.finite(probability$prediction_horizon) &&
    probability$prediction_horizon > 0 &&
    is.character(probability$time_unit) &&
    length(probability$time_unit) == 1L &&
    !is.na(probability$time_unit) && nzchar(probability$time_unit) &&
    identical(probability$source, sb_rlt_direct_survival_source()) &&
    identical(probability$ties, "not_applicable")
  if (!isTRUE(valid)) {
    sb_bad_learner_state(
      "The saved RLT direct survival-probability state is invalid."
    )
  }
  probability
}

sb_predict_rlt_state <- function(state, x) {
  learner <- if (is.list(state)) state$learner else NULL
  if (!is.character(learner) || length(learner) != 1L || is.na(learner) ||
    !learner %in% sb_rlt_learner_ids()) {
    sb_bad_learner_state("The saved RLT learner state is invalid.")
  }
  x <- sb_validate_learner_state_x(state, x, learner)
  sb_require_rlt_dependency()
  model <- sb_restore_rlt_model(state, learner)
  if (learner %in% sb_rlt_binary_learner_ids()) {
    return(list(probability = sb_predict_rlt_binary(model, x)))
  }
  curves <- sb_predict_rlt_survival_curves(model, x)
  result <- list(risk = curves$risk)
  if (!is.null(state$survival_probability)) {
    probability <- sb_validate_rlt_probability_state(state)
    eligible <- which(curves$time <= probability$prediction_horizon)
    result$prediction_horizon <- probability$prediction_horizon
    result$time_unit <- probability$time_unit
    result$survival_probability <- if (length(eligible) == 0L) {
      rep(1, length(curves$risk))
    } else {
      unname(curves$survival[, max(eligible)])
    }
  }
  result
}

sb_rlt_learner_adapter <- function(
  learner,
  params = sb_rlt_learner_params(learner)
) {
  if (!learner %in% sb_rlt_learner_ids()) {
    sb_abort(
      paste0("Unknown RLT learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  }
  adapter <- list(
    id = learner,
    params = params,
    fit = if (learner %in% sb_rlt_binary_learner_ids()) {
      function(x, outcome) {
        sb_fit_rlt_learner(x, list(outcome = outcome), learner, params)
      }
    } else {
      function(x, time, event) {
        sb_fit_rlt_learner(
          x,
          list(time = time, event = event),
          learner,
          params
        )
      }
    },
    predict = sb_predict_rlt_learner,
    extract_state = sb_extract_rlt_state,
    predict_state = sb_predict_rlt_state
  )
  if (learner %in% sb_rlt_survival_learner_ids()) {
    adapter$predict_survival <- sb_predict_rlt_survival
  }
  adapter
}
