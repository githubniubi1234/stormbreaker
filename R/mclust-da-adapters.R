sb_mclust_da_learner_ids <- function() {
  c("edda_binary_bic", "mclustda_binary_bic")
}

sb_mclust_da_multivariate_models <- function() {
  c(
    "EII", "VII", "EEI", "VEI", "EVI", "VVI", "EEE",
    "VEE", "EVE", "VVE", "EEV", "VEV", "EVV", "VVV"
  )
}

sb_mclust_da_params <- function(learner) {
  if (!learner %in% sb_mclust_da_learner_ids()) {
    sb_abort(
      paste0("Unknown mclust discriminant learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  }
  list(
    model_type = if (identical(learner, "edda_binary_bic")) {
      "EDDA"
    } else {
      "MclustDA"
    },
    components = if (identical(learner, "edda_binary_bic")) 1L else 1:3,
    univariate_covariance_models = c("E", "V"),
    multivariate_covariance_models = sb_mclust_da_multivariate_models(),
    criterion = "BIC",
    class_prior = "analysis_fold_empirical",
    parameter_prior = "none",
    scale = "original",
    em_control = list(
      eps = .Machine$double.eps,
      tol = c(1e-5, sqrt(.Machine$double.eps)),
      itmax = rep.int(.Machine$integer.max, 2L),
      equalPro = FALSE
    ),
    initialization = if (identical(learner, "edda_binary_bic")) {
      list(method = "known_class_labels")
    } else {
      list(
        method = "native_hierarchical",
        subset_max = 2000L,
        hc_model_name = "VVV",
        hc_use = "SVD",
        seed = 1L
      )
    },
    native_candidate_warnings = FALSE,
    all_model_failure = "typed"
  )
}

sb_mclust_available <- function() {
  requireNamespace("mclust", quietly = TRUE)
}

sb_require_mclust <- function() {
  if (!sb_mclust_available()) {
    sb_abort(
      "Package `mclust` is required for model-based discriminant learners.",
      "stormbreaker_missing_mclust"
    )
  }
  invisible(TRUE)
}

sb_mclust_da_active_models <- function(n_features, params) {
  if (identical(n_features, 1L)) {
    params$univariate_covariance_models
  } else {
    params$multivariate_covariance_models
  }
}

sb_validate_mclust_da_inputs <- function(x, outcome, learner) {
  x <- sb_validate_predictor_frame(x, "x")
  outcome <- sb_validate_binary_outcome_vector(outcome, nrow(x))
  if (!all(c(0L, 1L) %in% outcome)) {
    sb_abort(
      paste0("A `", learner, "` analysis fold must contain both classes."),
      "stormbreaker_mclust_da_single_class"
    )
  }
  list(x = x, outcome = outcome)
}

sb_mclust_da_fit_call <- function(matrix_x, grouping, params, learner) {
  # MclustDA resolves native fitting helpers in its caller environment.
  call_environment <- new.env(parent = asNamespace("mclust"))
  call_environment$data <- matrix_x
  call_environment$grouping <- grouping
  call_environment$components <- params$components
  call_environment$covariance_models <- params$active_covariance_models
  call_environment$model_type <- params$model_type
  call_environment$control <- do.call(mclust::emControl, params$em_control)
  call_environment$native_candidate_warnings <-
    params$native_candidate_warnings
  if (identical(params$model_type, "MclustDA")) {
    initialization <- params$initialization
    fixed_options <- list(
      subset = initialization$subset_max,
      hcModelName = initialization$hc_model_name,
      hcUse = initialization$hc_use
    )
    call_environment$mclust.options <- function(option) {
      if (option %in% names(fixed_options)) {
        return(fixed_options[[option]])
      }
      mclust::mclust.options(option)
    }
    local_mclust_bic <- mclust::mclustBIC
    environment(local_mclust_bic) <- call_environment
    call_environment$mclustBIC <- local_mclust_bic
  }

  tryCatch(
    withCallingHandlers(
      {
        fit_call <- quote(MclustDA(
          data = data,
          class = grouping,
          G = components,
          modelNames = covariance_models,
          modelType = model_type,
          prior = NULL,
          control = control,
          warn = native_candidate_warnings,
          verbose = FALSE
        ))
        if (identical(params$model_type, "MclustDA")) {
          withr::with_seed(
            params$initialization$seed,
            eval(fit_call, envir = call_environment)
          )
        } else {
          eval(fit_call, envir = call_environment)
        }
      },
      warning = function(warning) {
        sb_abort(
          paste0(
            "Learner `", learner, "` raised a fitting warning: ",
            conditionMessage(warning)
          ),
          "stormbreaker_mclust_da_fit_warning",
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
          "Learner `", learner, "` could not fit its mclust model: ",
          conditionMessage(error)
        ),
        "stormbreaker_mclust_da_fit_error",
        parent = error
      )
    }
  )
}

sb_validate_mclust_da_search_grid <- function(engine, params) {
  covariance_models <- if (!is.null(params$active_covariance_models)) {
    params$active_covariance_models
  } else {
    sb_mclust_da_active_models(as.integer(engine$d), params)
  }
  if (identical(params$model_type, "EDDA")) {
    small_class_models <- c("EII", "VII", "EEI", "VEI", "EVI", "VVI")
    native_models <- unique(unlist(lapply(engine$models, function(model) {
      if (model$n <= engine$d) {
        intersect(covariance_models, small_class_models)
      } else {
        covariance_models
      }
    }), use.names = FALSE))
    valid <- length(native_models) == length(covariance_models) &&
      setequal(native_models, covariance_models)
  } else {
    expected_components <- as.character(params$components)
    valid <- inherits(engine, "MclustDA") && is.list(engine$models) &&
      length(engine$models) == 2L
    if (valid) {
      valid <- all(vapply(engine$models, function(model) {
        bic <- model$BIC
        is.matrix(bic) &&
          identical(rownames(bic), expected_components) &&
          identical(colnames(bic), covariance_models) &&
          any(is.finite(bic))
      }, logical(1)))
    }
  }
  if (!valid) {
    sb_abort(
      paste0(
        "mclust reduced or did not return the predeclared BIC search grid; ",
        "the learner will not silently fit a narrower specification."
      ),
      "stormbreaker_mclust_da_search_reduction"
    )
  }
  invisible(TRUE)
}

sb_validate_mclust_da_engine <- function(engine, params) {
  expected_classes <- c("0", "1")
  valid <- inherits(engine, "MclustDA") &&
    identical(engine$type, params$model_type) &&
    is.numeric(engine$prop) && length(engine$prop) == 2L &&
    identical(names(engine$prop), expected_classes) &&
    !anyNA(engine$prop) && all(is.finite(engine$prop)) &&
    all(engine$prop > 0) && abs(sum(engine$prop) - 1) < 1e-12 &&
    is.list(engine$models) &&
    identical(names(engine$models), expected_classes) &&
    is.numeric(engine$bic) && length(engine$bic) == 1L &&
    !is.na(engine$bic) && is.finite(engine$bic)
  if (valid) {
    allowed_names <- c(
      params$active_covariance_models,
      "X", "XII", "XXI", "XXX"
    )
    valid <- all(vapply(engine$models, function(model) {
      is.list(model) &&
        is.character(model$modelName) && length(model$modelName) == 1L &&
        !is.na(model$modelName) && model$modelName %in% allowed_names &&
        is.numeric(model$G) && length(model$G) == 1L &&
        !is.na(model$G) && model$G %in% params$components &&
        is.list(model$parameters) &&
        is.numeric(model$parameters$pro) &&
        length(model$parameters$pro) == model$G &&
        !anyNA(model$parameters$pro) &&
        all(is.finite(model$parameters$pro)) &&
        all(model$parameters$pro > 0)
    }, logical(1)))
  }
  if (valid && identical(params$model_type, "EDDA")) {
    valid <- length(unique(vapply(
      engine$models,
      `[[`,
      character(1),
      "modelName"
    ))) == 1L
  }
  if (!valid) {
    sb_abort(
      "The fitted mclust discriminant model is incomplete or invalid.",
      "stormbreaker_incomplete_mclust_da_fit"
    )
  }
  sb_validate_mclust_da_search_grid(engine, params)
  engine
}

sb_mclust_da_audit <- function(engine, params) {
  selected <- data.frame(
    class = names(engine$models),
    model_name = vapply(engine$models, `[[`, character(1), "modelName"),
    components = as.integer(vapply(engine$models, `[[`, numeric(1), "G")),
    stringsAsFactors = FALSE
  )
  candidate_summary <- if (identical(params$model_type, "MclustDA")) {
    data.frame(
      class = names(engine$models),
      requested = vapply(
        engine$models,
        function(model) length(model$BIC),
        integer(1)
      ),
      estimable = vapply(
        engine$models,
        function(model) sum(is.finite(model$BIC)),
        integer(1)
      ),
      initialization_subset_size = vapply(
        engine$models,
        function(model) {
          initialization <- attr(model$BIC, "initialization", exact = TRUE)
          length(initialization$subset)
        },
        integer(1)
      ),
      stringsAsFactors = FALSE
    )
  } else {
    classes <- names(engine$models)
    data.frame(
      class = classes,
      requested = rep.int(
        length(params$active_covariance_models),
        length(classes)
      ),
      estimable = rep.int(NA_integer_, length(classes)),
      initialization_subset_size = rep.int(NA_integer_, length(classes)),
      stringsAsFactors = FALSE
    )
  }
  list(
    class_prior = engine$prop,
    selected = selected,
    candidate_summary = candidate_summary
  )
}

sb_mclust_da_prediction_engine <- function(engine) {
  models <- lapply(engine$models, function(model) {
    model[c("modelName", "parameters")]
  })
  names(models) <- names(engine$models)
  structure(
    list(type = engine$type, prop = engine$prop, models = models),
    class = "MclustDA"
  )
}

sb_validate_mclust_da_prediction_engine <- function(engine, params) {
  expected_classes <- c("0", "1")
  valid <- inherits(engine, "MclustDA") &&
    identical(engine$type, params$model_type) &&
    is.numeric(engine$prop) && length(engine$prop) == 2L &&
    identical(names(engine$prop), expected_classes) &&
    !anyNA(engine$prop) && all(is.finite(engine$prop)) &&
    all(engine$prop > 0) && abs(sum(engine$prop) - 1) < 1e-12 &&
    is.list(engine$models) &&
    identical(names(engine$models), expected_classes) &&
    all(vapply(engine$models, function(model) {
      is.list(model) &&
        is.character(model$modelName) && length(model$modelName) == 1L &&
        !is.na(model$modelName) && is.list(model$parameters)
    }, logical(1)))
  if (!valid) {
    sb_abort(
      "The saved mclust discriminant engine is incomplete or invalid.",
      "stormbreaker_bad_mclust_da_model"
    )
  }
  engine
}

sb_new_mclust_da_model <- function(learner, features, params, engine, audit) {
  structure(
    list(
      learner = learner,
      features = features,
      params = params,
      engine = engine,
      audit = audit
    ),
    class = "stormbreaker_mclust_da"
  )
}

sb_fit_mclust_da <- function(x, outcome, learner, params) {
  sb_require_mclust()
  inputs <- sb_validate_mclust_da_inputs(x, outcome, learner)
  active_params <- params
  active_params$active_covariance_models <- sb_mclust_da_active_models(
    ncol(inputs$x),
    params
  )
  engine <- sb_mclust_da_fit_call(
    data.matrix(inputs$x),
    factor(inputs$outcome, levels = c(0L, 1L)),
    active_params,
    learner
  )
  engine <- sb_validate_mclust_da_engine(engine, active_params)
  audit <- sb_mclust_da_audit(engine, active_params)
  engine <- sb_mclust_da_prediction_engine(engine)
  sb_validate_mclust_da_prediction_engine(engine, active_params)
  sb_new_mclust_da_model(
    learner,
    names(inputs$x),
    active_params,
    engine,
    audit
  )
}

sb_validate_mclust_da_model_x <- function(model, x) {
  valid <- inherits(model, "stormbreaker_mclust_da") &&
    is.character(model$learner) && length(model$learner) == 1L &&
    model$learner %in% sb_mclust_da_learner_ids() &&
    is.character(model$features) && length(model$features) > 0L &&
    !anyNA(model$features) && !anyDuplicated(model$features) &&
    is.list(model$params) && is.list(model$audit)
  if (!valid) {
    sb_abort(
      "The fitted mclust discriminant model is incomplete or invalid.",
      "stormbreaker_bad_mclust_da_model"
    )
  }
  sb_validate_mclust_da_prediction_engine(model$engine, model$params)
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), model$features)) {
    sb_abort(
      "Prediction columns do not match the fitted mclust model.",
      "stormbreaker_mclust_da_schema_mismatch"
    )
  }
  x
}

sb_predict_mclust_da <- function(model, x) {
  sb_require_mclust()
  x <- sb_validate_mclust_da_model_x(model, x)
  prediction <- tryCatch(
    withCallingHandlers(
      stats::predict(model$engine, newdata = data.matrix(x)),
      warning = function(warning) {
        sb_abort(
          paste0(
            "Learner `", model$learner, "` raised a prediction warning: ",
            conditionMessage(warning)
          ),
          "stormbreaker_mclust_da_prediction_warning",
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
          "Learner `", model$learner, "` could not produce probabilities: ",
          conditionMessage(error)
        ),
        "stormbreaker_mclust_da_prediction_error",
        parent = error
      )
    }
  )
  probability <- if (is.list(prediction) && is.matrix(prediction$z) &&
      "1" %in% colnames(prediction$z)) {
    unname(prediction$z[, "1"])
  } else {
    numeric()
  }
  valid <- length(probability) == nrow(x) && !anyNA(probability) &&
    all(is.finite(probability)) && all(probability >= 0 & probability <= 1)
  if (!valid) {
    sb_abort(
      paste0("Learner `", model$learner, "` produced invalid probabilities."),
      "stormbreaker_invalid_binary_prediction"
    )
  }
  probability
}

sb_extract_mclust_da_learner_state <- function(model, context = list()) {
  if (!inherits(model, "stormbreaker_mclust_da")) {
    sb_bad_learner_state("The fitted mclust discriminant state is invalid.")
  }
  structure(
    list(
      learner = model$learner,
      features = model$features,
      params = model$params,
      engine = model$engine,
      audit = model$audit
    ),
    class = c(
      "stormbreaker_mclust_da_learner_state",
      "stormbreaker_learner_state"
    )
  )
}

sb_predict_mclust_da_learner_state <- function(state, x) {
  learner <- if (is.list(state)) state[["learner"]] else NULL
  if (!is.character(learner) || length(learner) != 1L || is.na(learner) ||
      !learner %in% sb_mclust_da_learner_ids()) {
    sb_bad_learner_state("The saved mclust discriminant state is invalid.")
  }
  x <- sb_validate_learner_state_x(state, x, learner)
  model <- sb_new_mclust_da_model(
    learner,
    state$features,
    state$params,
    state$engine,
    state$audit
  )
  list(probability = sb_predict_mclust_da(model, x))
}

sb_mclust_da_learner_adapter <- function(
    learner,
    params = sb_mclust_da_params(learner)) {
  if (!learner %in% sb_mclust_da_learner_ids()) {
    sb_abort(
      paste0("Unknown mclust discriminant learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  }
  list(
    id = learner,
    params = params,
    fit = function(x, outcome) {
      sb_fit_mclust_da(x, outcome, learner, params)
    },
    predict = sb_predict_mclust_da,
    extract_state = sb_extract_mclust_da_learner_state,
    predict_state = sb_predict_mclust_da_learner_state
  )
}
