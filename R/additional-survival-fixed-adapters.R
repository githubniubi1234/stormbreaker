sb_survreg_learner_ids <- function() {
  c(
    "survreg_weibull",
    "survreg_exponential",
    "survreg_lognormal",
    "survreg_loglogistic"
  )
}

sb_additional_survival_learner_ids <- function() {
  c(
    sb_survreg_learner_ids(),
    "rpart_survival",
    "rsf_survival",
    "ctree_survival",
    "cforest_survival",
    "bart_survival",
    "aorsf_survival"
  )
}

sb_additional_survival_params <- function(learner) {
  switch(learner,
    survreg_weibull = list(
      distribution = "weibull",
      scale = "estimated"
    ),
    survreg_exponential = list(
      distribution = "exponential",
      scale = "fixed_1"
    ),
    survreg_lognormal = list(
      distribution = "lognormal",
      scale = "estimated"
    ),
    survreg_loglogistic = list(
      distribution = "loglogistic",
      scale = "estimated"
    ),
    rpart_survival = list(
      method = "exp",
      cp = 0.01,
      minsplit = 20L,
      maxdepth = 30L,
      xval = 0L
    ),
    rsf_survival = list(
      num_trees = 1000L,
      mtry_rule = "randomForestSRC_default",
      min_node_size = 5L,
      splitrule = "logrank",
      seed = 1L
    ),
    ctree_survival = list(
      mincriterion = 0.95,
      minsplit = 20L,
      minbucket = 7L,
      maxdepth = Inf
    ),
    cforest_survival = list(
      num_trees = 500L,
      mtry_rule = "ceiling_sqrt_p",
      replace = TRUE,
      sample_fraction = 1,
      mincriterion = 0,
      minsplit = 20L,
      minbucket = 7L,
      maxdepth = Inf,
      seed = 1L
    ),
    bart_survival = list(
      type = "pbart",
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
    aorsf_survival = list(
      num_trees = 500L,
      num_split = 5L,
      num_retry = 3L,
      mtry_rule = "ceiling_sqrt_p",
      sample_with_replacement = TRUE,
      sample_fraction = 0.632,
      control_method = "glm",
      control_scale_x = TRUE,
      control_ties = "efron",
      control_max_iter = 1L,
      seed = 1L,
      threads = 1L,
      importance = "none",
      attach_data = FALSE
    ),
    sb_abort(
      paste0("Unknown fixed survival learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  )
}

sb_additional_survival_dependency <- function(learner) {
  packages <- switch(learner,
    survreg_weibull = "survival",
    survreg_exponential = "survival",
    survreg_lognormal = "survival",
    survreg_loglogistic = "survival",
    rpart_survival = "rpart",
    rsf_survival = "randomForestSRC",
    ctree_survival = c("partykit", "coin"),
    cforest_survival = c("partykit", "coin"),
    bart_survival = "BART",
    aorsf_survival = "aorsf"
  )
  list(
    package = packages[[1L]],
    packages = packages,
    available = all(vapply(
      packages,
      requireNamespace,
      logical(1),
      quietly = TRUE
    ))
  )
}

sb_require_additional_survival_dependency <- function(learner) {
  dependency <- sb_additional_survival_dependency(learner)
  if (!dependency$available) {
    missing <- dependency$packages[!vapply(
      dependency$packages,
      requireNamespace,
      logical(1),
      quietly = TRUE
    )]
    sb_abort(
      paste0(
        "Package", if (length(missing) > 1L) "s" else "", " `",
        paste(missing, collapse = "`, `"), "` ",
        if (length(missing) > 1L) "are" else "is",
        " required for learner `", learner, "`."
      ),
      "stormbreaker_missing_additional_survival_dependency"
    )
  }
  invisible(dependency)
}

sb_validate_additional_survival_inputs <- function(x, time, event, learner) {
  x <- sb_validate_predictor_frame(x, "x")
  time <- sb_validate_endpoint_vector(time, "time", nrow(x))
  event <- sb_validate_endpoint_vector(event, "event", nrow(x))
  if (any(time <= 0)) {
    sb_abort(
      paste0("Learner `", learner, "` requires strictly positive times."),
      "stormbreaker_additional_survival_nonpositive_time"
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
      paste0("Learner `", learner, "` requires at least one event."),
      "stormbreaker_additional_survival_no_event"
    )
  }
  list(x = x, time = time, event = as.integer(event))
}

sb_additional_survival_formula_data <- function(x, time, event,
                                                random_forest_src = FALSE) {
  analysis <- x
  time_name <- ".stormbreaker_time"
  while (time_name %in% names(analysis)) {
    time_name <- paste0(time_name, "_")
  }
  event_name <- ".stormbreaker_event"
  while (event_name %in% c(names(analysis), time_name)) {
    event_name <- paste0(event_name, "_")
  }
  analysis[[time_name]] <- time
  analysis[[event_name]] <- event
  rhs <- paste(vapply(names(x), sb_quote_name, ""), collapse = " + ")
  response <- if (random_forest_src) {
    paste0(
      "Surv(", sb_quote_name(time_name), ", ",
      sb_quote_name(event_name), ")"
    )
  } else {
    paste0(
      "survival::Surv(", sb_quote_name(time_name), ", ",
      sb_quote_name(event_name), ")"
    )
  }
  formula <- stats::as.formula(
    paste(response, "~", rhs),
    env = if (random_forest_src) asNamespace("survival") else baseenv()
  )
  list(data = analysis, formula = formula)
}

sb_new_additional_survival_model <- function(
  learner, features, params, engine = NULL, prediction = NULL,
  effective_mtry = NULL
) {
  structure(
    list(
      learner = learner,
      features = features,
      params = params,
      engine = engine,
      prediction = prediction,
      effective_mtry = effective_mtry
    ),
    class = "stormbreaker_additional_survival_fixed"
  )
}

sb_bart_survival_prediction_engine <- function(engine) {
  keep <- c(
    "treedraws", "rm.const", "binaryOffset", "offset", "type", "times", "K"
  )
  engine[setdiff(names(engine), keep)] <- NULL
  engine
}

sb_fit_bart_survival <- function(x, time, event, params) {
  engine <- NULL
  invisible(utils::capture.output(
    engine <- withr::with_seed(
      params$seed,
      withCallingHandlers(
        BART::surv.bart(
          x.train = as.matrix(x),
          times = time,
          delta = event,
          K = params$time_grid_size,
          sparse = params$sparse,
          type = params$type,
          k = params$prior_k,
          power = params$tree_power,
          base = params$tree_base,
          ntree = params$num_trees,
          numcut = params$num_cutpoints,
          ndpost = params$posterior_draws,
          nskip = params$burn_in,
          keepevery = params$keep_every,
          printevery = .Machine$integer.max,
          seed = params$seed,
          mc.cores = params$cores
        ),
        warning = function(warning) {
          sb_abort(
            paste0(
              "Fixed survival BART raised a fitting warning: ",
              conditionMessage(warning)
            ),
            "stormbreaker_bart_survival_fit_warning",
            parent = warning
          )
        }
      )
    ),
    type = "output"
  ))
  valid <- inherits(engine, "survbart") && !is.null(engine$treedraws) &&
    is.numeric(engine$times) && length(engine$times) > 0L &&
    !anyNA(engine$times) && all(is.finite(engine$times)) &&
    all(engine$times > 0) && identical(as.integer(engine$K),
      as.integer(length(engine$times)))
  if (!valid) {
    sb_abort(
      "Fixed survival BART did not return posterior trees and a time grid.",
      "stormbreaker_incomplete_bart_survival_fit"
    )
  }
  sb_bart_survival_prediction_engine(engine)
}

sb_fit_survreg <- function(formula_data, learner, params) {
  arguments <- list(
    formula = formula_data$formula,
    data = formula_data$data,
    dist = params$distribution,
    model = FALSE,
    x = FALSE,
    y = FALSE
  )
  if (!identical(params$distribution, "exponential")) {
    arguments$scale <- 0
  }
  withCallingHandlers(
    do.call(survival::survreg, arguments),
    warning = function(warning) {
      if (identical(
        conditionMessage(warning),
        "Ran out of iterations and did not converge"
      )) {
        sb_abort(
          paste0(
            "Fixed ", params$distribution, " SurvReg did not converge: ",
            conditionMessage(warning)
          ),
          "stormbreaker_survreg_nonconvergence",
          parent = warning
        )
      }
    }
  )
}

# The native default callback closes over rfsrc's full training frame. Keep the
# same default SWOR sample-size rule in the package namespace instead.
sb_rsf_native_sample_size <- function(x) x * 0.632

sb_fit_additional_survival <- function(x, time, event, learner, params) {
  inputs <- sb_validate_additional_survival_inputs(x, time, event, learner)
  sb_require_additional_survival_dependency(learner)
  formula_x <- inputs$x
  if (identical(learner, "aorsf_survival")) {
    names(formula_x) <- make.names(names(formula_x), unique = TRUE)
  }
  formula_data <- sb_additional_survival_formula_data(
    formula_x,
    inputs$time,
    inputs$event,
    random_forest_src = identical(learner, "rsf_survival")
  )

  engine <- if (learner %in% sb_survreg_learner_ids()) {
    sb_fit_survreg(formula_data, learner, params)
  } else {
    switch(learner,
    rpart_survival = rpart::rpart(
      formula_data$formula,
      data = formula_data$data,
      method = params$method,
      model = FALSE,
      x = FALSE,
      y = FALSE,
      control = rpart::rpart.control(
        cp = params$cp,
        minsplit = params$minsplit,
        maxdepth = params$maxdepth,
        xval = params$xval
      )
    ),
    rsf_survival = withr::with_options(
      list(rf.cores = 1L, mc.cores = 1L),
      randomForestSRC::rfsrc(
        formula_data$formula,
        data = formula_data$data,
        ntree = params$num_trees,
        nodesize = params$min_node_size,
        sampsize = sb_rsf_native_sample_size,
        splitrule = params$splitrule,
        importance = FALSE,
        proximity = FALSE,
        forest = TRUE,
        seed = params$seed
      )
    ),
    ctree_survival = partykit::ctree(
      formula_data$formula,
      data = formula_data$data,
      control = partykit::ctree_control(
        mincriterion = params$mincriterion,
        minsplit = params$minsplit,
        minbucket = params$minbucket,
        maxdepth = params$maxdepth
      )
    ),
    cforest_survival = withr::with_seed(
      params$seed,
      partykit::cforest(
        formula_data$formula,
        data = formula_data$data,
        ntree = params$num_trees,
        mtry = max(1L, ceiling(sqrt(ncol(inputs$x)))),
        perturb = list(
          replace = params$replace,
          fraction = params$sample_fraction
        ),
        control = partykit::ctree_control(
          teststat = "quad",
          testtype = "Univ",
          mincriterion = params$mincriterion,
          saveinfo = FALSE,
          minsplit = params$minsplit,
          minbucket = params$minbucket,
          maxdepth = params$maxdepth
        )
      )
    ),
    bart_survival = sb_fit_bart_survival(
      inputs$x,
      inputs$time,
      inputs$event,
      params
    ),
    aorsf_survival = {
      tree_seeds <- withr::with_seed(
        params$seed,
        sample.int(.Machine$integer.max, params$num_trees)
      )
      aorsf::orsf(
        data = formula_data$data,
        formula = formula_data$formula,
        control = aorsf::orsf_control_survival(
          method = params$control_method,
          scale_x = params$control_scale_x,
          ties = params$control_ties,
          max_iter = params$control_max_iter
        ),
        n_tree = params$num_trees,
        n_split = params$num_split,
        n_retry = params$num_retry,
        n_thread = params$threads,
        mtry = max(1L, ceiling(sqrt(ncol(inputs$x)))),
        sample_with_replacement = params$sample_with_replacement,
        sample_fraction = params$sample_fraction,
        importance = params$importance,
        tree_seeds = tree_seeds,
        attach_data = params$attach_data,
        na_action = "fail",
        verbose_progress = FALSE
      )
    }
    )
  }

  prediction <- NULL
  effective_mtry <- NULL
  if (learner %in% sb_survreg_learner_ids()) {
    coefficients <- stats::coef(engine)
    coefficient_names <- colnames(stats::model.matrix(
      formula_data$formula, data = formula_data$data
    ))
    valid <- is.numeric(coefficients) &&
      identical(names(coefficients), coefficient_names) &&
      !anyNA(coefficients) && all(is.finite(coefficients)) &&
      is.numeric(engine$scale) && length(engine$scale) == 1L &&
      !is.na(engine$scale) && is.finite(engine$scale) && engine$scale > 0
    if (!valid) {
      sb_abort(
        paste0(
          "Fixed ", params$distribution,
          " SurvReg did not return finite prediction parameters."
        ),
        "stormbreaker_incomplete_survreg_fit"
      )
    }
    prediction <- list(
      intercept = unname(coefficients[[1L]]),
      coefficient = unname(coefficients[-1L]),
      scale = unname(engine$scale)
    )
  } else if (identical(learner, "rpart_survival")) {
    if (!inherits(engine, "rpart")) {
      sb_abort(
        "Fixed survival Rpart did not return a fitted tree.",
        "stormbreaker_incomplete_rpart_survival_fit"
      )
    }
  } else if (identical(learner, "rsf_survival")) {
    if (!inherits(engine, "rfsrc") || is.null(engine$forest)) {
      sb_abort(
        "Fixed RSF did not return a prediction forest.",
        "stormbreaker_incomplete_rsf_fit"
      )
    }
    effective_mtry <- as.integer(engine$mtry)
  } else if (identical(learner, "ctree_survival")) {
    if (!inherits(engine, "constparty")) {
      sb_abort(
        "Fixed survival Ctree did not return a fitted tree.",
        "stormbreaker_incomplete_ctree_survival_fit"
      )
    }
  } else if (identical(learner, "bart_survival")) {
    if (!inherits(engine, "survbart") || is.null(engine$treedraws)) {
      sb_abort(
        "Fixed survival BART did not return a prediction engine.",
        "stormbreaker_incomplete_bart_survival_fit"
      )
    }
  } else if (identical(learner, "aorsf_survival")) {
    if (!inherits(engine, "ObliqueForestSurvival")) {
      sb_abort(
        "Fixed aorsf survival did not return an oblique survival forest.",
        "stormbreaker_incomplete_aorsf_survival_fit"
      )
    }
    effective_mtry <- max(1L, ceiling(sqrt(ncol(inputs$x))))
  } else if (!inherits(engine, "cforest")) {
    sb_abort(
      "Fixed survival CForest did not return a fitted forest.",
      "stormbreaker_incomplete_cforest_survival_fit"
    )
  } else {
    effective_mtry <- max(1L, ceiling(sqrt(ncol(inputs$x))))
  }

  model <- sb_new_additional_survival_model(
    learner,
    names(inputs$x),
    params,
    engine,
    prediction,
    effective_mtry
  )
  risk <- sb_predict_additional_survival(model, inputs$x)
  if (length(risk) != nrow(inputs$x)) {
    sb_abort(
      paste0("Learner `", learner, "` did not produce one fitted risk per row."),
      "stormbreaker_invalid_additional_survival_prediction"
    )
  }
  model
}

sb_validate_additional_survival_model_x <- function(model, x) {
  if (!inherits(model, "stormbreaker_additional_survival_fixed") ||
    !model$learner %in% sb_additional_survival_learner_ids()) {
    sb_abort(
      "The supplied model is not a fitted additional survival learner.",
      "stormbreaker_bad_additional_survival_fit"
    )
  }
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), model$features)) {
    sb_abort(
      "Prediction columns do not match the fitted survival learner schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  sb_require_additional_survival_dependency(model$learner)
  x
}

sb_survfit_last_cumulative_hazard <- function(fit) {
  value <- fit$cumhaz
  if (is.null(value) || length(value) == 0L) {
    return(0)
  }
  unname(as.numeric(value[[length(value)]]))
}

sb_partykit_survival_distributions <- function(model, x) {
  batch_size <- if (inherits(model$engine, "cforest")) {
    nrow(model$engine$data)
  } else {
    nrow(x)
  }
  # Bound native forest-weight allocation when assessment exceeds training size.
  # Each observation uses the same fitted forest and native probability method.
  if (nrow(x) > batch_size) {
    batches <- split(seq_len(nrow(x)), ceiling(seq_len(nrow(x)) / batch_size))
    distributions <- unlist(lapply(batches, function(ii) {
      stats::predict(model$engine, newdata = x[ii, , drop = FALSE], type = "prob")
    }), recursive = FALSE, use.names = FALSE)
  } else {
    distributions <- stats::predict(model$engine, newdata = x, type = "prob")
  }
  valid <- is.list(distributions) && length(distributions) == nrow(x) &&
    all(vapply(distributions, inherits, logical(1), what = "survfit"))
  if (!valid) {
    sb_abort(
      "The partykit survival learner produced invalid distributions.",
      "stormbreaker_invalid_additional_survival_prediction"
    )
  }
  distributions
}

sb_rsf_survival_curves <- function(model, x) {
  prediction <- withr::with_options(
    list(rf.cores = 1L, mc.cores = 1L),
    stats::predict(model$engine, newdata = x)
  )
  risk <- unname(as.numeric(prediction$predicted))
  survival <- as.matrix(prediction$survival)
  time <- as.numeric(prediction$time.interest)
  valid <- length(risk) == nrow(x) && !anyNA(risk) && all(is.finite(risk)) &&
    nrow(survival) == nrow(x) && ncol(survival) == length(time) &&
    length(time) > 0L && !anyNA(time) && all(is.finite(time)) &&
    !anyNA(survival) && all(is.finite(survival)) &&
    all(survival >= 0 & survival <= 1)
  if (!valid) {
    sb_abort(
      "Fixed RSF produced invalid mortality or survival predictions.",
      "stormbreaker_invalid_additional_survival_prediction"
    )
  }
  list(risk = risk, time = time, survival = survival)
}

sb_bart_survival_curves <- function(model, x) {
  x <- sb_validate_additional_survival_model_x(model, x)
  grid <- unname(as.numeric(model$engine$times))
  n <- nrow(x)
  k <- length(grid)
  matrix_x <- as.matrix(x)
  augmented <- cbind(
    t = rep(grid, times = n),
    matrix_x[rep(seq_len(n), each = k), , drop = FALSE]
  )
  removed <- model$engine$rm.const
  if (length(removed) > 0L) {
    augmented <- augmented[, removed, drop = FALSE]
  }
  prediction <- NULL
  predict_method <- get(
    "predict.survbart",
    envir = asNamespace("BART"),
    inherits = FALSE
  )
  invisible(utils::capture.output(
    prediction <- predict_method(
      model$engine,
      augmented,
      mc.cores = 1L,
      openmp = FALSE
    ),
    type = "output"
  ))
  survival <- matrix(
    unname(as.numeric(prediction$surv.test.mean)),
    nrow = n,
    ncol = k,
    byrow = TRUE
  )
  valid <- length(grid) > 0L && !anyNA(grid) && all(is.finite(grid)) &&
    all(grid > 0) && !is.unsorted(grid) &&
    nrow(survival) == n && ncol(survival) == k &&
    !anyNA(survival) && all(is.finite(survival)) &&
    all(survival >= 0 & survival <= 1)
  if (!valid) {
    sb_abort(
      "Fixed survival BART produced invalid posterior survival curves.",
      "stormbreaker_invalid_additional_survival_prediction"
    )
  }
  list(time = grid, survival = survival)
}

sb_bart_survival_negative_rmst <- function(curves) {
  interval <- diff(c(0, curves$time))
  previous <- cbind(
    1,
    curves$survival[, -ncol(curves$survival), drop = FALSE]
  )
  -as.vector(previous %*% interval)
}

sb_additional_survival_prediction_details <- function(model, x) {
  x <- sb_validate_additional_survival_model_x(model, x)
  details <- if (identical(model$learner, "bart_survival")) {
    list(curves = sb_bart_survival_curves(model, x))
  } else {
    list(distributions = sb_partykit_survival_distributions(model, x))
  }
  details$risk <- sb_predict_additional_survival(model, x, details)
  details
}

sb_predict_additional_survival <- function(model, x, prediction_details = NULL) {
  x <- sb_validate_additional_survival_model_x(model, x)
  learner <- model$learner
  risk <- if (learner %in% sb_survreg_learner_ids()) {
    location <- model$prediction$intercept +
      as.vector(as.matrix(x) %*% model$prediction$coefficient)
    -location
  } else {
    switch(learner,
    rpart_survival = {
      rate <- as.numeric(stats::predict(model$engine, newdata = x, type = "vector"))
      if (anyNA(rate) || any(!is.finite(rate)) || any(rate <= 0)) {
        sb_abort(
          "Fixed survival Rpart produced invalid relative rates.",
          "stormbreaker_invalid_additional_survival_prediction"
        )
      }
      log(rate)
    },
    rsf_survival = sb_rsf_survival_curves(model, x)$risk,
    ctree_survival = vapply(
      sb_partykit_survival_distributions(model, x),
      sb_survfit_last_cumulative_hazard,
      numeric(1)
    ),
    cforest_survival = vapply(
      if (is.null(prediction_details)) sb_partykit_survival_distributions(model, x) else
        prediction_details$distributions,
      sb_survfit_last_cumulative_hazard,
      numeric(1)
    ),
    bart_survival = sb_bart_survival_negative_rmst(
      if (is.null(prediction_details)) sb_bart_survival_curves(model, x) else
        prediction_details$curves
    ),
    aorsf_survival = {
      names(x) <- model$engine$get_names_x()
      stats::predict(
        model$engine,
        new_data = x,
        pred_type = "mort",
        n_thread = model$params$threads
      )
    }
    )
  }
  risk <- unname(as.numeric(risk))
  if (length(risk) != nrow(x) || anyNA(risk) || any(!is.finite(risk))) {
    sb_abort(
      paste0("Learner `", learner, "` produced invalid risk predictions."),
      "stormbreaker_invalid_additional_survival_prediction"
    )
  }
  risk
}

sb_survival_curve_probability_at_horizon <- function(time, survival,
                                                     horizon) {
  eligible <- which(time <= horizon)
  if (length(eligible) == 0L) {
    return(1)
  }
  unname(as.numeric(survival[[max(eligible)]]))
}

sb_additional_direct_survival_probability <- function(model, x, horizon,
                                                      prediction_details = NULL) {
  x <- sb_validate_additional_survival_model_x(model, x)
  learner <- model$learner
  probability <- if (learner %in% sb_survreg_learner_ids()) {
    location <- model$prediction$intercept +
      as.vector(as.matrix(x) %*% model$prediction$coefficient)
    transformed <- (log(horizon) - location) / model$prediction$scale
    switch(model$params$distribution,
      weibull = exp(-exp(transformed)),
      exponential = exp(-exp(transformed)),
      lognormal = stats::pnorm(transformed, lower.tail = FALSE),
      loglogistic = stats::plogis(transformed, lower.tail = FALSE)
    )
  } else {
    switch(learner,
    rsf_survival = {
      curves <- sb_rsf_survival_curves(model, x)
      eligible <- which(curves$time <= horizon)
      if (length(eligible) == 0L) {
        rep(1, nrow(x))
      } else {
        unname(curves$survival[, max(eligible)])
      }
    },
    ctree_survival = vapply(
      sb_partykit_survival_distributions(model, x),
      function(distribution) {
        sb_survival_curve_probability_at_horizon(
          distribution$time,
          distribution$surv,
          horizon
        )
      },
      numeric(1)
    ),
    bart_survival = {
      curves <- if (is.null(prediction_details)) sb_bart_survival_curves(model, x) else
        prediction_details$curves
      eligible <- which(curves$time <= horizon)
      if (length(eligible) == 0L) {
        rep(1, nrow(x))
      } else {
        unname(curves$survival[, max(eligible)])
      }
    },
    cforest_survival = vapply(
      if (is.null(prediction_details)) sb_partykit_survival_distributions(model, x) else
        prediction_details$distributions,
      function(distribution) {
        sb_survival_curve_probability_at_horizon(
          distribution$time,
          distribution$surv,
          horizon
        )
      },
      numeric(1)
    ),
    aorsf_survival = {
      names(x) <- model$engine$get_names_x()
      stats::predict(
        model$engine,
        new_data = x,
        pred_type = "surv",
        pred_horizon = horizon,
        n_thread = model$params$threads
      )
    },
    sb_abort(
      "This learner requires Breslow calibration rather than direct survival.",
      "stormbreaker_bad_additional_survival_fit"
    )
    )
  }
  probability <- unname(as.numeric(probability))
  valid <- length(probability) == nrow(x) && !anyNA(probability) &&
    all(is.finite(probability)) &&
    all(probability >= 0 & probability <= 1)
  if (!valid) {
    sb_abort(
      paste0("Learner `", learner, "` produced invalid horizon survival."),
      "stormbreaker_invalid_additional_survival_prediction"
    )
  }
  probability
}

sb_additional_direct_survival_source <- function(learner) {
  if (learner %in% sb_survreg_learner_ids()) {
    params <- sb_additional_survival_params(learner)
    return(paste0(
      "survreg_", params$distribution, "_distribution_direct"
    ))
  }
  switch(learner,
    rsf_survival = "randomForestSRC_survival_forest_direct",
    ctree_survival = "partykit_ctree_distribution_direct",
    cforest_survival = "partykit_cforest_distribution_direct",
    bart_survival = "BART_discrete_time_posterior_survival_direct",
    aorsf_survival = "aorsf_oblique_survival_forest_direct",
    NULL
  )
}

sb_predict_additional_survival_survival <- function(
  model, x_analysis, time, event, x_assessment, prediction_horizon,
  prediction_details = NULL
) {
  inputs <- sb_validate_additional_survival_inputs(
    x_analysis,
    time,
    event,
    model$learner
  )
  if (!identical(names(inputs$x), model$features)) {
    sb_abort(
      "Analysis columns do not match the fitted survival learner schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  prediction_horizon <- sb_validate_prediction_horizon(
    prediction_horizon,
    inputs$time
  )
  if (identical(model$learner, "rpart_survival")) {
    analysis_risk <- sb_predict_additional_survival(model, inputs$x)
    assessment_risk <- sb_predict_additional_survival(
      model,
      x_assessment
    )
    baseline <- sb_breslow_baseline_hazard(
      inputs$time,
      inputs$event,
      analysis_risk,
      prediction_horizon
    )
    return(sb_survival_prediction_result(
      prediction_horizon,
      baseline,
      source = "rpart_survival_breslow_calibration",
      ties = "breslow",
      risk = assessment_risk
    ))
  }
  list(
    prediction_horizon = prediction_horizon,
    source = sb_additional_direct_survival_source(model$learner),
    ties = "not_applicable",
    probability = sb_additional_direct_survival_probability(
      model,
      x_assessment,
      prediction_horizon,
      prediction_details = prediction_details
    )
  )
}

sb_additional_survival_prediction_engine <- function(learner, engine) {
  if (identical(learner, "rpart_survival")) {
    engine[intersect(c("model", "x", "y", "where"), names(engine))] <- NULL
    return(engine)
  }
  if (identical(learner, "bart_survival")) {
    return(sb_bart_survival_prediction_engine(engine))
  }
  if (identical(learner, "aorsf_survival")) {
    return(engine)
  }
  if (!identical(learner, "rsf_survival")) {
    return(engine)
  }
  drop_fields <- c(
    "call", "yvar", "xvar", "subj", "proximity", "forest.wt",
    "case.depth", "distance", "membership", "tdc.membership", "inbag",
    "var.used", "imputed.indv", "imputed.data", "split.depth",
    "node.stats", "ensemble", "holdout.array", "empr.risk",
    "oob.empr.risk", "chf", "chf.oob", "predicted", "predicted.oob",
    "hazard", "hazard.oob", "survival", "survival.oob", "cif", "cif.oob",
    "err.rate"
  )
  engine[intersect(drop_fields, names(engine))] <- NULL
  engine
}

sb_extract_additional_survival_learner_state <- function(
  model, context = list()
) {
  if (!inherits(model, "stormbreaker_additional_survival_fixed")) {
    sb_bad_learner_state(
      "The fitted additional survival learner state is invalid."
    )
  }
  baseline_hazard <- NULL
  survival_probability <- NULL
  if ("prediction_horizon" %in% names(context)) {
    probability <- sb_predict_additional_survival_survival(
      model,
      context$x,
      context$time,
      context$event,
      context$x[1L, , drop = FALSE],
      context$prediction_horizon
    )
    if (identical(model$learner, "rpart_survival")) {
      baseline_hazard <- list(
        prediction_horizon = probability$prediction_horizon,
        time_unit = if ("time_unit" %in% names(context)) {
          context$time_unit
        } else {
          NA_character_
        },
        cumulative_baseline_hazard = probability$cumulative_baseline_hazard,
        baseline_survival = probability$baseline_survival,
        source = probability$source,
        ties = probability$ties
      )
    } else {
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
  }

  engine_raw <- NULL
  is_survreg <- model$learner %in% sb_survreg_learner_ids()
  if (!is_survreg) {
    engine <- sb_additional_survival_prediction_engine(
      model$learner,
      model$engine
    )
    engine_raw <- serialize(
      engine,
      connection = NULL,
      ascii = FALSE,
      xdr = TRUE,
      version = 3L
    )
  }
  state <- list(
    learner = model$learner,
    features = model$features,
    params = model$params,
    effective_mtry = model$effective_mtry,
    risk_reference = if (is_survreg) {
      paste0(model$params$distribution, "_aft_native")
    } else switch(model$learner,
      rpart_survival = "rpart_native",
      rsf_survival = "randomForestSRC_native",
      ctree_survival = "partykit_native",
      cforest_survival = "partykit_native",
      bart_survival = "BART_discrete_time_posterior",
      aorsf_survival = "aorsf_native"
    ),
    risk_formula = if (is_survreg) {
      "negative_log_time_location"
    } else switch(model$learner,
      rpart_survival = "log_relative_event_rate",
      rsf_survival = "mortality_sum_cumulative_hazard",
      ctree_survival = "terminal_distribution_cumulative_hazard",
      cforest_survival = "ensemble_distribution_cumulative_hazard",
      bart_survival = "negative_step_RMST_over_analysis_time_grid",
      aorsf_survival = "native_mortality_over_event_times"
    ),
    baseline_hazard = baseline_hazard,
    survival_probability = survival_probability
  )
  if (is_survreg) {
    state$intercept <- model$prediction$intercept
    state$coefficient <- model$prediction$coefficient
    state$scale <- model$prediction$scale
  } else {
    state$engine_raw <- engine_raw
  }
  if (model$learner %in% c(
    "ctree_survival",
    "cforest_survival"
  )) {
    state$training_data_retention <-
      "required_by_partykit_conditional_distribution_state"
  } else if (identical(model$learner, "rsf_survival")) {
    state$training_data_retention <-
      "event_time_metadata_required_by_randomForestSRC_forest"
  }
  structure(
    state,
    class = c(
      "stormbreaker_additional_survival_fixed_learner_state",
      "stormbreaker_learner_state"
    )
  )
}

sb_restore_additional_survival_model <- function(state, learner) {
  if (learner %in% sb_survreg_learner_ids()) {
    valid <- is.numeric(state$intercept) && length(state$intercept) == 1L &&
      !is.na(state$intercept) && is.finite(state$intercept) &&
      is.numeric(state$coefficient) &&
      length(state$coefficient) == length(state$features) &&
      !anyNA(state$coefficient) && all(is.finite(state$coefficient)) &&
      is.numeric(state$scale) && length(state$scale) == 1L &&
      !is.na(state$scale) && is.finite(state$scale) && state$scale > 0
    if (!valid) {
      sb_bad_learner_state("The saved fixed SurvReg state is invalid.")
    }
    return(sb_new_additional_survival_model(
      learner,
      state$features,
      state$params,
      prediction = list(
        intercept = state$intercept,
        coefficient = state$coefficient,
        scale = state$scale
      )
    ))
  }
  if (!is.raw(state$engine_raw) || length(state$engine_raw) == 0L) {
    sb_bad_learner_state("The saved survival learner engine is invalid.")
  }
  engine <- tryCatch(
    unserialize(state$engine_raw),
    error = function(error) {
      sb_abort(
        paste0(
          "The saved `", learner, "` engine could not be loaded: ",
          conditionMessage(error)
        ),
        "stormbreaker_invalid_additional_survival_state",
        parent = error
      )
    }
  )
  valid_engine <- switch(learner,
    rpart_survival = inherits(engine, "rpart"),
    rsf_survival = inherits(engine, "rfsrc"),
    ctree_survival = inherits(engine, "constparty"),
    cforest_survival = inherits(engine, "cforest"),
    bart_survival = inherits(engine, "survbart"),
    aorsf_survival = inherits(engine, "ObliqueForestSurvival")
  )
  if (!valid_engine) {
    sb_bad_learner_state("The saved survival learner engine is invalid.")
  }
  sb_new_additional_survival_model(
    learner,
    state$features,
    state$params,
    engine = engine,
    effective_mtry = state$effective_mtry
  )
}

sb_validate_additional_survival_probability_state <- function(state, learner) {
  probability <- state$survival_probability
  expected_source <- sb_additional_direct_survival_source(learner)
  valid <- is.list(probability) &&
    is.numeric(probability$prediction_horizon) &&
    length(probability$prediction_horizon) == 1L &&
    !is.na(probability$prediction_horizon) &&
    is.finite(probability$prediction_horizon) &&
    probability$prediction_horizon > 0 &&
    is.character(probability$time_unit) &&
    length(probability$time_unit) == 1L &&
    !is.na(probability$time_unit) && nzchar(probability$time_unit) &&
    identical(probability$source, expected_source) &&
    identical(probability$ties, "not_applicable")
  if (!valid) {
    sb_bad_learner_state(
      "The saved direct survival-probability state is invalid."
    )
  }
  probability
}

sb_predict_additional_survival_learner_state <- function(state, x) {
  learner <- if (is.list(state)) state[["learner"]] else NULL
  if (!is.character(learner) || length(learner) != 1L || is.na(learner) ||
    !learner %in% sb_additional_survival_learner_ids()) {
    sb_bad_learner_state(
      "The saved additional survival learner state is invalid."
    )
  }
  x <- sb_validate_learner_state_x(state, x, learner)
  sb_require_additional_survival_dependency(learner)
  model <- sb_restore_additional_survival_model(state, learner)
  prediction_details <- NULL
  if (!is.null(state$survival_probability) && learner %in%
      c("bart_survival", "cforest_survival")) {
    warned <- FALSE
    prediction_details <- withCallingHandlers(
      sb_additional_survival_prediction_details(model, x),
      warning = function(w) warned <<- TRUE)
    risk <- prediction_details$risk
    if (warned) prediction_details <- NULL
  } else {
    risk <- sb_predict_additional_survival(model, x)
  }
  result <- list(risk = risk)
  if (!is.null(state$baseline_hazard)) {
    baseline <- state$baseline_hazard
    valid <- identical(learner, "rpart_survival") &&
      is.list(baseline) &&
      is.numeric(baseline$prediction_horizon) &&
      length(baseline$prediction_horizon) == 1L &&
      !is.na(baseline$prediction_horizon) &&
      is.finite(baseline$prediction_horizon) &&
      baseline$prediction_horizon > 0 &&
      is.character(baseline$time_unit) && length(baseline$time_unit) == 1L &&
      !is.na(baseline$time_unit) && nzchar(baseline$time_unit)
    if (!valid) {
      sb_bad_learner_state(
        "The saved Rpart survival baseline-hazard state is invalid."
      )
    }
    result$prediction_horizon <- baseline$prediction_horizon
    result$time_unit <- baseline$time_unit
    result$survival_probability <- sb_survival_probability_from_hazard(
      risk,
      baseline$cumulative_baseline_hazard
    )
  } else if (!is.null(state$survival_probability)) {
    probability <- sb_validate_additional_survival_probability_state(
      state,
      learner
    )
    result$prediction_horizon <- probability$prediction_horizon
    result$time_unit <- probability$time_unit
    result$survival_probability <- sb_additional_direct_survival_probability(
      model,
      x,
      probability$prediction_horizon,
      prediction_details = prediction_details
    )
  }
  result
}

sb_additional_survival_learner_adapter <- function(
  learner,
  params = sb_additional_survival_params(learner)
) {
  if (!learner %in% sb_additional_survival_learner_ids()) {
    sb_abort(
      paste0("Unknown fixed survival learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  }
  list(
    id = learner,
    params = params,
    fit = function(x, time, event) {
      sb_fit_additional_survival(x, time, event, learner, params)
    },
    predict = sb_predict_additional_survival,
    predict_details = if (learner %in%
      c("bart_survival", "cforest_survival")) {
      sb_additional_survival_prediction_details
    } else NULL,
    predict_survival = sb_predict_additional_survival_survival,
    extract_state = function(model, context = list()) {
      sb_extract_additional_survival_learner_state(model, context)
    },
    predict_state = sb_predict_additional_survival_learner_state
  )
}
