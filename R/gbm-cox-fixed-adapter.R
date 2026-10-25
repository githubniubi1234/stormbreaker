sb_gbm_available <- function() {
  requireNamespace("gbm", quietly = TRUE)
}

sb_validate_gbm_shrinkage <- function(shrinkage) {
  valid <- is.numeric(shrinkage) && length(shrinkage) == 1L &&
    !is.na(shrinkage) && is.finite(shrinkage) && shrinkage > 0
  if (!valid) {
    sb_abort(
      "`gbm_shrinkage` must be one finite positive number.",
      "stormbreaker_bad_gbm_shrinkage"
    )
  }
  as.numeric(shrinkage)
}

sb_validate_gbm_bag_fraction <- function(bag_fraction) {
  valid <- is.numeric(bag_fraction) && length(bag_fraction) == 1L &&
    !is.na(bag_fraction) && is.finite(bag_fraction) &&
    bag_fraction > 0 && bag_fraction <= 1
  if (!valid) {
    sb_abort(
      "`gbm_bag_fraction` must be one finite number in (0, 1].",
      "stormbreaker_bad_gbm_bag_fraction"
    )
  }
  as.numeric(bag_fraction)
}

sb_gbm_cox_params <- function(
    n_trees = 100L,
    interaction_depth = 1L,
    min_obs_node = 10L,
    shrinkage = 0.1,
    bag_fraction = 1,
    seed = 1L) {
  list(
    n_trees = sb_check_positive_integer(
      n_trees,
      "gbm_n_trees",
      minimum = 1L
    ),
    interaction_depth = sb_check_positive_integer(
      interaction_depth,
      "gbm_interaction_depth",
      minimum = 1L
    ),
    min_obs_node = sb_check_positive_integer(
      min_obs_node,
      "gbm_min_obs_node",
      minimum = 1L
    ),
    shrinkage = sb_validate_gbm_shrinkage(shrinkage),
    bag_fraction = sb_validate_gbm_bag_fraction(bag_fraction),
    seed = sb_check_positive_integer(seed, "gbm_seed", minimum = 1L),
    distribution = "coxph",
    keep_data = FALSE
  )
}

sb_gbm_cox_cv_params <- function(inner_folds = 10L, seed = 1L,
                                 maximum_trees = 10000L) {
  list(
    inner_folds = sb_check_positive_integer(
      inner_folds,
      "learner_tuning_folds",
      minimum = 3L
    ),
    seed = sb_check_positive_integer(seed, "learner_seed", minimum = 1L),
    maximum_trees = sb_check_positive_integer(
      maximum_trees,
      "gbm_cv_max_trees",
      minimum = 1L
    ),
    interaction_depth = 3L,
    min_obs_node = 10L,
    shrinkage = 0.001,
    bag_fraction = 0.5,
    distribution = "coxph",
    keep_data = FALSE,
    fit_interface = "formula"
  )
}

sb_fit_gbm_cox <- function(x, time, event, params) {
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
      "A fixed GBM-Cox analysis fold must contain at least one event.",
      "stormbreaker_gbm_cox_no_event"
    )
  }
  if (!sb_gbm_available()) {
    sb_abort(
      "Package `gbm` is required for the fixed GBM-Cox learner.",
      "stormbreaker_missing_gbm"
    )
  }

  booster <- withr::with_seed(
    params$seed,
    gbm::gbm.fit(
      x = as.matrix(x),
      y = survival::Surv(time, event),
      distribution = params$distribution,
      n.trees = params$n_trees,
      interaction.depth = params$interaction_depth,
      n.minobsinnode = params$min_obs_node,
      shrinkage = params$shrinkage,
      bag.fraction = params$bag_fraction,
      nTrain = nrow(x),
      keep.data = params$keep_data,
      verbose = FALSE,
      var.names = names(x),
      response.name = "Surv(time, event)"
    )
  )
  if (!inherits(booster, "gbm") ||
      !identical(as.integer(booster$n.trees), params$n_trees)) {
    sb_abort(
      "Fixed GBM-Cox did not fit the declared number of trees.",
      "stormbreaker_incomplete_gbm_fit"
    )
  }

  structure(
    list(
      features = names(x),
      params = params,
      n_trees = params$n_trees,
      booster = booster
    ),
    class = "stormbreaker_gbm_cox"
  )
}

sb_fit_gbm_cox_cv <- function(x, time, event, params) {
  x <- sb_validate_predictor_frame(x, "x")
  time <- sb_validate_endpoint_vector(time, "time", nrow(x))
  event <- sb_validate_endpoint_vector(event, "event", nrow(x))
  if (any(time < 0)) {
    sb_abort(
      "`time` must contain non-negative values.",
      "stormbreaker_bad_adapter_time"
    )
  }
  if (!all(event %in% c(0, 1)) || !any(event == 1L)) {
    sb_abort(
      "A CV-tuned GBM-Cox analysis fold requires 0/1 events and an event.",
      "stormbreaker_gbm_cox_no_event"
    )
  }
  if (!sb_gbm_available()) {
    sb_abort(
      "Package `gbm` is required for the CV-tuned GBM-Cox learner.",
      "stormbreaker_missing_gbm"
    )
  }
  formula_data <- sb_additional_survival_formula_data(x, time, event)
  fit_gbm <- function(n_trees, cv_folds = params$inner_folds, search = NULL) {
    native <- gbm::gbm
    if (!is.null(search)) {
      reuse_cv <- function(...) {
        args <- list(...)
        ns <- asNamespace("gbm")
        i.train <- seq_len(args$nTrain)
        cv.group <- get("getCVgroup", ns)(args$distribution,
          args$class.stratify.cv, args$y, i.train, args$cv.folds, args$group)
        seeds <- as.integer(stats::runif(args$cv.folds, -(2^31 - 1), 2^31))
        # Serial native CV resets the seed for each fold. Its last fold's
        # terminal RNG is needed by the subsequent full-analysis booster.
        do_fold <- get("gbmDoFold", ns)
        last <- args[intersect(names(args), names(formals(do_fold)))]
        last$X <- args$cv.folds
        last$i.train <- i.train
        last$cv.group <- cv.group
        last$s <- seeds
        invisible(do.call(do_fold, last))
        list(error = search$cv.error[seq_len(args$n.trees)],
          predictions = search$cv.fitted)
      }
      environment(native) <- list2env(list(gbmCrossVal = reuse_cv),
        parent = environment(native))
    }
    booster <- NULL
    invisible(utils::capture.output(
      booster <- withr::with_seed(
        params$seed,
        native(
          formula = formula_data$formula,
          data = formula_data$data,
          distribution = params$distribution,
          n.trees = n_trees,
          interaction.depth = params$interaction_depth,
          n.minobsinnode = params$min_obs_node,
          shrinkage = params$shrinkage,
          bag.fraction = params$bag_fraction,
          cv.folds = cv_folds,
          n.cores = 1L,
          keep.data = params$keep_data,
          verbose = FALSE
        )
      )
    ))
    # Keep the same recorded native call for full and reused CV execution.
    booster$call[[1L]] <- quote(gbm::gbm)
    booster$m[[1L]] <- quote(gbm::gbm)
    booster
  }
  single_candidate <- params$maximum_trees == 1L
  # There is no tree-count choice with one candidate. Native gbm's one-row
  # CV reduction also fails; fit that candidate and record that no search ran.
  search_warned <- FALSE
  search <- if (single_candidate) NULL else withCallingHandlers(
    fit_gbm(params$maximum_trees),
    warning = function(w) search_warned <<- TRUE
  )
  best_trees <- if (single_candidate) 1L else as.integer(which.min(search$cv.error))
  if (length(best_trees) != 1L || is.na(best_trees) || best_trees < 1L ||
    best_trees > params$maximum_trees) {
    sb_abort(
      "CV-tuned GBM-Cox returned an invalid tree count.",
      "stormbreaker_invalid_gbm_cv"
    )
  }
  # Native gbm CV drops its error matrix when n.trees is one. The tree count
  # is already settled above; the final one-tree fit needs no further CV.
  booster <- if (best_trees == 1L) {
    fit_gbm(best_trees, cv_folds = 0L)
  } else if (search_warned) {
    fit_gbm(best_trees)
  } else if (best_trees == params$maximum_trees) {
    search
  } else {
    fit_gbm(best_trees, search = search)
  }
  if (!inherits(booster, "gbm") ||
    !identical(as.integer(booster$n.trees), best_trees)) {
    sb_abort(
      "CV-tuned GBM-Cox did not fit the selected number of trees.",
      "stormbreaker_incomplete_gbm_fit"
    )
  }
  tuning <- list(selected_trees = best_trees, cv_error = as.numeric(search$cv.error))
  if (single_candidate) tuning$search_status <- "single_candidate_no_search"
  structure(
    list(
      learner = "gbm_cox_cv",
      features = names(x),
      params = params,
      n_trees = best_trees,
      tuning = tuning,
      booster = booster
    ),
    class = "stormbreaker_gbm_cox"
  )
}

sb_predict_gbm_cox <- function(model, x) {
  if (!inherits(model, "stormbreaker_gbm_cox")) {
    sb_abort(
      "The supplied model is not a fitted fixed GBM-Cox learner.",
      "stormbreaker_bad_gbm_fit"
    )
  }
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), model$features)) {
    sb_abort(
      "Prediction columns do not match the fixed GBM-Cox schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  if (!sb_gbm_available()) {
    sb_abort(
      "Package `gbm` is required for fixed GBM-Cox prediction.",
      "stormbreaker_missing_gbm"
    )
  }

  risk <- as.numeric(stats::predict(
    model$booster,
    newdata = if (identical(model$params$fit_interface, "formula")) {
      x
    } else {
      as.matrix(x)
    },
    n.trees = model$n_trees,
    type = "link"
  ))
  if (length(risk) != nrow(x) || anyNA(risk) || any(!is.finite(risk))) {
    sb_abort(
      "Fixed GBM-Cox produced invalid log-hazard predictions.",
      "stormbreaker_invalid_gbm_prediction"
    )
  }
  unname(risk)
}

sb_predict_gbm_cox_survival <- function(
    model, x_analysis, time, event, x_assessment, prediction_horizon) {
  x_analysis <- sb_validate_predictor_frame(x_analysis, "x_analysis")
  time <- sb_validate_endpoint_vector(time, "time", nrow(x_analysis))
  event <- sb_validate_endpoint_vector(event, "event", nrow(x_analysis))
  if (!all(event %in% c(0, 1))) {
    sb_abort(
      "`event` must contain only 0/1 values.",
      "stormbreaker_bad_adapter_event"
    )
  }
  prediction_horizon <- sb_validate_prediction_horizon(
    prediction_horizon,
    time
  )
  analysis_risk <- sb_predict_gbm_cox(model, x_analysis)
  assessment_risk <- sb_predict_gbm_cox(model, x_assessment)
  cumulative_hazard <- sb_breslow_baseline_hazard(
    time,
    as.integer(event),
    analysis_risk,
    prediction_horizon
  )
  sb_survival_prediction_result(
    prediction_horizon,
    cumulative_hazard,
    source = paste0(
      if ("learner" %in% names(model)) model$learner else "gbm_cox",
      "_breslow_log_hazard"
    ),
    ties = "breslow",
    risk = assessment_risk
  )
}

sb_extract_gbm_cox_tree_model <- function(model) {
  if (!inherits(model, "stormbreaker_gbm_cox")) {
    sb_abort(
      "The supplied model is not a fitted fixed GBM-Cox learner.",
      "stormbreaker_bad_gbm_fit"
    )
  }
  if (!sb_gbm_available()) {
    sb_abort(
      "Package `gbm` is required to audit a GBM-Cox model.",
      "stormbreaker_missing_gbm"
    )
  }

  importance <- as.numeric(gbm::relative.influence(
    model$booster,
    n.trees = model$n_trees,
    scale. = FALSE,
    sort. = FALSE
  ))
  if (length(importance) != length(model$features) ||
      anyNA(importance) || any(!is.finite(importance))) {
    sb_abort(
      "Fixed GBM-Cox did not return one finite importance per feature.",
      "stormbreaker_invalid_gbm_importance"
    )
  }

  list(
    n_trees = model$n_trees,
    features = model$features,
    params = model$params,
    tuning = if ("tuning" %in% names(model)) model$tuning else NULL,
    importance_type = "relative_influence",
    model_serialization = "R_xdr_v3_gbm",
    model_raw = serialize(
      model$booster,
      connection = NULL,
      ascii = FALSE,
      xdr = TRUE,
      version = 3L
    ),
    feature_importance = data.frame(
      feature = model$features,
      importance = importance,
      stringsAsFactors = FALSE
    )
  )
}

sb_gbm_cox_learner_adapter <- function(
    n_trees = 100L,
    interaction_depth = 1L,
    min_obs_node = 10L,
    shrinkage = 0.1,
    bag_fraction = 1,
    seed = 1L) {
  params <- sb_gbm_cox_params(
    n_trees = n_trees,
    interaction_depth = interaction_depth,
    min_obs_node = min_obs_node,
    shrinkage = shrinkage,
    bag_fraction = bag_fraction,
    seed = seed
  )
  list(
    id = "gbm_cox",
    params = params,
    fit = function(x, time, event) {
      sb_fit_gbm_cox(x, time, event, params)
    },
    predict = sb_predict_gbm_cox,
    predict_survival = sb_predict_gbm_cox_survival,
    extract_tree_model = sb_extract_gbm_cox_tree_model,
    extract_state = function(model, context = list()) {
      sb_extract_gbm_cox_learner_state(model, context)
    },
    predict_state = sb_predict_gbm_cox_learner_state
  )
}

sb_gbm_cox_cv_learner_adapter <- function(inner_folds = 10L,
                                           seed = 1L,
                                           maximum_trees = 10000L) {
  params <- sb_gbm_cox_cv_params(
    inner_folds = inner_folds,
    seed = seed,
    maximum_trees = maximum_trees
  )
  list(
    id = "gbm_cox_cv",
    params = params,
    fit = function(x, time, event) {
      sb_fit_gbm_cox_cv(x, time, event, params)
    },
    predict = sb_predict_gbm_cox,
    predict_survival = sb_predict_gbm_cox_survival,
    extract_tree_model = sb_extract_gbm_cox_tree_model,
    extract_state = function(model, context = list()) {
      sb_extract_gbm_cox_learner_state(model, context)
    },
    predict_state = sb_predict_gbm_cox_learner_state
  )
}
