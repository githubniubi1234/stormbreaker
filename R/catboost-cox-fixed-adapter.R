sb_catboost_cox_params <- function(
    iterations = 100L,
    learning_rate = 0.1,
    depth = 6L,
    l2_leaf_reg = 3,
    seed = 1L) {
  list(
    iterations = sb_check_positive_integer(
      iterations,
      "catboost_iterations",
      minimum = 1L
    ),
    learning_rate = sb_validate_catboost_learning_rate(learning_rate),
    depth = sb_validate_catboost_depth(depth),
    l2_leaf_reg = sb_validate_catboost_l2_leaf_reg(l2_leaf_reg),
    random_seed = sb_check_positive_integer(
      seed,
      "catboost_seed",
      minimum = 1L
    ),
    loss_function = "Cox",
    eval_metric = "Cox",
    task_type = "CPU",
    thread_count = 1L,
    bootstrap_type = "Bayesian",
    bagging_temperature = 1,
    random_strength = 1,
    rsm = 1,
    grow_policy = "SymmetricTree",
    boosting_type = "Plain",
    use_best_model = FALSE,
    allow_writing_files = FALSE,
    save_snapshot = FALSE,
    logging_level = "Silent"
  )
}

sb_fit_catboost_cox <- function(x, time, event, params) {
  x <- sb_validate_predictor_frame(x, "x")
  time <- sb_validate_endpoint_vector(time, "time", nrow(x))
  event <- sb_validate_endpoint_vector(event, "event", nrow(x))
  if (any(time <= 0)) {
    sb_abort(
      paste0(
        "The fixed CatBoost Cox learner requires strictly positive ",
        "survival times."
      ),
      "stormbreaker_catboost_cox_nonpositive_time"
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
      "A fixed CatBoost Cox analysis fold must contain at least one event.",
      "stormbreaker_catboost_cox_no_event"
    )
  }
  if (!sb_catboost_available()) {
    sb_abort(
      "Package `catboost` is required for the fixed CatBoost Cox learner.",
      "stormbreaker_missing_catboost"
    )
  }

  training_pool <- catboost::catboost.load_pool(
    data = as.matrix(x),
    label = ifelse(event == 1L, time, -time)
  )
  booster <- catboost::catboost.train(
    learn_pool = training_pool,
    test_pool = NULL,
    params = params
  )
  tree_count <- as.integer(booster$tree_count)
  if (!identical(tree_count, params$iterations)) {
    sb_abort(
      "Fixed CatBoost Cox did not fit the declared number of iterations.",
      "stormbreaker_incomplete_catboost_fit"
    )
  }

  structure(
    list(
      features = names(x),
      params = params,
      iterations = params$iterations,
      booster = booster
    ),
    class = "stormbreaker_catboost_cox"
  )
}

sb_predict_catboost_cox <- function(model, x) {
  if (!inherits(model, "stormbreaker_catboost_cox")) {
    sb_abort(
      "The supplied model is not a fitted fixed CatBoost Cox learner.",
      "stormbreaker_bad_catboost_fit"
    )
  }
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), model$features)) {
    sb_abort(
      "Prediction columns do not match the fixed CatBoost Cox schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  if (!sb_catboost_available()) {
    sb_abort(
      "Package `catboost` is required for fixed CatBoost Cox prediction.",
      "stormbreaker_missing_catboost"
    )
  }

  prediction_pool <- catboost::catboost.load_pool(as.matrix(x))
  risk <- as.numeric(catboost::catboost.predict(
    model$booster,
    prediction_pool,
    prediction_type = "RawFormulaVal",
    ntree_start = 0L,
    ntree_end = model$iterations,
    thread_count = 1L
  ))
  if (length(risk) != nrow(x) || anyNA(risk) || any(!is.finite(risk))) {
    sb_abort(
      "Fixed CatBoost Cox produced invalid risk predictions.",
      "stormbreaker_invalid_catboost_prediction"
    )
  }
  unname(risk)
}

sb_predict_catboost_cox_survival <- function(
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
  analysis_risk <- sb_predict_catboost_cox(model, x_analysis)
  assessment_risk <- sb_predict_catboost_cox(model, x_assessment)
  cumulative_hazard <- sb_breslow_baseline_hazard(
    time,
    as.integer(event),
    analysis_risk,
    prediction_horizon
  )
  sb_survival_prediction_result(
    prediction_horizon,
    cumulative_hazard,
    source = "catboost_cox_breslow_raw_margin",
    ties = "breslow",
    risk = assessment_risk
  )
}

sb_catboost_cox_learner_adapter <- function(
    iterations = 100L,
    learning_rate = 0.1,
    depth = 6L,
    l2_leaf_reg = 3,
    seed = 1L) {
  params <- sb_catboost_cox_params(
    iterations = iterations,
    learning_rate = learning_rate,
    depth = depth,
    l2_leaf_reg = l2_leaf_reg,
    seed = seed
  )
  list(
    id = "catboost_cox",
    params = params,
    fit = function(x, time, event) {
      sb_fit_catboost_cox(x, time, event, params)
    },
    predict = sb_predict_catboost_cox,
    predict_survival = sb_predict_catboost_cox_survival,
    extract_tree_model = sb_extract_catboost_tree_model,
    extract_state = function(model, context = list()) {
      sb_extract_catboost_cox_learner_state(model, context)
    },
    predict_state = sb_predict_catboost_cox_learner_state
  )
}
