sb_wsrf_binary_learner_ids <- function() {
  "wsrf_binary"
}

sb_wsrf_available <- function() {
  requireNamespace("wsrf", quietly = TRUE)
}

sb_require_wsrf_dependency <- function() {
  if (!sb_wsrf_available()) {
    sb_abort(
      paste0(
        "Package `wsrf` is required for the weighted-subspace binary forest."
      ),
      "stormbreaker_missing_wsrf"
    )
  }
  invisible(TRUE)
}

sb_wsrf_binary_params <- function(
  learner = sb_wsrf_binary_learner_ids()
) {
  if (!learner %in% sb_wsrf_binary_learner_ids()) {
    sb_abort(
      paste0("Unknown weighted-subspace forest learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  }
  list(
    num_trees = 500L,
    mtry_rule = "floor(log2(p)+1)",
    weighted_subspace = TRUE,
    parallel = FALSE,
    importance = FALSE,
    leaf_node_size = 2L,
    seed = 1L,
    probability = "native_vote_proportion"
  )
}

sb_wsrf_binary_parameter_semantics <- function(learner) {
  params <- sb_wsrf_binary_params(learner)
  paste0(
    "weighted-subspace C4.5 classification forest with ",
    params$num_trees, " trees, mtry=floor(log2(p)+1), leaf node size ",
    params$leaf_node_size, ", native vote probabilities, serial execution, ",
    "and seed ", params$seed, " fixed before resampling"
  )
}



sb_wsrf_internal_names <- function(p) {
  paste0(".stormbreaker_wsrf_x", seq_len(p))
}

sb_wsrf_internal_data <- function(x) {
  internal <- x
  names(internal) <- sb_wsrf_internal_names(ncol(x))
  uses_dummy <- ncol(x) == 1L
  if (uses_dummy) {
    internal$.stormbreaker_wsrf_zero <- 0
  }
  list(
    data = internal,
    internal_features = names(internal),
    uses_zero_variance_bridge = uses_dummy
  )
}

sb_validate_wsrf_engine <- function(engine, internal_features, num_trees) {
  required <- c(
    "meta", "targetData", "trees", "treeOOBErrorRates", "OOBSets",
    "OOBPredictLabels", "treeIgrImportance"
  )
  valid <- inherits(engine, "wsrf") && is.list(engine) &&
    all(required %in% names(engine)) && is.list(engine$trees) &&
    length(engine$trees) == num_trees &&
    is.list(engine$meta) &&
    is.character(internal_features) && length(internal_features) > 0L
  if (!isTRUE(valid)) {
    sb_abort(
      "The weighted-subspace forest returned an invalid native fit.",
      "stormbreaker_invalid_wsrf_fit"
    )
  }
  engine
}

sb_wsrf_compact_prediction_engine <- function(engine, internal_features,
                                               num_trees) {
  engine <- sb_validate_wsrf_engine(engine, internal_features, num_trees)
  engine["targetData"] <- list(NULL)
  engine$treeOOBErrorRates <- rep(0, num_trees)
  engine$OOBSets <- replicate(num_trees, integer(), simplify = FALSE)
  engine$OOBPredictLabels <- replicate(
    num_trees,
    integer(),
    simplify = FALSE
  )
  engine$treeIgrImportance <- replicate(
    num_trees,
    numeric(),
    simplify = FALSE
  )
  removable <- c(
    "predicted", "oob.times", "confusion", "importance", "importanceSD",
    "RFOOBErrorRate", "strength", "correlation", "c_s2", "useweights",
    "mtry", "nodesize"
  )
  for (name in intersect(removable, names(engine))) {
    engine[name] <- list(NULL)
  }
  class(engine) <- "wsrf"
  sb_validate_wsrf_engine(engine, internal_features, num_trees)
}

sb_validate_wsrf_model <- function(model) {
  valid <- inherits(model, "stormbreaker_wsrf_binary_model") &&
    is.list(model) && identical(model$learner, sb_wsrf_binary_learner_ids()) &&
    is.character(model$features) && length(model$features) > 0L &&
    !anyNA(model$features) && !anyDuplicated(model$features) &&
    is.character(model$internal_features) &&
    length(model$internal_features) %in% c(
      length(model$features),
      length(model$features) + 1L
    ) &&
    is.logical(model$uses_zero_variance_bridge) &&
    length(model$uses_zero_variance_bridge) == 1L &&
    is.numeric(model$effective_mtry) && length(model$effective_mtry) == 1L &&
    !is.na(model$effective_mtry) && model$effective_mtry >= 1L &&
    is.list(model$params)
  if (!isTRUE(valid)) {
    sb_abort(
      "The supplied model is not a valid fitted weighted-subspace forest.",
      "stormbreaker_bad_wsrf_fit"
    )
  }
  sb_validate_wsrf_engine(
    model$engine,
    model$internal_features,
    model$params$num_trees
  )
  model
}

sb_fit_wsrf_binary <- function(x, outcome, learner, params) {
  sb_require_wsrf_dependency()
  inputs <- sb_validate_classical_binary_inputs(x, outcome, learner)
  prepared <- sb_wsrf_internal_data(inputs$x)
  effective_mtry <- as.integer(floor(log2(ncol(inputs$x)) + 1))
  grouping <- factor(
    inputs$outcome,
    levels = c(0L, 1L),
    labels = c("control", "case")
  )
  engine <- withr::with_seed(
    params$seed,
    wsrf::wsrf(
      x = prepared$data,
      y = grouping,
      mtry = effective_mtry,
      ntree = params$num_trees,
      weights = params$weighted_subspace,
      parallel = params$parallel,
      importance = params$importance,
      nodesize = params$leaf_node_size
    )
  )
  diagnostic_metrics <- list(importance = engine$importance,
    tree_oob = engine$treeOOBErrorRates, forest_oob = engine$RFOOBErrorRate)
  engine <- sb_wsrf_compact_prediction_engine(
    engine,
    prepared$internal_features,
    params$num_trees
  )
  model <- structure(
    list(
      learner = learner,
      features = names(inputs$x),
      internal_features = prepared$internal_features,
      uses_zero_variance_bridge = prepared$uses_zero_variance_bridge,
      effective_mtry = effective_mtry,
      diagnostic_metrics = diagnostic_metrics,
      params = params,
      engine = engine
    ),
    class = "stormbreaker_wsrf_binary_model"
  )
  sb_validate_wsrf_model(model)
}

sb_wsrf_prediction_data <- function(model, x) {
  model <- sb_validate_wsrf_model(model)
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), model$features)) {
    sb_abort(
      paste0(
        "Prediction columns do not match the fitted weighted-subspace forest ",
        "schema."
      ),
      "stormbreaker_learner_schema_mismatch"
    )
  }
  prepared <- sb_wsrf_internal_data(x)
  if (!identical(prepared$internal_features, model$internal_features) ||
    !identical(
      prepared$uses_zero_variance_bridge,
      model$uses_zero_variance_bridge
    )) {
    sb_abort(
      "The internal weighted-subspace forest predictor schema is invalid.",
      "stormbreaker_bad_wsrf_fit"
    )
  }
  list(model = model, x = x, internal = prepared$data)
}

sb_predict_wsrf_binary <- function(model, x) {
  prepared <- sb_wsrf_prediction_data(model, x)
  prediction <- stats::predict(
    prepared$model$engine,
    newdata = prepared$internal,
    type = "prob"
  )
  probability_matrix <- prediction$prob
  valid <- is.matrix(probability_matrix) &&
    nrow(probability_matrix) == nrow(prepared$x) &&
    identical(colnames(probability_matrix), c("control", "case")) &&
    !anyNA(probability_matrix) && all(is.finite(probability_matrix)) &&
    all(probability_matrix >= 0 & probability_matrix <= 1) &&
    all(abs(rowSums(probability_matrix) - 1) <= 1e-10)
  if (!isTRUE(valid)) {
    sb_abort(
      "The weighted-subspace forest produced invalid class probabilities.",
      "stormbreaker_invalid_wsrf_prediction"
    )
  }
  unname(probability_matrix[, "case"])
}

sb_extract_wsrf_binary_state <- function(model, context = list()) {
  model <- sb_validate_wsrf_model(model)
  structure(
    list(
      learner = model$learner,
      features = model$features,
      internal_features = model$internal_features,
      uses_zero_variance_bridge = model$uses_zero_variance_bridge,
      effective_mtry = model$effective_mtry,
      params = model$params,
      positive_class = 1L,
      native_positive_label = "case",
      probability = "native_vote_proportion",
      engine_serialization = "R_xdr_v3_wsrf_compact_prediction_engine",
      engine_raw = serialize(
        model$engine,
        connection = NULL,
        ascii = FALSE,
        xdr = TRUE,
        version = 3L
      ),
      training_outcome_retained = FALSE,
      oob_patient_indices_retained = FALSE
    ),
    class = c(
      "stormbreaker_wsrf_binary_state",
      "stormbreaker_learner_state"
    )
  )
}

sb_restore_wsrf_binary_model <- function(state) {
  if (!is.raw(state$engine_raw) || length(state$engine_raw) == 0L) {
    sb_bad_learner_state(
      "The saved weighted-subspace forest prediction engine is invalid."
    )
  }
  engine <- tryCatch(
    unserialize(state$engine_raw),
    error = function(error) {
      sb_abort(
        paste0(
          "The saved weighted-subspace forest engine could not be loaded: ",
          conditionMessage(error)
        ),
        "stormbreaker_invalid_wsrf_state",
        parent = error
      )
    }
  )
  model <- structure(
    list(
      learner = state$learner,
      features = state$features,
      internal_features = state$internal_features,
      uses_zero_variance_bridge = state$uses_zero_variance_bridge,
      effective_mtry = state$effective_mtry,
      params = state$params,
      engine = engine
    ),
    class = "stormbreaker_wsrf_binary_model"
  )
  sb_validate_wsrf_model(model)
}

sb_predict_wsrf_binary_state <- function(state, x) {
  valid <- inherits(state, "stormbreaker_wsrf_binary_state") &&
    is.list(state) && identical(state$learner, sb_wsrf_binary_learner_ids()) &&
    identical(state$positive_class, 1L) &&
    identical(state$native_positive_label, "case") &&
    identical(state$probability, "native_vote_proportion") &&
    identical(state$training_outcome_retained, FALSE) &&
    identical(state$oob_patient_indices_retained, FALSE)
  if (!isTRUE(valid)) {
    sb_bad_learner_state(
      "The saved weighted-subspace forest state is invalid."
    )
  }
  x <- sb_validate_learner_state_x(state, x, state$learner)
  sb_require_wsrf_dependency()
  model <- sb_restore_wsrf_binary_model(state)
  list(probability = sb_predict_wsrf_binary(model, x))
}

sb_wsrf_binary_learner_adapter <- function(
  learner,
  params = sb_wsrf_binary_params(learner)
) {
  list(
    id = learner,
    params = params,
    fit = function(x, outcome) {
      sb_fit_wsrf_binary(x, outcome, learner, params)
    },
    predict = sb_predict_wsrf_binary,
    extract_state = sb_extract_wsrf_binary_state,
    predict_state = sb_predict_wsrf_binary_state
  )
}
