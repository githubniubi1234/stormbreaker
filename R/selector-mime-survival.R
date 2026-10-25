sb_mime_rsf_selector_ids <- function() {
  c("rsf_minimal_depth", "rsf_then_stepcox")
}

sb_mime_tuned_selector_ids <- function() {
  c(
    "coxboost_nonzero",
    "lasso_nonzero",
    "enet_cox_nonzero",
    "coxboost_then_stepcox",
    "lasso_then_stepcox"
  )
}

sb_mime_composite_selector_ids <- function() {
  c(
    "rsf_then_stepcox",
    "coxboost_then_stepcox",
    "lasso_then_stepcox"
  )
}

sb_mime_survival_selector_ids <- function() {
  c(
    "rsf_minimal_depth",
    "coxboost_nonzero",
    "lasso_nonzero",
    "enet_cox_nonzero",
    sb_mime_composite_selector_ids()
  )
}

sb_validate_mime_survival_selector_inputs <- function(x, time, event,
                                                       label) {
  inputs <- sb_validate_classical_survival_inputs(x, time, event, label)
  if (any(inputs$time <= 0)) {
    sb_abort(
      paste0(label, " selection requires strictly positive survival times."),
      "stormbreaker_mime_selector_nonpositive_time"
    )
  }
  inputs
}

sb_validate_mime_selector_features <- function(selected_features,
                                                source_features,
                                                label) {
  valid <- is.character(selected_features) &&
    length(selected_features) > 0L &&
    !anyNA(selected_features) &&
    anyDuplicated(selected_features) == 0L &&
    all(selected_features %in% source_features)
  if (!valid) {
    sb_abort(
      paste0(
        label, " did not return a nonempty set of valid predictors (returned ",
        length(selected_features), " of ", length(source_features),
        " candidate predictors)."
      ),
      "stormbreaker_mime_selector_too_few_features"
    )
  }
  unname(selected_features)
}

sb_apply_mime_survival_selector <- function(fitted, x) {
  valid_fit <- inherits(fitted, "stormbreaker_selector_fit") &&
    fitted$id %in% c(
      "rsf_minimal_depth",
      "coxboost_nonzero",
      "lasso_nonzero",
      "enet_cox_nonzero",
      "rsf_then_stepcox",
      "coxboost_then_stepcox",
      "lasso_then_stepcox"
    ) &&
    is.character(fitted$source_features) &&
    is.character(fitted$selected_features) &&
    length(fitted$selected_features) > 0L &&
    anyDuplicated(fitted$selected_features) == 0L &&
    all(fitted$selected_features %in% fitted$source_features)
  if (!valid_fit) {
    sb_abort(
      "The fitted selector is not a valid MIME survival selector state.",
      "stormbreaker_bad_selector_fit"
    )
  }
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), fitted$source_features)) {
    sb_abort(
      paste0(
        "Predictor columns must exactly match the fitted selector schema: ",
        paste(fitted$source_features, collapse = ", "),
        "."
      ),
      "stormbreaker_selector_schema_mismatch"
    )
  }
  x[fitted$selected_features]
}

sb_rsf_parse_tree_in_place <- function(recursiveObject, ...) {
  # Keep the native traversal, arithmetic and diagnostics. Reference containers
  # avoid copying the accumulated arrays at every recursive node visit.
  state <- list2env(recursiveObject, parent = emptyenv())
  state$subtree <- list2env(state$subtree, parent = emptyenv())
  native <- get("rfsrcParseTree", asNamespace("randomForestSRC"))
  result <- native(state, ...)
  result$subtree <- as.list(result$subtree)
  as.list(result)
}

sb_rsf_minimal_depth_selection <- function(engine) {
  # Rebind only this invocation; never alter the dependency's namespace.
  native <- asNamespace("randomForestSRC")
  scope <- new.env(parent = native)
  select <- get("var.select.rfsrc", native)
  maximum <- randomForestSRC::max.subtree
  environment(select) <- environment(maximum) <- scope
  scope$max.subtree <- maximum
  scope$rfsrcParseTree <- sb_rsf_parse_tree_in_place
  select(object = engine, conservative = "high", verbose = FALSE)
}

sb_fit_rsf_minimal_depth <- function(x, time, event, nodesize, seed) {
  inputs <- sb_validate_mime_survival_selector_inputs(
    x,
    time,
    event,
    "RSF minimal-depth"
  )
  if (!requireNamespace("randomForestSRC", quietly = TRUE)) {
    sb_abort(
      "Package `randomForestSRC` is required for RSF feature selection.",
      "stormbreaker_missing_randomforestsrc"
    )
  }
  formula_data <- sb_additional_survival_formula_data(
    inputs$x,
    inputs$time,
    inputs$event,
    random_forest_src = TRUE
  )
  engine <- sb_classical_survival_stop_on_warning(
    withr::with_options(
      list(rf.cores = 1L, mc.cores = 1L),
      randomForestSRC::rfsrc(
        formula_data$formula,
        data = formula_data$data,
        ntree = 1000L,
        nodesize = nodesize,
        splitrule = "logrank",
        importance = TRUE,
        proximity = FALSE,
        forest = TRUE,
        seed = seed
      )
    ),
    "RSF minimal-depth selector fit",
    "stormbreaker_rsf_selector_warning"
  )
  selection <- sb_classical_survival_stop_on_warning(
    withr::with_options(
      list(rf.cores = 1L, mc.cores = 1L),
      sb_rsf_minimal_depth_selection(engine)
    ),
    "RSF minimal-depth variable selection",
    "stormbreaker_rsf_selector_warning"
  )
  selected_features <- sb_validate_mime_selector_features(
    selection$topvars,
    names(inputs$x),
    "RSF minimal-depth"
  )
  native_scores <- selection$varselect
  positions <- match(names(inputs$x), rownames(native_scores))
  if (anyNA(positions)) {
    sb_abort(
      "RSF minimal-depth selection returned an incomplete score table.",
      "stormbreaker_invalid_rsf_selector"
    )
  }
  scores <- data.frame(
    feature = names(inputs$x),
    minimal_depth = as.numeric(native_scores$depth[positions]),
    importance = as.numeric(native_scores$vimp[positions]),
    selected = names(inputs$x) %in% selected_features,
    stringsAsFactors = FALSE
  )

  structure(
    list(
      id = "rsf_minimal_depth",
      source_features = names(inputs$x),
      selected_features = selected_features,
      scores = scores,
      conservative = "high",
      num_trees = 1000L,
      nodesize = nodesize,
      splitrule = "logrank",
      seed = seed,
      model_size = as.integer(length(selected_features)),
      error_rate = as.numeric(selection$err.rate)
    ),
    class = "stormbreaker_selector_fit"
  )
}

sb_rsf_minimal_depth_selector_adapter <- function(nodesize = 5L,
                                                   seed = 1L) {
  nodesize <- sb_check_positive_integer(
    nodesize,
    "selector_rsf_nodesize",
    minimum = 1L
  )
  seed <- sb_check_positive_integer(seed, "selector_seed", minimum = 1L)
  list(
    id = "rsf_minimal_depth",
    params = list(nodesize = nodesize, seed = seed),
    needs_outcome = TRUE,
    fit = function(x, time, event) {
      sb_fit_rsf_minimal_depth(x, time, event, nodesize, seed)
    },
    apply = sb_apply_mime_survival_selector
  )
}

sb_fit_coxboost_nonzero <- function(x, time, event, inner_folds, seed, capture_history = FALSE) {
  params <- sb_coxboost_cv_params(inner_folds = inner_folds, seed = seed)
  native <- sb_fit_coxboost_cv_native(
    x,
    time,
    event,
    params,
    label = "CoxBoost selector"
  )
  inputs <- native$inputs
  coefficient <- native$native_coefficient
  selected_features <- sb_validate_mime_selector_features(
    names(inputs$x)[coefficient != 0],
    names(inputs$x),
    "CoxBoost"
  )
  scores <- data.frame(
    feature = names(inputs$x),
    coefficient = unname(coefficient),
    selected = coefficient != 0,
    stringsAsFactors = FALSE
  )
  structure(
    list(
      id = "coxboost_nonzero",
      source_features = names(inputs$x),
      selected_features = selected_features,
      scores = scores,
      penalty = native$penalty,
      start_penalty = params$start_penalty,
      optimal_step = native$optimal_step,
      maximum_steps = params$maximum_steps,
      inner_folds = inner_folds,
      fold_id = native$fold_id,
      seed = seed,
      tuning_type = params$tuning_type,
      coefficient_path = if (capture_history) sweep(as.matrix(native$engine$coefficients),
        2, native$engine$sdx, "/") else NULL,
      tuning_curve = native$curve
    ),
    class = "stormbreaker_selector_fit"
  )
}

sb_coxboost_nonzero_selector_adapter <- function(inner_folds = 10L,
                                                  seed = 1L, capture_history = FALSE) {
  inner_folds <- sb_check_positive_integer(
    inner_folds,
    "selector_inner_folds",
    minimum = 3L
  )
  seed <- sb_check_positive_integer(seed, "selector_seed", minimum = 1L)
  list(
    id = "coxboost_nonzero",
    params = c(list(inner_folds = inner_folds, seed = seed),
      if (capture_history) list(capture_history = TRUE)),
    needs_outcome = TRUE,
    fit = function(x, time, event) {
      sb_fit_coxboost_nonzero(x, time, event, inner_folds, seed, capture_history)
    },
    apply = sb_apply_mime_survival_selector
  )
}

sb_fit_lasso_nonzero <- function(x, time, event, inner_folds, seed, group_ids = NULL,
                                  alpha = 1, id = "lasso_nonzero", capture_history = FALSE) {
  inputs <- sb_validate_mime_survival_selector_inputs(
    x,
    time,
    event,
    "Lasso"
  )
  tuning <- sb_select_glmnet_cox_lambda(
    inputs$x,
    inputs$time,
    inputs$event,
    inner_folds = inner_folds,
    seed = seed,
    alpha = alpha,
    standardize = TRUE,
    selection_rule = "lambda.min",
    measure = "deviance",
    group_ids = group_ids,
    retain_analysis_fit = TRUE
  )
  learner <- sb_glmnet_cox_learner_adapter(
    alpha = alpha,
    lambda = tuning$lambda,
    standardize = TRUE
  )
  # cv.glmnet has already fitted this exact full-analysis lambda path.
  engine <- tuning$analysis_fit
  attr(engine, "stormbreaker_features") <- names(inputs$x)
  attr(engine, "stormbreaker_lambda") <- tuning$lambda
  scores <- learner$extract_coefficients(engine)
  selected_features <- sb_validate_mime_selector_features(
    scores$feature[scores$nonzero],
    names(inputs$x),
    "Lasso"
  )

  structure(
    list(
      id = id,
      source_features = names(inputs$x),
      selected_features = selected_features,
      scores = data.frame(
        feature = scores$feature,
        coefficient = scores$coefficient,
        selected = scores$nonzero,
        stringsAsFactors = FALSE
      ),
      lambda = tuning$lambda,
      alpha = alpha,
      standardize = TRUE,
      selection_rule = "lambda.min",
      measure = "deviance",
      inner_folds = inner_folds,
      fold_id = tuning$fold_id,
      seed = seed,
      tuning_curve = tuning$curve,
      coefficient_path = if (capture_history) list(lambda = engine$lambda,
        beta = engine$beta) else NULL,
      tuning_warnings = tuning$warnings
    ),
    class = "stormbreaker_selector_fit"
  )
}

sb_lasso_nonzero_selector_adapter <- function(inner_folds = 10L,
                                               seed = 1L, id = "lasso_nonzero", capture_history = FALSE) {
  inner_folds <- sb_check_positive_integer(
    inner_folds,
    "selector_inner_folds",
    minimum = 3L
  )
  seed <- sb_check_positive_integer(seed, "selector_seed", minimum = 1L)
  list(
    id = id,
    params = c(list(inner_folds = inner_folds, seed = seed),
      if (capture_history) list(capture_history = TRUE)),
    needs_outcome = TRUE,
    fit = function(x, time, event, group_ids = NULL) {
      sb_fit_lasso_nonzero(x, time, event, inner_folds, seed, group_ids,
        alpha = if (id == "enet_cox_nonzero") .5 else 1, id = id,
        capture_history = capture_history)
    },
    apply = sb_apply_mime_survival_selector
  )
}

sb_sequential_survival_selector_adapter <- function(first, second, id) {
  valid_adapter <- function(adapter) {
    is.list(adapter) && is.character(adapter$id) &&
      length(adapter$id) == 1L && is.list(adapter$params) &&
      isTRUE(adapter$needs_outcome) && is.function(adapter$fit) &&
      is.function(adapter$apply)
  }
  valid_id <- is.character(id) && length(id) == 1L && !is.na(id) &&
    id %in% c(
      "rsf_then_stepcox",
      "coxboost_then_stepcox",
      "lasso_then_stepcox"
    )
  if (!valid_adapter(first) || !valid_adapter(second) || !valid_id) {
    sb_abort(
      "Sequential survival selection requires two valid selector adapters.",
      "stormbreaker_bad_sequential_selector"
    )
  }

  list(
    id = id,
    params = list(
      first = list(id = first$id, params = first$params),
      second = list(id = second$id, params = second$params)
    ),
    needs_outcome = TRUE,
    fit = function(x, time, event) {
      inputs <- sb_validate_mime_survival_selector_inputs(
        x,
        time,
        event,
        "Sequential survival"
      )
      endpoint <- list(task_type = "survival", time = inputs$time,
        event = inputs$event)
      fitted_first <- sb_fit_runtime_selector(
        first,
        inputs$x,
        endpoint
      )
      selected_first <- first$apply(fitted_first, inputs$x)
      fitted_second <- sb_fit_runtime_selector(
        second,
        selected_first,
        endpoint
      )
      selected_features <- fitted_second$selected_features
      valid_selection <- is.character(selected_features) &&
        length(selected_features) > 0L && !anyNA(selected_features) &&
        anyDuplicated(selected_features) == 0L &&
        all(selected_features %in% fitted_first$selected_features)
      if (!valid_selection) {
        sb_abort(
          "Sequential survival selection returned invalid final features.",
          "stormbreaker_invalid_sequential_selector"
        )
      }
      scores <- data.frame(
        feature = names(inputs$x),
        first_stage_selected = names(inputs$x) %in%
          fitted_first$selected_features,
        selected = names(inputs$x) %in% selected_features,
        stringsAsFactors = FALSE
      )
      structure(
        list(
          id = id,
          source_features = names(inputs$x),
          selected_features = selected_features,
          scores = scores,
          stages = list(first = fitted_first, second = fitted_second)
        ),
        class = "stormbreaker_selector_fit"
      )
    },
    apply = sb_apply_mime_survival_selector
  )
}

sb_rsf_then_stepcox_selector_adapter <- function(nodesize = 5L,
                                                  seed = 1L,
                                                  direction) {
  adapter <- sb_sequential_survival_selector_adapter(
    sb_rsf_minimal_depth_selector_adapter(nodesize, seed),
    sb_stepcox_selector_adapter(direction),
    "rsf_then_stepcox"
  )
  adapter$params <- list(
    nodesize = sb_check_positive_integer(
      nodesize,
      "selector_rsf_nodesize",
      minimum = 1L
    ),
    seed = sb_check_positive_integer(
      seed,
      "selector_seed",
      minimum = 1L
    ),
    direction = sb_validate_stepcox_direction(direction)
  )
  adapter
}

sb_coxboost_then_stepcox_selector_adapter <- function(inner_folds = 10L,
                                                       seed = 1L,
                                                       direction) {
  adapter <- sb_sequential_survival_selector_adapter(
    sb_coxboost_nonzero_selector_adapter(inner_folds, seed),
    sb_stepcox_selector_adapter(direction),
    "coxboost_then_stepcox"
  )
  adapter$params <- list(
    inner_folds = sb_check_positive_integer(
      inner_folds,
      "selector_inner_folds",
      minimum = 3L
    ),
    seed = sb_check_positive_integer(
      seed,
      "selector_seed",
      minimum = 1L
    ),
    direction = sb_validate_stepcox_direction(direction)
  )
  adapter
}

sb_lasso_then_stepcox_selector_adapter <- function(inner_folds = 10L,
                                                    seed = 1L,
                                                    direction) {
  adapter <- sb_sequential_survival_selector_adapter(
    sb_lasso_nonzero_selector_adapter(inner_folds, seed),
    sb_stepcox_selector_adapter(direction),
    "lasso_then_stepcox"
  )
  adapter$params <- list(
    inner_folds = sb_check_positive_integer(
      inner_folds,
      "selector_inner_folds",
      minimum = 3L
    ),
    seed = sb_check_positive_integer(
      seed,
      "selector_seed",
      minimum = 1L
    ),
    direction = sb_validate_stepcox_direction(direction)
  )
  adapter
}
