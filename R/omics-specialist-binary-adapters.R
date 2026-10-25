sb_omics_specialist_binary_learner_ids <- function() {
  c(
    "brglm2_logistic",
    "bayesglm_logistic_cauchy",
    "bst_hinge_binary",
    "hda_binary",
    "dwd_binary",
    "ktsp_binary"
  )
}

sb_omics_specialist_binary_params <- function(learner) {
  switch(learner,
    brglm2_logistic = list(
      type = "AS_mean",
      link = "logit"
    ),
    bayesglm_logistic_cauchy = list(
      link = "logit",
      prior_mean = 0,
      prior_scale = 2.5,
      prior_df = 1,
      intercept_prior_mean = 0,
      intercept_prior_scale = 10,
      intercept_prior_df = 1,
      scaled = TRUE,
      maximum_iterations = 100L
    ),
    bst_hinge_binary = list(
      family = "hinge",
      cost = 0.5,
      base_learner = "ls",
      iterations = 50L,
      learning_rate = 0.1,
      center = TRUE,
      fit_intercept = FALSE,
      probability_mapping = "binary_softmax_native_margin"
    ),
    hda_binary = list(
      reduced_dimension = 1L,
      classification_rule = "naive_bayes",
      regularization_lambda = NULL,
      regularization_gamma = NULL,
      initial_loadings = "identity",
      outer_iterations = 7L,
      inner_iterations = 10L
    ),
    dwd_binary = list(
      exponent = 1,
      penalty_rule = "DWDLargeR_penaltyParameter",
      remove_zero_norm_features = TRUE,
      scale_features = TRUE,
      tolerance = 1e-5,
      maximum_iterations = 2000L,
      probability_mapping = "binary_softmax_native_margin"
    ),
    ktsp_binary = list(
      maximum_pairs = 10L,
      feature_filter = "wilcoxon",
      feature_cap = 100L,
      disjoint_pairs = TRUE,
      handle_ties = FALSE,
      probability_mapping = "positive_pair_vote_fraction"
    ),
    sb_abort(
      paste0("Unknown omics-specialist binary learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  )
}

sb_omics_specialist_binary_dependency <- function(learner) {
  package <- switch(learner,
    brglm2_logistic = "brglm2",
    bayesglm_logistic_cauchy = "arm",
    bst_hinge_binary = "bst",
    hda_binary = "hda",
    dwd_binary = "DWDLargeR",
    ktsp_binary = "switchBox",
    sb_abort(
      paste0("Unknown omics-specialist binary learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  )
  list(
    package = package,
    packages = package,
    available = requireNamespace(package, quietly = TRUE)
  )
}

sb_require_omics_specialist_binary_dependency <- function(learner) {
  dependency <- sb_omics_specialist_binary_dependency(learner)
  if (!dependency$available) {
    sb_abort(
      paste0(
        "Package `", dependency$package, "` is required for learner `",
        learner, "`."
      ),
      "stormbreaker_missing_omics_specialist_binary_dependency"
    )
  }
  invisible(dependency)
}

sb_omics_specialist_binary_parameter_semantics <- function(learner) {
  switch(learner,
    brglm2_logistic = paste0(
      "mean bias-reduced logistic adjusted-score fit with logit link fixed ",
      "before resampling"
    ),
    bayesglm_logistic_cauchy = paste0(
      "Bayesian logistic posterior-mode fit with scaled Cauchy coefficient ",
      "priors and a Cauchy intercept prior fixed before resampling"
    ),
    bst_hinge_binary = paste0(
      "50-step componentwise linear hinge-loss boosting with shrinkage 0.1, ",
      "centering, and a fixed softmax mapping of the native margin"
    ),
    hda_binary = paste0(
      "unregularized one-dimensional heteroscedastic discriminant ",
      "transformation with identity initialization, nonsingular class ",
      "covariances, and a native Gaussian naive-Bayes rule"
    ),
    dwd_binary = paste0(
      "distance-weighted discrimination exponent 1 with native ",
      "analysis-fold penalty rule, scaling, and all-zero feature removal"
    ),
    ktsp_binary = paste0(
      "disjoint K-TSP with analysis-fold Wilcoxon filtering, at most ten ",
      "pairs, and native K selection"
    )
  )
}

sb_new_omics_specialist_binary_model <- function(
  learner, features, params, engine
) {
  structure(
    list(
      learner = learner,
      features = features,
      params = params,
      engine = engine
    ),
    class = "stormbreaker_omics_specialist_binary_model"
  )
}

sb_fit_brglm2_logistic <- function(x, outcome, params) {
  design <- cbind("(Intercept)" = 1, data.matrix(x))
  fit <- brglm2::brglmFit(
    x = design,
    y = outcome,
    family = stats::binomial(link = params$link),
    control = brglm2::brglmControl(type = params$type),
    intercept = TRUE,
    singular.ok = FALSE
  )
  coefficients <- unname(fit$coefficients)
  valid <- isTRUE(fit$converged) && identical(fit$type, params$type) &&
    identical(as.integer(fit$rank), as.integer(ncol(design))) &&
    length(coefficients) == ncol(design) &&
    !anyNA(coefficients) && all(is.finite(coefficients))
  if (!valid) {
    sb_abort(
      paste0(
        "Mean bias-reduced logistic regression did not return a converged, ",
        "full-rank finite fit."
      ),
      "stormbreaker_incomplete_brglm2_logistic_fit"
    )
  }
  list(
    intercept = coefficients[[1L]],
    coefficient = stats::setNames(coefficients[-1L], names(x)),
    iterations = as.integer(fit$iter),
    type = fit$type
  )
}

sb_fit_bayesglm_logistic_cauchy <- function(x, outcome, params) {
  design <- cbind("(Intercept)" = 1, data.matrix(x))
  fit <- suppressWarnings(arm::bayesglm.fit(
    x = design,
    y = outcome,
    family = stats::binomial(link = params$link),
    control = list(maxit = params$maximum_iterations),
    intercept = TRUE,
    prior.mean = params$prior_mean,
    prior.scale = params$prior_scale,
    prior.df = params$prior_df,
    prior.mean.for.intercept = params$intercept_prior_mean,
    prior.scale.for.intercept = params$intercept_prior_scale,
    prior.df.for.intercept = params$intercept_prior_df,
    scaled = params$scaled,
    Warning = FALSE
  ))
  coefficients <- unname(fit$coefficients)
  valid <- isTRUE(fit$converged) &&
    identical(as.integer(fit$rank), as.integer(ncol(design))) &&
    length(coefficients) == ncol(design) &&
    !anyNA(coefficients) && all(is.finite(coefficients))
  if (!valid) {
    sb_abort(
      "Bayesian Cauchy logistic regression did not return a finite fit.",
      "stormbreaker_incomplete_bayesglm_logistic_fit"
    )
  }
  list(
    intercept = coefficients[[1L]],
    coefficient = stats::setNames(coefficients[-1L], names(x)),
    iterations = as.integer(fit$iter)
  )
}

sb_fit_bst_hinge_binary <- function(x, outcome, params) {
  matrix_x <- data.matrix(x)
  labels <- ifelse(outcome == 1L, 1, -1)
  fit <- bst::bst(
    matrix_x,
    labels,
    cost = params$cost,
    family = params$family,
    learner = params$base_learner,
    ctrl = bst::bst_control(
      mstop = params$iterations,
      nu = params$learning_rate,
      center = params$center,
      intercept = params$fit_intercept,
      trace = FALSE
    )
  )
  coefficient <- numeric(ncol(matrix_x))
  for (iteration in seq_len(params$iterations)) {
    feature <- fit$ensemble[[iteration]]
    coefficient[[feature]] <- coefficient[[feature]] +
      params$learning_rate * fit$coef[[iteration]]
  }
  intercept <- fit$offset - sum(coefficient * fit$meanx)
  valid <- inherits(fit, "bst") && identical(fit$family, params$family) &&
    identical(fit$learner, params$base_learner) &&
    length(coefficient) == ncol(matrix_x) &&
    !anyNA(coefficient) && all(is.finite(coefficient)) &&
    length(intercept) == 1L && is.finite(intercept)
  if (!valid) {
    sb_abort(
      "Hinge-loss boosting did not return a finite linear margin.",
      "stormbreaker_incomplete_bst_hinge_fit"
    )
  }
  list(
    intercept = unname(intercept),
    coefficient = stats::setNames(coefficient, names(x)),
    selected_features = names(x)[coefficient != 0],
    final_training_risk = unname(fit$risk[[params$iterations]]),
    diagnostic_path = list(risk = fit$risk, selected = fit$ensemble,
      updates = fit$coef, meanx = fit$meanx, offset = fit$offset,
      learning_rate = params$learning_rate, family = fit$family,
      learner = fit$learner, intercept_updates = fit$int)
  )
}

sb_hda_binary_prediction_engine <- function(engine) {
  keep <- c("hda.loadings", "reduced.dimension", "naivebayes")
  lean <- engine[keep]
  class(lean) <- class(engine)
  lean
}

sb_predict_hda_binary <- function(engine, matrix_x) {
  prediction <- suppressWarnings(stats::predict(
    engine,
    newdata = matrix_x,
    task = "c"
  ))
  positive <- match("1", engine$naivebayes$levels)
  posterior <- prediction$posteriors
  if (is.na(positive) || !is.matrix(posterior) ||
    nrow(posterior) != nrow(matrix_x) || ncol(posterior) < positive) {
    sb_abort(
      "HDA did not return a class-one posterior.",
      "stormbreaker_invalid_hda_binary_prediction"
    )
  }
  unname(posterior[, positive])
}

sb_fit_hda_binary <- function(x, outcome, params) {
  engine <- NULL
  invisible(utils::capture.output(
    engine <- hda::hda(
      data.matrix(x),
      grouping = factor(outcome, levels = c(0L, 1L)),
      newdim = params$reduced_dimension,
      crule = TRUE,
      reg.lamb = params$regularization_lambda,
      reg.gamm = params$regularization_gamma,
      initial.loadings = NULL,
      noutit = params$outer_iterations,
      ninit = params$inner_iterations,
      verbose = FALSE
    ),
    type = "output"
  ))
  engine <- sb_hda_binary_prediction_engine(engine)
  valid <- inherits(engine, "hda") &&
    is.matrix(engine$hda.loadings) &&
    nrow(engine$hda.loadings) == ncol(x) &&
    identical(engine$reduced.dimension, params$reduced_dimension) &&
    inherits(engine$naivebayes, "naiveBayes") &&
    identical(engine$naivebayes$levels, c("0", "1"))
  if (!valid) {
    sb_abort(
      "HDA did not return a two-class discriminant prediction engine.",
      "stormbreaker_incomplete_hda_binary_fit"
    )
  }
  sb_predict_hda_binary(engine, data.matrix(x))
  engine
}

sb_dwd_native_feature_map <- function(matrix_x, params) {
  feature_norm <- sqrt(colSums(matrix_x * matrix_x))
  kept <- seq_len(ncol(matrix_x))
  if (isTRUE(params$remove_zero_norm_features) && any(feature_norm == 0)) {
    kept <- which(feature_norm > 0)
    feature_norm <- c(feature_norm[kept], 0)
  }

  multiplier <- rep(1, length(feature_norm))
  dimension <- length(feature_norm)
  n <- nrow(matrix_x)
  if (isTRUE(params$scale_features) && dimension > 0.5 * n &&
    max(feature_norm) > 2 * min(feature_norm)) {
    denominator <- if (dimension > 3 * n) {
      sqrt(feature_norm)
    } else {
      feature_norm
    }
    multiplier <- 1 / pmax(1, denominator)
  }
  list(
    kept = kept,
    multiplier = multiplier
  )
}

sb_fit_dwd_binary <- function(x, outcome, params) {
  matrix_x <- data.matrix(x)
  feature_map <- sb_dwd_native_feature_map(matrix_x, params)
  labels <- ifelse(outcome == 1L, 1, -1)
  sparse_x <- methods::as(t(matrix_x), "dgCMatrix")
  fit <- NULL
  penalty <- NULL
  withr::with_preserve_seed(
    invisible(utils::capture.output(
      {
        penalty <- DWDLargeR::penaltyParameter(
          sparse_x,
          labels,
          expon = params$exponent
        )
        fit <- DWDLargeR::genDWD(
          sparse_x,
          labels,
          C = penalty,
          expon = params$exponent,
          tol = params$tolerance,
          maxIter = params$maximum_iterations,
          printDetails = 0,
          rmzeroFea = as.integer(params$remove_zero_norm_features),
          scaleFea = as.integer(params$scale_features)
        )
      },
      type = "output"
    ))
  )
  native_coefficient <- if (is.list(fit)) as.numeric(fit$w) else numeric()
  intercept <- if (is.list(fit)) as.numeric(fit$beta) else numeric()
  valid <- is.numeric(penalty) && length(penalty) == 1L &&
    !is.na(penalty) && is.finite(penalty) && penalty > 0 &&
    length(intercept) == 1L && !is.na(intercept) && is.finite(intercept) &&
    length(native_coefficient) == length(feature_map$multiplier) &&
    !anyNA(native_coefficient) && all(is.finite(native_coefficient))
  if (!valid) {
    sb_abort(
      "Distance-weighted discrimination did not return a finite linear fit.",
      "stormbreaker_incomplete_dwd_binary_fit"
    )
  }
  coefficient <- numeric(ncol(matrix_x))
  if (length(feature_map$kept) > 0L) {
    kept_index <- seq_along(feature_map$kept)
    coefficient[feature_map$kept] <-
      native_coefficient[kept_index] * feature_map$multiplier[kept_index]
  }
  list(
    intercept = unname(intercept),
    coefficient = stats::setNames(coefficient, names(x)),
    penalty = unname(as.numeric(penalty)),
    iterations = if (is.list(fit$info)) as.integer(fit$info$iter) else NA_integer_,
    training_error = if (is.list(fit$info)) {
      as.numeric(fit$info$trainerr)
    } else {
      NA_real_
    }
  )
}

sb_fit_ktsp_binary <- function(x, outcome, params) {
  if (ncol(x) < 4L) {
    sb_abort(
      "K-TSP requires at least four features in the current analysis fold.",
      "stormbreaker_ktsp_too_few_features"
    )
  }
  maximum_pairs <- min(params$maximum_pairs, floor(ncol(x) / 2L))
  classifier <- switchBox::SWAP.Train.KTSP(
    inputMat = t(data.matrix(x)),
    phenoGroup = factor(outcome, levels = c(0L, 1L)),
    krange = seq_len(maximum_pairs),
    FilterFunc = switchBox::SWAP.Filter.Wilcoxon,
    handleTies = params$handle_ties,
    disjoint = params$disjoint_pairs,
    featureNo = min(params$feature_cap, ncol(x)),
    verbose = FALSE
  )
  pairs <- classifier$TSPs
  valid <- is.matrix(pairs) && ncol(pairs) == 2L && nrow(pairs) >= 1L &&
    nrow(pairs) <= maximum_pairs && all(pairs %in% names(x)) &&
    identical(as.character(classifier$labels), c("0", "1")) &&
    length(classifier$tieVote) == nrow(pairs) &&
    all(as.character(classifier$tieVote) == "both")
  if (!valid) {
    sb_abort(
      "K-TSP did not return an oriented set of class-one feature pairs.",
      "stormbreaker_incomplete_ktsp_binary_fit"
    )
  }
  list(
    pairs = unname(pairs),
    pair_scores = unname(as.numeric(classifier$score)),
    selected_pairs = as.integer(nrow(pairs)),
    labels = c("0", "1")
  )
}

sb_fit_omics_specialist_binary <- function(x, outcome, learner, params) {
  inputs <- sb_validate_classical_binary_inputs(x, outcome, learner)
  sb_require_omics_specialist_binary_dependency(learner)
  engine <- switch(learner,
    brglm2_logistic = sb_fit_brglm2_logistic(
      inputs$x,
      inputs$outcome,
      params
    ),
    bayesglm_logistic_cauchy =
      sb_fit_bayesglm_logistic_cauchy(
        inputs$x,
        inputs$outcome,
        params
      ),
    bst_hinge_binary = sb_fit_bst_hinge_binary(
      inputs$x,
      inputs$outcome,
      params
    ),
    hda_binary = sb_fit_hda_binary(
      inputs$x,
      inputs$outcome,
      params
    ),
    dwd_binary = sb_fit_dwd_binary(
      inputs$x,
      inputs$outcome,
      params
    ),
    ktsp_binary = sb_fit_ktsp_binary(
      inputs$x,
      inputs$outcome,
      params
    )
  )
  diagnostic_path <- NULL
  if (learner == "bst_hinge_binary") {
    diagnostic_path <- engine$diagnostic_path
    engine$diagnostic_path <- NULL
  }
  model <- sb_new_omics_specialist_binary_model(
    learner,
    names(inputs$x),
    params,
    engine
  )
  if (!is.null(diagnostic_path)) model$diagnostic_path <- diagnostic_path
  model
}

sb_predict_omics_specialist_binary <- function(model, x) {
  if (!inherits(model, "stormbreaker_omics_specialist_binary_model")) {
    sb_abort(
      "The supplied model is not a fitted omics-specialist binary learner.",
      "stormbreaker_bad_omics_specialist_binary_fit"
    )
  }
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), model$features)) {
    sb_abort(
      paste0(
        "Prediction columns do not match the omics-specialist binary ",
        "learner schema."
      ),
      "stormbreaker_learner_schema_mismatch"
    )
  }
  matrix_x <- data.matrix(x)
  probability <- switch(model$learner,
    brglm2_logistic = stats::plogis(
      model$engine$intercept +
        as.vector(matrix_x %*% model$engine$coefficient)
    ),
    bayesglm_logistic_cauchy = stats::plogis(
      model$engine$intercept +
        as.vector(matrix_x %*% model$engine$coefficient)
    ),
    bst_hinge_binary = {
      margin <- model$engine$intercept +
        as.vector(matrix_x %*% model$engine$coefficient)
      stats::plogis(2 * margin)
    },
    hda_binary = sb_predict_hda_binary(
      model$engine,
      matrix_x
    ),
    dwd_binary = {
      margin <- model$engine$intercept +
        as.vector(matrix_x %*% model$engine$coefficient)
      stats::plogis(2 * margin)
    },
    ktsp_binary = {
      first <- match(model$engine$pairs[, 1L], model$features)
      second <- match(model$engine$pairs[, 2L], model$features)
      first_value <- matrix_x[, first, drop = FALSE]
      second_value <- matrix_x[, second, drop = FALSE]
      rowMeans((first_value > second_value) +
        0.5 * (first_value == second_value))
    }
  )
  sb_class_one_probability(probability, nrow(x), model$learner)
}

sb_extract_omics_specialist_binary_coefficients <- function(model) {
  if (!model$learner %in% c(
    "brglm2_logistic",
    "bayesglm_logistic_cauchy",
    "bst_hinge_binary",
    "dwd_binary"
  )) {
    sb_abort(
      "This omics-specialist learner has no linear coefficient vector.",
      "stormbreaker_omics_specialist_coefficients_not_applicable"
    )
  }
  data.frame(
    term = c("(Intercept)", model$features),
    coefficient = unname(c(
      model$engine$intercept,
      model$engine$coefficient
    )),
    stringsAsFactors = FALSE
  )
}

sb_extract_omics_specialist_binary_state <- function(model) {
  if (!inherits(model, "stormbreaker_omics_specialist_binary_model")) {
    sb_bad_learner_state(
      "The fitted omics-specialist binary learner state is invalid."
    )
  }
  structure(
    list(
      learner = model$learner,
      features = model$features,
      params = model$params,
      engine = model$engine,
      training_data_retention = "none"
    ),
    class = c(
      "stormbreaker_omics_specialist_binary_state",
      "stormbreaker_learner_state"
    )
  )
}

sb_predict_omics_specialist_binary_state <- function(state, x) {
  valid <- inherits(state, "stormbreaker_omics_specialist_binary_state") &&
    is.list(state) &&
    state$learner %in% sb_omics_specialist_binary_learner_ids() &&
    is.list(state$params) && is.list(state$engine) &&
    identical(state$training_data_retention, "none")
  if (!isTRUE(valid)) {
    sb_bad_learner_state(
      "The saved omics-specialist binary learner state is invalid."
    )
  }
  x <- sb_validate_learner_state_x(state, x, state$learner)
  model <- sb_new_omics_specialist_binary_model(
    state$learner,
    state$features,
    state$params,
    state$engine
  )
  list(probability = sb_predict_omics_specialist_binary(model, x))
}

sb_omics_specialist_binary_learner_adapter <- function(
  learner,
  params = sb_omics_specialist_binary_params(learner)
) {
  adapter <- list(
    id = learner,
    params = params,
    fit = function(x, outcome) {
      sb_fit_omics_specialist_binary(x, outcome, learner, params)
    },
    predict = sb_predict_omics_specialist_binary,
    extract_state = function(model, context = list()) {
      sb_extract_omics_specialist_binary_state(model)
    },
    predict_state = sb_predict_omics_specialist_binary_state
  )
  if (learner %in% c(
    "brglm2_logistic",
    "bayesglm_logistic_cauchy",
    "bst_hinge_binary",
    "dwd_binary"
  )) {
    adapter$extract_coefficients <-
      sb_extract_omics_specialist_binary_coefficients
  }
  adapter
}
