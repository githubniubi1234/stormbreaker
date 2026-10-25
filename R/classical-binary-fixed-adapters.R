sb_lvq_binary_learner_ids <- function() {
  c(
    "lvq1_binary",
    "olvq1_binary",
    "lvq21_binary",
    "lvq3_binary"
  )
}

sb_liblinear_binary_learner_ids <- function() {
  c(
    "liblinear_l2_logistic_binary",
    "liblinear_l1_logistic_binary",
    "liblinear_l2_l2loss_svc_binary",
    "liblinear_l2_l1loss_svc_binary",
    "liblinear_l1_l2loss_svc_binary"
  )
}

sb_liblinear_binary_params <- function(learner) {
  switch(learner,
    liblinear_l2_logistic_binary = list(
      solver_type = 0L,
      cost = 1,
      tolerance = 0.01,
      bias = 1,
      standardize = TRUE,
      probability = TRUE
    ),
    liblinear_l1_logistic_binary = list(
      solver_type = 6L,
      penalty = "L1",
      loss = "logistic",
      optimizer = "primal",
      cost = 1,
      tolerance = 0.01,
      bias = 1,
      standardize = TRUE,
      probability = TRUE,
      score_mapping = "native_class_one_probability",
      seed = 1L
    ),
    liblinear_l2_l2loss_svc_binary = list(
      solver_type = 1L,
      penalty = "L2",
      loss = "L2-loss SVC",
      optimizer = "dual",
      cost = 1,
      tolerance = 0.1,
      bias = 1,
      standardize = TRUE,
      probability = FALSE,
      score_mapping = "binary_softmax_native_decision_values",
      seed = 1L
    ),
    liblinear_l2_l1loss_svc_binary = list(
      solver_type = 3L,
      penalty = "L2",
      loss = "L1-loss SVC",
      optimizer = "dual",
      cost = 1,
      tolerance = 0.1,
      bias = 1,
      standardize = TRUE,
      probability = FALSE,
      score_mapping = "binary_softmax_native_decision_values",
      seed = 1L
    ),
    liblinear_l1_l2loss_svc_binary = list(
      solver_type = 5L,
      penalty = "L1",
      loss = "L2-loss SVC",
      optimizer = "primal",
      cost = 1,
      tolerance = 0.01,
      bias = 1,
      standardize = TRUE,
      probability = FALSE,
      score_mapping = "binary_softmax_native_decision_values",
      seed = 1L
    ),
    sb_abort(
      paste0("Unknown fixed LIBLINEAR learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  )
}

sb_liblinear_binary_parameter_semantics <- function(learner) {
  switch(learner,
    liblinear_l2_logistic_binary = paste0(
      "LIBLINEAR L2-logistic primal solver, cost, tolerance, bias, and ",
      "scaling fixed before resampling"
    ),
    liblinear_l1_logistic_binary = paste0(
      "LIBLINEAR L1-logistic primal solver, cost, tolerance, bias, scaling, ",
      "and seed fixed before resampling"
    ),
    liblinear_l2_l2loss_svc_binary = paste0(
      "LIBLINEAR L2-regularized L2-loss SVC dual solver, cost, tolerance, ",
      "bias, scaling, seed, and fixed decision-score mapping"
    ),
    liblinear_l2_l1loss_svc_binary = paste0(
      "LIBLINEAR L2-regularized L1-loss SVC dual solver, cost, tolerance, ",
      "bias, scaling, seed, and fixed decision-score mapping"
    ),
    liblinear_l1_l2loss_svc_binary = paste0(
      "LIBLINEAR L1-regularized L2-loss SVC primal solver, cost, tolerance, ",
      "bias, scaling, seed, and fixed decision-score mapping"
    ),
    sb_abort(
      paste0("Unknown fixed LIBLINEAR learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  )
}

sb_lvq_binary_params <- function(learner) {
  method <- switch(learner,
    lvq1_binary = "lvq1",
    olvq1_binary = "olvq1",
    lvq21_binary = "lvq2.1",
    lvq3_binary = "lvq3",
    sb_abort(
      paste0("Unknown fixed LVQ learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  )
  params <- list(
    method = method,
    iterations_per_prototype = if (identical(method, "olvq1")) 40L else 100L,
    alpha = if (identical(method, "olvq1")) 0.3 else 0.03,
    codebook_size_rule = "max_two_native_default",
    prototype_allocation = "equal_prior_target",
    initialization_k = 5L,
    standardize = TRUE,
    seed = 1L,
    score_mapping =
      "logistic_nearest_prototype_squared_distance_margin"
  )
  if (method %in% c("lvq2.1", "lvq3")) {
    params$window <- 0.3
  }
  if (identical(method, "lvq3")) {
    params$epsilon <- 0.1
  }
  params
}

sb_classical_binary_learner_ids <- function() {
  c(
    "nearest_centroid_binary", "bernoulli_naive_bayes_binary",
    "lda_binary",
    "qda_binary",
    "kernel_naive_bayes_binary",
    "rda_binary",
    "knn_binary",
    "nnet_binary",
    "decision_tree_binary",
    "random_forest_binary",
    "svm_linear_binary",
    "svm_radial_binary",
    "svm_radial_weighted_binary",
    "svm_polynomial_binary",
    "svm_sigmoid_binary",
    "gbm_binary",
    "bart_binary",
    "c50_binary",
    "kernel_weighted_knn_binary",
    "gaussian_process_binary",
    "pamr_binary",
    "neuralnet_multilayer_binary",
    "ctree_binary",
    "cforest_binary",
    "bagging_binary",
    "gamboost_binary",
    "blackboost_binary",
    "loclda_binary",
    "aorsf_binary",
    "earth_binary",
    "fda_binary",
    "mda_binary",
    sb_liblinear_binary_learner_ids(),
    "evtree_binary",
    "sda_binary",
    "rotation_forest_binary",
    "rrf_binary",
    "pda_binary",
    "slda_binary",
    "tree_binary",
    "c50_rules_binary",
    "random_ferns_binary",
    "robust_logistic_mqle",
    "lssvm_radial_binary",
    sb_lvq_binary_learner_ids()
  )
}

sb_classical_binary_params <- function(learner) {
  if (learner %in% sb_lvq_binary_learner_ids()) {
    return(sb_lvq_binary_params(learner))
  }
  if (learner %in% sb_liblinear_binary_learner_ids()) {
    return(sb_liblinear_binary_params(learner))
  }
  switch(learner,
    lda_binary = list(prior = "empirical"),
    nearest_centroid_binary = list(metric = "euclidean", prior = "uniform",
      standardize = FALSE, shrinkage = FALSE, score_mapping = "squared_distance_softmax"),
    bernoulli_naive_bayes_binary = list(binarize = 0, alpha = 1,
      prior = "empirical"),
    qda_binary = list(prior = "empirical"),
    kernel_naive_bayes_binary = list(
      prior = "empirical",
      usekernel = TRUE,
      fL = 0,
      adjust = 1
    ),
    rda_binary = list(
      prior = "empirical",
      gamma = 0.5,
      lambda = 0.5,
      estimate_error = FALSE
    ),
    knn_binary = list(
      k = 5L,
      standardize = TRUE,
      tie_rule = "distance_then_analysis_row"
    ),
    nnet_binary = list(
      size = 5L,
      decay = 0.01,
      maxit = 2000L,
      rang = 0.1,
      entropy = TRUE,
      standardize = TRUE,
      max_weights = 1000000L,
      seed = 1L
    ),
    decision_tree_binary = list(
      cp = 0.01,
      minsplit = 10L,
      maxdepth = 30L,
      xval = 0L
    ),
    random_forest_binary = list(
      num_trees = 500L,
      mtry = "floor_sqrt_p",
      min_node_size = 1L,
      seed = 1L
    ),
    svm_linear_binary = list(
      kernel = "linear",
      cost = 1,
      standardize = TRUE,
      probability = TRUE,
      seed = 1L
    ),
    svm_radial_binary = list(
      kernel = "radial",
      cost = 1,
      gamma = "one_over_p",
      standardize = TRUE,
      probability = TRUE,
      seed = 1L
    ),
    svm_radial_weighted_binary = list(
      kernel = "radial",
      cost = 1,
      gamma = "one_over_p",
      standardize = TRUE,
      probability = TRUE,
      class_weights = "inverse_frequency_normalized",
      seed = 1L
    ),
    svm_polynomial_binary = list(
      kernel = "polynomial",
      cost = 1,
      gamma = "one_over_p",
      degree = 3L,
      coef0 = 0,
      standardize = TRUE,
      probability = TRUE,
      seed = 1L
    ),
    svm_sigmoid_binary = list(
      kernel = "sigmoid",
      cost = 1,
      gamma = "one_over_p",
      coef0 = 0,
      standardize = TRUE,
      probability = TRUE,
      seed = 1L
    ),
    gbm_binary = list(
      n_trees = 100L,
      interaction_depth = 1L,
      min_obs_node = 10L,
      shrinkage = 0.1,
      bag_fraction = 1,
      seed = 1L
    ),
    bart_binary = list(
      num_trees = 50L,
      num_cutpoints = 100L,
      posterior_draws = 250L,
      burn_in = 100L,
      keep_every = 1L,
      prior_k = 2,
      tree_power = 2,
      tree_base = 0.95,
      sparse = FALSE,
      seed = 1L
    ),
    c50_binary = list(
      trials = 10L,
      rules = FALSE,
      winnow = FALSE,
      early_stopping = FALSE,
      seed = 1L
    ),
    kernel_weighted_knn_binary = list(
      k = 7L,
      distance = 2,
      kernel = "optimal",
      scale = TRUE
    ),
    gaussian_process_binary = list(
      kernel = "rbfdot",
      kernel_parameters = "automatic",
      scaled = TRUE,
      tolerance = 0.0005,
      cross_validation = 0L,
      keep_fitted = FALSE,
      seed = 1L
    ),
    pamr_binary = list(
      threshold = 1,
      scale_sd = TRUE,
      offset_percent = 50,
      remove_zero_thresholds = FALSE,
      prior = "empirical"
    ),
    neuralnet_multilayer_binary = list(
      hidden = c(3L, 2L),
      threshold = 0.01,
      stepmax = 500000L,
      algorithm = "rprop+",
      error_function = "sse",
      activation = "logistic",
      linear_output = FALSE,
      standardize = TRUE,
      seed = 1L
    ),
    ctree_binary = list(
      teststat = "quad",
      testtype = "Bonferroni",
      mincriterion = 0.95,
      minsplit = 20L,
      minbucket = 7L,
      maxdepth = Inf,
      saveinfo = FALSE
    ),
    cforest_binary = list(
      num_trees = 500L,
      mtry_rule = "ceiling_sqrt_p",
      replace = TRUE,
      sample_fraction = 1,
      teststat = "quad",
      testtype = "Univ",
      mincriterion = 0,
      minsplit = 20L,
      minbucket = 7L,
      maxdepth = Inf,
      saveinfo = FALSE,
      seed = 1L,
      cores = 1L
    ),
    bagging_binary = list(
      num_trees = 100L,
      cp = 0,
      minsplit = 20L,
      minbucket = 7L,
      maxdepth = 30L,
      xval = 0L,
      seed = 1L,
      coob = FALSE,
      keep_x = FALSE
    ),
    gamboost_binary = list(
      mstop = 100L,
      nu = 0.1,
      risk = "inbag",
      family = "Binomial_logit",
      baselearner = "bbs",
      dfbase = 4,
      seed = 1L
    ),
    blackboost_binary = list(
      mstop = 100L,
      nu = 0.1,
      risk = "inbag",
      family = "Binomial_logit",
      teststat = "quad",
      testtype = "Teststatistic",
      mincriterion = 0,
      minsplit = 20L,
      minbucket = 7L,
      maxdepth = 2L,
      saveinfo = FALSE,
      seed = 1L
    ),
    loclda_binary = list(
      k_max = 30L,
      k_rule = "dimension_aware",
      weighted_apriori = TRUE,
      prior = "local_empirical"
    ),
    aorsf_binary = list(
      num_trees = 500L,
      num_split = 5L,
      num_retry = 3L,
      mtry_rule = "ceiling_sqrt_p",
      sample_with_replacement = TRUE,
      sample_fraction = 0.632,
      control_method = "glm",
      control_scale_x = TRUE,
      control_max_iter = 1L,
      seed = 1L,
      threads = 1L,
      importance = "none",
      attach_data = FALSE
    ),
    earth_binary = list(
      degree = 1L,
      nk_rule = "max_21_2p_plus1",
      nprune = 15L,
      pmethod = "backward",
      glm_family = "binomial_logit",
      nfold = 0L,
      keep_xy = FALSE
    ),
    fda_binary = list(
      method = "polyreg",
      dimension = 1L,
      degree = 1L,
      monomial = FALSE,
      keep_fitted = FALSE
    ),
    mda_binary = list(
      subclasses_per_class = 2L,
      dimension = 3L,
      iterations = 10L,
      method = "polyreg",
      degree = 1L,
      keep_fitted = FALSE,
      seed = 1L
    ),
    evtree_binary = list(
      minbucket = 7L,
      minsplit = 20L,
      maxdepth = 5L,
      iterations = 1000L,
      population_size = 50L,
      complexity_weight = 1,
      seed = 1L
    ),
    sda_binary = list(
      covariance_shrinkage = "analytic",
      variance_shrinkage = "analytic",
      class_frequency_shrinkage = "analytic",
      diagonal = FALSE
    ),
    rotation_forest_binary = list(
      num_trees = 10L,
      subset_rule = "round_p_over_3_bounded_by_1_and_floor_p_over_2",
      cp = 0.01,
      minsplit = 10L,
      maxdepth = 10L,
      xval = 0L,
      seed = 1L
    ),
    rrf_binary = list(
      num_trees = 500L,
      mtry_rule = "floor_sqrt_p",
      min_node_size = 1L,
      replace = TRUE,
      regularization_mode = 1L,
      regularization_coefficient = 0.5,
      seed = 1L
    ),
    pda_binary = list(
      method = "gen.ridge",
      dimension = 1L,
      lambda = 1,
      prior = "empirical",
      keep_fitted = FALSE
    ),
    slda_binary = list(
      score_rule = "native_PC_q",
      q = NULL,
      prior = "empirical"
    ),
    tree_binary = list(
      split = "deviance",
      mincut = 5L,
      minsize = 10L,
      mindev = 0.01,
      keep_model = FALSE,
      keep_x = FALSE,
      keep_y = FALSE,
      keep_weights = FALSE
    ),
    c50_rules_binary = list(
      trials = 10L,
      rules = TRUE,
      winnow = FALSE,
      early_stopping = FALSE,
      seed = 1L
    ),
    random_ferns_binary = list(
      depth = 5L,
      num_ferns = 500L,
      importance = "none",
      save_forest = TRUE,
      threads = 1L,
      seed = 1L,
      score_mapping = "binary_softmax"
    ),
    robust_logistic_mqle = list(
      method = "Mqle",
      weights_on_x = "none",
      tuning_constant = 3.5,
      convergence_tolerance = 1e-4,
      convergence_test = "coef",
      max_iterations = 50L,
      family = "binomial_logit"
    ),
    lssvm_radial_binary = list(
      kernel = "rbfdot",
      kernel_parameters = "automatic",
      sigma_seed = 1L,
      tau = 0.01,
      scaled = TRUE,
      reduced = FALSE,
      cross_validation = 0L,
      keep_fitted = FALSE,
      score_mapping = "binary_softmax"
    ),
    sb_abort(
      paste0("Unknown fixed binary learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  )
}

sb_classical_binary_dependency <- function(learner) {
  package <- if (learner %in% sb_liblinear_binary_learner_ids()) {
    "LiblineaR"
  } else {
    switch(learner,
      lda_binary = "MASS",
      nearest_centroid_binary = NA_character_,
      bernoulli_naive_bayes_binary = NA_character_,
      qda_binary = "MASS",
      kernel_naive_bayes_binary = "klaR",
      rda_binary = "klaR",
      knn_binary = NA_character_,
      nnet_binary = "nnet",
      decision_tree_binary = "rpart",
      random_forest_binary = "randomForest",
      svm_linear_binary = "e1071",
      svm_radial_binary = "e1071",
      svm_radial_weighted_binary = "e1071",
      svm_polynomial_binary = "e1071",
      svm_sigmoid_binary = "e1071",
      gbm_binary = "gbm",
      bart_binary = "BART",
      c50_binary = "C50",
      kernel_weighted_knn_binary = "kknn",
      gaussian_process_binary = "kernlab",
      pamr_binary = "pamr",
      neuralnet_multilayer_binary = "neuralnet",
      ctree_binary = "partykit",
      cforest_binary = "partykit",
      bagging_binary = "ipred",
      gamboost_binary = "mboost",
      blackboost_binary = "mboost",
      loclda_binary = "klaR",
      aorsf_binary = "aorsf",
      earth_binary = "earth",
      fda_binary = "mda",
      mda_binary = "mda",
      evtree_binary = "evtree",
      sda_binary = "sda",
      rotation_forest_binary = "rotationForest",
      rrf_binary = "RRF",
      pda_binary = "mda",
      slda_binary = "ipred",
      tree_binary = "tree",
      c50_rules_binary = "C50",
      random_ferns_binary = "rFerns",
      robust_logistic_mqle = "robustbase",
      lssvm_radial_binary = "kernlab",
      lvq1_binary = "class",
      olvq1_binary = "class",
      lvq21_binary = "class",
      lvq3_binary = "class"
    )
  }
  list(
    package = package,
    available = is.na(package) || requireNamespace(package, quietly = TRUE)
  )
}

sb_require_classical_binary_dependency <- function(learner) {
  dependency <- sb_classical_binary_dependency(learner)
  if (!dependency$available) {
    sb_abort(
      paste0(
        "Package `", dependency$package, "` is required for learner `",
        learner, "`."
      ),
      "stormbreaker_missing_classical_binary_dependency"
    )
  }
  invisible(dependency)
}

sb_validate_classical_binary_inputs <- function(x, outcome, learner) {
  x <- sb_validate_predictor_frame(x, "x")
  outcome <- sb_validate_binary_outcome_vector(outcome, nrow(x))
  if (!all(c(0L, 1L) %in% outcome)) {
    sb_abort(
      paste0("A `", learner, "` analysis fold must contain both classes."),
      "stormbreaker_classical_binary_single_class"
    )
  }
  list(x = x, outcome = outcome)
}

sb_classical_binary_formula_data <- function(x, outcome) {
  analysis <- x
  outcome_name <- ".stormbreaker_outcome"
  while (outcome_name %in% names(analysis)) {
    outcome_name <- paste0(outcome_name, "_")
  }
  analysis[[outcome_name]] <- factor(outcome, levels = c(0L, 1L))
  rhs <- paste(vapply(names(x), sb_quote_name, ""), collapse = " + ")
  formula <- stats::as.formula(
    paste(sb_quote_name(outcome_name), "~", rhs),
    env = baseenv()
  )
  list(data = analysis, formula = formula)
}

sb_tree_binary_internal_names <- function(p) {
  paste0(".stormbreaker_tree_x", seq_len(p))
}

sb_tree_binary_formula_data <- function(x, outcome) {
  analysis <- x
  names(analysis) <- sb_tree_binary_internal_names(ncol(x))
  outcome_name <- ".stormbreaker_tree_outcome"
  analysis[[outcome_name]] <- factor(outcome, levels = c(0L, 1L))
  list(
    data = analysis,
    formula = stats::reformulate(
      names(analysis)[names(analysis) != outcome_name],
      response = outcome_name,
      env = baseenv()
    )
  )
}

sb_new_classical_binary_model <- function(learner, features, params,
                                                engine) {
  structure(
    list(
      learner = learner,
      features = features,
      params = params,
      engine = engine
    ),
    class = "stormbreaker_classical_binary_fixed"
  )
}

sb_fit_lda_binary <- function(matrix_x, grouping, priors) {
  withCallingHandlers(
    MASS::lda(
      x = matrix_x,
      grouping = grouping,
      prior = priors
    ),
    warning = function(warning) {
      if (identical(conditionMessage(warning), "variables are collinear")) {
        sb_abort(
          "LDA cannot fit collinear analysis-fold predictors.",
          "stormbreaker_lda_collinear_predictors",
          parent = warning
        )
      }
    }
  )
}

sb_bart_binary_prediction_engine <- function(engine) {
  keep <- c("treedraws", "rm.const", "binaryOffset", "offset")
  engine[setdiff(names(engine), keep)] <- NULL
  engine
}

sb_fit_bart_binary <- function(matrix_x, outcome, learner, params) {
  engine <- NULL
  invisible(utils::capture.output(
    engine <- withr::with_seed(
      params$seed,
      BART::pbart(
        x.train = matrix_x,
        y.train = outcome,
        sparse = params$sparse,
        k = params$prior_k,
        power = params$tree_power,
        base = params$tree_base,
        ntree = params$num_trees,
        numcut = params$num_cutpoints,
        ndpost = params$posterior_draws,
        nskip = params$burn_in,
        keepevery = params$keep_every,
        nkeeptrain = 0L,
        nkeeptest = 0L,
        nkeeptreedraws = params$posterior_draws,
        printevery = .Machine$integer.max
      )
    ),
    type = "output"
  ))
  if (!inherits(engine, "pbart") || is.null(engine$treedraws)) {
    sb_abort(
      "Fixed binary BART did not return posterior tree draws.",
      "stormbreaker_incomplete_bart_binary_fit"
    )
  }
  sb_bart_binary_prediction_engine(engine)
}

sb_fit_c50_binary <- function(x, grouping, learner, params) {
  C50::C5.0(
    x = x,
    y = grouping,
    trials = params$trials,
    rules = params$rules,
    control = C50::C5.0Control(
      winnow = params$winnow,
      seed = params$seed,
      earlyStopping = params$early_stopping
    )
  )
}

sb_fit_kernel_weighted_knn_binary <- function(
  x, grouping, learner, params
) {
  if (nrow(x) < params$k) {
    sb_abort(
      "The kernel-weighted KNN analysis fold has fewer rows than its k.",
      "stormbreaker_kernel_weighted_knn_too_few_rows"
    )
  }
  structure(
    list(
      training = x,
      outcome = grouping,
      k = params$k,
      distance = params$distance,
      kernel = params$kernel,
      scale = params$scale
    ),
    class = "stormbreaker_kernel_weighted_knn_binary"
  )
}

sb_fit_gaussian_process_binary <- function(
  matrix_x, grouping, learner, params
) {
  engine <- NULL
  invisible(utils::capture.output(
    engine <- withr::with_seed(
      params$seed,
      suppressMessages(kernlab::gausspr(
        x = matrix_x,
        y = grouping,
        scaled = params$scaled,
        type = "classification",
        kernel = params$kernel,
        kpar = params$kernel_parameters,
        tol = params$tolerance,
        cross = params$cross_validation,
        fit = params$keep_fitted
      ))
    ),
    type = "output"
  ))
  engine
}

sb_pamr_binary_prediction_engine <- function(engine) {
  keep <- c(
    "centroids", "centroid.overall", "sd", "threshold.scale", "se.scale",
    "sign.contrast", "prior", "norm.cent"
  )
  lean <- engine[intersect(keep, names(engine))]
  class(lean) <- class(engine)
  lean
}

sb_fit_pamr_binary <- function(x, grouping, priors, learner, params) {
  data <- list(
    x = t(as.matrix(x)),
    y = grouping,
    geneid = names(x),
    genenames = names(x)
  )
  engine <- NULL
  invisible(utils::capture.output(
    engine <- pamr::pamr.train(
      data,
      threshold = params$threshold,
      scale.sd = params$scale_sd,
      offset.percent = params$offset_percent,
      prior = priors,
      remove.zeros = params$remove_zero_thresholds,
      sign.contrast = "both"
    ),
    type = "output"
  ))
  if (!inherits(engine, "pamrtrained")) {
    sb_abort(
      "Fixed PAM did not return a nearest shrunken-centroid model.",
      "stormbreaker_incomplete_pamr_fit"
    )
  }
  sb_pamr_binary_prediction_engine(engine)
}

sb_fit_neuralnet_multilayer_binary <- function(
  x, outcome, matrix_x, learner, params
) {
  center <- colMeans(matrix_x)
  scale <- apply(matrix_x, 2L, stats::sd)
  if (any(!is.finite(scale)) || any(scale <= 0)) {
    sb_abort(
      "The fixed multilayer neural network requires non-constant predictors.",
      "stormbreaker_neuralnet_constant_predictor"
    )
  }
  standardized <- sweep(
    sweep(matrix_x, 2L, center, "-"),
    2L,
    scale,
    "/"
  )
  analysis <- as.data.frame(standardized)
  names(analysis) <- names(x)
  outcome_name <- ".stormbreaker_outcome"
  while (outcome_name %in% names(analysis)) {
    outcome_name <- paste0(outcome_name, "_")
  }
  analysis[[outcome_name]] <- outcome
  rhs <- paste(vapply(names(x), sb_quote_name, ""), collapse = " + ")
  formula <- stats::as.formula(
    paste(sb_quote_name(outcome_name), "~", rhs),
    env = baseenv()
  )
  engine <- withr::with_seed(
    params$seed,
    withCallingHandlers(
      neuralnet::neuralnet(
        formula,
        data = analysis,
        hidden = params$hidden,
        threshold = params$threshold,
        stepmax = params$stepmax,
        rep = 1L,
        algorithm = params$algorithm,
        err.fct = params$error_function,
        act.fct = params$activation,
        linear.output = params$linear_output,
        likelihood = identical(params$error_function, "ce"),
        lifesign = "none"
      ),
      warning = function(warning) {
        message <- conditionMessage(warning)
        if (grepl("Algorithm did not converge", message, fixed = TRUE)) {
          sb_abort(
            "The fixed multilayer neural-network learner did not converge.",
            "stormbreaker_neuralnet_nonconvergence",
            parent = warning
          )
        }
        sb_abort(
          paste0(
            "Learner `", learner, "` raised a fitting warning: ", message
          ),
          "stormbreaker_classical_binary_fit_warning",
          parent = warning
        )
      }
    )
  )
  reached <- if (is.matrix(engine$result.matrix) &&
    "reached.threshold" %in% rownames(engine$result.matrix)) {
    unname(engine$result.matrix["reached.threshold", 1L])
  } else {
    NA_real_
  }
  if (!inherits(engine, "nn") || length(engine$weights) != 1L ||
    !is.finite(reached) || reached > params$threshold) {
    sb_abort(
      "The fixed multilayer neural-network learner did not converge.",
      "stormbreaker_neuralnet_nonconvergence"
    )
  }
  lean <- list(
    model.list = engine$model.list,
    weights = engine$weights,
    act.fct = engine$act.fct,
    linear.output = engine$linear.output
  )
  class(lean) <- "nn"
  list(model = lean, center = center, scale = scale)
}

sb_fit_earth_binary <- function(x, outcome, learner, params) {
  engine <- earth::earth(
    x = x,
    y = outcome,
    pmethod = params$pmethod,
    keepxy = params$keep_xy,
    trace = 0,
    glm = list(family = stats::binomial(link = "logit")),
    degree = params$degree,
    nprune = params$nprune,
    nfold = params$nfold,
    nk = max(21L, 2L * ncol(x) + 1L)
  )
  valid <- inherits(engine, "earth") &&
    is.matrix(engine$glm.coefficients) &&
    ncol(engine$glm.coefficients) == 1L &&
    !anyNA(engine$glm.coefficients) &&
    all(is.finite(engine$glm.coefficients))
  if (!valid) {
    sb_abort(
      paste0("Learner `", learner, "` did not return a finite Earth GLM."),
      "stormbreaker_incomplete_earth_fit"
    )
  }
  engine
}

sb_fit_pda_binary <- function(formula_data, params) {
  withCallingHandlers(
    mda::fda(
      formula_data$formula,
      data = formula_data$data,
      dimension = params$dimension,
      method = mda::gen.ridge,
      keep.fitted = params$keep_fitted,
      lambda = params$lambda
    ),
    warning = function(warning) {
      if (identical(
        conditionMessage(warning),
        "degenerate problem; no discrimination"
      )) {
        sb_abort(
          "PDA found no discriminant dimension in the analysis fold.",
          "stormbreaker_pda_no_discrimination",
          parent = warning
        )
      }
    }
  )
}

sb_fit_liblinear_binary <- function(matrix_x, outcome, learner, params) {
  center <- colMeans(matrix_x)
  scale <- apply(matrix_x, 2L, stats::sd)
  if (any(!is.finite(scale)) || any(scale <= 0)) {
    sb_abort(
      paste0(
        "Learner `", learner,
        "` requires non-constant predictors in each analysis fold."
      ),
      "stormbreaker_liblinear_constant_predictor"
    )
  }
  standardized <- sweep(
    sweep(matrix_x, 2L, center, "-"),
    2L,
    scale,
    "/"
  )
  fit_engine <- function() {
    LiblineaR::LiblineaR(
      data = standardized,
      target = outcome,
      type = params$solver_type,
      cost = params$cost,
      epsilon = params$tolerance,
      bias = params$bias,
      cross = 0L,
      verbose = FALSE,
      findC = FALSE
    )
  }
  engine <- if ("seed" %in% names(params)) {
    withr::with_seed(params$seed, fit_engine())
  } else {
    fit_engine()
  }
  if (!inherits(engine, "LiblineaR") ||
    !identical(sort(as.integer(engine$ClassNames)), c(0L, 1L))) {
    sb_abort(
      paste0("Learner `", learner, "` did not return a binary LiblineaR fit."),
      "stormbreaker_incomplete_liblinear_fit"
    )
  }
  list(model = engine, center = center, scale = scale)
}

sb_fit_rotation_forest_binary <- function(x, grouping, learner, params) {
  predictor_count <- ncol(x)
  if (predictor_count < 2L) {
    sb_abort(
      paste0("Learner `", learner, "` requires at least two predictors."),
      "stormbreaker_rotation_forest_too_few_predictors"
    )
  }
  max_subsets <- floor(predictor_count / 2L)
  subsets <- max(
    1L,
    min(as.integer(round(predictor_count / 3L)), max_subsets)
  )
  engine <- withr::with_seed(
    params$seed,
    rotationForest::rotationForest(
      x = x,
      y = grouping,
      K = subsets,
      L = params$num_trees,
      control = rpart::rpart.control(
        cp = params$cp,
        minsplit = params$minsplit,
        maxdepth = params$maxdepth,
        xval = params$xval
      ),
      model = FALSE
    )
  )
  if (!inherits(engine, "rotationForest") ||
    length(engine$models) != params$num_trees) {
    sb_abort(
      paste0("Learner `", learner, "` returned an incomplete rotation forest."),
      "stormbreaker_incomplete_rotation_forest_fit"
    )
  }
  engine
}

sb_fit_random_ferns_binary <- function(x, grouping, learner, params) {
  engine <- withr::with_seed(
    params$seed,
    rFerns::rFerns(
      x = x,
      y = grouping,
      depth = params$depth,
      ferns = params$num_ferns,
      importance = params$importance,
      saveForest = params$save_forest,
      threads = params$threads
    )
  )
  valid <- inherits(engine, "rFerns") &&
    identical(engine$classLabels, c("0", "1")) &&
    is.list(engine$model) && length(engine$model) > 0L &&
    identical(unname(engine$parameters[["classes"]]), 2L) &&
    identical(unname(engine$parameters[["depth"]]), params$depth) &&
    identical(unname(engine$parameters[["ferns"]]), params$num_ferns)
  if (!valid) {
    sb_abort(
      "Fixed Random Ferns returned an incomplete native forest.",
      "stormbreaker_incomplete_random_ferns_fit"
    )
  }
  engine[c(
    "oobScores", "importance", "oobPreds", "oobConfusionMatrix",
    "oobErr", "timeTaken"
  )] <- NULL
  engine
}

sb_fit_robust_logistic_mqle <- function(x, grouping, learner, params) {
  internal_names <- paste0(".stormbreaker_robust_x", seq_len(ncol(x)))
  analysis <- x
  names(analysis) <- internal_names
  outcome_name <- ".stormbreaker_robust_outcome"
  analysis[[outcome_name]] <- grouping
  formula <- stats::reformulate(
    internal_names,
    response = outcome_name,
    env = baseenv()
  )
  native <- robustbase::glmrob(
    formula = formula,
    family = stats::binomial(link = "logit"),
    data = analysis,
    method = params$method,
    weights.on.x = params$weights_on_x,
    control = robustbase::glmrobMqle.control(
      acc = params$convergence_tolerance,
      test.acc = params$convergence_test,
      maxit = params$max_iterations,
      tcc = params$tuning_constant
    ),
    model = FALSE,
    x = FALSE
  )
  coefficient <- stats::coef(native)
  valid <- isTRUE(native$converged) &&
    length(coefficient) == ncol(x) + 1L &&
    identical(names(coefficient), c("(Intercept)", internal_names)) &&
    !anyNA(coefficient) && all(is.finite(coefficient)) &&
    is.numeric(native$iter) && length(native$iter) == 1L &&
    is.finite(native$iter) && native$iter <= params$max_iterations
  if (!valid) {
    sb_abort(
      paste0(
        "Fixed robust Mqle logistic regression did not return a converged ",
        "finite full-rank fit."
      ),
      "stormbreaker_incomplete_robust_logistic_fit"
    )
  }
  list(
    intercept = unname(coefficient[[1L]]),
    coefficient = stats::setNames(
      unname(coefficient[-1L]),
      names(x)
    ),
    iterations = as.integer(native$iter)
  )
}

sb_fit_lssvm_radial_binary <- function(matrix_x, grouping, learner,
                                             params) {
  native <- NULL
  invisible(utils::capture.output(
    native <- withr::with_seed(
      params$sigma_seed,
      kernlab::lssvm(
        x = matrix_x,
        y = grouping,
        scaled = params$scaled,
        kernel = params$kernel,
        kpar = params$kernel_parameters,
        type = "classification",
        tau = params$tau,
        reduced = params$reduced,
        cross = params$cross_validation,
        fit = params$keep_fitted
      )
    ),
    type = "output"
  ))
  scaling <- kernlab::scaling(native)
  scaled <- scaling$scaled
  scaled_training <- kernlab::xmatrix(native)
  alpha <- kernlab::alpha(native)
  intercept <- kernlab::b(native)
  labels <- kernlab::lev(native)
  sigma <- unname(kernlab::kpar(kernlab::kernelf(native))[["sigma"]])
  valid <- inherits(native, "lssvm") &&
    identical(kernlab::type(native), "classification") &&
    is.logical(scaled) && length(scaled) == ncol(matrix_x) &&
    is.matrix(scaled_training) &&
    identical(dim(scaled_training), dim(matrix_x)) &&
    is.matrix(alpha) &&
    identical(dim(alpha), c(nrow(matrix_x), 2L)) &&
    is.numeric(intercept) && length(intercept) == 2L &&
    is.numeric(sigma) && length(sigma) == 1L &&
    is.finite(sigma) && sigma > 0 &&
    identical(labels, c("0", "1")) &&
    !anyNA(scaled_training) && all(is.finite(scaled_training)) &&
    !anyNA(alpha) && all(is.finite(alpha)) &&
    !anyNA(intercept) && all(is.finite(intercept))
  center <- stats::setNames(rep(0, ncol(matrix_x)), colnames(matrix_x))
  scale <- stats::setNames(rep(1, ncol(matrix_x)), colnames(matrix_x))
  if (valid && any(scaled)) {
    native_center <- scaling$x.scale[["scaled:center"]]
    native_scale <- scaling$x.scale[["scaled:scale"]]
    valid <- is.numeric(native_center) && is.numeric(native_scale) &&
      length(native_center) == sum(scaled) &&
      length(native_scale) == sum(scaled) &&
      !anyNA(native_center) && all(is.finite(native_center)) &&
      !anyNA(native_scale) && all(is.finite(native_scale)) &&
      all(native_scale > 0)
    if (valid) {
      center[scaled] <- unname(native_center)
      scale[scaled] <- unname(native_scale)
    }
  }
  if (!valid) {
    sb_abort(
      "Fixed radial LS-SVM returned an incomplete native fit.",
      "stormbreaker_incomplete_lssvm_fit"
    )
  }
  list(
    center = center,
    scale = scale,
    scaled_training = unname(scaled_training),
    alpha = unname(alpha),
    intercept = unname(intercept),
    sigma = sigma,
    class_labels = labels
  )
}

sb_fit_lvq_binary <- function(matrix_x, grouping, learner, params) {
  center <- colMeans(matrix_x)
  scale <- apply(matrix_x, 2L, stats::sd)
  scale[!is.finite(scale) | scale <= 0] <- 1
  standardized <- sweep(
    sweep(matrix_x, 2L, center, "-"),
    2L,
    scale,
    "/"
  )
  native_size <- min(
    round(0.4 * 2 * (2 - 1 + ncol(standardized) / 2), 0),
    nrow(standardized)
  )
  requested_size <- max(2L, as.integer(native_size))

  fitted <- withr::with_seed(params$seed, {
    initial <- class::lvqinit(
      x = standardized,
      cl = grouping,
      size = requested_size,
      prior = c(0.5, 0.5),
      k = params$initialization_k
    )
    initial_classes <- as.character(initial$cl)
    valid_initial <- is.matrix(initial$x) && ncol(initial$x) == ncol(matrix_x) &&
      nrow(initial$x) >= 2L && !anyNA(initial$x) &&
      all(is.finite(initial$x)) &&
      identical(sort(unique(initial_classes)), c("0", "1"))
    if (!valid_initial) {
      sb_abort(
        paste0(
          "Learner `", learner,
          "` could not initialize at least one eligible prototype per class."
        ),
        "stormbreaker_incomplete_lvq_initialization"
      )
    }
    iterations <- params$iterations_per_prototype * nrow(initial$x)
    native <- switch(params$method,
      lvq1 = class::lvq1(
        standardized,
        grouping,
        initial,
        niter = iterations,
        alpha = params$alpha
      ),
      olvq1 = class::olvq1(
        standardized,
        grouping,
        initial,
        niter = iterations,
        alpha = params$alpha
      ),
      `lvq2.1` = class::lvq2(
        standardized,
        grouping,
        initial,
        niter = iterations,
        alpha = params$alpha,
        win = params$window
      ),
      lvq3 = class::lvq3(
        standardized,
        grouping,
        initial,
        niter = iterations,
        alpha = params$alpha,
        win = params$window,
        epsilon = params$epsilon
      )
    )
    list(native = native, iterations = as.integer(iterations))
  })

  prototypes <- as.matrix(fitted$native$x)
  prototype_class <- as.character(fitted$native$cl)
  valid <- ncol(prototypes) == ncol(matrix_x) && nrow(prototypes) >= 2L &&
    !anyNA(prototypes) && all(is.finite(prototypes)) &&
    identical(sort(unique(prototype_class)), c("0", "1"))
  if (!valid) {
    sb_abort(
      paste0("Learner `", learner, "` returned an incomplete LVQ codebook."),
      "stormbreaker_incomplete_lvq_fit"
    )
  }
  rownames(prototypes) <- NULL
  colnames(prototypes) <- NULL
  list(
    center = unname(center),
    scale = unname(scale),
    prototypes = unname(prototypes),
    prototype_class = unname(prototype_class),
    codebook_size = as.integer(nrow(prototypes)),
    iterations = fitted$iterations
  )
}

sb_fit_classical_binary <- function(x, outcome, learner, params) {
  inputs <- sb_validate_classical_binary_inputs(x, outcome, learner)
  sb_require_classical_binary_dependency(learner)
  matrix_x <- as.matrix(inputs$x)
  grouping <- factor(inputs$outcome, levels = c(0L, 1L))
  priors <- as.numeric(table(grouping) / length(grouping))
  formula_data <- if (identical(learner, "tree_binary")) {
    sb_tree_binary_formula_data(inputs$x, inputs$outcome)
  } else if (learner %in% c(
    "ctree_binary",
    "cforest_binary",
    "bagging_binary",
    "gamboost_binary",
    "blackboost_binary",
    "aorsf_binary",
    "fda_binary",
    "pda_binary",
    "mda_binary",
    "evtree_binary",
    "tree_binary"
  )) {
    sb_classical_binary_formula_data(inputs$x, inputs$outcome)
  } else {
    NULL
  }

  engine <- switch(learner,
    nearest_centroid_binary = list(centroids = rbind(
      colMeans(matrix_x[grouping == "0", , drop = FALSE]),
      colMeans(matrix_x[grouping == "1", , drop = FALSE]))),
    bernoulli_naive_bayes_binary = sb_fit_bernoulli_binary(
      matrix_x, inputs$outcome, params),
    lda_binary = sb_fit_lda_binary(matrix_x, grouping, priors),
    qda_binary = MASS::qda(
      x = matrix_x,
      grouping = grouping,
      prior = priors
    ),
    kernel_naive_bayes_binary = klaR::NaiveBayes(
      x = inputs$x,
      grouping = grouping,
      prior = priors,
      usekernel = params$usekernel,
      fL = params$fL,
      adjust = params$adjust
    ),
    rda_binary = klaR::rda(
      x = inputs$x,
      grouping = grouping,
      prior = priors,
      gamma = params$gamma,
      lambda = params$lambda,
      estimate.error = params$estimate_error,
      output = FALSE
    ),
    knn_binary = {
      if (nrow(matrix_x) < params$k) {
        sb_abort(
          "The fixed KNN analysis fold has fewer rows than its declared k.",
          "stormbreaker_knn_too_few_rows"
        )
      }
      center <- colMeans(matrix_x)
      scale <- apply(matrix_x, 2L, stats::sd)
      if (any(!is.finite(scale)) || any(scale <= 0)) {
        sb_abort(
          "Fixed KNN requires non-constant predictors in each analysis fold.",
          "stormbreaker_knn_constant_predictor"
        )
      }
      list(
        center = center,
        scale = scale,
        training = sweep(sweep(matrix_x, 2L, center, "-"), 2L, scale, "/"),
        outcome = inputs$outcome,
        k = params$k
      )
    },
    nnet_binary = {
      center <- colMeans(matrix_x)
      scale <- apply(matrix_x, 2L, stats::sd)
      if (any(!is.finite(scale)) || any(scale <= 0)) {
        sb_abort(
          paste0(
            "Learner `", learner,
            "` requires non-constant predictors in each analysis fold."
          ),
          "stormbreaker_nnet_constant_predictor"
        )
      }
      standardized <- sweep(
        sweep(matrix_x, 2L, center, "-"),
        2L,
        scale,
        "/"
      )
      network <- withr::with_seed(
        params$seed,
        nnet::nnet(
          x = standardized,
          y = inputs$outcome,
          size = params$size,
          decay = params$decay,
          maxit = params$maxit,
          rang = params$rang,
          entropy = params$entropy,
          trace = FALSE,
          MaxNWts = params$max_weights
        )
      )
      if (!identical(network$convergence, 0L)) {
        sb_abort(
          "The fixed neural-network learner did not converge.",
          "stormbreaker_nnet_nonconvergence"
        )
      }
      list(model = network, center = center, scale = scale)
    },
    decision_tree_binary = {
      analysis <- inputs$x
      outcome_name <- ".stormbreaker_outcome"
      while (outcome_name %in% names(analysis)) {
        outcome_name <- paste0(outcome_name, "_")
      }
      analysis[[outcome_name]] <- grouping
      rhs <- paste(vapply(names(inputs$x), sb_quote_name, ""), collapse = " + ")
      formula <- stats::as.formula(
        paste(sb_quote_name(outcome_name), "~", rhs),
        env = baseenv()
      )
      rpart::rpart(
        formula,
        data = analysis,
        method = "class",
        model = FALSE,
        x = FALSE,
        y = FALSE,
        control = rpart::rpart.control(
          cp = params$cp,
          minsplit = params$minsplit,
          maxdepth = params$maxdepth,
          xval = params$xval
        )
      )
    },
    random_forest_binary = withr::with_seed(
      params$seed,
      randomForest::randomForest(
        x = matrix_x,
        y = grouping,
        ntree = params$num_trees,
        mtry = max(1L, floor(sqrt(ncol(matrix_x)))),
        nodesize = params$min_node_size,
        importance = FALSE,
        keep.forest = TRUE
      )
    ),
    svm_linear_binary = withr::with_seed(
      params$seed,
      e1071::svm(
        x = matrix_x,
        y = grouping,
        type = "C-classification",
        kernel = params$kernel,
        cost = params$cost,
        scale = params$standardize,
        probability = params$probability,
        fitted = FALSE
      )
    ),
    svm_radial_binary = withr::with_seed(
      params$seed,
      e1071::svm(
        x = matrix_x,
        y = grouping,
        type = "C-classification",
        kernel = params$kernel,
        cost = params$cost,
        gamma = 1 / ncol(matrix_x),
        scale = params$standardize,
        probability = params$probability,
        fitted = FALSE
      )
    ),
    svm_radial_weighted_binary = withr::with_seed(
      params$seed,
      {
        counts <- table(grouping)
        class_weights <- stats::setNames(
          length(grouping) / (2 * as.numeric(counts)),
          names(counts)
        )
        e1071::svm(
          x = matrix_x,
          y = grouping,
          type = "C-classification",
          kernel = params$kernel,
          cost = params$cost,
          gamma = 1 / ncol(matrix_x),
          scale = params$standardize,
          probability = params$probability,
          class.weights = class_weights,
          fitted = FALSE
        )
      }
    ),
    svm_polynomial_binary = withr::with_seed(
      params$seed,
      e1071::svm(
        x = matrix_x,
        y = grouping,
        type = "C-classification",
        kernel = params$kernel,
        cost = params$cost,
        gamma = 1 / ncol(matrix_x),
        degree = params$degree,
        coef0 = params$coef0,
        scale = params$standardize,
        probability = params$probability,
        fitted = FALSE
      )
    ),
    svm_sigmoid_binary = withr::with_seed(
      params$seed,
      e1071::svm(
        x = matrix_x,
        y = grouping,
        type = "C-classification",
        kernel = params$kernel,
        cost = params$cost,
        gamma = 1 / ncol(matrix_x),
        coef0 = params$coef0,
        scale = params$standardize,
        probability = params$probability,
        fitted = FALSE
      )
    ),
    gbm_binary = withr::with_seed(
      params$seed,
      gbm::gbm.fit(
        x = matrix_x,
        y = inputs$outcome,
        distribution = "bernoulli",
        n.trees = params$n_trees,
        interaction.depth = params$interaction_depth,
        n.minobsinnode = params$min_obs_node,
        shrinkage = params$shrinkage,
        bag.fraction = params$bag_fraction,
        nTrain = nrow(matrix_x),
        keep.data = FALSE,
        verbose = FALSE,
        var.names = names(inputs$x),
        response.name = "outcome"
      )
    ),
    bart_binary = sb_fit_bart_binary(
      matrix_x,
      inputs$outcome,
      learner,
      params
    ),
    c50_binary = sb_fit_c50_binary(
      inputs$x,
      grouping,
      learner,
      params
    ),
    kernel_weighted_knn_binary =
      sb_fit_kernel_weighted_knn_binary(
        inputs$x,
        grouping,
        learner,
        params
      ),
    gaussian_process_binary = sb_fit_gaussian_process_binary(
      matrix_x,
      grouping,
      learner,
      params
    ),
    pamr_binary = sb_fit_pamr_binary(
      inputs$x,
      grouping,
      priors,
      learner,
      params
    ),
    neuralnet_multilayer_binary =
      sb_fit_neuralnet_multilayer_binary(
        inputs$x,
        inputs$outcome,
        matrix_x,
        learner,
        params
      ),
    ctree_binary = partykit::ctree(
      formula_data$formula,
      data = formula_data$data,
      control = partykit::ctree_control(
        teststat = params$teststat,
        testtype = params$testtype,
        mincriterion = params$mincriterion,
        minsplit = params$minsplit,
        minbucket = params$minbucket,
        maxdepth = params$maxdepth,
        saveinfo = params$saveinfo
      )
    ),
    cforest_binary = withr::with_seed(
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
          teststat = params$teststat,
          testtype = params$testtype,
          mincriterion = params$mincriterion,
          minsplit = params$minsplit,
          minbucket = params$minbucket,
          maxdepth = params$maxdepth,
          saveinfo = params$saveinfo
        ),
        cores = params$cores,
        trace = FALSE
      )
    ),
    bagging_binary = withr::with_seed(
      params$seed,
      ipred::bagging(
        formula_data$formula,
        data = formula_data$data,
        nbagg = params$num_trees,
        coob = params$coob,
        keepX = params$keep_x,
        control = rpart::rpart.control(
          cp = params$cp,
          minsplit = params$minsplit,
          minbucket = params$minbucket,
          maxdepth = params$maxdepth,
          xval = params$xval
        )
      )
    ),
    gamboost_binary = withr::with_seed(
      params$seed,
      mboost::gamboost(
        formula_data$formula,
        data = formula_data$data,
        family = mboost::Binomial(type = "glm", link = "logit"),
        control = mboost::boost_control(
          mstop = params$mstop,
          nu = params$nu,
          risk = params$risk
        ),
        baselearner = mboost::bbs,
        dfbase = params$dfbase
      )
    ),
    blackboost_binary = withr::with_seed(
      params$seed,
      mboost::blackboost(
        formula_data$formula,
        data = formula_data$data,
        family = mboost::Binomial(type = "glm", link = "logit"),
        control = mboost::boost_control(
          mstop = params$mstop,
          nu = params$nu,
          risk = params$risk
        ),
        tree_controls = partykit::ctree_control(
          teststat = params$teststat,
          testtype = params$testtype,
          mincriterion = params$mincriterion,
          minsplit = params$minsplit,
          minbucket = params$minbucket,
          maxdepth = params$maxdepth,
          saveinfo = params$saveinfo
        )
      )
    ),
    loclda_binary = klaR::loclda(
      x = matrix_x,
      grouping = grouping,
      k = sb_loclda_neighbors(params, ncol(matrix_x), nrow(matrix_x), nlevels(grouping)),
      weighted.apriori = params$weighted_apriori
    ),
    aorsf_binary = {
      tree_seeds <- withr::with_seed(
        params$seed,
        sample.int(.Machine$integer.max, params$num_trees)
      )
      aorsf::orsf(
        data = formula_data$data,
        formula = formula_data$formula,
        control = aorsf::orsf_control_classification(
          method = params$control_method,
          scale_x = params$control_scale_x,
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
    },
    earth_binary = sb_fit_earth_binary(
      inputs$x,
      inputs$outcome,
      learner,
      params
    ),
    fda_binary = mda::fda(
      formula_data$formula,
      data = formula_data$data,
      dimension = params$dimension,
      method = mda::polyreg,
      keep.fitted = params$keep_fitted,
      degree = params$degree,
      monomial = params$monomial
    ),
    pda_binary = sb_fit_pda_binary(formula_data, params),
    mda_binary = withr::with_seed(
      params$seed,
      mda::mda(
        formula_data$formula,
        data = formula_data$data,
        subclasses = params$subclasses_per_class,
        dimension = params$dimension,
        iter = params$iterations,
        method = mda::polyreg,
        keep.fitted = params$keep_fitted,
        trace = FALSE,
        degree = params$degree
      )
    ),
    liblinear_l2_logistic_binary = sb_fit_liblinear_binary(
      matrix_x,
      inputs$outcome,
      learner,
      params
    ),
    liblinear_l1_logistic_binary = sb_fit_liblinear_binary(
      matrix_x,
      inputs$outcome,
      learner,
      params
    ),
    liblinear_l2_l2loss_svc_binary = sb_fit_liblinear_binary(
      matrix_x,
      inputs$outcome,
      learner,
      params
    ),
    liblinear_l2_l1loss_svc_binary = sb_fit_liblinear_binary(
      matrix_x,
      inputs$outcome,
      learner,
      params
    ),
    liblinear_l1_l2loss_svc_binary = sb_fit_liblinear_binary(
      matrix_x,
      inputs$outcome,
      learner,
      params
    ),
    evtree_binary = withr::with_seed(
      params$seed,
      evtree::evtree(
        formula_data$formula,
        data = formula_data$data,
        control = evtree::evtree.control(
          minbucket = params$minbucket,
          minsplit = params$minsplit,
          maxdepth = params$maxdepth,
          niterations = params$iterations,
          ntrees = params$population_size,
          alpha = params$complexity_weight,
          seed = params$seed
        )
      )
    ),
    sda_binary = sda::sda(
      Xtrain = matrix_x,
      L = grouping,
      diagonal = params$diagonal,
      verbose = FALSE
    ),
    slda_binary = ipred::slda(
      grouping,
      X = matrix_x,
      q = params$q,
      prior = priors
    ),
    rotation_forest_binary = sb_fit_rotation_forest_binary(
      inputs$x,
      grouping,
      learner,
      params
    ),
    rrf_binary = withr::with_seed(
      params$seed,
      RRF::RRF(
        x = matrix_x,
        y = grouping,
        ntree = params$num_trees,
        mtry = max(1L, floor(sqrt(ncol(matrix_x)))),
        replace = params$replace,
        nodesize = params$min_node_size,
        importance = FALSE,
        localImp = FALSE,
        proximity = FALSE,
        keep.forest = TRUE,
        keep.inbag = FALSE,
        coefReg = params$regularization_coefficient,
        flagReg = params$regularization_mode
      )
    ),
    tree_binary = tree::tree(
      formula_data$formula,
      data = formula_data$data,
      control = tree::tree.control(
        nobs = nrow(inputs$x),
        mincut = params$mincut,
        minsize = params$minsize,
        mindev = params$mindev
      ),
      split = params$split,
      model = params$keep_model,
      x = params$keep_x,
      y = params$keep_y,
      wts = params$keep_weights
    ),
    c50_rules_binary = sb_fit_c50_binary(
      inputs$x,
      grouping,
      learner,
      params
    ),
    random_ferns_binary = sb_fit_random_ferns_binary(
      inputs$x,
      grouping,
      learner,
      params
    ),
    robust_logistic_mqle = sb_fit_robust_logistic_mqle(
      inputs$x,
      grouping,
      learner,
      params
    ),
    lssvm_radial_binary = sb_fit_lssvm_radial_binary(
      matrix_x,
      grouping,
      learner,
      params
    ),
    lvq1_binary = sb_fit_lvq_binary(
      matrix_x,
      grouping,
      learner,
      params
    ),
    olvq1_binary = sb_fit_lvq_binary(
      matrix_x,
      grouping,
      learner,
      params
    ),
    lvq21_binary = sb_fit_lvq_binary(
      matrix_x,
      grouping,
      learner,
      params
    ),
    lvq3_binary = sb_fit_lvq_binary(
      matrix_x,
      grouping,
      learner,
      params
    )
  )

  sb_new_classical_binary_model(
    learner,
    names(inputs$x),
    params,
    engine
  )
}

sb_class_one_probability <- function(probability, n, learner) {
  probability <- unname(as.numeric(probability))
  valid <- length(probability) == n && !anyNA(probability) &&
    all(is.finite(probability)) &&
    all(probability >= 0 & probability <= 1)
  if (!valid) {
    sb_abort(
      paste0("Learner `", learner, "` produced invalid probabilities."),
      "stormbreaker_invalid_classical_binary_prediction"
    )
  }
  probability
}

sb_predict_knn_binary <- function(engine, matrix_x) {
  standardized <- sweep(
    sweep(matrix_x, 2L, engine$center, "-"),
    2L,
    engine$scale,
    "/"
  )
  vapply(seq_len(nrow(standardized)), function(index) {
    difference <- sweep(engine$training, 2L, standardized[index, ], "-")
    distance <- rowSums(difference^2)
    neighbors <- order(distance, seq_along(distance))[seq_len(engine$k)]
    mean(engine$outcome[neighbors])
  }, numeric(1))
}

sb_predict_klar_rda_binary <- function(engine, x) {
  predict_method <- get(
    "predict.rda",
    envir = asNamespace("klaR"),
    inherits = FALSE
  )
  predict_method(engine, newdata = x)$posterior[, "1"]
}

sb_bart_prediction_matrix <- function(engine, matrix_x) {
  removed <- engine$rm.const
  if (length(removed) > 0L) {
    matrix_x <- matrix_x[, removed, drop = FALSE]
  }
  matrix_x
}

sb_predict_bart_binary <- function(engine, matrix_x) {
  prediction <- NULL
  predict_method <- get(
    "predict.pbart",
    envir = asNamespace("BART"),
    inherits = FALSE
  )
  invisible(utils::capture.output(
    prediction <- predict_method(
      engine,
      sb_bart_prediction_matrix(engine, matrix_x),
      mc.cores = 1L,
      openmp = FALSE
    ),
    type = "output"
  ))
  prediction$prob.test.mean
}

sb_predict_c50_binary <- function(engine, x) {
  predict_method <- get(
    "predict.C5.0",
    envir = asNamespace("C50"),
    inherits = FALSE
  )
  predict_method(engine, newdata = x, type = "prob")[, "1"]
}

sb_predict_tree_binary <- function(engine, x) {
  internal_x <- x
  names(internal_x) <- sb_tree_binary_internal_names(ncol(x))
  stats::predict(engine, newdata = internal_x, type = "vector")[, "1"]
}

sb_predict_kernel_weighted_knn_binary <- function(engine, x) {
  training <- engine$training
  outcome_name <- ".stormbreaker_outcome"
  while (outcome_name %in% names(training)) {
    outcome_name <- paste0(outcome_name, "_")
  }
  training[[outcome_name]] <- engine$outcome
  rhs <- paste(vapply(names(x), sb_quote_name, ""), collapse = " + ")
  formula <- stats::as.formula(
    paste(sb_quote_name(outcome_name), "~", rhs),
    env = baseenv()
  )
  prediction <- kknn::kknn(
    formula,
    train = training,
    test = x,
    k = engine$k,
    distance = engine$distance,
    kernel = engine$kernel,
    scale = engine$scale
  )
  prediction$prob[, "1"]
}

sb_predict_gaussian_process_binary <- function(engine, matrix_x) {
  probability <- kernlab::predict(
    engine,
    matrix_x,
    type = "probabilities"
  )
  probability[, "1"]
}

sb_predict_cforest_binary <- function(engine, x) {
  n_training <- nrow(engine$data)
  if (nrow(x) <= n_training) {
    return(stats::predict(engine, newdata = x, type = "prob")[, "1"])
  }

  # partykit 1.2-25's native rfweights reads past training weights when a
  # prediction batch is larger. Row-wise batches preserve its forest weights.
  starts <- seq.int(1L, nrow(x), by = n_training)
  unlist(lapply(starts, function(start) {
    rows <- seq.int(start, min(nrow(x), start + n_training - 1L))
    stats::predict(engine, newdata = x[rows, , drop = FALSE], type = "prob")[, "1"]
  }), use.names = FALSE)
}

sb_predict_pamr_binary <- function(engine, matrix_x, threshold) {
  posterior <- pamr::pamr.predict(
    engine,
    newx = t(matrix_x),
    threshold = threshold,
    type = "posterior"
  )
  posterior[, "1"]
}

sb_predict_neuralnet_multilayer_binary <- function(engine, x) {
  matrix_x <- as.matrix(x)
  standardized <- sweep(
    sweep(matrix_x, 2L, engine$center, "-"),
    2L,
    engine$scale,
    "/"
  )
  standardized <- as.data.frame(standardized)
  names(standardized) <- names(x)
  predict_method <- get(
    "predict.nn",
    envir = asNamespace("neuralnet"),
    inherits = FALSE
  )
  as.numeric(predict_method(engine$model, newdata = standardized))
}

sb_predict_loclda_binary <- function(engine, x) {
  predict_method <- get(
    "predict.loclda",
    envir = asNamespace("klaR"),
    inherits = FALSE
  )
  posterior <- predict_method(engine, newdata = x)$posterior
  if (is.null(colnames(posterior))) {
    # predict.loclda constructs columns in sorted observed-group order;
    # its uniform-probability branch can drop the column names.
    groups <- as.character(sort(unique(engine$grouping)))
    if (!is.matrix(posterior) || ncol(posterior) != length(groups) ||
        !setequal(groups, c("0", "1"))) {
      sb_abort(
        "locLDA did not retain identifiable binary posterior columns.",
        "stormbreaker_invalid_classical_binary_prediction"
      )
    }
    colnames(posterior) <- groups
  }
  posterior[, "1"]
}

sb_bagging_binary_prediction_engine <- function(engine) {
  engine$y <- factor(character(), levels = levels(engine$y))
  engine[c("X", "OOB", "call")] <- NULL
  for (index in seq_along(engine$mtrees)) {
    engine$mtrees[[index]]$bindx <- NULL
    engine$mtrees[[index]]$btree[c("y", "where", "call")] <- NULL
    environment(engine$mtrees[[index]]$btree$terms) <- baseenv()
  }
  engine
}

sb_predict_earth_binary <- function(engine, x) {
  model_matrix_method <- get(
    "model.matrix.earth",
    envir = asNamespace("earth"),
    inherits = FALSE
  )
  basis <- model_matrix_method(engine, x)
  stats::plogis(drop(basis %*% engine$glm.coefficients))
}

sb_predict_liblinear_binary <- function(engine, matrix_x) {
  standardized <- sweep(
    sweep(matrix_x, 2L, engine$center, "-"),
    2L,
    engine$scale,
    "/"
  )
  if (engine$model$Type %in% c(0L, 6L, 7L)) {
    prediction <- stats::predict(
      engine$model,
      standardized,
      proba = TRUE,
      decisionValues = FALSE
    )
    probability <- prediction$probabilities
    if (!is.matrix(probability) || !"1" %in% colnames(probability)) {
      return(numeric())
    }
    return(probability[, "1"])
  }

  prediction <- stats::predict(
    engine$model,
    standardized,
    proba = FALSE,
    decisionValues = TRUE
  )
  decision <- prediction$decisionValues
  if (!is.matrix(decision) ||
    !all(c("0", "1") %in% colnames(decision))) {
    return(numeric())
  }
  stats::plogis(decision[, "1"] - decision[, "0"])
}

sb_predict_sda_binary <- function(engine, matrix_x) {
  predict_method <- get(
    "predict.sda",
    envir = asNamespace("sda"),
    inherits = FALSE
  )
  predict_method(engine, matrix_x, verbose = FALSE)$posterior[, "1"]
}

sb_predict_slda_binary <- function(engine, matrix_x) {
  scores <- matrix_x %*% engine$scores
  colnames(scores) <- paste0("V", seq_len(ncol(scores)))
  stats::predict(
    engine$mylda,
    newdata = as.data.frame(scores)
  )$posterior[, "1"]
}

sb_predict_random_ferns_binary <- function(engine, x) {
  scores <- stats::predict(engine, x = x, scores = TRUE)
  scores <- as.matrix(scores)
  valid <- nrow(scores) == nrow(x) && ncol(scores) == 2L &&
    identical(colnames(scores), c("0", "1")) &&
    !anyNA(scores) && all(is.finite(scores))
  if (!valid) {
    sb_abort(
      "Fixed Random Ferns returned invalid native class scores.",
      "stormbreaker_invalid_random_ferns_scores"
    )
  }
  stats::plogis(scores[, "1"] - scores[, "0"])
}

sb_predict_robust_logistic_mqle <- function(engine, matrix_x) {
  linear_predictor <- unname(
    engine$intercept + as.vector(matrix_x %*% engine$coefficient)
  )
  stats::plogis(linear_predictor)
}

sb_predict_lssvm_radial_binary <- function(engine, matrix_x) {
  standardized <- sweep(
    sweep(matrix_x, 2L, engine$center, "-"),
    2L,
    engine$scale,
    "/"
  )
  squared_distance <- outer(
    rowSums(standardized^2),
    rowSums(engine$scaled_training^2),
    "+"
  ) - 2 * tcrossprod(standardized, engine$scaled_training)
  kernel <- exp(-engine$sigma * pmax(squared_distance, 0))
  scores <- kernel %*% engine$alpha
  scores <- sweep(scores, 2L, engine$intercept, "+")
  stats::plogis(scores[, 2L] - scores[, 1L])
}

sb_predict_lvq_binary <- function(engine, matrix_x) {
  standardized <- sweep(
    sweep(matrix_x, 2L, engine$center, "-"),
    2L,
    engine$scale,
    "/"
  )
  squared_distance <- outer(
    rowSums(standardized^2),
    rowSums(engine$prototypes^2),
    "+"
  ) - 2 * tcrossprod(standardized, engine$prototypes)
  squared_distance <- pmax(squared_distance, 0)
  distance_zero <- apply(
    squared_distance[, engine$prototype_class == "0", drop = FALSE],
    1L,
    min
  )
  distance_one <- apply(
    squared_distance[, engine$prototype_class == "1", drop = FALSE],
    1L,
    min
  )
  stats::plogis(distance_zero - distance_one)
}

sb_predict_classical_binary <- function(model, x) {
  if (!inherits(model, "stormbreaker_classical_binary_fixed")) {
    sb_abort(
      "The supplied model is not a fitted fixed binary learner.",
      "stormbreaker_bad_classical_binary_fit"
    )
  }
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), model$features)) {
    sb_abort(
      "Prediction columns do not match the fitted binary learner schema.",
      "stormbreaker_learner_schema_mismatch"
    )
  }
  learner <- model$learner
  sb_require_classical_binary_dependency(learner)
  matrix_x <- as.matrix(x)

  probability <- switch(learner,
    nearest_centroid_binary = sb_predict_centroid_binary(model$engine, matrix_x),
    bernoulli_naive_bayes_binary = sb_predict_bernoulli_binary(
      model$engine, matrix_x, model$params),
    lda_binary = stats::predict(model$engine, newdata = matrix_x)$posterior[, "1"],
    qda_binary = stats::predict(model$engine, newdata = matrix_x)$posterior[, "1"],
    kernel_naive_bayes_binary = stats::predict(
      model$engine,
      newdata = x,
      type = "raw"
    )$posterior[, "1"],
    rda_binary = sb_predict_klar_rda_binary(model$engine, x),
    knn_binary = sb_predict_knn_binary(model$engine, matrix_x),
    nnet_binary = {
      standardized <- sweep(
        sweep(matrix_x, 2L, model$engine$center, "-"),
        2L,
        model$engine$scale,
        "/"
      )
      as.numeric(stats::predict(
        model$engine$model,
        newdata = standardized,
        type = "raw"
      ))
    },
    decision_tree_binary = stats::predict(
      model$engine,
      newdata = x,
      type = "prob"
    )[, "1"],
    random_forest_binary = stats::predict(
      model$engine,
      newdata = matrix_x,
      type = "prob"
    )[, "1"],
    svm_linear_binary = {
      prediction <- stats::predict(
        model$engine,
        newdata = matrix_x,
        probability = TRUE
      )
      attr(prediction, "probabilities", exact = TRUE)[, "1"]
    },
    svm_radial_binary = {
      prediction <- stats::predict(
        model$engine,
        newdata = matrix_x,
        probability = TRUE
      )
      attr(prediction, "probabilities", exact = TRUE)[, "1"]
    },
    svm_radial_weighted_binary = {
      prediction <- stats::predict(
        model$engine,
        newdata = matrix_x,
        probability = TRUE
      )
      attr(prediction, "probabilities", exact = TRUE)[, "1"]
    },
    svm_polynomial_binary = {
      prediction <- stats::predict(
        model$engine,
        newdata = matrix_x,
        probability = TRUE
      )
      attr(prediction, "probabilities", exact = TRUE)[, "1"]
    },
    svm_sigmoid_binary = {
      prediction <- stats::predict(
        model$engine,
        newdata = matrix_x,
        probability = TRUE
      )
      attr(prediction, "probabilities", exact = TRUE)[, "1"]
    },
    gbm_binary = stats::predict(
      model$engine,
      newdata = matrix_x,
      n.trees = model$params$n_trees,
      type = "response"
    ),
    bart_binary = sb_predict_bart_binary(
      model$engine,
      matrix_x
    ),
    c50_binary = sb_predict_c50_binary(model$engine, x),
    kernel_weighted_knn_binary =
      sb_predict_kernel_weighted_knn_binary(model$engine, x),
    gaussian_process_binary =
      sb_predict_gaussian_process_binary(model$engine, matrix_x),
    pamr_binary = sb_predict_pamr_binary(
      model$engine,
      matrix_x,
      model$params$threshold
    ),
    neuralnet_multilayer_binary =
      sb_predict_neuralnet_multilayer_binary(model$engine, x),
    ctree_binary = stats::predict(
      model$engine,
      newdata = x,
      type = "prob"
    )[, "1"],
    cforest_binary = sb_predict_cforest_binary(model$engine, x),
    bagging_binary = stats::predict(
      model$engine,
      newdata = x,
      type = "prob",
      aggregation = "average"
    )[, "1"],
    gamboost_binary = stats::predict(
      model$engine,
      newdata = x,
      type = "response"
    ),
    blackboost_binary = stats::predict(
      model$engine,
      newdata = x,
      type = "response"
    ),
    loclda_binary = sb_predict_loclda_binary(model$engine, x),
    aorsf_binary = stats::predict(
      model$engine,
      new_data = x,
      pred_type = "prob",
      n_thread = model$params$threads
    )[, "1"],
    earth_binary = sb_predict_earth_binary(model$engine, x),
    fda_binary = sb_predict_fda_probability(model$engine, x),
    pda_binary = sb_predict_fda_probability(model$engine, x),
    mda_binary = stats::predict(
      model$engine,
      newdata = x,
      type = "posterior"
    )[, "1"],
    liblinear_l2_logistic_binary =
      sb_predict_liblinear_binary(model$engine, matrix_x),
    liblinear_l1_logistic_binary =
      sb_predict_liblinear_binary(model$engine, matrix_x),
    liblinear_l2_l2loss_svc_binary =
      sb_predict_liblinear_binary(model$engine, matrix_x),
    liblinear_l2_l1loss_svc_binary =
      sb_predict_liblinear_binary(model$engine, matrix_x),
    liblinear_l1_l2loss_svc_binary =
      sb_predict_liblinear_binary(model$engine, matrix_x),
    evtree_binary = stats::predict(
      model$engine,
      newdata = x,
      type = "prob"
    )[, "1"],
    sda_binary = sb_predict_sda_binary(model$engine, matrix_x),
    slda_binary = sb_predict_slda_binary(model$engine, matrix_x),
    rotation_forest_binary = stats::predict(
      model$engine,
      newdata = x
    ),
    rrf_binary = stats::predict(
      model$engine,
      newdata = matrix_x,
      type = "prob"
    )[, "1"],
    tree_binary = sb_predict_tree_binary(model$engine, x),
    c50_rules_binary = sb_predict_c50_binary(model$engine, x),
    random_ferns_binary =
      sb_predict_random_ferns_binary(model$engine, x),
    robust_logistic_mqle =
      sb_predict_robust_logistic_mqle(model$engine, matrix_x),
    lssvm_radial_binary =
      sb_predict_lssvm_radial_binary(model$engine, matrix_x),
    lvq1_binary = sb_predict_lvq_binary(model$engine, matrix_x),
    olvq1_binary = sb_predict_lvq_binary(model$engine, matrix_x),
    lvq21_binary = sb_predict_lvq_binary(model$engine, matrix_x),
    lvq3_binary = sb_predict_lvq_binary(model$engine, matrix_x)
  )
  sb_class_one_probability(probability, nrow(x), learner)
}

sb_earth_binary_prediction_engine <- function(engine) {
  drop_fields <- c(
    "rss", "rsq", "gcv", "grsq", "bx", "prune.terms",
    "fitted.values", "residuals", "rss.per.response", "rsq.per.response",
    "gcv.per.response", "grsq.per.response", "rss.per.subset",
    "gcv.per.subset", "leverages", "weights", "glm.list", "glm.stats",
    "call"
  )
  engine[intersect(drop_fields, names(engine))] <- NULL
  engine
}

sb_fda_binary_prediction_engine <- function(engine) {
  engine[intersect(c("call", "confusion"), names(engine))] <- NULL
  if (!is.null(engine$terms)) {
    environment(engine$terms) <- baseenv()
  }
  engine$fit[intersect(
    c("call", "residuals", "x"),
    names(engine$fit)
  )] <- NULL
  engine
}

sb_slda_binary_prediction_engine <- function(engine) {
  engine$mylda$call <- NULL
  engine
}

sb_mda_binary_prediction_engine <- function(engine) {
  engine[intersect(c("call", "confusion", "weights"), names(engine))] <- NULL
  if (!is.null(engine$terms)) {
    environment(engine$terms) <- baseenv()
  }
  engine
}

sb_evtree_binary_prediction_engine <- function(engine) {
  engine$data <- engine$data[0L, , drop = FALSE]
  environment(engine$terms) <- baseenv()
  engine$info$call <- NULL
  engine
}

sb_rotation_forest_binary_prediction_engine <- function(engine) {
  for (index in seq_along(engine$models)) {
    tree <- engine$models[[index]]
    tree[intersect(c("model", "x", "y", "where", "call"), names(tree))] <-
      NULL
    environment(tree$terms) <- baseenv()
    engine$models[[index]] <- tree
  }
  engine
}

sb_rrf_binary_prediction_engine <- function(engine) {
  drop_fields <- c(
    "call", "predicted", "err.rate", "confusion", "votes", "oob.times",
    "y", "test", "inbag", "localImportance", "proximity", "importanceSD"
  )
  engine[intersect(drop_fields, names(engine))] <- NULL
  engine
}

sb_tree_binary_prediction_engine <- function(engine) {
  engine[intersect(c("call", "where"), names(engine))] <- NULL
  environment(engine$terms) <- baseenv()
  engine
}

sb_classical_binary_prediction_engine <- function(learner, engine) {
  if (identical(learner, "bagging_binary")) {
    return(sb_bagging_binary_prediction_engine(engine))
  }
  if (identical(learner, "earth_binary")) {
    return(sb_earth_binary_prediction_engine(engine))
  }
  if (learner %in% c("fda_binary", "pda_binary")) {
    return(sb_fda_binary_prediction_engine(engine))
  }
  if (identical(learner, "slda_binary")) {
    return(sb_slda_binary_prediction_engine(engine))
  }
  if (identical(learner, "mda_binary")) {
    return(sb_mda_binary_prediction_engine(engine))
  }
  if (identical(learner, "evtree_binary")) {
    return(sb_evtree_binary_prediction_engine(engine))
  }
  if (identical(learner, "rotation_forest_binary")) {
    return(sb_rotation_forest_binary_prediction_engine(engine))
  }
  if (identical(learner, "rrf_binary")) {
    return(sb_rrf_binary_prediction_engine(engine))
  }
  if (identical(learner, "tree_binary")) {
    return(sb_tree_binary_prediction_engine(engine))
  }
  if (learner %in% c(
    "gaussian_process_binary",
    "aorsf_binary",
    "random_ferns_binary",
    "robust_logistic_mqle",
    "lssvm_radial_binary",
    sb_liblinear_binary_learner_ids(),
    sb_lvq_binary_learner_ids()
  )) {
    return(engine)
  }
  if (identical(learner, "kernel_naive_bayes_binary")) {
    engine$x <- engine$x[0L, , drop = FALSE]
  }
  if (identical(learner, "nnet_binary")) {
    engine$model[c("call", "fitted.values", "residuals")] <- NULL
  }
  drop_fields <- switch(learner,
    lda_binary = "call",
    qda_binary = "call",
    kernel_naive_bayes_binary = "call",
    rda_binary = c("call", "error.rate", "converged", "iter"),
    knn_binary = character(),
    nnet_binary = character(),
    decision_tree_binary = c("model", "x", "y", "where"),
    random_forest_binary = c(
      "call", "predicted", "err.rate", "confusion", "votes", "oob.times",
      "y", "test", "inbag", "localImportance", "proximity", "importanceSD"
    ),
    svm_linear_binary = c(
      "call", "index", "na.action", "xlevels", "fitted", "decision.values"
    ),
    svm_radial_binary = c(
      "call", "index", "na.action", "xlevels", "fitted", "decision.values"
    ),
    svm_radial_weighted_binary = c(
      "call", "index", "na.action", "xlevels", "fitted", "decision.values"
    ),
    svm_polynomial_binary = c(
      "call", "index", "na.action", "xlevels", "fitted", "decision.values"
    ),
    svm_sigmoid_binary = c(
      "call", "index", "na.action", "xlevels", "fitted", "decision.values"
    ),
    gbm_binary = c(
      "fit", "train.error", "valid.error", "oobag.improve"
    ),
    bart_binary = character(),
    c50_binary = c("call", "output", "boostResults", "size"),
    c50_rules_binary = c("call", "output", "boostResults", "size"),
    kernel_weighted_knn_binary = character(),
    pamr_binary = character(),
    neuralnet_multilayer_binary = character(),
    ctree_binary = character(),
    cforest_binary = character(),
    bagging_binary = character(),
    gamboost_binary = character(),
    blackboost_binary = character(),
    loclda_binary = character(),
    sda_binary = character()
  )
  engine[intersect(drop_fields, names(engine))] <- NULL
  engine
}

sb_extract_classical_binary_learner_state <- function(model) {
  if (!inherits(model, "stormbreaker_classical_binary_fixed")) {
    sb_bad_learner_state("The fitted fixed binary learner state is invalid.")
  }
  engine <- sb_classical_binary_prediction_engine(model$learner, model$engine)
  state <- list(
    learner = model$learner,
    features = model$features,
    params = model$params,
    engine_raw = serialize(
      engine,
      connection = NULL,
      ascii = FALSE,
      xdr = TRUE,
      version = 3L
    )
  )
  state$training_data_retention <- switch(model$learner,
    ctree_binary =
      "required_by_partykit_class_probability_state",
    cforest_binary =
      "required_by_partykit_class_probability_state",
    gamboost_binary =
      "required_by_mboost_prediction_state",
    blackboost_binary =
      "required_by_mboost_prediction_state",
    loclda_binary =
      "required_by_klar_localized_lda_prediction_state",
    evtree_binary =
      "required_by_evtree_terminal_class_probability_state",
    lssvm_radial_binary =
      "required_scaled_analysis_predictors_for_radial_kernel_prediction",
    NULL
  )
  structure(
    state,
    class = c(
      "stormbreaker_classical_binary_fixed_learner_state",
      "stormbreaker_learner_state"
    )
  )
}

sb_predict_classical_binary_learner_state <- function(state, x) {
  learner <- if (is.list(state)) state[["learner"]] else NULL
  if (!is.character(learner) || length(learner) != 1L || is.na(learner) ||
    !learner %in% sb_classical_binary_learner_ids()) {
    sb_bad_learner_state("The saved fixed binary learner state is invalid.")
  }
  x <- sb_validate_learner_state_x(state, x, learner)
  if (!is.raw(state[["engine_raw"]]) || length(state[["engine_raw"]]) == 0L) {
    sb_bad_learner_state("The saved fixed binary learner state is invalid.")
  }
  sb_require_classical_binary_dependency(learner)
  engine <- tryCatch(
    unserialize(state[["engine_raw"]]),
    error = function(error) {
      sb_abort(
        paste0("The saved `", learner, "` model could not be loaded."),
        "stormbreaker_invalid_classical_binary_state",
        parent = error
      )
    }
  )
  model <- sb_new_classical_binary_model(
    learner,
    state[["features"]],
    state[["params"]],
    engine
  )
  list(probability = sb_predict_classical_binary(model, x))
}

sb_classical_binary_learner_adapter <- function(
  learner,
  params = sb_classical_binary_params(learner)
) {
  if (!learner %in% sb_classical_binary_learner_ids()) {
    sb_abort(
      paste0("Unknown fixed binary learner `", learner, "`."),
      "stormbreaker_unsupported_learner"
    )
  }
  list(
    id = learner,
    params = params,
    fit = function(x, outcome) {
      sb_fit_classical_binary(x, outcome, learner, params)
    },
    predict = sb_predict_classical_binary,
    extract_state = function(model, context = list()) {
      sb_extract_classical_binary_learner_state(model)
    },
    predict_state = sb_predict_classical_binary_learner_state
  )
}
