sb_multiclass_specialist_catalog <- function() {
  data.frame(method = c("edda", "mclustda", "hdda", "hda", "rlt", "aorsf", "wsrf"),
    label = c("EDDA", "Mixture discriminant analysis (mclust)", "High-dimensional discriminant analysis",
      "Heteroscedastic discriminant analysis", "Reinforcement-learning forest",
      "Oblique random forest", "Weighted-subspace forest"), supported = TRUE,
    engine = c("mclust", "mclust", "HDclassif", "hda", "RLT", "aorsf", "wsrf"),
    feature_result = "model_features", selection_rule = "No standalone gene-selection rule")
}

sb_fit_multiclass_specialist <- function(fit, x, y) {
  method <- fit$method
  catalog <- sb_multiclass_specialist_catalog()
  fit$package <- catalog$engine[match(method, catalog$method)]
  if (!requireNamespace(fit$package, quietly = TRUE)) sb_abort(
    paste("Install", fit$package, "for this method."), "stormbreaker_missing_dependency")
  fit$feature_codes <- paste0("feature", seq_along(fit$features))
  names(x) <- fit$feature_codes
  matrix_x <- as.matrix(x)
  params <- list()
  record_warning <- function(w) fit$warnings <<- c(fit$warnings, conditionMessage(w))
  fit$engine <- withCallingHandlers(withr::with_seed(fit$seed, {
    if (method %in% c("edda", "mclustda")) {
      params <- sb_mclust_da_params(paste0(method, "_binary_bic"))
      params$active_covariance_models <- sb_mclust_da_active_models(ncol(x), params)
      if (method == "mclustda") params$initialization$seed <- fit$seed
      sb_mclust_da_fit_call(matrix_x, y, params, method)
    } else if (method == "hdda") {
      params <- sb_hdda_binary_params()
      params$probability_mapping <- "native_multiclass_posterior"
      HDclassif::hdda(matrix_x, cls = y, model = params$model, graph = FALSE,
        d_select = params$dimension_selection, threshold = params$threshold, show = FALSE,
        scaling = params$scaling, noise.ctrl = params$noise_control)
    } else if (method == "hda") {
      params <- list(reduced_dimension = min(nlevels(y) - 1L, ncol(x)),
        outer_iterations = 7L, inner_iterations = 10L, classification_rule = "naive_bayes")
      native <- NULL
      invisible(utils::capture.output(native <- hda::hda(matrix_x, grouping = y,
        newdim = params$reduced_dimension, crule = TRUE, reg.lamb = NULL, reg.gamm = NULL,
        initial.loadings = NULL, noutit = params$outer_iterations, ninit = params$inner_iterations,
        verbose = FALSE)))
      native
    } else if (method == "rlt") {
      params <- sb_rlt_learner_params("rlt_reinforcement_binary")
      RLT::RLT(x, y, model = "classification", print.summary = 0, use.cores = 1,
        ntrees = params$num_trees, nmin = params$native_nmin, reinforcement = TRUE,
        combsplit = 1L, embed.nmin = params$embedded_stop_size, importance = FALSE, track.obs = FALSE)
    } else if (method == "aorsf") {
      params <- sb_classical_binary_params("aorsf_binary")
      d <- x
      d$.outcome <- y
      aorsf::orsf(data = d, formula = stats::as.formula(".outcome ~ .", env = baseenv()),
        control = aorsf::orsf_control_classification(method = params$control_method,
          scale_x = params$control_scale_x, max_iter = params$control_max_iter),
        n_tree = params$num_trees, n_split = params$num_split, n_retry = params$num_retry,
        n_thread = 1L, mtry = max(1L, ceiling(sqrt(ncol(x)))),
        sample_with_replacement = params$sample_with_replacement, sample_fraction = params$sample_fraction,
        importance = "none", tree_seeds = sample.int(.Machine$integer.max, params$num_trees),
        attach_data = FALSE, na_action = "fail", verbose_progress = FALSE)
    } else {
      params <- sb_wsrf_binary_params()
      wsrf::wsrf(x, y, ntree = params$num_trees, mtry = min(ncol(x), floor(log2(ncol(x)) + 1)),
        weights = TRUE, parallel = FALSE, importance = FALSE, nodesize = params$leaf_node_size)
    }
  }), warning = record_warning)
  params$seed <- fit$seed
  fit$parameters <- params
  structure(fit, class = "stormbreaker_multiclass_fit")
}

sb_predict_multiclass_specialist <- function(fit, x) {
  if (!requireNamespace(fit$package, quietly = TRUE)) sb_abort(
    paste("Install", fit$package, "to predict with this saved fit."), "stormbreaker_missing_dependency")
  names(x) <- fit$feature_codes
  p <- switch(fit$method,
    edda = , mclustda = stats::predict(fit$engine, newdata = as.matrix(x))$z,
    hdda = stats::predict(fit$engine, data = as.matrix(x))$posterior,
    hda = stats::predict(fit$engine, newdata = as.matrix(x), task = "c")$posteriors,
    rlt = {
      native <- stats::predict(fit$engine, x)$ProbPrediction
      colnames(native) <- fit$engine$ylevels
      native
    },
    aorsf = stats::predict(fit$engine, new_data = x, pred_type = "prob", n_thread = 1L),
    wsrf = stats::predict(fit$engine, newdata = x, type = "prob")$prob)
  p[, fit$codes, drop = FALSE]
}
