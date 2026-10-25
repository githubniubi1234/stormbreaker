sb_multiclass_multinomial_catalog <- function() {
  data.frame(method = c("logistic", "stepglm", "brglm2", "gam", "pre", "slope"),
    label = c("Multinomial logistic", "Stepwise multinomial logistic", "Bias-reduced multinomial",
      "Multinomial GAM", "Prediction rule ensemble", "Multinomial SLOPE"),
    supported = c(rep(TRUE, 5L), FALSE),
    engine = c("nnet", "nnet", "brglm2", "mgcv", "pre", "SLOPE"),
    feature_result = "model_features", selection_rule = c(rep(
      "Predictive model; no standalone gene-selection interface", 5L),
      "Experimental adapter; multiclass validation incomplete; not available"))
}

sb_fit_multiclass_multinomial <- function(fit, x, y) {
  method <- fit$method
  catalog <- sb_multiclass_multinomial_catalog()
  fit$package <- catalog$engine[match(method, catalog$method)]
  if (!requireNamespace(fit$package, quietly = TRUE)) sb_abort(
    paste("Install", fit$package, "for this method."), "stormbreaker_missing_dependency")
  fit$feature_codes <- paste0("feature", seq_along(fit$features))
  names(x) <- fit$feature_codes
  d <- x
  d$.outcome <- y
  formula <- stats::as.formula(".outcome ~ .", env = environment())
  params <- list()
  record_warning <- function(w) fit$warnings <<- c(fit$warnings, conditionMessage(w))
  fit$engine <- withCallingHandlers(withr::with_seed(fit$seed, {
    if (method %in% c("logistic", "stepglm")) {
      params <- list(maxit = 200L, decay = 0, reference_class = fit$classes[1L])
      native <- nnet::multinom(formula, data = d, trace = FALSE, maxit = params$maxit,
        MaxNWts = 1000000L, decay = 0, model = FALSE)
      if (method == "stepglm") {
        params$selection <- "both_direction_AIC"
        native <- MASS::stepAIC(native, direction = "both", trace = FALSE)
      }
      native
    } else if (method == "brglm2") {
      params <- list(type = "AS_mean", reference_class = fit$classes[1L])
      brglm2::brmultinom(formula, data = d, ref = 1L, type = "AS_mean")
    } else if (method == "gam") {
      params <- list(family = "multinom", basis = "cs", basis_dimension = 5L, method = "REML")
      d$.outcome <- as.integer(y) - 1L
      rhs <- paste0("s(", names(x), ",bs='cs',k=5)", collapse = " + ")
      formulas <- lapply(seq_len(nlevels(y) - 1L), function(k)
        stats::as.formula(paste(if (k == 1L) ".outcome" else "", "~", rhs),
          env = asNamespace("mgcv")))
      mgcv::gam(formulas, data = d, family = mgcv::multinom(K = nlevels(y) - 1L),
        method = "REML", control = mgcv::gam.control(nthreads = 1L))
    } else if (method == "pre") {
      params <- sb_rule_aft_params("pre_binary_rule_ensemble")
      params$family <- "multinomial"
      params$internal_folds <- length(unique(fit$fold_id))
      pre::pre(formula, data = d, family = "multinomial", type = params$ensemble_type,
        sampfrac = params$sample_fraction, maxdepth = params$maximum_depth,
        learnrate = params$learning_rate, mtry = max(1L, ceiling(sqrt(ncol(x)))),
        ntrees = params$num_trees, normalize = TRUE, standardize = FALSE,
        nfolds = params$internal_folds, foldid = fit$fold_id, tree.unbiased = TRUE,
        removecomplements = TRUE, removeduplicates = TRUE, verbose = FALSE)
    } else {
      params <- sb_slope_randomglm_binary_params("slope_binary")
      params$family <- "multinomial"
      # The native gamma=0 relaxed multinomial refit can silently produce
      # NaN intercepts. Use ordinary penalized SLOPE, not relaxed refitting.
      params$gamma <- 1
      params$internal_resampling <- "task_stratified_or_grouped_folds"
      params$inner_folds <- length(unique(fit$fold_id))
      params$probability_mapping <- "multinomial_inverse_link_last_class_reference"
      native <- sb_multiclass_slope_cv(as.matrix(x), y, fit$fold_id, params)
      params$optimum <- native$optima
      fit$tuning <- native$summary
      if (any(!is.finite(unlist(native$model$intercepts))) ||
          any(!is.finite(as.matrix(native$model$coefficients[[1L]])))) sb_abort(
        "SLOPE returned non-finite multinomial coefficients.", "stormbreaker_bad_multiclass_fit")
      native$model
    }
  }), warning = record_warning)
  fit$parameters <- params
  fit$parameters$seed <- fit$seed
  structure(fit, class = "stormbreaker_multiclass_fit")
}

sb_predict_multiclass_multinomial <- function(fit, x) {
  if (!requireNamespace(fit$package, quietly = TRUE)) sb_abort(
    paste("Install", fit$package, "to predict with this saved fit."), "stormbreaker_missing_dependency")
  names(x) <- fit$feature_codes
  if (fit$method == "slope") {
    # SLOPE stores K-1 logits against the last class. Its native response
    # method omits that class and uses sum(1 + exp(eta)) as denominator.
    eta <- stats::predict(fit$engine, as.matrix(x), type = "link", simplify = FALSE)
    if (!identical(dim(eta), c(nrow(x), length(fit$codes) - 1L, 1L))) sb_abort(
      "Expected one selected SLOPE model with K-1 class logits.", "stormbreaker_bad_multiclass_fit")
    eta <- cbind(matrix(eta, nrow(x), length(fit$codes) - 1L), 0)
    p <- exp(eta - apply(eta, 1L, max))
    p <- p / rowSums(p)
    colnames(p) <- fit$engine$class_names
    return(p[, fit$codes, drop = FALSE])
  }
  p <- switch(fit$method,
    logistic = , stepglm = stats::predict(fit$engine, x, type = "probs"),
    brglm2 = stats::predict(fit$engine, x, type = "probs"),
    gam = stats::predict(fit$engine, x, type = "response"),
    pre = stats::predict(fit$engine, x, type = "response", penalty.par.val = "lambda.1se"))
  if (fit$method == "gam") return(p)
  if (is.null(dim(p)) && nrow(x) == 1L) p <- matrix(p, nrow = 1L, dimnames = list(NULL, names(p)))
  p[, fit$codes, drop = FALSE]
}

sb_multiclass_slope_cv <- function(x, y, folds, params) {
  # The installed cvSLOPE multinomial path gives NaN fold scores even when
  # SLOPE fits on those training rows are finite. Evaluate native fits here.
  train <- function(rows, q, alpha = "path") SLOPE::SLOPE(
    x[rows, , drop = FALSE], y[rows], family = "multinomial", q = q,
    alpha = alpha, gamma = params$gamma, center = params$center, scale = params$scale,
    path_length = params$path_length, threads = 1L,
    tol_dev_change = 0, tol_dev_ratio = 1)
  summaries <- lapply(params$q, function(q) {
    path <- train(seq_len(nrow(x)), q)
    alpha <- path$alpha
    loss <- matrix(NA_real_, nrow(x), length(alpha))
    for (fold in sort(unique(folds))) {
      test <- which(folds == fold)
      model <- train(which(folds != fold), q, alpha)
      eta <- stats::predict(model, x[test, , drop = FALSE], type = "link", simplify = FALSE)
      if (!identical(dim(eta), c(length(test), nlevels(y) - 1L, length(alpha))) ||
          any(!is.finite(eta))) sb_abort(
        "SLOPE returned an incomplete or non-finite inner-CV path.", "stormbreaker_bad_multiclass_fit")
      for (j in seq_along(alpha)) {
        logits <- cbind(matrix(eta[, , j], length(test), nlevels(y) - 1L), 0)
        maximum <- apply(logits, 1L, max)
        normalizer <- maximum + log(rowSums(exp(logits - maximum)))
        loss[test, j] <- 2 * (normalizer - logits[cbind(seq_along(test), as.integer(y[test]))])
      }
    }
    data.frame(q = q, gamma = params$gamma, alpha = alpha, measure = "deviance",
      mean = colMeans(loss))
  })
  summary <- do.call(rbind, summaries)
  optimum <- summary[which.min(summary$mean), , drop = FALSE]
  list(model = train(seq_len(nrow(x)), optimum$q, optimum$alpha),
    optima = optimum, summary = summary)
}
