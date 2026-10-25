sb_multiclass_discriminant_catalog <- function() {
  data.frame(method = c("fda", "mda", "pda", "slda", "sda", "loclda", "pamr",
    "gaussian_process", "kernel_weighted_knn", "liblinear_l2_logistic", "liblinear_l1_logistic"),
    label = c("Flexible discriminant analysis", "Mixture discriminant analysis",
      "Penalized discriminant analysis", "Stabilized LDA", "Shrinkage discriminant analysis",
      "Local LDA", "Nearest shrunken centroids", "Gaussian process", "Kernel-weighted KNN",
      "Liblinear L2 logistic", "Liblinear L1 logistic"), supported = TRUE,
    engine = c(rep("mda", 3), "ipred", "sda", "klaR", "pamr", "kernlab", "kknn", rep("LiblineaR", 2)),
    feature_result = "model_features", selection_rule = "No standalone gene-selection rule")
}

sb_fit_multiclass_discriminant <- function(fit, x, y) {
  method <- fit$method
  catalog <- sb_multiclass_discriminant_catalog()
  fit$package <- catalog$engine[match(method, catalog$method)]
  if (!requireNamespace(fit$package, quietly = TRUE)) sb_abort(
    paste("Install", fit$package, "for this method."), "stormbreaker_missing_dependency")
  params <- sb_classical_binary_params(paste0(method, "_binary"))
  params$seed <- fit$seed
  fit$feature_codes <- paste0("feature", seq_along(fit$features))
  names(x) <- fit$feature_codes
  matrix_x <- as.matrix(x)
  d <- x
  d$.outcome <- y
  formula <- stats::as.formula(".outcome ~ .", env = baseenv())
  priors <- as.numeric(table(y)) / length(y)
  if (method %in% c("fda", "pda")) params$dimension <- min(nlevels(y) - 1L, ncol(x))
  record_warning <- function(w) fit$warnings <<- c(fit$warnings, conditionMessage(w))
  fit$engine <- withCallingHandlers(withr::with_seed(fit$seed, {
    switch(method,
      fda = mda::fda(formula, data = d, dimension = params$dimension, method = mda::polyreg,
        degree = params$degree, monomial = params$monomial, keep.fitted = FALSE),
      pda = mda::fda(formula, data = d, dimension = params$dimension, method = mda::gen.ridge,
        lambda = params$lambda, keep.fitted = FALSE),
      mda = mda::mda(formula, data = d, subclasses = params$subclasses_per_class,
        dimension = params$dimension, iter = params$iterations, method = mda::polyreg,
        degree = params$degree, keep.fitted = FALSE, trace = FALSE),
      slda = ipred::slda(formula, data = d, q = params$q, prior = priors),
      sda = sda::sda(matrix_x, y, diagonal = FALSE, verbose = FALSE),
      loclda = klaR::loclda(matrix_x, y,
        k = sb_loclda_neighbors(params, ncol(x), nrow(x), nlevels(y)),
        weighted.apriori = params$weighted_apriori),
      pamr = {
        native <- NULL
        invisible(utils::capture.output(native <- pamr::pamr.train(
          list(x = t(matrix_x), y = y, geneid = fit$feature_codes, genenames = fit$feature_codes),
          threshold = params$threshold, scale.sd = params$scale_sd,
          offset.percent = params$offset_percent, prior = priors,
          remove.zeros = params$remove_zero_thresholds, sign.contrast = "both")))
        native
      },
      gaussian_process = kernlab::gausspr(matrix_x, y, type = "classification",
        kernel = params$kernel, kpar = "automatic", scaled = TRUE,
        tol = params$tolerance, cross = 0L, fit = FALSE),
      kernel_weighted_knn = list(training = d),
      {
        fit$center <- colMeans(matrix_x)
        fit$scale <- apply(matrix_x, 2L, stats::sd)
        fit$scale[fit$scale == 0] <- 1
        z <- sweep(sweep(matrix_x, 2L, fit$center, "-"), 2L, fit$scale, "/")
        params$score_mapping <- "native_multiclass_OVR_probabilities"
        LiblineaR::LiblineaR(z, y, type = params$solver_type, cost = params$cost,
          epsilon = params$tolerance, bias = params$bias, verbose = FALSE)
      })
  }), warning = record_warning)
  fit$parameters <- params
  structure(fit, class = "stormbreaker_multiclass_fit")
}

sb_predict_multiclass_discriminant <- function(fit, x) {
  if (!requireNamespace(fit$package, quietly = TRUE)) sb_abort(
    paste("Install", fit$package, "to predict with this saved fit."), "stormbreaker_missing_dependency")
  names(x) <- fit$feature_codes
  method <- fit$method
  p <- if (method %in% c("fda", "mda", "pda")) {
    # mda's quadratic polybasis drops dimensions for one query row.
    query <- if (method == "fda" && nrow(x) == 1L) x[c(1L, 1L), , drop = FALSE] else x
    native <- stats::predict(fit$engine, newdata = query, type = "posterior")
    native[seq_len(nrow(x)), , drop = FALSE]
  } else if (method == "sda") {
    stats::predict(fit$engine, as.matrix(x), verbose = FALSE)$posterior
  } else if (method == "slda") {
    stats::predict(fit$engine, x)$posterior
  } else if (method == "loclda") {
    native <- stats::predict(fit$engine, as.matrix(x))$posterior
    colnames(native) <- fit$engine$lev
    native
  } else if (method == "pamr") {
    pamr::pamr.predict(fit$engine, t(as.matrix(x)), threshold = fit$parameters$threshold,
      type = "posterior")
  } else if (method == "gaussian_process") {
    kernlab::predict(fit$engine, as.matrix(x), type = "probabilities")
  } else if (method == "kernel_weighted_knn") {
    kknn::kknn(stats::as.formula(".outcome ~ .", env = baseenv()),
      train = fit$engine$training, test = x, k = fit$parameters$k,
      distance = fit$parameters$distance, kernel = fit$parameters$kernel,
      scale = fit$parameters$scale)$prob
  } else {
    z <- sweep(sweep(as.matrix(x), 2L, fit$center, "-"), 2L, fit$scale, "/")
    stats::predict(fit$engine, z, proba = TRUE)$probabilities
  }
  p[, fit$codes, drop = FALSE]
}
