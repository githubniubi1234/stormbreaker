sb_multiclass_native_catalog <- function() {
  basic <- data.frame(method = c("svm_linear", "svm_radial", "svm_radial_weighted",
    "svm_polynomial", "svm_sigmoid", "kernel_naive_bayes", "rda", "nnet", "knn"),
    label = c("Linear SVM", "Radial SVM", "Weighted radial SVM", "Polynomial SVM",
      "Sigmoid SVM", "Kernel naive Bayes", "Regularized discriminant analysis",
      "Neural network", "K-nearest neighbors"), supported = TRUE,
    engine = c(rep("e1071_svm", 5), "klaR_NaiveBayes", "klaR_rda", "nnet_softmax", "knn_votes"),
    feature_result = "model_features", selection_rule = "No standalone gene-selection rule")
  rbind(basic, sb_multiclass_tree_catalog(), sb_multiclass_discriminant_catalog(),
    sb_multiclass_boosting_catalog(), sb_multiclass_specialist_catalog(),
    sb_multiclass_multinomial_catalog())
}

sb_fit_multiclass_native <- function(fit, x, y) {
  if (fit$method %in% sb_multiclass_multinomial_catalog()$method) return(sb_fit_multiclass_multinomial(fit, x, y))
  if (fit$method %in% sb_multiclass_specialist_catalog()$method) return(sb_fit_multiclass_specialist(fit, x, y))
  if (fit$method %in% sb_multiclass_boosting_catalog()$method) return(sb_fit_multiclass_boosting(fit, x, y))
  if (fit$method %in% sb_multiclass_discriminant_catalog()$method) return(
    sb_fit_multiclass_discriminant(fit, x, y))
  if (fit$method %in% sb_multiclass_tree_catalog()$method) return(sb_fit_multiclass_tree(fit, x, y))
  method <- fit$method
  package <- if (startsWith(method, "svm_")) "e1071" else
    switch(method, kernel_naive_bayes = , rda = "klaR", nnet = "nnet", knn = NA_character_)
  if (!is.na(package) && !requireNamespace(package, quietly = TRUE)) sb_abort(
    paste("Install", package, "for this method."), "stormbreaker_missing_dependency")
  fit$package <- package
  params <- sb_classical_binary_params(paste0(method, "_binary"))
  fit$parameters <- params
  record_warning <- function(w) fit$warnings <<- c(fit$warnings, conditionMessage(w))
  if (startsWith(method, "svm_")) {
    weights <- if (method == "svm_radial_weighted") {
      counts <- table(y)
      stats::setNames(length(y) / (length(counts) * as.numeric(counts)), names(counts))
    } else NULL
    fit$engine <- withCallingHandlers(withr::with_seed(fit$seed,
      e1071::svm(x, y, type = "C-classification", kernel = params$kernel,
        cost = params$cost, gamma = 1 / ncol(x), degree = 3L, coef0 = 0,
        scale = TRUE, probability = TRUE, class.weights = weights)), warning = record_warning)
    fit$parameters$gamma <- 1 / ncol(x)
    fit$parameters$class_weights <- weights
    fit$parameters$multiclass_strategy <- "native_one_vs_one_probability_coupling"
    fit$parameters$probability_calibration <- "libsvm_internal_row_based"
  } else if (method == "kernel_naive_bayes") {
    fit$engine <- withCallingHandlers(klaR::NaiveBayes(x, y, usekernel = TRUE,
      fL = params$fL, adjust = params$adjust), warning = record_warning)
  } else if (method == "rda") {
    fit$engine <- withCallingHandlers(klaR::rda(x, y, gamma = params$gamma,
      lambda = params$lambda, estimate.error = FALSE), warning = record_warning)
  } else {
    fit$center <- colMeans(x)
    fit$scale <- vapply(x, stats::sd, numeric(1))
    fit$scale[fit$scale == 0] <- 1
    z <- sweep(sweep(as.matrix(x), 2L, fit$center, "-"), 2L, fit$scale, "/")
    if (method == "nnet") {
      target <- outer(as.character(y), fit$codes, `==`) * 1
      colnames(target) <- fit$codes
      fit$engine <- withCallingHandlers(withr::with_seed(fit$seed,
        nnet::nnet(z, target, size = params$size, decay = params$decay,
          maxit = params$maxit, rang = params$rang, softmax = TRUE,
          MaxNWts = params$max_weights, trace = FALSE)), warning = record_warning)
      fit$parameters$entropy <- NULL
      fit$parameters$softmax <- TRUE
      fit$parameters$convergence <- fit$engine$convergence
      if (fit$engine$convergence != 0L) fit$warnings <- c(fit$warnings,
        "nnet reached its iteration limit; convergence was not established.")
    } else {
      if (nrow(x) < params$k) sb_abort("KNN requires at least k training rows.",
        "stormbreaker_bad_multiclass_fit")
      fit$engine <- list(training = z, outcome = y, k = params$k)
    }
  }
  fit$parameters$seed <- fit$seed
  structure(fit, class = "stormbreaker_multiclass_fit")
}

sb_predict_multiclass_native <- function(fit, x) {
  if (fit$method %in% sb_multiclass_multinomial_catalog()$method) return(sb_predict_multiclass_multinomial(fit, x))
  if (fit$method %in% sb_multiclass_specialist_catalog()$method) return(sb_predict_multiclass_specialist(fit, x))
  if (fit$method %in% sb_multiclass_boosting_catalog()$method) return(sb_predict_multiclass_boosting(fit, x))
  if (fit$method %in% sb_multiclass_discriminant_catalog()$method) return(
    sb_predict_multiclass_discriminant(fit, x))
  if (fit$method %in% sb_multiclass_tree_catalog()$method) return(sb_predict_multiclass_tree(fit, x))
  if (!is.na(fit$package) && !requireNamespace(fit$package, quietly = TRUE)) sb_abort(
    paste("Install", fit$package, "to predict with this saved fit."), "stormbreaker_missing_dependency")
  if (startsWith(fit$method, "svm_")) return(
    attr(stats::predict(fit$engine, x, probability = TRUE), "probabilities")[, fit$codes, drop = FALSE])
  if (fit$method == "rda") {
    predict_method <- get("predict.rda", envir = asNamespace("klaR"), inherits = FALSE)
    return(predict_method(fit$engine, newdata = x)$posterior[, fit$codes, drop = FALSE])
  }
  if (fit$method == "kernel_naive_bayes") return(
    stats::predict(fit$engine, x)$posterior[, fit$codes, drop = FALSE])
  z <- sweep(sweep(as.matrix(x), 2L, fit$center, "-"), 2L, fit$scale, "/")
  if (fit$method == "nnet") return(stats::predict(fit$engine, z, type = "raw")[, fit$codes, drop = FALSE])
  t(vapply(seq_len(nrow(z)), function(i) {
    distance <- rowSums(sweep(fit$engine$training, 2L, z[i, ], "-")^2)
    index <- order(distance, seq_along(distance))[seq_len(fit$engine$k)]
    tabulate(as.integer(fit$engine$outcome[index]), nbins = length(fit$codes)) / fit$engine$k
  }, numeric(length(fit$codes))))
}

sb_multiclass_native_evidence <- function(fit) {
  if (fit$method %in% sb_multiclass_boosting_catalog()$method) return(sb_multiclass_boosting_evidence(fit))
  if (fit$method %in% sb_multiclass_tree_catalog()$method) return(sb_multiclass_tree_evidence(fit))
  list(features = data.frame(feature = fit$features, selected = NA), by_class = NULL)
}

sb_multiclass_native_plots <- function(fit, new_data = NULL) {
  if (fit$method == "loclda") return(sb_multiclass_loclda_plots(fit, new_data))
  if (fit$method %in% c("knn", "kernel_weighted_knn")) return(sb_multiclass_neighbor_plots(fit, new_data))
  if (fit$method == "pamr") return(sb_multiclass_pamr_plots(fit))
  if (fit$method == "gam") return(sb_multiclass_gam_plots(fit))
  if (fit$method == "pre") return(sb_multiclass_pre_plots(fit))
  if (fit$method == "gaussian_process") return(sb_multiclass_gp_plots(fit))
  if (fit$method %in% c("cforest", "bagging", "rlt", "aorsf", "wsrf")) return(sb_multiclass_forest_structure_plots(fit))
  if (fit$method %in% c("edda", "mclustda", "hdda", "hda")) return(sb_multiclass_specialist_plots(fit))
  if (fit$method %in% c("liblinear_l2_logistic", "liblinear_l1_logistic", "xgboost_linear")) return(sb_multiclass_linear_weight_plots(fit))
  if (fit$method %in% c("fda", "mda", "pda", "slda", "sda")) return(sb_multiclass_projection_plots(fit))
  if (fit$method %in% c("kernel_naive_bayes", "rda")) return(sb_multiclass_distribution_plots(fit))
  if (fit$method %in% c("logistic", "stepglm", "brglm2")) return(sb_multiclass_logistic_plots(fit))
  if (fit$method == "nnet") return(sb_multiclass_neural_plots(fit))
  if (fit$method %in% sb_multiclass_boosting_catalog()$method) return(sb_multiclass_boosting_plots(fit))
  if (fit$method %in% sb_multiclass_tree_catalog()$method) return(sb_multiclass_tree_plots(fit))
  if (startsWith(fit$method, "svm_")) {
    codes <- fit$engine$levels[fit$engine$labels]
    tab <- data.frame(class = fit$classes[match(codes, fit$codes)], count = fit$engine$nSV)
    tab$class <- factor(tab$class, levels = fit$classes)
    class <- count <- NULL
    plot <- ggplot2::ggplot(tab, ggplot2::aes(class, count)) +
      ggplot2::geom_col(fill = "#2878A0") + sb_viz_theme() +
      ggplot2::labs(title = "SVM: Support vectors", x = NULL, y = "Support vectors")
    decisions <- fit$engine$decision.values
    if (is.null(decisions)) return(list(plots = list(support_vectors = plot), tables = list(support_vectors = tab)))
    pairs <- colnames(decisions)
    for (i in order(nchar(fit$codes), decreasing = TRUE))
      pairs <- gsub(fit$codes[i], fit$classes[i], pairs, fixed = TRUE)
    scores <- data.frame(sample = rep(seq_len(nrow(decisions)), ncol(decisions)),
      pair = rep(pairs, each = nrow(decisions)), decision = as.vector(decisions))
    decision <- NULL
    distribution <- ggplot2::ggplot(scores, ggplot2::aes(decision)) +
      ggplot2::stat_ecdf(geom = "step", colour = "#2878A0") +
      ggplot2::geom_vline(xintercept = 0, linetype = 2, colour = "grey60") +
      ggplot2::facet_wrap(~pair) + sb_viz_theme() +
      ggplot2::labs(title = "SVM: Pairwise decision scores", x = "Decision score", y = "Cumulative fraction",
        caption = "Original training-row one-vs-one decision scores, not probabilities or held-out performance. Every training row is retained for every native class pair; zero is the pairwise decision boundary.")
    return(list(plots = list(support_vectors = plot, decision_distribution = distribution),
      tables = list(support_vectors = tab, decision_distribution = scores)))
  }
  list(plots = list(), tables = list())
}
