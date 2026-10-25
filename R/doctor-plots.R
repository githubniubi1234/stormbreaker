sb_compact_plot <- function(p) {
  notes <- unique(c(attr(p, "stormbreaker_notes"), attr(p, "stormbreaker_figure_notes"),
    p$labels$subtitle, p$labels$caption))
  title <- p$labels$title
  if (is.character(title)) {
    title <- sub("^(Original|Native|Saved) ", "", title)
    title <- paste0(toupper(substr(title, 1, 1)), substring(title, 2))
  }
  p <- p + ggplot2::labs(title = title, subtitle = NULL, caption = NULL)
  attr(p, "stormbreaker_notes") <- notes
  p
}

sb_small_model_plots <- function(fit) {
  id <- fit$pipeline$learner
  if (!id %in% c("nearest_centroid_binary", "bernoulli_naive_bayes_binary")) return(list())
  model <- fit$diagnostic_artifacts$learner
  if (is.null(model)) stop("Refit with capture_diagnostics=TRUE to retain algorithm evidence.")
  m <- if (id == "nearest_centroid_binary") model$engine$centroids else model$engine$feature_probability
  d <- data.frame(feature = rep(colnames(m), each = 2), class = rep(c("0", "1"), ncol(m)), value = as.vector(m))
  feature <- class <- value <- NULL
  list(class_profiles = ggplot2::ggplot(d, ggplot2::aes(value, feature, colour = class)) +
    ggplot2::geom_point(position = ggplot2::position_dodge(width = .4), size = 2.4) +
    ggplot2::labs(x = if (id == "nearest_centroid_binary") "Centroid" else "P(feature = 1 | class)", y = NULL, colour = "Class") + sb_viz_theme())
}

#' Draw available algorithm and selector evidence through one entrypoint
#' @param fit Successful final fit (capture_diagnostics=TRUE), or a saved
#'   multiclass fit from sb_analyze_methods(). Multiclass SHAP is not supported.
#' @param new_data Observations for local diagnostics or SHAP. Required for
#'   multiclass prediction neighborhoods; training rows are not substituted.
#' @param id_col Observation ID column for SHAP.
#' @param background_data Optional development reference data. Supplying this
#'   and new_data explicitly requests SHAP (approximate for supported trees,
#'   LogitBoost and RDA).
#' @return List with plots, diagnostics/status and optional SHAP evidence.
#'   Nothing is refitted. Missing evidence is recorded, never invented.
#'   Explanatory metadata stays in returned results, not tiny plot footnotes.
#' @export
sb_plots <- function(fit, new_data = NULL, id_col = NULL, background_data = NULL) {
  if (inherits(fit, "stormbreaker_multiclass_fit")) {
    if (!is.null(background_data)) stop("Multiclass SHAP is not implemented through sb_plots.", call. = FALSE)
    drawn <- sb_multiclass_plots(fit, sb_multiclass_evidence(fit), new_data)
    plots <- lapply(drawn$plots, sb_compact_plot)
    return(list(plots = plots, pages = sb_plot_pages(plots), diagnostics = drawn$tables,
      notes = lapply(plots, attr, which = "stormbreaker_notes")))
  }
  fit <- sb_validate_final_fit_for_prediction(fit)
  sb_require_curve_plotting()
  plots <- list()
  status <- list()
  add <- function(name, expression) {
    p <- tryCatch(force(expression), error = identity)
    if (inherits(p, "error")) status[[name]] <<- conditionMessage(p) else plots[[name]] <<- sb_compact_plot(p)
  }
  if (fit$pipeline$learner %in% c("nearest_centroid_binary", "bernoulli_naive_bayes_binary")) {
    small <- tryCatch(sb_small_model_plots(fit), error = identity)
    if (inherits(small, "error")) status$native <- conditionMessage(small) else plots <- small
    native <- NULL
  } else {
    diagnostic_data <- if (is.null(new_data)) fit$diagnostic_artifacts$predictors else new_data
    native <- sb_native_diagnostics(fit, complete = TRUE, newdata = diagnostic_data)
    missing <- native$status[native$status$status != "available", , drop = FALSE]
    if (nrow(missing)) status$native <- missing
    for (name in names(native$diagnostics)) add(name, sb_plot_algorithm_diagnostic(native$diagnostics[[name]]))
  }
  state <- fit$selection$state
  stages <- if (!is.null(state$stages)) state$stages else list(state)
  for (i in seq_along(stages)) {
    s <- stages[[i]]
    if (is.null(s) || s$id == "none") next
    if (s$id %in% sb_abess_selector_ids()) {
      criterion <- s$criterion_path
      criterion$selected <- criterion$support_size == s$selected_support_size
      diagnostic <- structure(list(type = "support_criterion", table = criterion,
        selected_support_size = s$selected_support_size), class = "stormbreaker_algorithm_diagnostic")
      add(paste0("selector_", i, "_criterion"), sb_plot_algorithm_diagnostic(diagnostic))
      if (is.null(s$coefficient_path)) status[[paste0("selector_", i, "_coefficient_path")]] <- "ABESS path was not retained; refit with capture_diagnostics=TRUE."
      else {
        b <- as.matrix(s$coefficient_path$beta)
        diagnostic$type <- "abess_coefficient_path"
        diagnostic$table <- data.frame(support_size = rep(s$coefficient_path$support_size, each = nrow(b)),
          feature = rep(s$source_features, ncol(b)), coefficient = as.vector(b))
        add(paste0("selector_", i, "_coefficient_path"), sb_plot_algorithm_diagnostic(diagnostic))
      }
    } else if (s$id %in% sb_survival_weight_selector_ids()) {
      d <- s$scores
      d$selected <- ifelse(d$selected, "Yes", "No")
      feature <- score <- selected <- NULL
      add(paste0("selector_", i), ggplot2::ggplot(d,
        ggplot2::aes(score, stats::reorder(feature, score), colour = selected)) +
        ggplot2::geom_point(size = 2.3) +
        ggplot2::scale_colour_manual(values = c(Yes = "#2878A0", No = "grey65")) +
        ggplot2::labs(x = "Screening weight", y = NULL, colour = "Selected") + sb_viz_theme())
      if (s$id == "ridge_cox_weights") {
        tab <- s$evidence$tuning_curve
        tab$mean_value <- tab$mean_deviance
        diagnostic <- list(type = "cv_curve", table = tab, measure = "CV deviance",
          selected_lambda = s$evidence$lambda, selection_rule = "lambda.min", alpha = 0)
        add(paste0("selector_", i, "_cv"), sb_plot_glmnet_tuning_diagnostic(diagnostic))
        add(paste0("selector_", i, "_coefficient_path"),
          sb_selector_path_plot(s$evidence$coefficient_path, s$evidence$lambda))
      }
      if (s$id == "plsrcox_weights" && !is.null(s$evidence$loadings)) {
        m <- s$evidence$loadings
        tab <- data.frame(feature = rep(s$source_features, ncol(m)),
          component = rep(seq_len(ncol(m)), each = nrow(m)), loading = as.vector(m))
        add(paste0("selector_", i, "_loadings"), sb_plot_pls_diagnostic(list(
          type = "pls_loadings", table = tab, features = s$source_features,
          loading_definition = "PLS loadings", interpretation = NULL, learner = "PLS-Cox")))
      }
    } else if (s$id %in% sb_survival_boosting_selector_ids()) {
      d <- s$scores
      d$selected <- ifelse(d$selected, "Yes", "No")
      feature <- score <- selected <- step <- loss <- NULL
      add(paste0("selector_", i), ggplot2::ggplot(d,
        ggplot2::aes(score, stats::reorder(feature, score), colour = selected)) +
        ggplot2::geom_point(size = 2.3, show.legend = length(unique(d$selected)) > 1L) +
        ggplot2::scale_colour_manual(values = c(Yes = "#2878A0", No = "grey65")) +
        ggplot2::labs(y = NULL, colour = "Selected",
          x = if (s$id == "xgboost_cox_importance") "Gain" else if (s$id == "gbm_cox_importance") "Importance" else
            if (s$id == "superpc_survival_scores") "SuperPC score" else
              if (s$id == "gamboost_cox_nonzero") "Spline coefficient norm" else "Coefficient") + sb_viz_theme())
      if (s$id == "superpc_survival_scores") {
        name <- paste0("selector_", i)
        plots[[name]] <- plots[[name]] + ggplot2::geom_vline(xintercept = c(-1, 1) * s$params$threshold, linetype = 2)
      } else if (nrow(s$evidence$training_loss)) add(paste0("selector_", i, "_training_loss"), ggplot2::ggplot(s$evidence$training_loss,
        ggplot2::aes(step, loss)) + ggplot2::geom_line() +
        ggplot2::labs(x = "Boosting step", y = "Training loss") + sb_viz_theme())
    } else if (s$id == "svm_rfe_binary") {
      d <- s$evidence$curve
      n_features <- value <- NULL
      for (measure in c("accuracy", "error")) {
        d$value <- d[[measure]]
        add(paste0("selector_", i, "_cv_", measure),
          ggplot2::ggplot(d, ggplot2::aes(n_features, value)) +
            ggplot2::geom_line(colour = "#2878A0") + ggplot2::geom_point() +
            ggplot2::geom_vline(xintercept = d$n_features[d$selected], linetype = 2) +
            ggplot2::labs(x = "Number of features", y = paste("CV", measure)) + sb_viz_theme())
      }
    } else if (s$id %in% c(sb_binary_model_selector_ids(), sb_pamr_selector_ids(), "lasso_nonzero", "enet_cox_nonzero", "coxboost_nonzero", "stepcox")) {
      d <- s$scores
      if (!"score" %in% names(d)) d$score <- d$coefficient
      d <- d[is.finite(d$score), ]
      d$selected <- ifelse(d$selected, "Yes", "No")
      feature <- score <- selected <- NULL
      add(paste0("selector_", i), ggplot2::ggplot(d, ggplot2::aes(score, stats::reorder(feature, score), colour = selected)) +
        ggplot2::geom_point(size = 2.3) + ggplot2::labs(x = if (s$id %in% sb_pamr_selector_ids()) "Shrunken-centroid score" else if (s$id == "plsrglm_binary_weights") "|Coefficient| \u00d7 SD" else if (s$id == "xgboost_binary_importance") "Gain" else if (s$id %in% c("rf_binary_importance", "gbm_binary_importance")) "Importance" else if (s$id == "gamboost_binary_nonzero") "Spline coefficient norm" else "Coefficient", y = NULL, colour = "Selected") + sb_viz_theme())
      if (s$id == "gbm_binary_importance") {
        trees <- loss <- NULL
        add(paste0("selector_", i, "_training_loss"), ggplot2::ggplot(s$evidence$training_loss, ggplot2::aes(trees, loss)) +
          ggplot2::geom_line() + ggplot2::labs(x = "Trees", y = "Training loss") + sb_viz_theme())
      }
      if (s$id %in% c("stepglm_binary_aic", "glmboost_binary_nonzero", "gamboost_binary_nonzero")) {
        history <- if (s$id == "stepglm_binary_aic") s$evidence$selection_path else s$evidence$training_loss
        if (is.null(history)) status[[paste0("selector_", i, "_history")]] <- "History was not retained; refit with capture_diagnostics=TRUE."
        else {
          d <- data.frame(step = seq_len(nrow(history)) - 1L,
            value = if (s$id == "stepglm_binary_aic") history$AIC else history$loss)
          step <- value <- NULL
          add(paste0("selector_", i, "_history"), ggplot2::ggplot(d, ggplot2::aes(step, value)) +
            ggplot2::geom_line() + ggplot2::geom_point() +
            ggplot2::labs(x = "Step", y = if (s$id == "stepglm_binary_aic") "AIC" else "Training loss") + sb_viz_theme())
        }
      }
      curve <- if (!is.null(s$tuning_curve)) s$tuning_curve else s$evidence$tuning_curve
      if (s$id == "coxboost_nonzero") {
        diagnostic <- structure(list(type = "coxboost_cv_curve", table = curve,
          stopping_iteration = s$optimal_step, learner = s$id, measure = "CoxBoost CV"),
          class = "stormbreaker_algorithm_diagnostic")
        add(paste0("selector_", i, "_cv"), sb_plot_coxboost_diagnostic(diagnostic))
        if (is.null(s$coefficient_path)) status[[paste0("selector_", i, "_coefficient_path")]] <- "CoxBoost path was not retained; refit with capture_diagnostics=TRUE."
        else {
          b <- s$coefficient_path
          diagnostic$type <- "coxboost_coefficient_path"
          diagnostic$measure <- "CoxBoost coefficient path"
          diagnostic$table <- data.frame(step = rep(seq_len(nrow(b)) - 1L, ncol(b)),
            feature = rep(s$source_features, each = nrow(b)), coefficient = as.vector(b))
          add(paste0("selector_", i, "_coefficient_path"), sb_plot_coxboost_diagnostic(diagnostic))
        }
      }
      selected_lambda <- if (!is.null(s$lambda)) s$lambda else s$evidence$lambda
      if (s$id %in% c("lasso_binary_nonzero", "enet_binary_nonzero", "lasso_nonzero", "enet_cox_nonzero")) {
        path <- if (!is.null(s$coefficient_path)) s$coefficient_path else s$evidence$coefficient_path
        add(paste0("selector_", i, "_coefficient_path"), sb_selector_path_plot(path, selected_lambda))
      }
      if (s$id == "rf_binary_importance") add(paste0("selector_", i, "_oob_error"),
        sb_selector_oob_plot(s$evidence$oob_error))
      if (is.data.frame(curve) && "lambda" %in% names(curve)) {
        y <- intersect(c("mean_deviance", "mean_auc", "mean_value"), names(curve))[1L]
        if (!is.na(y)) {
          cvdata <- data.frame(lambda = curve$lambda, mean = curve[[y]], se = curve$standard_error)
          lambda <- mean <- se <- NULL
          ticks <- unique(round(seq(1, nrow(curve), length.out = min(7L, nrow(curve)))))
          axis <- if ("nonzero" %in% names(curve)) ggplot2::scale_x_log10(
            sec.axis = ggplot2::dup_axis(breaks = curve$lambda[ticks],
              labels = curve$nonzero[ticks], name = "Nonzero coefficients")) else ggplot2::scale_x_log10()
          add(paste0("selector_", i, "_cv"), ggplot2::ggplot(cvdata, ggplot2::aes(lambda, mean)) +
            ggplot2::geom_errorbar(ggplot2::aes(ymin = mean - se, ymax = mean + se), width = 0, colour = "grey70") +
            ggplot2::geom_point(colour = "#2878B5", size = 1.5) + axis +
            ggplot2::geom_vline(xintercept = selected_lambda, linetype = 2, colour = "grey50") +
            ggplot2::labs(x = "Lambda", y = if (y == "mean_auc") "CV AUC" else "CV loss") + sb_viz_theme())
        }
      }
      if (s$id == "stepcox" && is.data.frame(s$selection_path)) {
        path <- data.frame(step = seq_len(nrow(s$selection_path)) - 1L, AIC = s$selection_path$AIC)
        step <- AIC <- NULL
        add(paste0("selector_", i, "_aic"), ggplot2::ggplot(path, ggplot2::aes(step, AIC)) +
          ggplot2::geom_line() + ggplot2::geom_point() + ggplot2::labs(x = "Step", y = "AIC") + sb_viz_theme())
      }
    } else {
      views <- sb_selector_score_views(s$id)
      for (v in names(views)) add(paste0("selector_", i, "_", v), sb_plot_selector_scores(s, v))
      types <- if (s$id %in% sb_boruta_selector_ids()) c("boruta_summary", "boruta_history", "boruta_boxplot") else
        if (s$id == "varselrf_binary_oob_1se") "selection_history" else
          if (s$id == "vsurf_binary_prediction") c("vsurf_importance", "vsurf_stages") else character()
      for (t in types) add(paste0("selector_", i, "_", t), sb_plot_selector_diagnostic(sb_selector_diagnostics(s, t)))
      if (!length(types) && !length(views)) status[[paste0("selector_", i)]] <- "No specialized selector plot route; selected genes remain in fit$selection."
    }
  }
  shap <- NULL
  if (!is.null(background_data)) {
    shap <- tryCatch(if (fit$pipeline$learner %in% sb_tree_shap_ids()) {
      sb_tree_shap(fit, new_data, id_col, background_data)
    } else if (fit$pipeline$learner %in% c("logitboost_binary", "rda_binary")) {
      sb_model_shap(fit, new_data, id_col, background_data)
    } else sb_additive_shap(fit, new_data, id_col, background_data), error = identity)
    if (inherits(shap, "error")) { status$shap <- conditionMessage(shap); shap <- NULL } else {
      add("shap_importance", sb_plot_shap_importance(shap))
      add("shap_beeswarm", sb_plot_shap_beeswarm(shap, feature_values = new_data, id_col = id_col))
      add("shap_waterfall", sb_plot_shap_waterfall(shap))
      plots$shap_importance <- plots$shap_importance + ggplot2::labs(title = "SHAP importance", x = "Mean |SHAP|")
      plots$shap_beeswarm <- plots$shap_beeswarm + ggplot2::labs(title = "SHAP summary", x = "SHAP value", colour = "Feature value")
      plots$shap_waterfall <- plots$shap_waterfall + ggplot2::labs(title = "Individual prediction",
        x = if (shap$score_scale == "class_1_probability") "Predicted probability" else "Model score")
      top <- utils::head(names(sort(colMeans(abs(shap$contributions)), decreasing = TRUE)), 6L)
      for (gene in top) add(paste0("shap_dependence_", gene), sb_plot_shap_dependence(shap, gene, new_data, id_col))
    }
  }
  notes <- lapply(plots, function(p) attr(p, "stormbreaker_notes"))
  if (!is.null(native)) for (name in names(native$diagnostics)) {
    notes[[name]] <- unique(c(notes[[name]], native$diagnostics[[name]]$interpretation))
    if (name %in% names(plots)) attr(plots[[name]], "stormbreaker_notes") <- notes[[name]]
  }
  list(plots = plots, pages = sb_plot_pages(plots), notes = notes,
    unavailable = status, native = native, shap = shap)
}
