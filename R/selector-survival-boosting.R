sb_survival_boosting_selector_ids <- function() c("gbm_cox_importance", "glmboost_cox_nonzero", "superpc_survival_scores", "gamboost_cox_nonzero", "xgboost_cox_importance")

sb_mboost_spline_selector_scores <- function(model, features) {
  beta <- stats::coef(model)
  index <- match(sub("^bbs\\(([^,)]+).*", "\\1", names(beta)), features)
  keep <- !is.na(index)
  score <- numeric(length(features))
  score[index[keep]] <- vapply(beta[keep], function(b) sqrt(sum(b^2)), numeric(1))
  score
}

sb_survival_boosting_selector_params <- function(id) {
  if (id == "xgboost_cox_importance") return(sb_xgboost_params("xgboost_cox"))
  if (id == "superpc_survival_scores") return(sb_superpc_params())
  if (id == "gamboost_cox_nonzero") return(sb_mboost_cox_params("gamboost_cox"))
  if (id == "gbm_cox_importance") sb_gbm_cox_params() else
    sb_mboost_cox_params("glmboost_cox")
}

sb_survival_boosting_selector_adapter <- function(id, params) {
  list(id = id, params = params, needs_outcome = TRUE,
    fit = function(x, time, event) {
      x <- sb_validate_predictor_frame(x, "x")
      work <- x
      names(work) <- paste0("feature_", seq_len(ncol(x)))
      if (id == "superpc_survival_scores") {
        inputs <- sb_validate_classical_survival_inputs(work, time, event, "SuperPC selector")
        if (!requireNamespace("superpc", quietly = TRUE)) stop("Install superpc for SuperPC selection.")
        if (ncol(work) < 2L) stop("SuperPC selection needs at least two input features.")
        model <- superpc::superpc.train(list(x = t(as.matrix(inputs$x)),
          y = inputs$time, censoring.status = inputs$event, featurenames = names(work)),
          type = "survival", s0.perc = params$s0_perc)
        score <- as.numeric(model$feature.scores)
        selected <- abs(score) >= params$threshold
        history <- numeric()
      } else if (id == "xgboost_cox_importance") {
        model <- sb_fit_xgboost(work, list(time = time, event = event), "xgboost_cox", params)
        score <- sb_xgboost_selector_gain(model)
        selected <- score > 0
        history <- numeric()
      } else if (id == "gbm_cox_importance") {
        model <- sb_fit_gbm_cox(work, time, event, params)$booster
        importance <- summary(model, n.trees = params$n_trees, plotit = FALSE)
        score <- importance$rel.inf[match(names(work), importance$var)]
        selected <- score > 0
        history <- model$train.error
      } else {
        spline <- id == "gamboost_cox_nonzero"
        model <- sb_fit_mboost_cox(work, time, event,
          if (spline) "gamboost_cox" else "glmboost_cox", params)$engine
        if (spline) score <- sb_mboost_spline_selector_scores(model, names(work)) else {
          beta <- stats::coef(model, off2int = FALSE)
          score <- numeric(ncol(work))
          index <- match(names(beta), names(work))
          keep <- !is.na(index)
          score[index[keep]] <- as.numeric(beta[keep])
        }
        selected <- score != 0
        history <- mboost::risk(model)
      }
      if (length(score) != ncol(x) || any(!is.finite(score))) stop("Invalid native selector scores.")
      if (!any(selected)) stop("No features met the selection rule; no fallback genes were added.")
      structure(list(id = id, params = params, source_features = names(x),
        selected_features = names(x)[selected],
        scores = data.frame(feature = names(x), score = score, selected = selected),
        evidence = list(training_loss = data.frame(
          step = seq_along(history) - as.integer(id %in% c("glmboost_cox_nonzero", "gamboost_cox_nonzero")), loss = as.numeric(history)))),
        class = "stormbreaker_selector_fit")
    },
    apply = function(fitted, x) {
      x <- sb_validate_predictor_frame(x, "x")
      if (!identical(names(x), fitted$source_features)) stop("Predictor columns differ from the fitted selector.")
      x[fitted$selected_features]
    })
}
