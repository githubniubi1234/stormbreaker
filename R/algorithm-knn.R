#' Explain original nearest-neighbor predictions
#' @param fit A successful ordinary or kernel-weighted binary KNN final fit.
#' @param newdata Explicit query rows with the original fitted features.
#'   NULL reports unavailable rather than silently querying the training set.
#' @return A diagnostic with one row per query and neighbor, original training
#'   row index, distance, class, native weight and normalized voting weight.
#'   Prediction probabilities are returned separately. Ordinary KNN uses its
#'   original standardization, Euclidean distance and training-row tie rule.
#'   Weighted KNN uses the original kknn distance/kernel/scaling settings.
#'   Training-row indices and classes are sensitive. Query rows that duplicate
#'   training rows retain self-neighbors, exactly as prediction does; this is
#'   not leave-one-out validation. No tuning or patient card is generated.
#' @export
sb_knn_neighborhood <- function(fit, newdata = NULL) {
  fit <- sb_validate_final_fit_for_prediction(fit)
  s <- fit$learner_state
  if (!s$learner %in% c("knn_binary", "kernel_weighted_knn_binary")) sb_abort(
    "This neighborhood view requires binary KNN.", "stormbreaker_unsupported_diagnostic")
  if (is.null(newdata)) sb_abort("Supply explicit query data for a nearest-neighbor explanation.",
    "stormbreaker_unavailable_diagnostic")
  x <- sb_select_final_feature_frame(newdata, s$features, "Query data", "query fitted features")
  if (!nrow(x)) sb_abort("Query data must contain at least one row.", "stormbreaker_invalid_diagnostic")
  e <- unserialize(s$engine_raw)
  n <- nrow(x)
  k <- e$k
  if (s$learner == "knn_binary") {
    z <- sweep(sweep(as.matrix(x), 2L, e$center, "-"), 2L, e$scale, "/")
    native <- lapply(seq_len(n), function(i) {
      distance <- rowSums(sweep(e$training, 2L, z[i, ], "-")^2)
      index <- order(distance, seq_along(distance))[seq_len(k)]
      list(index = index, distance = sqrt(distance[index]))
    })
    index <- do.call(rbind, lapply(native, `[[`, "index"))
    distance <- do.call(rbind, lapply(native, `[[`, "distance"))
    weight <- matrix(1 / k, n, k)
    classes <- matrix(as.character(e$outcome[index]), n, k)
    probability <- rowMeans(classes == "1")
    metric <- "Euclidean distance in original training-standardized feature space"
  } else {
    training <- e$training
    outcome_name <- ".stormbreaker_outcome"
    while (outcome_name %in% names(training)) outcome_name <- paste0(outcome_name, "_")
    training[[outcome_name]] <- e$outcome
    rhs <- paste(vapply(names(x), sb_quote_name, ""), collapse = " + ")
    formula <- stats::as.formula(paste(sb_quote_name(outcome_name), "~", rhs), env = baseenv())
    native <- kknn::kknn(formula, train = training, test = x, k = k,
      distance = e$distance, kernel = e$kernel, scale = e$scale)
    index <- native$C
    distance <- native$D
    weight <- native$W
    classes <- native$CL
    probability <- native$prob[, "1"]
    metric <- "Native kknn distance under original distance/kernel/scaling settings"
  }
  if (any(vapply(list(index, distance, weight, classes), function(z)
      !is.matrix(z) || !identical(dim(z), c(n, as.integer(k))), logical(1))) ||
      any(!is.finite(index)) || any(index != floor(index)) || any(index < 1 | index > nrow(e$training)) ||
      any(!is.finite(distance)) || any(distance < 0) || any(!is.finite(weight)) ||
      any(weight < 0) || any(rowSums(weight) <= 0) || any(!classes %in% c("0", "1"))) sb_abort(
    "Original nearest-neighbor result is invalid.", "stormbreaker_invalid_diagnostic")
  normalized <- weight / rowSums(weight)
  tab <- data.frame(query = rep(seq_len(n), each = k), rank = rep(seq_len(k), n),
    training_row = as.integer(t(index)), distance = as.numeric(t(distance)),
    class = as.character(t(classes)), native_weight = as.numeric(t(weight)),
    voting_weight = as.numeric(t(normalized)))
  structure(list(type = "knn_neighborhood", table = tab,
    predictions = data.frame(query = seq_len(n), probability = as.numeric(probability)),
    settings = s$params,
    learner = s$learner, pipeline_id = fit$pipeline$id, sensitive_training_data = TRUE,
    source = metric, interpretation = paste("Explicit query rows; original neighbor selection and voting.",
      "Native training row indices and classes are sensitive. Self-neighbors are not removed.",
      "Not validation, a tuning curve, causal effects or a patient prediction card.")),
    class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_knn_neighborhood <- function(x) {
  distance <- rank <- class <- voting_weight <- query <- NULL
  p <- ggplot2::ggplot(x$table, ggplot2::aes(rank, distance, colour = class, size = voting_weight)) +
    ggplot2::geom_point(alpha = 0.85) +
    ggplot2::facet_wrap(~query, labeller = ggplot2::label_both) +
    ggplot2::scale_colour_manual(values = c("0" = "#2878A0", "1" = "#B86236")) +
    ggplot2::scale_x_continuous(breaks = function(z) unique(round(pretty(z)))) +
    ggplot2::labs(title = "Prediction neighborhoods", x = "Neighbor rank", y = "Distance",
      colour = "Class", size = "Voting weight") + sb_viz_theme()
  attr(p, "stormbreaker_suggested_size") <- c(width = 11, height = 5)
  p
}
