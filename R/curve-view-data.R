sb_require_curve_plotting <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    sb_abort(
      "Plotting requires ggplot2; install it with install.packages('ggplot2').",
      "stormbreaker_missing_plot_dependency"
    )
  }
}

sb_curve_view_data <- function(x, kind, cohort_label = NULL) {
  x <- sb_viz_result_input(x)
  if (is.list(x) && !is.object(x) && length(x) && !is.null(names(x))) {
    parts <- lapply(sb_viz_leaves(x, cohort_label), function(leaf) {
      value <- leaf$value
      attr(value, "stormbreaker_view_context") <- leaf$context
      sb_curve_view_data(value, kind, leaf$dataset)
    })
    data <- sb_view_bind(lapply(parts, `[[`, "data"))
    status <- sb_view_bind(lapply(parts, `[[`, "status"))
    data$panel <- factor(data$panel, levels = unique(status$panel))
    return(list(data = data, status = status))
  }
  if (sb_is_composite_view(x) || !is.null(attr(x, "stormbreaker_view_context"))) {
    return(sb_composite_curves(x, kind, cohort_label))
  }
  is_cv <- function(value) {
    inherits(value, c("stormbreaker_cv_result", "stormbreaker_binary_cv_result"))
  }
  is_external <- function(value) inherits(value, "stormbreaker_external_evaluation")
  if (inherits(x, "stormbreaker_cv_comparison")) {
    results <- x$results
    external <- FALSE
  } else if (inherits(x, "stormbreaker_external_comparison")) {
    results <- x$evaluations
    external <- TRUE
  } else if (is_cv(x)) {
    results <- list(x)
    external <- FALSE
  } else if (is_external(x)) {
    results <- list(x)
    external <- TRUE
  } else {
    sb_abort(
      "Supply a CV result/comparison or an external evaluation/comparison.",
      "stormbreaker_bad_curve_view"
    )
  }
  predicate <- if (external) is_external else is_cv
  if (!is.list(results) || length(results) == 0L ||
      !all(vapply(results, predicate, logical(1)))) {
    sb_abort("The comparison contains no supported results.",
             "stormbreaker_bad_curve_view")
  }
  if (!is.null(cohort_label) &&
      (!is.character(cohort_label) || length(cohort_label) != 1L ||
       is.na(cohort_label) || !nzchar(trimws(cohort_label)))) {
    sb_abort("`cohort_label` must be NULL or one non-empty string.",
             "stormbreaker_bad_cohort_label")
  }
  boundary <- if (external) "External evaluation" else "Development CV"
  cohort <- if (is.null(cohort_label)) boundary else cohort_label
  value_or <- function(value, default) if (is.null(value)) default else value

  parts <- lapply(results, function(result) {
    task_type <- if (external) {
      result$task_type
    } else if (inherits(result, "stormbreaker_binary_cv_result")) {
      "binary"
    } else {
      "survival"
    }
    pipeline_id <- if (external) result$pipeline_id else result$pipeline$id
    descriptor <- if (external) result$endpoint else result$task
    key <- if (identical(kind, "dca")) {
      "decision_curve"
    } else if (identical(task_type, "binary")) {
      "moderate_calibration"
    } else {
      "calibration"
    }
    metric <- result$metrics[[key]]
    horizon <- value_or(metric$prediction_horizon,
                        value_or(descriptor$prediction_horizon, NA_real_))
    unit <- if (identical(task_type, "survival")) {
      value_or(metric$time_unit, value_or(descriptor$time_unit, NA_character_))
    } else {
      NA_character_
    }
    status <- if (!is.null(metric)) {
      metric$status
    } else if (identical(result$status, "success") &&
               (identical(kind, "dca") ||
                (identical(task_type, "survival") && is.na(horizon)))) {
      "not_requested"
    } else {
      "not_saved"
    }
    scope <- value_or(metric$scope, NA_character_)
    target <- if (identical(task_type, "survival")) {
      if (is.na(horizon)) "horizon not declared" else paste(horizon, unit)
    } else {
      "class 1"
    }
    panel <- paste0(
      pipeline_id, "\n", cohort,
      if (!identical(cohort, boundary)) paste0(" | ", boundary) else "",
      " | ", target, "\ncurve: ", status, " | run: ", result$status
    )
    curve <- metric$curve
    n_rows <- if (is.data.frame(curve)) nrow(curve) else 0L
    if (n_rows == 0L) {
      coords <- data.frame(
        source_row = NA_integer_, strategy = "model", x = NA_real_, y = NA_real_,
        row_usable = FALSE, stringsAsFactors = FALSE
      )
    } else if (identical(kind, "dca")) {
      columns <- c("threshold_probability", "model_net_benefit",
                   "treat_all_net_benefit", "treat_none_net_benefit")
      if (!all(columns %in% names(curve))) {
        sb_abort("The saved DCA curve lacks its documented coordinate columns.",
                 "stormbreaker_bad_saved_curve")
      }
      row_usable <- if ("status" %in% names(curve)) {
        !is.na(curve$status) & curve$status == "success"
      } else {
        rep(TRUE, n_rows)
      }
      coords <- do.call(rbind, lapply(c("model", "treat_all", "treat_none"),
                                     function(strategy) {
        data.frame(
          source_row = seq_len(n_rows), strategy = strategy,
          x = curve$threshold_probability,
          y = curve[[paste0(strategy, "_net_benefit")]],
          row_usable = row_usable, stringsAsFactors = FALSE
        )
      }))
    } else {
      columns <- if (identical(task_type, "binary")) {
        c("probability", "calibrated_probability")
      } else {
        c("predicted_event_probability", "calibrated_event_probability")
      }
      if (!all(columns %in% names(curve))) {
        sb_abort("The saved calibration curve lacks its documented coordinate columns.",
                 "stormbreaker_bad_saved_curve")
      }
      coords <- data.frame(
        source_row = seq_len(n_rows), strategy = "model",
        x = curve[[columns[[1]]]], y = curve[[columns[[2]]]],
        row_usable = TRUE, stringsAsFactors = FALSE
      )
    }
    coords$draw <- identical(status, "success") &&
      identical(result$status, "success")
    coords$draw <- coords$draw & coords$row_usable &
      is.finite(coords$x) & is.finite(coords$y)
    coords$segment <- cumsum(!coords$draw)
    info <- data.frame(
      pipeline_id = pipeline_id, task_type = task_type,
      validation = if (external) "external_evaluation" else "development_cv",
      cohort = cohort, curve_key = key, metric_status = status,
      run_status = result$status, scope = scope,
      prediction_horizon = horizon, time_unit = unit,
      error_class = value_or(metric$error_class,
                             value_or(result$error_class, NA_character_)),
      error_message = value_or(metric$error_message,
                               value_or(result$error_message, NA_character_)),
      n_curve_rows = as.integer(n_rows), n_drawn_coordinates = as.integer(sum(coords$draw)),
      panel = panel, stringsAsFactors = FALSE
    )
    data <- cbind(info[rep(1L, nrow(coords)), , drop = FALSE], coords)
    list(data = data, status = info)
  })
  data <- do.call(rbind, lapply(parts, `[[`, "data"))
  status <- do.call(rbind, lapply(parts, `[[`, "status"))
  rownames(data) <- rownames(status) <- NULL
  data$panel <- factor(data$panel, levels = status$panel)
  data$line_group <- interaction(data$pipeline_id, data$strategy, data$segment,
                                 drop = TRUE)
  list(data = data, status = status)
}

sb_curve_view_lines <- function(data) {
  rows <- data[data$draw, , drop = FALSE]
  counts <- table(rows$line_group)
  rows[as.integer(counts[as.character(rows$line_group)]) >= 2L, , drop = FALSE]
}

sb_curve_view_annotations <- function(view) {
  info <- view$status
  info <- info[info$n_drawn_coordinates == 0L, , drop = FALSE]
  info$label <- if (nrow(info) > 0L) {
    paste0("No curve\n", info$metric_status, " / run: ", info$run_status)
  } else {
    character()
  }
  info$panel <- factor(info$panel, levels = levels(view$data$panel))
  info
}

sb_curve_view_scope_label <- function(view) {
  scopes <- unique(view$status$scope[!is.na(view$status$scope)])
  if (length(scopes) == 0L) "No saved curve scope" else paste(scopes, collapse = " | ")
}
