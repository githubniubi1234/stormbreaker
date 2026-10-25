sb_is_composite_view <- function(x) {
  inherits(x, c("stormbreaker_repeated_cv_comparison",
                "stormbreaker_internal_external_cv", "stormbreaker_stack_cv_result"))
}

sb_view_bind <- function(rows) {
  columns <- unique(unlist(lapply(rows, names), use.names = FALSE))
  rows <- lapply(rows, function(row) {
    for (name in setdiff(columns, names(row))) row[[name]] <- rep(NA, nrow(row))
    row[columns]
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

sb_view_context <- function() {
  list(repeat_id = NA_integer_, heldout_cohort_id = NA_character_,
       model_role = "pipeline", validation_scheme = "development_cv",
       view_id = NA_character_)
}

sb_view_leaves <- function(x, context = sb_view_context(), path = "result") {
  if (inherits(x, "stormbreaker_repeated_cv_comparison")) {
    if (!length(x$repetitions)) sb_abort("No saved repetitions.", "stormbreaker_bad_result_view")
    return(unlist(lapply(seq_along(x$repetitions), function(i) {
      ctx <- context
      ctx$repeat_id <- as.integer(i)
      ctx$validation_scheme <- "repeated_cv"
      sb_view_leaves(x$repetitions[[i]], ctx, paste0(path, "/repeat_", i))
    }), recursive = FALSE))
  }
  if (inherits(x, "stormbreaker_internal_external_cv")) {
    if (!length(x$folds)) sb_abort("No saved IECV folds.", "stormbreaker_bad_result_view")
    return(unlist(lapply(seq_along(x$folds), function(i) {
      fold <- x$folds[[i]]
      ctx <- context
      ctx$heldout_cohort_id <- fold$heldout_cohort_id
      ctx$validation_scheme <- "internal_external_cv"
      if (identical(x$strategy_type, "stack")) ctx$model_role <- "stack"
      sb_view_leaves(fold$comparison, ctx, paste0(path, "/cohort_", i))
    }), recursive = FALSE))
  }
  if (inherits(x, "stormbreaker_stack_cv_result")) {
    ctx <- context
    ctx$model_role <- "member"
    members <- sb_view_leaves(x$members, ctx, paste0(path, "/members"))
    # This local adapter only exposes the already saved stack outputs to readers.
    stack <- x[c("status", "task_type", "predictions", "metrics", "runs", "task")]
    stack$pipeline <- list(id = x$stack$id)
    # Stack descriptors are compact; retain unambiguous metadata saved per row.
    for (field in c("prediction_horizon", "time_unit")) {
      saved <- unique(stack$predictions[[field]])
      saved <- saved[!is.na(saved)]
      if (is.null(stack$task[[field]]) && length(saved) == 1L) {
        stack$task[[field]] <- saved[[1L]]
      }
    }
    for (component in c("predictions", "runs")) {
      if (!"pipeline_id" %in% names(stack[[component]])) {
        stack[[component]]$pipeline_id <- rep(x$stack$id, nrow(stack[[component]]))
      }
    }
    class(stack) <- if (identical(x$task_type, "binary")) {
      "stormbreaker_binary_cv_result"
    } else "stormbreaker_cv_result"
    ctx$model_role <- "stack"
    return(c(members, sb_view_leaves(stack, ctx, paste0(path, "/stack"))))
  }
  if (inherits(x, c("stormbreaker_cv_comparison", "stormbreaker_external_comparison"))) {
    values <- if (inherits(x, "stormbreaker_cv_comparison")) x$results else x$evaluations
    if (!length(values)) sb_abort("No saved comparison results.", "stormbreaker_bad_result_view")
    return(unlist(lapply(seq_along(values), function(i) {
      sb_view_leaves(values[[i]], context, paste0(path, "/pipeline_", i))
    }), recursive = FALSE))
  }
  if (!inherits(x, c("stormbreaker_binary_cv_result", "stormbreaker_cv_result",
                    "stormbreaker_external_evaluation"))) {
    sb_abort("Unsupported saved result in composite object.", "stormbreaker_bad_result_view")
  }
  if (inherits(x, "stormbreaker_external_evaluation") &&
      identical(context$validation_scheme, "development_cv")) {
    context$validation_scheme <- "external_evaluation"
  }
  context$view_id <- path
  list(list(value = x, context = context))
}

sb_view_annotate <- function(tab, context) {
  for (name in names(context)) tab[[name]] <- rep(context[[name]], nrow(tab))
  tab
}

sb_view_boundary <- function(context) {
  switch(context$validation_scheme,
    internal_external_cv = "Development IECV (not final independent external validation)",
    repeated_cv = "Development repeated CV", external_evaluation = "External evaluation",
    "Development CV")
}

sb_view_label <- function(context) {
  paste(c(if (!is.na(context$repeat_id)) paste("Repeat", context$repeat_id),
    if (!is.na(context$heldout_cohort_id)) paste("Held-out cohort", context$heldout_cohort_id),
    paste("Role:", context$model_role)), collapse = " | ")
}

sb_view_wrap_label <- function(labels) {
  vapply(strsplit(as.character(labels), "\n", fixed = TRUE), function(lines) {
    paste(unlist(lapply(lines, function(line) {
      if (!nzchar(line)) return("")
      starts <- seq.int(1L, nchar(line), by = 50L)
      substring(line, starts, pmin(starts + 49L, nchar(line)))
    }), use.names = FALSE), collapse = "\n")
  }, character(1), USE.NAMES = FALSE)
}

sb_composite_table <- function(x, component, metric, cohort_label) {
  context <- attr(x, "stormbreaker_view_context")
  attr(x, "stormbreaker_view_context") <- NULL
  leaves <- if (is.null(context)) sb_view_leaves(x) else list(list(value = x, context = context))
  sb_view_bind(lapply(leaves, function(leaf) {
    sb_view_annotate(sb_result_table(leaf$value, component, metric, cohort_label), leaf$context)
  }))
}

sb_composite_curves <- function(x, kind, cohort_label) {
  context <- attr(x, "stormbreaker_view_context")
  attr(x, "stormbreaker_view_context") <- NULL
  leaves <- if (is.null(context)) sb_view_leaves(x) else list(list(value = x, context = context))
  parts <- lapply(leaves, function(leaf) {
    view <- sb_curve_view_data(leaf$value, kind, cohort_label)
    for (name in c("data", "status")) {
      tab <- sb_view_annotate(view[[name]], leaf$context)
      tab$validation <- rep(leaf$context$validation_scheme, nrow(tab))
      target <- ifelse(tab$task_type == "binary", "class 1",
        ifelse(is.na(tab$prediction_horizon), "horizon not declared",
               paste(tab$prediction_horizon, tab$time_unit)))
      tab$panel <- paste(paste(tab$pipeline_id, target, sep = " | "),
        paste(c(cohort_label, sb_view_label(leaf$context)), collapse = " | "),
        sb_view_boundary(leaf$context),
        paste("curve:", tab$metric_status, "run:", tab$run_status),
        paste("Saved entry:", leaf$context$view_id), sep = "\n")
      view[[name]] <- tab
    }
    view
  })
  data <- sb_view_bind(lapply(parts, `[[`, "data"))
  status <- sb_view_bind(lapply(parts, `[[`, "status"))
  data$panel <- factor(data$panel, levels = unique(status$panel))
  data$line_group <- interaction(data$view_id, data$strategy, data$segment, drop = TRUE)
  list(data = data, status = status)
}

sb_composite_failures <- function(x, label) {
  context <- attr(x, "stormbreaker_view_context")
  attr(x, "stormbreaker_view_context") <- NULL
  leaves <- if (is.null(context)) sb_view_leaves(x) else list(list(value = x, context = context))
  sb_view_bind(lapply(leaves, function(leaf) {
    tab <- sb_view_annotate(sb_failure_rows(leaf$value, label), leaf$context)
    tab$validation <- rep(leaf$context$validation_scheme, nrow(tab))
    tab
  }))
}
