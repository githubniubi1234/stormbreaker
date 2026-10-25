sb_page_scalar <- function(value, default = NULL) {
  if (is.null(value) || length(value) != 1L || is.na(value)) default else unname(value)
}

sb_page_json_string <- function(value) {
  codes <- utf8ToInt(enc2utf8(as.character(value)))
  paste0('"', paste(vapply(codes, function(code) {
    if (code == 34L) return('\\"')
    if (code == 92L) return("\\\\")
    if (code < 32L || code %in% c(60L, 62L, 38L, 8232L, 8233L)) {
      return(sprintf("\\u%04x", code))
    }
    intToUtf8(code)
  }, character(1)), collapse = ""), '"')
}

sb_page_json <- function(value) {
  if (is.null(value) || (!is.list(value) && length(value) == 1L && is.na(value))) return("null")
  if (is.list(value)) {
    if (!is.null(names(value))) {
      pairs <- vapply(seq_along(value), function(i) {
        paste0(sb_page_json_string(names(value)[i]), ":", sb_page_json(value[[i]]))
      }, character(1))
      return(paste0("{", paste(pairs, collapse = ","), "}"))
    }
    return(paste0("[", paste(vapply(value, sb_page_json, character(1)), collapse = ","), "]"))
  }
  if (length(value) != 1L) return(sb_page_json(unname(as.list(value))))
  if (is.logical(value)) return(if (value) "true" else "false")
  if (is.numeric(value)) return(if (is.finite(value)) sprintf("%.17g", value) else "null")
  sb_page_json_string(value)
}

sb_page_model <- function(result, task_type) {
  external <- inherits(result, "stormbreaker_external_evaluation")
  pipeline <- if (external) NULL else result$pipeline
  id <- if (external) result$pipeline_id else pipeline$id
  tokens <- strsplit(id, "__", fixed = TRUE)[[1L]]
  selector <- sb_page_scalar(pipeline$selector, tokens[1L])
  learner <- sb_page_scalar(pipeline$learner,
    sb_page_scalar(result$learner, if (length(tokens) >= 2L) tokens[2L] else "Model"))
  learner_params <- pipeline$learner_params
  selector_params <- pipeline$selector_params
  if (is.null(learner_params)) learner_params <- list()
  if (is.null(selector_params)) selector_params <- list()
  if (is.null(learner_params$alpha) && grepl("__alpha", id, fixed = TRUE)) {
    token <- tokens[grepl("^alpha[0-9]", tokens)]
    if (length(token)) learner_params$alpha <- suppressWarnings(as.numeric(sub("^alpha", "", token[1])))
  }
  if (is.null(selector_params$n_features)) {
    k <- tokens[grepl("^k[0-9]+$", tokens)]
    if (length(k)) selector_params$n_features <- as.integer(sub("^k", "", k[1]))
  }
  list(id = id, task = task_type, selector_id = selector, learner_id = learner,
    selector = sb_display_selector(selector, selector_params),
    learner = sb_display_learner(learner, learner_params),
    parameters = sb_display_parameters(learner_params),
    selection_parameters = sb_display_parameters(selector_params))
}

sb_page_metric_cell <- function(value, metric, task_type, dataset, context) {
  external <- inherits(value, "stormbreaker_external_evaluation")
  id <- if (external) value$pipeline_id else value$pipeline$id
  tab <- if (metric %in% names(value$metrics) || metric %in% sb_horizon_view_metrics()) {
    sb_result_table(value, metric = metric)
  } else {
    data.frame(pipeline_id = id, task_type = task_type, metric = metric,
      estimate = NA_real_, status = "not_saved", run_status = value$status,
      scope = NA_character_, aggregation = NA_character_)
  }
  if (nrow(tab) != 1L) sb_abort("A saved page cell must contain one estimate.",
    "stormbreaker_ambiguous_result_page")
  available <- identical(tab$status[1], "success") &&
    identical(tab$run_status[1], "success") && is.finite(tab$estimate[1])
  descriptor <- if (external) value$endpoint else value$task
  n <- sb_page_scalar(if (external) value$n else descriptor$n, NA_real_)
  list(model = id, task = task_type, dataset = dataset, metric = metric,
    estimate = sb_page_scalar(tab$estimate, NA_real_), available = available,
    metric_status = sb_page_scalar(tab$status, "not_saved"),
    run_status = sb_page_scalar(tab$run_status, "not_saved"),
    scope = sb_page_scalar(tab$scope, ""), aggregation = sb_page_scalar(tab$aggregation, ""),
    horizon = sb_page_scalar(tab$prediction_horizon, NA_real_),
    time_unit = sb_page_scalar(tab$time_unit, ""),
    n = n, context = context$validation_scheme,
    error_stage = sb_page_scalar(value$error_stage, ""),
    error_class = sb_page_scalar(value$error_class, ""))
}

sb_result_page_data <- function(x, roles = NULL, title = "Model comparison") {
  if (inherits(x, "stormbreaker_cohort_comparison")) {
    if (is.null(roles)) roles <- x$cohort_roles
    x <- x$results
  }
  if (inherits(x, "stormbreaker_external_collection")) x <- x$comparisons
  if (!is.list(x) || is.object(x)) {
    label <- if (inherits(x, c("stormbreaker_external_evaluation",
      "stormbreaker_external_comparison"))) "External" else "Development"
    x <- stats::setNames(list(x), label)
  }
  labels <- names(x)
  if (!length(x) || is.null(labels) || anyNA(labels) ||
      any(!nzchar(trimws(labels))) || anyDuplicated(labels)) {
    sb_abort("Use a named list with distinct dataset names.", "stormbreaker_bad_result_page")
  }
  if (!is.null(roles) && (!is.character(roles) || anyNA(roles) ||
      is.null(names(roles)) || anyDuplicated(names(roles)) ||
      !setequal(names(roles), labels) ||
      !all(roles %in% c("development_cv", "training", "validation", "test", "external")))) {
    sb_abort("roles must name each dataset as development_cv, training, validation, test or external.",
      "stormbreaker_bad_cohort_roles")
  }
  leaves <- sb_viz_leaves(x)
  models <- datasets <- cells <- list()
  for (leaf in leaves) {
    value <- leaf$value
    external <- inherits(value, "stormbreaker_external_evaluation")
    task_type <- if (inherits(value, "stormbreaker_binary_cv_result")) "binary" else
      if (external) value$task_type else "survival"
    ctx <- leaf$context
    suffix <- c(if (!is.na(ctx$repeat_id)) paste("CV", ctx$repeat_id),
      if (!is.na(ctx$heldout_cohort_id)) ctx$heldout_cohort_id)
    dataset <- paste(c(leaf$dataset, suffix), collapse = " / ")
    role <- if (!is.null(roles)) unname(roles[[leaf$dataset]]) else
      if (external && ctx$validation_scheme == "external_evaluation") "external" else "development_cv"
    if ((!external && role != "development_cv") ||
        (external && ctx$validation_scheme == "external_evaluation" && role == "development_cv")) {
      sb_abort("Dataset roles must agree with saved CV or evaluation results.",
        "stormbreaker_bad_cohort_roles")
    }
    if (ctx$validation_scheme == "internal_external_cv") role <- "internal_external_cv"
    if (ctx$validation_scheme == "repeated_cv") role <- "repeated_cv"
    if (is.null(datasets[[dataset]])) datasets[[dataset]] <- list(
      id = dataset, label = dataset, role = role,
      context = ctx$validation_scheme)
    model <- sb_page_model(value, task_type)
    key <- paste(task_type, model$id, sep = "::")
    if (is.null(models[[key]]) || !external) models[[key]] <- model
    metrics <- if (task_type == "binary") {
      c("roc_auc", "average_precision", "brier_score", "log_loss")
    } else c("harrell_c", "time_auc", "uno_c", "ipcw_brier")
    for (metric in metrics) cells[[length(cells) + 1L]] <-
      sb_page_metric_cell(value, metric, task_type, dataset, ctx)
  }
  keys <- vapply(cells, function(cell) {
    paste(cell$task, cell$model, cell$dataset, cell$metric, sep = "\r")
  }, character(1))
  if (anyDuplicated(keys)) {
    sb_abort("More than one saved result occupies a model/dataset cell; use distinct dataset labels.",
      "stormbreaker_ambiguous_result_page")
  }
  development <- vapply(datasets, function(dataset) dataset$role %in%
    c("development_cv", "repeated_cv", "internal_external_cv"), logical(1))
  default_sort <- if (any(development)) names(datasets)[which(development)[1]] else ""
  list(title = title, version = "0.3.6", datasets = unname(datasets),
    models = unname(models), cells = cells, default_sort = default_sort)
}

sb_check_page_request <- function(view, path) {
  if (!is.logical(view) || length(view) != 1L || is.na(view)) {
    sb_abort("view/open must be TRUE or FALSE.", "stormbreaker_bad_result_page")
  }
  if (!is.null(path) && (!is.character(path) || length(path) != 1L ||
      is.na(path) || !nzchar(trimws(path)))) {
    sb_abort("Supply one HTML file path.", "stormbreaker_bad_result_page")
  }
  invisible(NULL)
}

sb_open_result_page <- function(path) {
  utils::browseURL(path)
}

#' Open a model-by-cohort results matrix
#'
#' Create a standalone HTML page from saved results. Rows show feature selectors
#' and learner names; columns show development, validation and test cohorts.
#' Metrics are read from the supplied results without refitting models.
#'
#' @param x A saved CV/evaluation result, comparison, external collection,
#'   [sb_compare_cohorts()] result, or uniquely named list of saved results.
#' @param roles Optional named dataset roles: `"development_cv"`, `"training"`,
#'   `"validation"`, `"test"` or `"external"`. CV objects retain their CV role.
#'   Unlabeled external evaluations remain external, not assumed validation.
#' @param path New HTML file path. Defaults to a temporary report.
#' @param open Open the page in the configured browser; defaults to `interactive()`.
#' @param title Page title.
#' @param overwrite Replace an existing report only when explicitly TRUE.
#' @return Invisibly, the normalized path of the saved HTML file.
#' @details The page contains summary metrics and model metadata, not patient
#'   rows or raw error messages. Original model IDs and parameter settings
#'   remain available in details. Default sorting uses development CV when
#'   present; it never selects or refits a winning model. Optional validation
#'   means exclude training/test cohorts and require all selected validation
#'   cells to be available for the same metric and horizon.
#' @export
#' @examples
#' d <- sb_example_binary(50, seed = 1)
#' task <- sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b"))
#' result <- sb_compare_models(task, models = "logistic", folds = 3, view = FALSE)
#' page <- sb_view_results(result, open = FALSE)
#' file.exists(page)
sb_view_results <- function(x, roles = NULL, path = tempfile(fileext = ".html"),
                            open = interactive(), title = "Model comparison",
                            overwrite = FALSE) {
  sb_check_page_request(open, path)
  if (!is.logical(overwrite) || length(overwrite) != 1L || is.na(overwrite) ||
      !is.character(title) || length(title) != 1L || is.na(title)) {
    sb_abort("Supply one title and a TRUE/FALSE overwrite setting.", "stormbreaker_bad_result_page")
  }
  if (!dir.exists(dirname(path))) sb_abort("The report folder does not exist.",
    "stormbreaker_bad_result_page")
  if (file.exists(path) && !overwrite) sb_abort("The report exists; use a new path or overwrite = TRUE.",
    "stormbreaker_result_page_exists")
  data <- sb_result_page_data(x, roles, title)
  template <- system.file("leaderboard", "result-page-template.html", package = "stormbreaker")
  if (!nzchar(template)) sb_abort("The results-page template is missing.", "stormbreaker_bad_result_page")
  text <- paste(readLines(template, encoding = "UTF-8", warn = FALSE), collapse = "\n")
  parts <- strsplit(text, "/*SB_REPORT_DATA*/null", fixed = TRUE)[[1L]]
  if (length(parts) != 2L) sb_abort("The results-page template is invalid.", "stormbreaker_bad_result_page")
  html <- paste0(parts[1L], sb_page_json(data), parts[2L])
  writeLines(enc2utf8(html), path, useBytes = TRUE)
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  if (open) {
    tryCatch(sb_open_result_page(path), error = function(error) {
      warning(paste0("The page was saved, but the browser could not open it: ", path),
        call. = FALSE)
    })
  }
  invisible(path)
}

sb_attach_result_page <- function(result, view, report_file) {
  if (!view && is.null(report_file)) return(result)
  if (is.null(report_file)) report_file <- tempfile("stormbreaker-results-", fileext = ".html")
  attempt <- tryCatch(
    list(path = sb_view_results(result, path = report_file, open = view), error = NULL),
    error = function(error) list(path = NULL, error = conditionMessage(error))
  )
  if (is.null(attempt$error)) {
    result$report_file <- attempt$path
  } else {
    result$report_error <- attempt$error
    warning(paste("The comparison is complete; the results page could not be saved:", attempt$error),
      call. = FALSE)
  }
  result
}
