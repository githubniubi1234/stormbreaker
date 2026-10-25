#' Rank candidate genes by saved cross-validation selection counts
#'
#' @param x A saved CV result, CV comparison, repeated CV comparison, or a
#'   cohort comparison (only its development CV is used).
#' @param source Learned selector counts or sparse glmnet nonzero counts.
#' @param candidate_features Optional complete input gene list. Evaluated genes
#'   with no selections receive zero; genes absent from the saved candidate
#'   sets receive NA, not zero. Default: saved task predictors in input order.
#' @param cohort_label Optional name for this development cohort.
#' @return An object with `ranking`, `by_pipeline`, `groups` and `status`.
#'   Ranking columns include selected_count, eligible_count, selection_frequency
#'   and selected_methods. All supplied genes remain present. Tied counts share
#'   a rank and retain input order. No model is fitted or selected.
#' @details Each successful outer selection fold is one opportunity. Copies of
#'   the same selector configuration on the same saved split assignment are
#'   counted once, even if several downstream learners use them. Conflicting
#'   saved counts for the same configuration are marked unavailable. Missing
#'   split identity prevents deduplication and is recorded in groups.
#'   A failed learner does not erase a separately successful selector record.
#'   None/fixed feature panels are not learned selections. Glmnet active sets
#'   require alpha > 0 and remain separate from selector counts. Repeated CV
#'   repetitions remain distinct counting units. External outcomes are not read.
#'   One call covers one development cohort and task; do not pool cohorts.
#' @export
sb_gene_selection_ranking <- function(x, source = c("selector", "glmnet_active_set"),
                                      candidate_features = NULL, cohort_label = NULL) {
  source <- match.arg(source)
  if (inherits(x, "stormbreaker_cohort_comparison")) x <- x$development
  supported <- c(
    "stormbreaker_binary_cv_result", "stormbreaker_cv_result",
    "stormbreaker_cv_comparison", "stormbreaker_repeated_cv_comparison"
  )
  if (!inherits(x, supported)) {
    sb_abort("Supply one development CV result or comparison.", "stormbreaker_bad_gene_ranking")
  }
  leaves <- sb_viz_leaves(x, cohort_label)
  tasks <- unique(vapply(leaves, function(leaf) {
    if (inherits(leaf$value, "stormbreaker_binary_cv_result")) "binary" else "survival"
  }, character(1)))
  if (length(tasks) != 1L) sb_abort("Rank binary and survival tasks separately.", "stormbreaker_bad_gene_ranking")
  valid_genes <- function(z) {
    is.character(z) && length(z) > 0L &&
      !anyNA(z) && all(nzchar(trimws(z))) && !anyDuplicated(z)
  }
  saved_candidates <- unique(unlist(lapply(leaves, function(leaf) leaf$value$task$predictors), use.names = FALSE))
  if (is.null(candidate_features)) candidate_features <- saved_candidates
  if (!valid_genes(candidate_features)) sb_abort("Supply distinct non-empty candidate gene names.", "stormbreaker_bad_gene_ranking")
  groups <- list()
  pipeline_rows <- list()
  for (i in seq_along(leaves)) {
    leaf <- leaves[[i]]
    result <- leaf$value
    pipeline <- result$pipeline
    summary <- if (source == "selector") result$selection_stability else result$glmnet_active_set_stability
    applicable <- if (source == "selector") {
      !is.null(pipeline$selector) && !pipeline$selector %in% c("none", "fixed_features")
    } else {
      pipeline$learner %in% c("glmnet_cox", "glmnet_logistic") &&
        isTRUE(pipeline$learner_params$alpha > 0)
    }
    status <- if (!applicable) "not_applicable" else if (is.null(summary)) "not_saved" else summary$status
    tab <- if (applicable) summary$feature_frequency else NULL
    n <- if (applicable) summary$n_successful_folds else 0L
    valid <- is.data.frame(tab) && all(c("feature", "selected_count") %in% names(tab)) &&
      valid_genes(tab$feature) && is.numeric(n) && length(n) == 1L && is.finite(n) && n >= 0 && n == floor(n) &&
      is.numeric(tab$selected_count) && all(is.finite(tab$selected_count)) &&
      all(tab$selected_count >= 0 & tab$selected_count <= n & tab$selected_count == floor(tab$selected_count))
    if (applicable && !is.null(tab) && !valid) {
      sb_abort("Saved selection counts are malformed.", "stormbreaker_bad_gene_ranking")
    }
    usable <- applicable && valid && n > 0
    group_id <- NA_integer_
    duplicate <- FALSE
    if (usable) {
      tab <- tab[order(tab$feature), c("feature", "selected_count"), drop = FALSE]
      rownames(tab) <- NULL
      assignment <- result$provenance$fold_assignment
      if (is.data.frame(assignment) && all(c("patient_id", "fold") %in% names(assignment))) {
        assignment <- assignment[order(assignment$patient_id), c("patient_id", "fold"), drop = FALSE]
        rownames(assignment) <- NULL
      } else {
        assignment <- result$folds
      }
      identity_known <- !is.null(assignment) && length(assignment) > 0L
      key <- list(
        source = source, selector = pipeline$selector,
        selector_params = pipeline$selector_params,
        learner = if (source == "glmnet_active_set") pipeline$learner else NULL,
        learner_params = if (source == "glmnet_active_set") pipeline$learner_params else NULL,
        candidates = sort(result$task$predictors), assignment = assignment,
        repeat_id = leaf$context$repeat_id,
        separate_entry = if (!identity_known) i else NULL
      )
      matches <- which(vapply(groups, function(group) identical(group$key, key), logical(1)))
      counts <- list(
        table = tab, denominator = as.integer(n),
        per_fold = summary$per_fold[intersect(c("fold", "status", "n_selected"), names(summary$per_fold))]
      )
      if (length(matches)) {
        group_id <- matches[1]
        duplicate <- TRUE
        groups[[group_id]]$pipeline_ids <- c(groups[[group_id]]$pipeline_ids, pipeline$id)
        if (!identical(groups[[group_id]]$counts, counts)) groups[[group_id]]$status <- "conflicting_counts"
      } else {
        group_id <- length(groups) + 1L
        groups[[group_id]] <- list(
          key = key, counts = counts, status = "success",
          selector = pipeline$selector, learner = pipeline$learner,
          pipeline_ids = pipeline$id, split_identity_saved = identity_known
        )
      }
    }
    j <- if (valid) match(candidate_features, tab$feature) else rep(NA_integer_, length(candidate_features))
    pipeline_rows[[i]] <- data.frame(
      pipeline_id = pipeline$id, source = source,
      feature = candidate_features, selected_count = if (usable) tab$selected_count[j] else NA_integer_,
      eligible_count = if (usable) ifelse(is.na(j), NA_integer_, as.integer(n)) else NA_integer_,
      status = if (is.null(status)) "not_saved" else status,
      group_id = group_id, duplicate_copy = duplicate, stringsAsFactors = FALSE
    )
  }
  successful_groups <- Filter(function(group) identical(group$status, "success"), groups)
  selected <- eligible <- methods <- integer(length(candidate_features))
  for (group in successful_groups) {
    j <- match(candidate_features, group$counts$table$feature)
    present <- !is.na(j)
    count <- group$counts$table$selected_count[j[present]]
    selected[present] <- selected[present] + count
    eligible[present] <- eligible[present] + group$counts$denominator
    methods[present] <- methods[present] + as.integer(count > 0)
  }
  selected[eligible == 0] <- NA_integer_
  methods[eligible == 0] <- NA_integer_
  ranking <- data.frame(
    rank = as.integer(rank(-selected, ties.method = "min", na.last = "keep")),
    feature = candidate_features, selected_count = selected, eligible_count = eligible,
    selection_frequency = ifelse(eligible > 0, selected / eligible, NA_real_),
    selected_methods = methods, stringsAsFactors = FALSE
  )
  ranking <- ranking[order(-ranking$selected_count, seq_len(nrow(ranking)), na.last = TRUE), , drop = FALSE]
  rownames(ranking) <- NULL
  group_table <- if (length(groups)) {
    sb_view_bind(lapply(seq_along(groups), function(i) {
      group <- groups[[i]]
      data.frame(
        group_id = i, source = source, selector = group$selector,
        learner = if (source == "glmnet_active_set") group$learner else NA_character_,
        status = group$status, successful_folds = group$counts$denominator,
        pipeline_copies = length(group$pipeline_ids),
        split_identity_saved = group$split_identity_saved,
        pipeline_ids = paste(unique(group$pipeline_ids), collapse = ";"), stringsAsFactors = FALSE
      )
    }))
  } else {
    data.frame(
      group_id = integer(), source = character(), selector = character(),
      learner = character(), status = character(), successful_folds = integer(),
      pipeline_copies = integer(), split_identity_saved = logical(), pipeline_ids = character()
    )
  }
  by_pipeline <- sb_view_bind(pipeline_rows)
  status <- if (!length(successful_groups)) {
    "not_available"
  } else if (
    any(group_table$status != "success") || any(by_pipeline$status %in% c("incomplete", "not_saved", "failed"))) {
    "incomplete"
  } else {
    "success"
  }
  structure(list(
    status = status, source = source, task_type = tasks,
    cohort = leaves[[1]]$dataset, counting_unit = "unique_configuration_successful_outer_folds",
    ranking = ranking, by_pipeline = by_pipeline, groups = group_table
  ), class = "stormbreaker_gene_selection_ranking")
}

#' Plot ranked gene selection counts
#' @param x A [sb_gene_selection_ranking()] object.
#' @param top_n Maximum genes to display; the ranking table retains all genes.
#' @return A ggplot with counts and successful-fold denominators.
#' @export
sb_plot_gene_selection_ranking <- function(x, top_n = 20L) {
  if (!inherits(x, "stormbreaker_gene_selection_ranking")) {
    sb_abort("Supply a gene-selection ranking.", "stormbreaker_bad_gene_ranking")
  }
  sb_require_curve_plotting()
  top_n <- sb_viz_positive_integer(top_n, "top_n")
  d <- utils::head(x$ranking[is.finite(x$ranking$selected_count), , drop = FALSE], top_n)
  if (!nrow(d)) {
    return(sb_empty_display_plot("Gene selection counts", "No saved learned-selection counts.", ""))
  }
  d$feature_axis <- factor(d$feature, levels = rev(d$feature))
  d$count_label <- paste(d$selected_count, "/", d$eligible_count)
  feature_axis <- selected_count <- count_label <- NULL
  ggplot2::ggplot(d, ggplot2::aes(x = selected_count, y = feature_axis)) +
    ggplot2::geom_col(fill = "#2878B5", width = 0.62) +
    ggplot2::geom_text(ggplot2::aes(label = count_label), hjust = -0.15, family = "Arial", size = 3.5) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.23)),
      breaks = function(limits) {
        b <- pretty(limits)
        b[b >= 0 & b == floor(b)]
      }
    ) +
    ggplot2::labs(
      title = if (x$source == "selector") "Gene selection counts" else "LASSO / Elastic Net nonzero counts",
      x = "Selected / eligible folds", y = NULL,
      subtitle = if (x$status == "incomplete") "Incomplete records" else NULL
    ) +
    sb_shap_theme()
}
