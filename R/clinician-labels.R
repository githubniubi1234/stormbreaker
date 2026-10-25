sb_clinician_models <- function(leaves) {
  rows <- lapply(leaves, function(leaf) {
    m <- sb_page_model(leaf$value, leaf$value$task_type)
    data.frame(id = m$id, model = m$learner, selector = m$selector,
      parameters = m$parameters, selection_parameters = m$selection_parameters)
  })
  models <- unique(do.call(rbind, rows))
  models <- models[!duplicated(models$id), , drop = FALSE]
  repeated <- duplicated(models$model) | duplicated(models$model, fromLast = TRUE)
  models$model[repeated] <- paste0(models$model[repeated], " (", seq_len(sum(repeated)), ")")
  models
}

sb_clinician_labeler <- function(models) {
  function(labels) vapply(as.character(labels), function(label) {
    for (i in order(nchar(models$id), decreasing = TRUE)) {
      label <- gsub(models$id[i], models$model[i], label, fixed = TRUE)
    }
    label <- gsub("Logistic regression", "Logistic", label, fixed = TRUE)
    label <- gsub("Gaussian naive Bayes", "Gaussian NB", label, fixed = TRUE)
    lines <- trimws(strsplit(label, "\n", fixed = TRUE)[[1L]])
    lines <- lines[!grepl("^(Unnamed cohort|Role: pipeline|Saved entry:|curve: success run: success)$", lines) &
      !grepl("^Saved entry:", lines)]
    lines <- gsub("Role: pipeline", "", lines, fixed = TRUE)
    lines <- gsub("Assessment fold", "Fold", lines, fixed = TRUE)
    lines <- gsub("Development CV", "CV", lines, fixed = TRUE)
    lines <- gsub("Cross-validation", "CV", lines, fixed = TRUE)
    lines <- gsub(" | class 1", "", lines, fixed = TRUE)
    lines <- gsub("External patients", "", lines, fixed = TRUE)
    lines <- gsub("curve: success", "", lines, fixed = TRUE)
    lines <- gsub("run: success", "", lines, fixed = TRUE)
    lines <- gsub("(^[ |]+|[ |]+$)", "", lines)
    if (any(lines == "CV")) lines <- lines[lines != "Out-of-fold"]
    paste(unique(lines[nzchar(lines)]), collapse = "\n")
  }, character(1), USE.NAMES = FALSE)
}

sb_clinician_clean_plot <- function(plot, models) {
  clean <- sb_clinician_labeler(models)
  if (length(intersect(all.vars(plot$mapping$x), c("pipeline_id", "display_label")))) {
    plot$scales <- plot$scales$clone()
    axis <- plot$scales$get_scales("x")
    if (!is.null(axis)) axis$labels <- clean
  }
  original <- plot$facet$params$labeller
  if (!is.null(original)) {
    facet <- plot$facet
    plot$facet <- ggplot2::ggproto(NULL, facet)
    plot$facet$params$labeller <- function(labels) {
      if ("panel_key" %in% names(labels)) return(lapply(original(labels), clean))
      lapply(labels, clean)
    }
  }
  if (identical(all.vars(plot$mapping$colour), "pipeline_id")) {
    if (is.null(plot$scales$get_scales("colour"))) {
      plot <- plot + ggplot2::scale_colour_discrete(name = NULL, labels = clean)
    } else {
      plot$scales <- plot$scales$clone()
      scale <- plot$scales$get_scales("colour")
      scale$labels <- clean
      scale$name <- NULL
    }
  }
  plot
}
