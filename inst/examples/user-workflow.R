# Source this file to define helpers; sourcing does not fit or select a model.
# All generated templates/demo results are SYNTHETIC, not clinical evidence.
# Requires an installed stormbreaker package. See guides/workflow.md.

workflow_new_directory <- function(path) {
  # Output files can contain patient data; never replace an existing run.
  if (file.exists(path)) stop("Choose a new output directory: ", path)
  if (!dir.create(path, recursive = TRUE)) stop("Cannot create: ", path)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

workflow_read_csv <- function(path) {
  # Preserve IDs such as 00123 and original feature column names.
  header <- names(utils::read.csv(path, nrows = 0, check.names = FALSE))
  classes <- rep(NA_character_, length(header))
  classes[header == "patient_id"] <- "character"
  utils::read.csv(path,
    colClasses = classes, check.names = FALSE,
    stringsAsFactors = FALSE, fileEncoding = "UTF-8"
  )
}

workflow_make_templates <- function(directory) {
  directory <- workflow_new_directory(directory)
  for (kind in c("binary", "survival")) {
    generator <- if (kind == "binary") {
      stormbreaker::bi_example
    } else {
      stormbreaker::surv_example
    }
    development <- generator(n = 100, seed = 20260819)
    external <- generator(n = 40, seed = 20260821)
    development$patient_id <- sprintf("%s_dev_%03d", kind, seq_len(100))
    external$patient_id <- sprintf("%s_ext_%03d", kind, seq_len(40))
    outcomes <- if (kind == "binary") "outcome" else c("time", "event")
    tables <- list(
      development = development,
      external_predictors = external[c("patient_id", "gene_a", "gene_b", "gene_c")],
      external_outcomes = external[c("patient_id", outcomes)]
    )
    for (role in names(tables)) {
      utils::write.csv(tables[[role]],
        file.path(directory, paste0(kind, "_", role, ".csv")),
        row.names = FALSE
      )
    }
  }
  invisible(directory)
}

workflow_develop <- function(development_csv, output_dir,
                             kind = c("binary", "survival"),
                             predictors = c("gene_a", "gene_b", "gene_c"),
                             time_unit = "month", prediction_horizon = 6,
                             folds = 5L, seed = 5201314L,
                             near_tie = 0.01) {
  kind <- match.arg(kind)
  if (!is.numeric(near_tie) || length(near_tie) != 1L ||
    !is.finite(near_tie) || near_tie < 0) {
    stop("near_tie must be nonnegative.")
  }
  directory <- workflow_new_directory(output_dir)
  data <- workflow_read_csv(development_csv)
  task <- if (kind == "binary") {
    stormbreaker::bi_task(
      data,
      "patient_id", "outcome", predictors
    )
  } else {
    stormbreaker::surv_task(data,
      "patient_id", "time", "event", predictors, time_unit,
      prediction_horizon = prediction_horizon
    )
  }
  # Declare BOTH candidates and the decision rule before running CV.
  # Here 'simpler' means the predeclared modelling procedure, not fitted size.
  pipelines <- if (kind == "binary") {
    list(
      logistic = stormbreaker::sb_pipeline(learner = "logistic"),
      naive_bayes = stormbreaker::sb_pipeline(learner = "gaussian_naive_bayes")
    )
  } else {
    list(
      all_features_cox = stormbreaker::sb_pipeline(learner = "coxph"),
      top2_cox = stormbreaker::sb_pipeline(
        learner = "coxph",
        selector = "univ_cox_top_k", selector_n_features = 2L
      )
    )
  }
  metric <- if (kind == "binary") "roc_auc" else "harrell_c"
  rule <- list(
    primary_metric = metric, direction = "higher",
    near_tie = near_tie, simplicity_order = names(pipelines),
    stability_review = paste(
      "Review fold failures, selection stability,",
      "calibration and secondary losses before confirming.",
      "No universal Jaccard cutoff; not_applicable is not failure."
    )
  )
  saveRDS(rule, file.path(directory, "selection_rule_before_cv.rds"))
  compared <- stormbreaker::compare_cv(task, pipelines, folds, seed)
  tab <- compared$comparison
  eligible <- tab$status == "success" & tab$complete_oof &
    tab[[paste0(metric, "_status")]] == "success" & is.finite(tab[[metric]])
  eligible[is.na(eligible)] <- FALSE
  candidate_ids <- vapply(pipelines, function(x) x$id, character(1))
  tab$candidate_name <- names(candidate_ids)[match(tab$pipeline_id, candidate_ids)]
  tab$eligible_for_primary_rule <- eligible
  tab$near_best <- FALSE
  suggested <- NA_character_
  if (any(eligible)) {
    tab$near_best <- eligible & tab[[metric]] >= max(tab[[metric]][eligible]) - near_tie
    tab$near_best[is.na(tab$near_best)] <- FALSE
    suggested <- names(pipelines)[which(names(pipelines) %in%
      tab$candidate_name[tab$near_best])[1L]]
  }
  utils::write.csv(tab, file.path(directory, "candidate_comparison.csv"), row.names = FALSE)
  utils::write.csv(compared$fold_assignment, file.path(directory, "shared_folds.csv"),
    row.names = FALSE
  )
  for (source in c("selector", "glmnet_active_set")) {
    ranking <- stormbreaker::gene_rank(compared, source,
      candidate_features = predictors
    )
    stem <- if (source == "selector") "gene_selection" else "gene_nonzero"
    saveRDS(ranking, file.path(directory, paste0(stem, "_ranking.rds")))
    utils::write.csv(ranking$ranking, file.path(directory, paste0(stem, "_ranking.csv")), row.names = FALSE)
    utils::write.csv(ranking$by_pipeline, file.path(directory, paste0(stem, "_by_pipeline.csv")), row.names = FALSE)
    utils::write.csv(ranking$groups, file.path(directory, paste0(stem, "_counting_groups.csv")), row.names = FALSE)
    if (ranking$status != "not_available" && requireNamespace("ggplot2", quietly = TRUE)) {
      ggplot2::ggsave(file.path(directory, paste0(stem, "_ranking.png")),
        stormbreaker::plot_gene_rank(ranking),
        width = 9, height = 7, dpi = 300
      )
    }
  }
  state <- list(
    kind = kind, task = task, pipelines = pipelines, cv = compared,
    comparison = tab, rule = rule, suggested_candidate = suggested,
    predictors = predictors, time_unit = time_unit,
    prediction_horizon = prediction_horizon
  )
  saveRDS(state, file.path(directory, "development.rds"))
  message(
    "Development only. Review candidate_comparison.csv. Suggested for review: ",
    suggested, ". Nothing has been selected or fitted as a final model."
  )
  invisible(state)
}

workflow_finalize <- function(development_rds, selected_candidate, selection_reason,
                              external_predictors_csv, external_outcomes_csv,
                              output_dir, render_report = TRUE) {
  if (missing(selected_candidate) || missing(selection_reason) ||
    !is.character(selection_reason) || length(selection_reason) != 1L ||
    is.na(selection_reason) || !nzchar(trimws(selection_reason))) {
    stop("Explicit selected_candidate and selection_reason are required.")
  }
  state <- readRDS(development_rds)
  if (!is.character(selected_candidate) || length(selected_candidate) != 1L ||
    is.na(selected_candidate) || !selected_candidate %in% names(state$pipelines)) {
    stop("selected_candidate must name a declared candidate.")
  }
  row <- state$comparison[state$comparison$candidate_name == selected_candidate, ]
  if (!isTRUE(row$eligible_for_primary_rule)) stop("Selected candidate did not complete primary evaluation.")
  directory <- workflow_new_directory(output_dir)
  decision <- list(
    selected_candidate = selected_candidate, reason = selection_reason,
    suggested_candidate = state$suggested_candidate, rule = state$rule,
    chosen_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
    package_version = as.character(utils::packageVersion("stormbreaker"))
  )
  saveRDS(decision, file.path(directory, "selection_decision.rds"))
  fit <- stormbreaker::fit_final(state$task, state$pipelines[[selected_candidate]])
  if (!identical(fit$status, "success")) stop("Final fit failed; no replacement model.")
  saveRDS(fit, file.path(directory, "final_model.rds"))
  # Reload the artifact to demonstrate that prediction does not require refitting.
  fit <- readRDS(file.path(directory, "final_model.rds"))
  new_data <- workflow_read_csv(external_predictors_csv)
  # Deliberately supply only the identifier and the declared predictors.
  new_data <- new_data[c("patient_id", state$predictors)]
  prediction <- stormbreaker::predict_final(fit, new_data, id_col = "patient_id")
  saveRDS(prediction, file.path(directory, "external_predictions.rds"))
  utils::write.csv(prediction$predictions,
    file.path(directory, "external_predictions.csv"),
    row.names = FALSE
  )
  # External outcomes enter for the first time AFTER the model and predictions exist.
  outcomes <- workflow_read_csv(external_outcomes_csv)
  endpoint <- if (state$kind == "binary") {
    stormbreaker::bi_endpoint(
      outcomes, "patient_id", "outcome"
    )
  } else {
    stormbreaker::surv_endpoint(
      outcomes, "patient_id", "time", "event", state$time_unit,
      prediction_horizon = state$prediction_horizon
    )
  }
  evaluation <- stormbreaker::eval_external(prediction, endpoint)
  saveRDS(evaluation, file.path(directory, "external_evaluation.rds"))
  capture.output(print(evaluation), file = file.path(directory, "evaluation_summary.txt"))
  if (isTRUE(render_report)) {
    stormbreaker::sb_render_report(
      list(development_cv = state$cv, external_validation = evaluation),
      file = file.path(directory, "report.html")
    )
  }
  invisible(list(
    decision = decision, fit = fit, prediction = prediction,
    evaluation = evaluation, directory = directory
  ))
}

workflow_demo <- function(output_dir) {
  directory <- workflow_new_directory(output_dir)
  input <- file.path(directory, "synthetic_inputs")
  workflow_make_templates(input)
  results <- list()
  for (kind in c("binary", "survival")) {
    development_dir <- file.path(directory, paste0(kind, "_development"))
    state <- workflow_develop(file.path(input, paste0(kind, "_development.csv")),
      development_dir,
      kind = kind
    )
    # Synthetic demonstration ONLY: enact the predeclared suggestion.
    # Real projects stop after workflow_develop for human review and confirmation.
    results[[kind]] <- workflow_finalize(file.path(development_dir, "development.rds"),
      selected_candidate = state$suggested_candidate,
      selection_reason = "SYNTHETIC DEMO ONLY: enact predeclared near-tie/simplicity rule.",
      external_predictors_csv = file.path(input, paste0(kind, "_external_predictors.csv")),
      external_outcomes_csv = file.path(input, paste0(kind, "_external_outcomes.csv")),
      output_dir = file.path(directory, paste0(kind, "_final"))
    )
  }
  invisible(results)
}
