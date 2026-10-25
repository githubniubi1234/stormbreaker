sb_repeated_cv_specifications <- function(folds, seeds) {
  if (is.list(folds)) {
    if (length(folds) == 0L) {
      sb_abort(
        "`folds` must contain at least one patient-ID fold assignment.",
        "stormbreaker_bad_repeated_folds"
      )
    }
    if (!is.null(seeds)) {
      sb_abort(
        paste0(
          "`seeds` must be `NULL` when `folds` is a list of explicit ",
          "patient-ID assignments."
        ),
        "stormbreaker_repeated_fold_seed_conflict"
      )
    }
    valid_plan <- vapply(folds, function(plan) {
      plan_names <- names(plan)
      is.numeric(plan) && is.null(dim(plan)) && length(plan) > 1L &&
        !is.null(plan_names) && !anyNA(plan_names) &&
        all(nzchar(plan_names)) && anyDuplicated(plan_names) == 0L
    }, logical(1L))
    if (!all(valid_plan)) {
      sb_abort(
        paste0(
          "Every `folds` list element must be one multi-patient numeric ",
          "vector named by unique, non-empty patient IDs."
        ),
        "stormbreaker_bad_repeated_folds"
      )
    }
    return(list(
      folds = folds,
      seeds = rep(NA_integer_, length(folds))
    ))
  }

  if (length(folds) != 1L) {
    sb_abort(
      paste0(
        "For repeated generated folds, `folds` must be one fold count. ",
        "Wrap explicit patient-ID assignments in a list."
      ),
      "stormbreaker_bad_repeated_folds"
    )
  }
  valid_seeds <- is.numeric(seeds) && length(seeds) > 0L &&
    !anyNA(seeds) && all(is.finite(seeds)) && all(seeds > 0) &&
    all(seeds <= .Machine$integer.max) && all(seeds == floor(seeds))
  if (!valid_seeds) {
    sb_abort(
      "`seeds` must contain one or more finite positive whole numbers.",
      "stormbreaker_bad_repeat_seeds"
    )
  }

  list(
    folds = rep(list(folds), length(seeds)),
    seeds = as.integer(seeds)
  )
}

sb_run_repeated_cv_iteration <- function(task, pipelines, folds, seed) {
  if (inherits(pipelines, "stormbreaker_stack")) {
    return(sb_run_stack_cv(task, pipelines, folds = folds, seed = seed))
  }
  sb_compare_cv(task, pipelines, folds = folds, seed = seed)
}

sb_repeated_cv_result_metadata <- function(result) {
  if (inherits(result, "stormbreaker_stack_cv_result")) {
    return(list(
      folds = result$provenance$folds,
      seed = result$provenance$seed,
      fold_source = result$provenance$fold_source,
      fold_assignment = result$provenance$fold_assignment,
      n_patients = result$task$n
    ))
  }
  list(
    folds = result$folds,
    seed = result$seed,
    fold_source = result$fold_source,
    fold_assignment = result$fold_assignment,
    n_patients = result$n_patients
  )
}

sb_learner_uses_caller_rng <- function(id) {
  id %in% c("blackboost_cox", "rsf_survival", "adaboost_binary")
}

sb_parallel_cv_repetitions <- function(task, pipelines, specifications) {
  members <- if (inherits(pipelines, "stormbreaker_stack")) pipelines$pipelines else pipelines
  if (is.list(members) && any(vapply(members, function(p) {
    inherits(p, "stormbreaker_pipeline") && sb_learner_uses_caller_rng(p$learner)
  }, logical(1)))) {
    sb_abort(paste0("These pipelines include a learner that consumes the caller's ",
      "sequential RNG stream; use workers = 1 to preserve its results."),
      "stormbreaker_parallel_sequential_rng")
  }
  threads <- c(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1", BLIS_NUM_THREADS = "1")
  cluster <- withr::with_envvar(threads, parallel::makePSOCKcluster(
    min(2L, length(specifications$folds))))
  on.exit(parallel::stopCluster(cluster), add = TRUE)
  initialize <- function(paths, package_path, rng, kind) {
    .libPaths(paths)
    options(rf.cores = 1L, mc.cores = 1L,
      stormbreaker.ranger_cindex_threads = 1L)
    do.call(RNGkind, as.list(kind))
    if (is.null(rng)) {
      if (exists(".Random.seed", .GlobalEnv, inherits = FALSE))
        rm(".Random.seed", envir = .GlobalEnv)
    } else assign(".Random.seed", rng, .GlobalEnv)
    if (file.exists(file.path(package_path, "Meta", "package.rds"))) {
      loadNamespace("stormbreaker", lib.loc = dirname(package_path))
    } else {
      pkgload::load_all(package_path, quiet = TRUE, helpers = FALSE, attach = FALSE)
    }
    invisible(NULL)
  }
  environment(initialize) <- baseenv()
  rng <- if (exists(".Random.seed", .GlobalEnv, inherits = FALSE))
    get(".Random.seed", .GlobalEnv) else NULL
  parallel::clusterCall(cluster, initialize, .libPaths(),
    getNamespaceInfo(asNamespace("stormbreaker"), "path"), rng, RNGkind())
  run <- function(i, task, pipelines, specifications, rng) {
    if (!is.null(rng)) {
      assign(".Random.seed", rng, .GlobalEnv)
    } else if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
    warnings <- list()
    value <- tryCatch(withCallingHandlers({
      seed <- specifications$seeds[[i]]
      get("sb_run_repeated_cv_iteration", asNamespace("stormbreaker"))(
        task, pipelines, specifications$folds[[i]], if (is.na(seed)) 1L else seed)
    }, warning = function(w) {
      warnings[[length(warnings) + 1L]] <<- w
      invokeRestart("muffleWarning")
    }), error = function(e) e)
    list(value = value, warnings = warnings, failed = inherits(value, "error"))
  }
  environment(run) <- baseenv()
  outputs <- parallel::parLapplyLB(cluster, seq_along(specifications$folds), run,
    task = task, pipelines = pipelines, specifications = specifications, rng = rng)
  lapply(outputs, function(output) {
    for (condition in output$warnings) warning(condition)
    if (output$failed) stop(output$value)
    output$value
  })
}

#' Compare pipelines across repeated outer-fold assignments
#'
#' `sb_compare_repeated_cv()` repeats either [sb_compare_cv()] for an ordered
#' pipeline list or [sb_run_stack_cv()] for one predeclared stack, without
#' changing either function's fitting, prediction, metric, or failure behavior.
#' Within a pipeline comparison every pipeline uses the same outer folds; within
#' a stack repetition every member and the stack do so. Each patient's
#' out-of-fold predictions remain inside that repetition's complete result.
#'
#' There are two explicit resampling modes. When `folds` is one fold count,
#' each value in `seeds` generates one event- or outcome-stratified assignment.
#' When `folds` is a list, each element must be one numeric vector named by all
#' task patient IDs and `seeds` must be `NULL`. Explicit assignments are passed
#' unchanged to [sb_compare_cv()], which aligns them by patient ID and neither
#' repairs nor rebalances them.
#'
#' The long `comparison` table contains one row per repetition and pipeline.
#' It copies the corresponding single-comparison rows and does not pool patient
#' predictions across repetitions. Repeated splits show sensitivity to outer-
#' fold assignment; they are not independent cohorts or a confidence interval.
#' This function calculates no cross-repeat mean, standard deviation, interval,
#' ranking, winner, or final-model choice. Duplicate seeds or fold assignments
#' are preserved and remain visible rather than silently replaced.
#'
#' @param task A task returned by [sb_surv_task()] or [sb_binary_task()].
#' @param pipelines Either a non-empty list of distinct pipeline specifications
#'   returned by [sb_pipeline()] or one stack returned by [sb_stack()]. Input
#'   and stack-member order are preserved within every repetition. A stack
#'   repetition retains its member-plus-stack comparison; it does not select
#'   among those rows.
#' @param folds Either one positive whole-number outer-fold count or a non-empty
#'   list of explicit patient-ID fold assignments. Each explicit assignment
#'   follows the `folds` vector contract of [sb_compare_cv()].
#' @param seeds For a scalar fold count, one or more positive whole-number seeds,
#'   one per repetition. Use `NULL` with a list of explicit fold assignments.
#' @param workers Either `1L` (the unchanged serial default) or `2L`. Two local
#'   processes execute independent repetitions, preserving within-repetition
#'   sharing and input order. Engines remain single-threaded; memory use increases.
#'   This does not parallelize pipelines within one repetition. Worker warnings
#'   are relayed in repetition order and failures are not replaced or retried.
#'   Native BlackBoost-Cox, RSF and AdaBoost learners consume a sequential caller
#'   RNG stream and currently require `workers = 1`; they are not silently reseeded.
#'
#' @return A `stormbreaker_repeated_cv_comparison` containing a long
#'   repeat-by-pipeline `comparison` table, every complete
#'   single-run comparison or stack result in `repetitions`, and a long
#'   `fold_assignments` table keyed by repetition and patient ID.
#' @export
#'
#' @examples
#' d <- sb_example_binary(n = 48, seed = 1)
#' task <- sb_binary_task(
#'   d, "patient_id", "outcome", c("gene_a", "gene_b", "gene_c")
#' )
#' repeated <- sb_compare_repeated_cv(
#'   task,
#'   list(
#'     sb_pipeline(learner = "logistic"),
#'     sb_pipeline(learner = "gaussian_naive_bayes")
#'   ),
#'   folds = 3,
#'   seeds = c(11, 12)
#' )
#' repeated
sb_compare_repeated_cv <- function(task, pipelines, folds = 5L,
                                   seeds = 1:5, workers = 1L) {
  if (inherits(pipelines, "stormbreaker_stack")) {
    pipelines <- sb_upgrade_legacy_names(pipelines)
  } else if (is.list(pipelines)) pipelines <- lapply(pipelines, sb_upgrade_legacy_names)

  if (!is.numeric(workers) || length(workers) != 1L || is.na(workers) ||
      !workers %in% c(1L, 2L)) {
    sb_abort("`workers` must be 1 or 2.", "stormbreaker_bad_workers")
  }
  specifications <- sb_repeated_cv_specifications(folds, seeds)
  repetitions <- if (workers == 2L) {
    sb_parallel_cv_repetitions(task, pipelines, specifications)
  } else lapply(
    seq_along(specifications$folds),
    function(repeat_id) {
      seed <- specifications$seeds[[repeat_id]]
      sb_run_repeated_cv_iteration(
        task,
        pipelines,
        folds = specifications$folds[[repeat_id]],
        seed = if (is.na(seed)) 1L else seed
      )
    }
  )
  names(repetitions) <- paste0("repeat_", seq_along(repetitions))
  metadata <- lapply(repetitions, sb_repeated_cv_result_metadata)

  comparison_rows <- lapply(seq_along(repetitions), function(repeat_id) {
    repetition <- repetitions[[repeat_id]]
    repeat_metadata <- metadata[[repeat_id]]
    data.frame(
      repeat_id = as.integer(repeat_id),
      folds = repeat_metadata$folds,
      seed = repeat_metadata$seed,
      fold_source = repeat_metadata$fold_source,
      repetition$comparison,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })
  comparison <- do.call(rbind, comparison_rows)
  rownames(comparison) <- NULL

  fold_rows <- lapply(seq_along(repetitions), function(repeat_id) {
    data.frame(
      repeat_id = as.integer(repeat_id),
      metadata[[repeat_id]]$fold_assignment,
      stringsAsFactors = FALSE
    )
  })
  fold_assignments <- do.call(rbind, fold_rows)
  rownames(fold_assignments) <- NULL

  complete <- all(vapply(
    repetitions,
    function(repetition) identical(repetition$status, "success"),
    logical(1L)
  ))
  result <- list(
    status = if (complete) "success" else "incomplete",
    task_type = repetitions[[1L]]$task_type,
    comparison = comparison,
    repetitions = repetitions,
    n_patients = metadata[[1L]]$n_patients,
    n_repeats = as.integer(length(repetitions)),
    fold_assignments = fold_assignments
  )
  if (inherits(pipelines, "stormbreaker_stack")) {
    result$strategy_type <- "stack"
  }
  structure(result, class = "stormbreaker_repeated_cv_comparison")
}

#' @export
print.stormbreaker_repeated_cv_comparison <- function(x, ...) {
  cat("<stormbreaker_repeated_cv_comparison:", x$status, ">\n")
  count_label <- if (identical(x$strategy_type, "stack")) {
    "Member/stack rows:"
  } else {
    "Pipelines:"
  }
  cat(
    count_label, length(unique(x$comparison$pipeline_id)),
    "| Patients:", x$n_patients,
    "| Outer repetitions:", x$n_repeats, "\n"
  )
  metric_columns <- if (identical(x$task_type, "binary")) {
    c("roc_auc", "average_precision", "brier_score", "log_loss")
  } else {
    columns <- "harrell_c"
    if ("ipcw_brier" %in% names(x$comparison)) {
      columns <- c(columns, "ipcw_brier", "time_auc", "uno_c")
    }
    columns
  }
  print(
    x$comparison[c(
      "repeat_id", "pipeline_id", "status", metric_columns
    )],
    row.names = FALSE,
    digits = 3
  )
  cat(
    paste0(
      "No cross-repeat mean, confidence interval, ranking, or winner ",
      "is computed.\n"
    )
  )
  invisible(x)
}
