sb_mime117_alpha_text <- function(alpha) {
  format(alpha, nsmall = 1L, trim = TRUE, scientific = FALSE)
}

sb_mime117_component_label <- function(family, alpha = NA_real_,
                                        direction = NA_character_,
                                        single = FALSE) {
  if (identical(family, "Enet")) {
    return(paste0("Enet[\u03b1=", sb_mime117_alpha_text(alpha), "]"))
  }
  if (identical(family, "StepCox")) {
    return(paste0("StepCox[", direction, "]"))
  }
  if (identical(family, "survival-SVM") && single) {
    return("survival - SVM")
  }
  family
}

sb_mime117_variant_text <- function(family, alpha, direction) {
  if (identical(family, "Enet")) {
    return(paste0("alpha=", sb_mime117_alpha_text(alpha)))
  }
  if (identical(family, "StepCox")) {
    return(paste0("direction=", direction))
  }
  ""
}

sb_mime117_source_grid <- function() {
  alphas <- seq(0.1, 0.9, by = 0.1)
  directions <- c("both", "backward", "forward")
  rows <- list()

  add <- function(block, stage1_family, stage2_family = "",
                  alpha = NA_real_, direction = NA_character_) {
    has_stage2 <- nzchar(stage2_family)
    stage1_alpha <- if (identical(stage1_family, "Enet")) {
      alpha
    } else {
      NA_real_
    }
    stage2_alpha <- if (identical(stage2_family, "Enet")) {
      alpha
    } else {
      NA_real_
    }
    stage1_direction <- if (identical(stage1_family, "StepCox")) {
      direction
    } else {
      NA_character_
    }
    stage2_direction <- if (identical(stage2_family, "StepCox")) {
      direction
    } else {
      NA_character_
    }
    stage1_label <- sb_mime117_component_label(
      stage1_family,
      alpha = stage1_alpha,
      direction = stage1_direction,
      single = !has_stage2
    )
    label <- if (has_stage2) {
      paste(
        stage1_label,
        sb_mime117_component_label(
          stage2_family,
          alpha = stage2_alpha,
          direction = stage2_direction
        ),
        sep = " + "
      )
    } else {
      stage1_label
    }

    rows[[length(rows) + 1L]] <<- data.frame(
      config_id = "",
      official_order = 0L,
      block = block,
      pipeline_type = if (has_stage2) "two_stage" else "single",
      official_model_label = label,
      stage1_family = stage1_family,
      stage1_variant = sb_mime117_variant_text(
        stage1_family,
        stage1_alpha,
        stage1_direction
      ),
      stage2_family = stage2_family,
      stage2_variant = if (has_stage2) {
        sb_mime117_variant_text(
          stage2_family,
          stage2_alpha,
          stage2_direction
        )
      } else {
        ""
      },
      enet_alpha = if (identical(stage1_family, "Enet") ||
        identical(stage2_family, "Enet")) {
        alpha
      } else {
        NA_real_
      },
      stepcox_direction = if (identical(stage1_family, "StepCox") ||
        identical(stage2_family, "StepCox")) {
        direction
      } else {
        NA_character_
      },
      stringsAsFactors = FALSE
    )
  }

  add("RSF", "RSF")
  add("RSF", "RSF", "CoxBoost")
  for (alpha in alphas) {
    add("RSF", "RSF", "Enet", alpha = alpha)
  }
  for (stage2 in c("GBM", "Lasso", "plsRcox", "Ridge")) {
    add("RSF", "RSF", stage2)
  }
  for (direction in directions) {
    add("RSF", "RSF", "StepCox", direction = direction)
  }
  add("RSF", "RSF", "SuperPC")
  add("RSF", "RSF", "survival-SVM")

  for (alpha in alphas) {
    add("Enet", "Enet", alpha = alpha)
  }

  for (direction in directions) {
    add("StepCox", "StepCox", direction = direction)
  }
  for (direction in directions) {
    add("StepCox", "StepCox", "CoxBoost", direction = direction)
    for (alpha in alphas) {
      add(
        "StepCox",
        "StepCox",
        "Enet",
        alpha = alpha,
        direction = direction
      )
    }
    for (stage2 in c(
      "GBM", "Lasso", "plsRcox", "Ridge", "RSF", "SuperPC",
      "survival-SVM"
    )) {
      add("StepCox", "StepCox", stage2, direction = direction)
    }
  }

  add("CoxBoost", "CoxBoost")
  for (alpha in alphas) {
    add("CoxBoost", "CoxBoost", "Enet", alpha = alpha)
  }
  for (stage2 in c("GBM", "Lasso", "plsRcox", "Ridge")) {
    add("CoxBoost", "CoxBoost", stage2)
  }
  for (direction in directions) {
    add("CoxBoost", "CoxBoost", "StepCox", direction = direction)
  }
  add("CoxBoost", "CoxBoost", "SuperPC")
  add("CoxBoost", "CoxBoost", "survival-SVM")

  add("plsRcox", "plsRcox")
  add("SuperPC", "SuperPC")
  add("GBM", "GBM")
  add("survival-SVM", "survival-SVM")
  add("Ridge", "Ridge")

  add("Lasso", "Lasso")
  for (stage2 in c("CoxBoost", "GBM", "plsRcox", "RSF")) {
    add("Lasso", "Lasso", stage2)
  }
  for (direction in directions) {
    add("Lasso", "Lasso", "StepCox", direction = direction)
  }
  add("Lasso", "Lasso", "SuperPC")
  add("Lasso", "Lasso", "survival-SVM")

  grid <- do.call(rbind, rows)
  rownames(grid) <- NULL
  grid$official_order <- seq_len(nrow(grid))
  grid$config_id <- sprintf("MIME117-%03d", grid$official_order)
  grid
}

sb_mime117_learner <- function(family) {
  switch(
    family,
    RSF = "rsf_survival",
    Enet = "glmnet_cox",
    Lasso = "glmnet_cox",
    Ridge = "glmnet_cox",
    CoxBoost = "coxboost_cv",
    GBM = "gbm_cox_cv",
    plsRcox = "plsrcox_cv",
    SuperPC = "superpc_cv",
    `survival-SVM` = "survivalsvm",
    StepCox = "coxph"
  )
}

sb_mime117_selector <- function(stage1_family, final_family,
                                has_stage2) {
  if (!has_stage2) {
    if (identical(stage1_family, "StepCox")) "stepcox" else "none"
  } else if (identical(final_family, "StepCox")) {
    switch(
      stage1_family,
      RSF = "rsf_then_stepcox",
      CoxBoost = "coxboost_then_stepcox",
      Lasso = "lasso_then_stepcox"
    )
  } else {
    switch(
      stage1_family,
      RSF = "rsf_minimal_depth",
      StepCox = "stepcox",
      CoxBoost = "coxboost_nonzero",
      Lasso = "lasso_nonzero"
    )
  }
}

sb_mime117_pipeline_from_row <- function(row, inner_folds, inner_seed,
                                         selector_inner_folds,
                                         selector_seed, rsf_nodesize,
                                         learner_tuning_folds,
                                         plsrcox_tuning_folds,
                                         learner_seed,
                                         gbm_cv_max_trees) {
  has_stage2 <- identical(row$pipeline_type[[1L]], "two_stage")
  final_family <- if (has_stage2) {
    row$stage2_family[[1L]]
  } else {
    row$stage1_family[[1L]]
  }
  selector <- sb_mime117_selector(
    row$stage1_family[[1L]],
    final_family,
    has_stage2
  )
  learner <- sb_mime117_learner(final_family)
  args <- list(selector = selector, learner = learner)

  if (selector %in% c("stepcox", sb_mime_composite_selector_ids())) {
    args$stepcox_direction <- row$stepcox_direction[[1L]]
  }
  if (selector %in% sb_mime_rsf_selector_ids()) {
    args$selector_rsf_nodesize <- rsf_nodesize
    args$selector_seed <- selector_seed
  } else if (selector %in% sb_mime_tuned_selector_ids()) {
    args$selector_inner_folds <- selector_inner_folds
    args$selector_seed <- selector_seed
  }
  if (identical(learner, "glmnet_cox")) {
    alpha <- switch(
      final_family,
      Enet = row$enet_alpha[[1L]],
      Lasso = 1,
      Ridge = 0
    )
    args$inner_folds <- inner_folds
    args$inner_seed <- inner_seed
    args$glmnet_alpha <- alpha
    args$glmnet_standardize <- TRUE
    args$glmnet_lambda_rule <- "lambda.min"
    args$glmnet_tuning_measure <- "deviance"
  }
  if (learner %in% sb_analysis_tuned_survival_learner_ids()) {
    args$learner_tuning_folds <- if (identical(learner, "plsrcox_cv")) {
      plsrcox_tuning_folds
    } else {
      learner_tuning_folds
    }
    args$learner_seed <- learner_seed
    if (identical(learner, "gbm_cox_cv")) {
      args$gbm_cv_max_trees <- gbm_cv_max_trees
    }
  }
  do.call(sb_pipeline, args)
}

#' Construct the runnable MIME117 survival pipeline suite
#'
#' Reconstructs the 117 method configurations generated by the published
#' Mime1 loop structure as real STORMBREAKER pipeline objects in the original
#' source order. Shared selectors and learners are implemented once and reused
#' across configurations. Calling this function constructs and validates
#' specifications only: it does not fit models, inspect outcomes, rank routes,
#' or choose a winner.
#'
#' Learned selectors are fold-local when the returned pipelines are evaluated.
#' The suite preserves the method structure but intentionally does not copy
#' same-cohort apparent evaluation, unbounded retries, parallel defaults, or a
#' known inconsistent Ridge-alpha path from the reference implementation.
#'
#' @param inner_folds Number of inner folds used by downstream Cox glmnet
#'   learners.
#' @param inner_seed Positive integer seed for downstream Cox glmnet folds.
#' @param selector_inner_folds Number of inner folds used by CoxBoost and Lasso
#'   first-stage selectors.
#' @param selector_seed Positive integer seed used by learned first-stage
#'   selectors.
#' @param rsf_nodesize Positive integer terminal-node size used by RSF
#'   first-stage selection.
#' @param learner_tuning_folds Number of internal folds for the reference
#'   CoxBoost, SuperPC, and GBM learner fits.
#' @param plsrcox_tuning_folds Number of internal folds for the reference
#'   plsRcox learner fit. Mime1 uses the plsRcox default of 5.
#' @param learner_seed Positive integer seed used by analysis-fold tuned
#'   downstream learners.
#' @param gbm_cv_max_trees Maximum number of trees searched by reference GBM.
#'   A value of one fits the sole candidate without an inner tree-count search
#'   or CV curve; outer performance evaluation is unchanged.
#'
#' @return A 117-row data frame in reference source order. The `pipeline`
#'   list-column contains validated `stormbreaker_pipeline` objects; the other
#'   columns describe the reference method label and its two stage roles.
#' @export
sb_mime117_pipelines <- function(inner_folds = 10L,
                                 inner_seed = 5201314L,
                                 selector_inner_folds = 10L,
                                 selector_seed = 5201314L,
                                 rsf_nodesize = 5L,
                                 learner_tuning_folds = 10L,
                                 plsrcox_tuning_folds = 5L,
                                 learner_seed = 5201314L,
                                 gbm_cv_max_trees = 10000L) {
  inner_folds <- sb_check_positive_integer(
    inner_folds,
    "inner_folds",
    minimum = 3L
  )
  inner_seed <- sb_check_positive_integer(inner_seed, "inner_seed")
  selector_inner_folds <- sb_check_positive_integer(
    selector_inner_folds,
    "selector_inner_folds",
    minimum = 3L
  )
  selector_seed <- sb_check_positive_integer(selector_seed, "selector_seed")
  rsf_nodesize <- sb_check_positive_integer(rsf_nodesize, "rsf_nodesize")
  learner_tuning_folds <- sb_check_positive_integer(
    learner_tuning_folds,
    "learner_tuning_folds",
    minimum = 3L
  )
  plsrcox_tuning_folds <- sb_check_positive_integer(
    plsrcox_tuning_folds,
    "plsrcox_tuning_folds",
    minimum = 3L
  )
  learner_seed <- sb_check_positive_integer(learner_seed, "learner_seed")
  gbm_cv_max_trees <- sb_check_positive_integer(
    gbm_cv_max_trees,
    "gbm_cv_max_trees"
  )

  suite <- sb_mime117_source_grid()
  pipelines <- lapply(seq_len(nrow(suite)), function(index) {
    sb_mime117_pipeline_from_row(
      suite[index, , drop = FALSE],
      inner_folds = inner_folds,
      inner_seed = inner_seed,
      selector_inner_folds = selector_inner_folds,
      selector_seed = selector_seed,
      rsf_nodesize = rsf_nodesize,
      learner_tuning_folds = learner_tuning_folds,
      plsrcox_tuning_folds = plsrcox_tuning_folds,
      learner_seed = learner_seed,
      gbm_cv_max_trees = gbm_cv_max_trees
    )
  })
  suite$operational_fidelity <- paste0(
    "fold_local_mime_training_reproduction_not_bug_for_bug"
  )
  suite[["pipeline"]] <- pipelines
  suite
}
