#' stormbreaker: auditable transcriptomic modeling and benchmarking
#'
#' STORMBREAKER stands for Systematic transcriptomic Risk Modeling and
#' Benchmarking for Reproducible, Evidence-Aware biomarker Evaluation in R.
#' It provides a method-focused workflow for numeric transcriptomic predictors
#' in patient-level binary and ordinary right-censored survival prediction.
#'
#' @section Core workflow:
#' * [sb_binary_task()] and [sb_surv_task()] define explicit patient IDs,
#'   outcomes, numeric predictors, and any study-declared prediction horizon or
#'   action thresholds.
#' * [sb_pipeline()] declares one admitted selector-and-learner route, including
#'   an optional outcome-independent fixed feature block. [sb_run_cv()] keeps
#'   supervised selection and tuning inside resampling and retains patient-level
#'   out-of-fold predictions, fitted-state audits, and failures.
#'   [sb_compare_cv()] and [sb_run_internal_external_cv()] reuse predeclared
#'   pipelines without choosing a winner.
#' * [sb_fit_final()] creates one deployment state from an explicitly supplied
#'   pipeline. [sb_predict_final()] reads only patient IDs and required numeric
#'   predictors; outcomes are attached later through [sb_binary_endpoint()] or
#'   [sb_surv_endpoint()] and evaluated with [sb_evaluate_external()].
#' * [sb_bootstrap_external()] quantifies external-patient sampling uncertainty
#'   conditional on the locked model and saved predictions. Separate collection
#'   and meta-analysis functions can retain or synthesize compatible results
#'   across named independent cohorts.
#'
#' @section Validation boundary:
#' Ordinary cross-validation and external bootstrap operate on patient rows;
#' internal-external validation instead holds out complete cohorts. Assessment
#' outcomes never enter selector fitting, hyperparameter tuning, final
#' prediction, feature mapping, or model choice. Failures remain visible and
#' never trigger a fallback learner. Fold-to-fold variability is not reported
#' as a confidence interval, and external outcomes are not used to choose a
#' model, feature, risk direction, horizon, or threshold.
#' Glmnet's optional internal predictor standardization is refitted inside each
#' relevant analysis fit; the package does not infer a general normalization or
#' imputation procedure.
#'
#' Discrimination, prediction error, calibration, and decision utility remain
#' separate estimands. Survival summaries retain their time unit, horizon, and
#' censoring assumptions. Bootstrap intervals keep development fitting and
#' selection uncertainty outside their stated external-sampling target.
#'
#' @section Formula and feature bridges:
#' [sb_linear_signature()] accepts only a completely stated published linear
#' formula and performs no fitting or inferred preprocessing.
#' [sb_apply_feature_map()] executes only a caller-reviewed one-to-one mapping
#' before prediction. It performs no annotation lookup, ambiguous-match
#' resolution, expression aggregation, normalization, imputation, or zero
#' filling.
#'
#' @section Scope:
#' The package supports numeric predictors, binary outcomes, and ordinary
#' right-censored survival outcomes. It does not currently claim support for
#' multiclass outcomes, competing risks, categorical predictors, automatic
#' model ranking, automatic threshold selection, or automatic biological
#' interpretation. Registry admission means that a built-in execution route
#' exists; it is not evidence that a method is appropriate for a particular
#' dataset or supported by a publication.
#'
#' @references McLernon DJ, Giardiello D, Van Calster B, et al. (2023).
#'   Assessing performance and clinical usefulness in prediction models with
#'   survival outcomes: practical guidance for Cox proportional hazards models.
#'   *Annals of Internal Medicine*, 176(1), 105-114.
#'   \doi{10.7326/M22-0844}.
#'
#' Uno H, Cai T, Pencina MJ, D'Agostino RB, Wei LJ (2011). On the C-statistics
#'   for evaluating overall adequacy of risk prediction procedures with
#'   censored survival data. *Statistics in Medicine*, 30(10), 1105-1117.
#'   \doi{10.1002/sim.4154}.
#'
#' @keywords internal
"_PACKAGE"
