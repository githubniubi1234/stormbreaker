sb_run_one_fold <- function(task, pipeline, runtime, split, seed) {
  id <- as.character(task$data[[task$id_col]])
  event <- as.integer(task$data[[task$event_col]])
  assessment_ids <- id[split$assessment_rows]
  prediction_horizon <- if ("prediction_horizon" %in% names(task)) {
    task[["prediction_horizon"]]
  } else {
    NULL
  }
  is_catboost <- pipeline$learner %in% c(
    "catboost_cox",
    "catboost_aft_survival"
  )
  is_gbm <- pipeline$learner %in% c("gbm_cox", "gbm_cox_cv")
  is_ranger <- pipeline$learner %in% sb_ranger_survival_learner_ids()
  empty_predictions <- data.frame(
    patient_id = assessment_ids,
    fold = split$fold,
    time = task$data[[task$time_col]][split$assessment_rows],
    event = as.integer(task$data[[task$event_col]][split$assessment_rows]),
    time_unit = task$time_unit,
    risk = NA_real_,
    status = "failed",
    pipeline_id = pipeline$id,
    learner = pipeline$learner,
    stringsAsFactors = FALSE
  )
  if (!is.null(prediction_horizon)) {
    probability_columns <- c(
      "prediction_horizon", "survival_probability",
      "survival_probability_status"
    )
    empty_predictions$prediction_horizon <- rep(
      as.numeric(prediction_horizon),
      length(assessment_ids)
    )
    empty_predictions$survival_probability <- rep(
      NA_real_,
      length(assessment_ids)
    )
    empty_predictions$survival_probability_status <- rep(
      "not_fitted",
      length(assessment_ids)
    )
    prediction_order <- append(
      setdiff(names(empty_predictions), probability_columns),
      probability_columns,
      after = match("risk", names(empty_predictions))
    )
    empty_predictions <- empty_predictions[prediction_order]
  }
  empty_selection <- list(
    fold = split$fold,
    selector = pipeline$selector,
    status = "failed",
    selected_features = character(),
    scores = NULL,
    error_class = NA_character_,
    error_message = NA_character_
  )
  empty_tuning <- list(
    fold = split$fold,
    learner = pipeline$learner,
    status = if (identical(pipeline$learner, "glmnet_cox")) {
      "failed"
    } else {
      "not_required"
    },
    lambda = NA_real_,
    alpha = if (identical(pipeline$learner, "glmnet_cox") ||
      pipeline$learner %in% sb_ncvreg_survival_learner_ids()) {
      pipeline$learner_params$alpha
    } else {
      NA_real_
    },
    selection_rule = if (identical(pipeline$learner, "glmnet_cox")) {
      pipeline$learner_params$selection_rule
    } else {
      NA_character_
    },
    measure = if (identical(pipeline$learner, "glmnet_cox")) {
      pipeline$learner_params$measure
    } else {
      NA_character_
    },
    inner_folds = if (identical(pipeline$learner, "glmnet_cox")) {
      pipeline$learner_params$inner_folds
    } else {
      NA_integer_
    },
    inner_seed = if (identical(pipeline$learner, "glmnet_cox")) {
      pipeline$learner_params$inner_seed
    } else {
      NA_integer_
    },
    fold_id = integer(),
    alignment = if (identical(pipeline$learner, "glmnet_cox")) {
      "lambda"
    } else {
      NA_character_
    },
    inner_selection_method = if (identical(
      pipeline$learner,
      "glmnet_cox"
    )) {
      sb_pipeline_inner_selection_method(pipeline)
    } else {
      NA_character_
    },
    inner_selections = list(),
    curve = NULL,
    warnings = character(),
    warning_validation = NULL,
    error_class = NA_character_,
    error_message = NA_character_,
    standardize = if (identical(pipeline$learner, "glmnet_cox")) {
      pipeline$learner_params$standardize
    } else {
      NA
    }
  )
  empty_coefficients <- list(
    fold = split$fold,
    learner = pipeline$learner,
    status = if (identical(pipeline$learner, "glmnet_cox") ||
      is.function(runtime$learner$extract_coefficients)) {
      "not_fitted"
    } else {
      "not_applicable"
    },
    lambda = NA_real_,
    alpha = if (identical(pipeline$learner, "glmnet_cox")) {
      pipeline$learner_params$alpha
    } else {
      NA_real_
    },
    n_nonzero = 0L,
    nonzero_features = character(),
    values = data.frame(
      feature = character(),
      coefficient = numeric(),
      nonzero = logical(),
      direction = character(),
      stringsAsFactors = FALSE
    )
  )
  empty_baseline_hazard <- if (is.null(prediction_horizon)) {
    NULL
  } else {
    list(
      fold = split$fold,
      learner = pipeline$learner,
      status = "not_fitted",
      prediction_horizon = as.numeric(prediction_horizon),
      time_unit = task$time_unit,
      cumulative_baseline_hazard = NA_real_,
      baseline_survival = NA_real_,
      calibration_coefficient = NA_real_,
      source = NA_character_,
      ties = NA_character_,
      error_class = NA_character_,
      error_message = NA_character_
    )
  }
  empty_tree_model <- if (is_catboost) {
    list(
      fold = split$fold,
      learner = pipeline$learner,
      status = "not_fitted",
      iterations = pipeline$learner_params$iterations,
      tree_count = NA_integer_,
      features = character(),
      params = pipeline$learner_params,
      model_serialization = NA_character_,
      model_raw = raw(),
      feature_importance = data.frame(
        feature = character(),
        importance = numeric(),
        stringsAsFactors = FALSE
      )
    )
  } else if (is_gbm) {
    list(
      fold = split$fold,
      learner = pipeline$learner,
      status = "not_fitted",
      n_trees = if ("n_trees" %in% names(pipeline$learner_params)) {
        pipeline$learner_params$n_trees
      } else {
        NA_integer_
      },
      features = character(),
      params = pipeline$learner_params,
      importance_type = "relative_influence",
      model_serialization = NA_character_,
      model_raw = raw(),
      feature_importance = data.frame(
        feature = character(),
        importance = numeric(),
        stringsAsFactors = FALSE
      )
    )
  } else if (is_ranger) {
    list(
      fold = split$fold,
      learner = pipeline$learner,
      status = "not_fitted",
      num_trees = pipeline$learner_params$num_trees,
      effective_mtry = NA_integer_,
      features = character(),
      params = pipeline$learner_params,
      importance_type = pipeline$learner_params$importance,
      model_serialization = NA_character_,
      model_raw = raw(),
      feature_importance = data.frame(
        feature = character(),
        importance = numeric(),
        stringsAsFactors = FALSE
      )
    )
  } else {
    NULL
  }
  selection <- empty_selection
  tuning <- empty_tuning
  coefficients <- empty_coefficients
  baseline_hazard <- empty_baseline_hazard
  tree_model <- empty_tree_model

  tryCatch(
    {
      runtime <- sb_validate_pipeline_runtime(runtime, pipeline)
      x_analysis <- task$data[
        split$analysis_rows,
        task$predictors,
        drop = FALSE
      ]
      x_assessment <- task$data[
        split$assessment_rows,
        task$predictors,
        drop = FALSE
      ]
      analysis_endpoint <- list(
        task_type = "survival",
        group_ids = sb_task_group_ids(task)[split$analysis_rows],
        time = task$data[[task$time_col]][split$analysis_rows],
        event = as.integer(
          task$data[[task$event_col]][split$analysis_rows]
        )
      )
      fitted_selector <- sb_fit_runtime_selector(
        runtime$selector,
        x_analysis,
        analysis_endpoint
      )
      selection <- list(
        fold = split$fold,
        selector = runtime$selector$id,
        status = "success",
        selected_features = fitted_selector$selected_features,
        scores = if (!is.null(fitted_selector$scores)) {
          fitted_selector$scores
        } else {
          NULL
        },
        error_class = NA_character_,
        error_message = NA_character_
      )
      selected_analysis <- runtime$selector$apply(
        fitted_selector,
        x_analysis
      )
      selected_assessment <- runtime$selector$apply(
        fitted_selector,
        x_assessment
      )
      prepared <- sb_prepare_runtime_learner(
        runtime$learner,
        pipeline,
        selected_analysis,
        analysis_endpoint,
        id[split$analysis_rows],
        group_ids = sb_task_group_ids(task, split$analysis_rows),
        selection_context = list(
          selector = runtime$selector,
          fitted_selector = fitted_selector,
          x = x_analysis
        )
      )
      learner <- prepared$learner
      if (!is.null(prepared$tuning)) {
        tuned <- prepared$tuning
        tuning <- list(
          fold = split$fold,
          learner = pipeline$learner,
          status = "success",
          lambda = tuned$lambda,
          alpha = tuned$alpha,
          selection_rule = tuned$selection_rule,
          measure = tuned$measure,
          inner_folds = tuned$inner_folds,
          inner_seed = tuned$seed,
          fold_id = tuned$fold_id,
          alignment = tuned$alignment,
          inner_selection_method = tuned$inner_selection_method,
          inner_selections = tuned$inner_selections,
          curve = tuned$curve,
          warnings = tuned$warnings,
          warning_validation = tuned$warning_validation,
          error_class = NA_character_,
          error_message = NA_character_,
          standardize = tuned$standardize
        )
      }
      model <- sb_fit_reusable_runtime_learner(
        learner,
        selected_analysis,
        analysis_endpoint,
        id[split$analysis_rows],
        pipeline$selector
      )
      fit_warnings <- attr(model, "stormbreaker_fit_warnings", exact = TRUE)
      if (length(fit_warnings) > 0L) {
        tuning$warnings <- c(
          tuning$warnings, paste0("Learner refit: ", fit_warnings)
        )
      }
      if (is.function(learner$extract_coefficients)) {
        coefficient_values <- learner$extract_coefficients(model)
        coefficients <- list(
          fold = split$fold,
          learner = pipeline$learner,
          status = "success",
          lambda = if (identical(pipeline$learner, "glmnet_cox")) {
            tuned$lambda
          } else if (pipeline$learner %in%
            sb_ncvreg_survival_learner_ids()) {
            model$effective$selected_lambda
          } else if (pipeline$learner %in%
            sb_ahaz_survival_learner_ids()) {
            model$audit$selected_lambda
          } else {
            NA_real_
          },
          alpha = if (identical(pipeline$learner, "glmnet_cox")) {
            tuned$alpha
          } else if (pipeline$learner %in%
            sb_ncvreg_survival_learner_ids()) {
            pipeline$learner_params$alpha
          } else if (identical(pipeline$learner, "ahaz_lasso_bic")) {
            pipeline$learner_params$alpha
          } else {
            NA_real_
          },
          n_nonzero = as.integer(sum(coefficient_values$nonzero)),
          nonzero_features = coefficient_values$feature[
            coefficient_values$nonzero
          ],
          values = coefficient_values
        )
      }
      if (is.function(learner$extract_tree_model)) {
        extracted <- learner$extract_tree_model(model)
        tree_model <- c(
          list(
            fold = split$fold,
            learner = pipeline$learner,
            status = "success"
          ),
          extracted
        )
      }
      prediction_details <- NULL
      if (!is.null(prediction_horizon) &&
          is.function(learner$predict_details)) {
        prediction_warned <- FALSE
        prediction_details <- withCallingHandlers(
          learner$predict_details(model, selected_assessment),
          warning = function(warning) prediction_warned <<- TRUE
        )
        risk <- prediction_details$risk
        # Retain ordinary repeated warning behavior rather than reusing a
        # warned prediction. Details never leave this assessment-fold call.
        if (prediction_warned) prediction_details <- NULL
      } else {
        risk <- learner$predict(model, selected_assessment)
      }
      successful_predictions <- empty_predictions
      successful_predictions$risk <- risk
      successful_predictions$status <- "success"
      if (!is.null(prediction_horizon)) {
        probability_result <- tryCatch({
          probability_args <- list(
            model,
            selected_analysis,
            task$data[[task$time_col]][split$analysis_rows],
            as.integer(task$data[[task$event_col]][split$analysis_rows]),
            selected_assessment,
            prediction_horizon
          )
          if (!is.null(prediction_details)) {
            probability_args$prediction_details <- prediction_details
          }
          do.call(learner$predict_survival, probability_args)
        },
          error = function(error) error
        )
        if (inherits(probability_result, "error")) {
          successful_predictions$survival_probability_status <- "failed"
          baseline_hazard$status <- "failed"
          baseline_hazard$error_class <- class(probability_result)[1]
          baseline_hazard$error_message <- conditionMessage(probability_result)
        } else {
          successful_predictions$survival_probability <-
            probability_result$probability
          successful_predictions$survival_probability_status <- "success"
          has_baseline_hazard <- all(c(
            "cumulative_baseline_hazard",
            "baseline_survival"
          ) %in% names(probability_result))
          baseline_hazard <- list(
            fold = split$fold,
            learner = pipeline$learner,
            status = if (has_baseline_hazard) {
              "success"
            } else {
              "not_applicable"
            },
            prediction_horizon = probability_result$prediction_horizon,
            time_unit = task$time_unit,
            cumulative_baseline_hazard = if (has_baseline_hazard) {
              probability_result$cumulative_baseline_hazard
            } else {
              NA_real_
            },
            baseline_survival = if (has_baseline_hazard) {
              probability_result$baseline_survival
            } else {
              NA_real_
            },
            calibration_coefficient = if (
              "calibration_coefficient" %in% names(probability_result)
            ) {
              probability_result$calibration_coefficient
            } else {
              NA_real_
            },
            source = probability_result$source,
            ties = probability_result$ties,
            error_class = NA_character_,
            error_message = NA_character_
          )
        }
      }
      list(
        predictions = successful_predictions,
        selection = selection,
        tuning = tuning,
        coefficients = coefficients,
        baseline_hazard = baseline_hazard,
        tree_model = tree_model,
        run = data.frame(
          fold = split$fold,
          seed = seed,
          pipeline_id = pipeline$id,
          selector_requested = pipeline$selector,
          learner_requested = pipeline$learner,
          learner_used = pipeline$learner,
          n_analysis = length(split$analysis_rows),
          n_analysis_events = sum(event[split$analysis_rows]),
          n_assessment = length(split$assessment_rows),
          n_assessment_events = sum(event[split$assessment_rows]),
          status = "success",
          error_class = NA_character_,
          error_message = NA_character_,
          stringsAsFactors = FALSE
        )
      )
    },
    error = function(error) {
      failed_selection <- selection
      if (identical(failed_selection$status, "failed")) {
        failed_selection$error_class <- class(error)[1]
        failed_selection$error_message <- conditionMessage(error)
      }
      failed_tuning <- tuning
      if (identical(pipeline$learner, "glmnet_cox") &&
        identical(failed_tuning$status, "failed")) {
        if (is.character(error$glmnet_warnings)) {
          failed_tuning$warnings <- error$glmnet_warnings
        }
        failed_tuning$error_class <- class(error)[1]
        failed_tuning$error_message <- conditionMessage(error)
      }
      list(
        predictions = empty_predictions,
        selection = failed_selection,
        tuning = failed_tuning,
        coefficients = coefficients,
        baseline_hazard = baseline_hazard,
        tree_model = tree_model,
        run = data.frame(
          fold = split$fold,
          seed = seed,
          pipeline_id = pipeline$id,
          selector_requested = pipeline$selector,
          learner_requested = pipeline$learner,
          learner_used = NA_character_,
          n_analysis = length(split$analysis_rows),
          n_analysis_events = sum(event[split$analysis_rows]),
          n_assessment = length(split$assessment_rows),
          n_assessment_events = sum(event[split$assessment_rows]),
          status = "failed",
          error_class = class(error)[1],
          error_message = conditionMessage(error),
          stringsAsFactors = FALSE
        )
      )
    }
  )
}

#' Run audited out-of-fold cross-validation
#'
#' `sb_run_cv()` dispatches from the supplied task to a survival or binary
#' pipeline. A scalar `folds` value creates deterministic folds stratified by
#' event or outcome. A multi-element numeric vector named by every task patient
#' ID instead supplies the exact outer-fold assignment; names, not vector
#' position, bind patients to folds. The package preserves a structurally valid
#' supplied plan without shuffling, stratifying, repairing, or rebalancing it.
#' Such acceptance does not certify design quality, independence,
#' representativeness, or transportability. A real study should declare or
#' import its plan independently of model performance.
#'
#' For each fold, the function fits the requested pipeline on analysis patients
#' and predicts assessment patients exactly once. A glmnet Cox or logistic
#' pipeline selects `lambda` by stratified inner CV using only that outer
#' fold's analysis patients, then refits on all of them. With a matching Top-K
#' selector, selection is repeated separately in every inner analysis fold;
#' it is never learned once from the corresponding inner assessment outcomes.
#' Ahaz Lasso and stepwise-SCAD instead select one lambda by BIC on the current
#' outer analysis fold only; they do not inspect its assessment patients.
#' A `fixed_features` selector instead applies one predeclared column block in
#' the same order in every fold. It never reads an outcome. A glmnet learner
#' tunes lambda within that already fixed block without treating it as
#' fold-local supervised selection.
#' Ordinary Cox, glmnet Cox, fixed CoxBoost, plsRcox, SuperPC, and the calibrated
#' linear survival-SVM return the zero-reference linear predictor `x_beta`;
#' additive and RBF survival-SVM return their negative native kernel scores times
#' their analysis-fold Cox calibration coefficients. Fixed CatBoost Cox returns
#' its native raw margin. Fixed CatBoost SurvivalAft and XGBoost AFT return the
#' negative native log-time margin; CatBoost uses a fixed Normal distribution
#' for direct horizon survival. Fixed GBM-Cox and XGBoost-Cox return their
#' native link-scale log hazards, while mboost AFT returns the negative native
#' distributional log-time link. Ahaz learners return their native additive-
#' hazard `x_beta`; their
#' declared-horizon probability is separately labelled analysis-fold
#' Cox/Breslow calibration. Bujar returns negative predicted log time.
#' Discrete-time survival BART summarizes its
#' posterior mean curve as negative stepwise restricted mean survival over the
#' analysis-fold time grid; its horizon survival remains a direct posterior
#' probability rather than a proportional-hazards calibration. The RLT
#' survival learner reconstructs native discrete hazard increments into a
#' survival curve and reports negative stepwise restricted mean survival over
#' its native event grid. Full Aalen
#' regression reports cumulative hazard at the last analysis-fold event time;
#' proportional-odds survival regression reports `x_beta`. Every binary
#' learner returns the
#' probability of class `1`. Gaussian Naive Bayes estimates empirical priors
#' and smoothed
#' class-feature Gaussian parameters separately in every outer analysis fold.
#' Cancerclass fits its native all-feature Welch/correlation score in each
#' analysis fold and learns its probability mapping only from five-fold inner
#' out-of-fold scores; it is not admitted after a learned upstream selector.
#' The four LVQ learners instead expose one declared logistic mapping of the
#' nearest class-specific squared-prototype distance margin. It preserves the
#' native nearest-prototype classification boundary but is not described as a
#' calibrated native posterior. Scaling, codebook initialization, and prototype
#' updates are repeated separately inside each analysis fold.
#' When glmnet predictor standardization is enabled in the pipeline, glmnet
#' learns it separately inside each inner analysis fit and repeats it for the
#' outer-analysis refit; it is not a package-wide preprocessing step.
#' The fixed tree learners use only predeclared settings and supply no
#' assessment patients for early stopping or best-model selection. A binary
#' AUC Top-K selector, when requested, is scored and fitted separately within
#' every outer analysis fold. Harrell
#' concordance, ROC AUC, and non-interpolated average precision use only saved
#' OOF predictions within the same outer assessment fold. Average precision is
#' aggregated by fold positive count without ranking patients scored by
#' different models. Binary Brier score averages patient-level squared errors,
#' while non-clipped log loss averages patient-level negative log likelihood
#' from the saved class-1 OOF probabilities. Each reference probability is
#' learned from the corresponding outer analysis fold's class-1 prevalence;
#' Brier and log-loss skill compare model and reference loss sums. Descriptive
#' calibration-in-the-large and calibration slope are fitted only to the pooled
#' saved OOF probabilities and outcomes. They assess weak calibration and never
#' alter an OOF prediction; boundary probabilities are not clipped, and a
#' separated recalibration fit is not reported as a finite slope. A separate
#' binary moderate-calibration curve fits one three-knot restricted natural
#' cubic spline logistic model to `qlogis()` of the pooled OOF probabilities.
#' Its ICI, E50, and E90 summarize the empirical distribution of absolute
#' differences between predicted and smoothed observed class-1 probability.
#' Exact zero or one is not clipped and leaves only that finite-logit curve
#' unavailable. A post-fit bounded signed-margin linear program also leaves a
#' complete or quasi-completely separated curve unavailable even if logistic
#' IRLS reports numerical convergence. The curve is not applied as recalibration
#' and has no universal acceptable cutoff. If the binary task declares decision
#' thresholds, the saved OOF probabilities also enter a patient-count decision
#' curve using `probability >= threshold`. Net benefit, treat-all, and
#' treat-none are reported at the declared thresholds. The same TP/FP/TN/FN
#' counts give sensitivity, specificity, positive predictive value, and
#' negative predictive value. A zero denominator remains `NA`; predictive
#' values describe the pooled OOF prevalence. The package does not derive a
#' grid, smooth or integrate the curve, or select a threshold. Counts are
#' additive over patients and do not compare predictions from different folds.
#' For a learned selector, the already
#' saved outer-analysis feature sets are also summarized by per-feature
#' selection frequency and every unique
#' fold-pair Jaccard similarity. The mean Jaccard is withheld unless every
#' requested selector fold succeeded. Because outer analysis folds overlap,
#' this is descriptive stability without uncertainty or chance correction.
#' A glmnet learner with `alpha > 0` separately summarizes the exact nonzero
#' predictor coefficients at each outer fold's selected `lambda`. This active-
#' set audit is distinct from any explicit Top-K selector. Ridge (`alpha = 0`)
#' and non-glmnet learners are not applicable, no coefficient cutoff is added,
#' and a fold pair with two empty active sets is left not evaluable rather than
#' called perfectly stable. This descriptive audit does not invalidate
#' otherwise complete OOF predictions.
#' Every fitted Cox-like linear survival learner and each ahaz learner also
#' retain their outer-fold coefficients so the declared linear predictor can be
#' reconstructed as `x_beta`. Additive and RBF survival-SVM outer-fold results
#' do not retain their kernel basis or native weights. A later `sb_fit_final()`
#' stores the scaled kernel basis and native weights required for outcome-blind
#' prediction, plus the actual RBF kernel parameter when applicable, without
#' outcomes or the fitted engine. Fixed CatBoost Cox and SurvivalAft, GBM-Cox,
#' and XGBoost survival learners instead retain their serialized models and
#' fixed prediction metadata.
#' Fixed ranger binary and
#' survival split-rule forests retain their serialized forest and descriptive impurity
#' importance. Fixed AdaBoost.M1 and LogitBoost retain their serialized native
#' models, requested and fitted iteration counts, and descriptive split-count
#' importance. Ranger survival risk is mortality, the row sum of the cumulative
#' hazard function. When the survival task declares one
#' `prediction_horizon`, ordinary Cox retains an Efron baseline. Glmnet Cox,
#' fixed CoxBoost, plsRcox, SuperPC, calibrated survival-SVM, CatBoost
#' Cox, GBM-Cox, XGBoost-Cox, and ahaz use a Breslow baseline matched to the
#' saved risk score; ahaz first fits its explicitly labelled nonnegative
#' analysis-fold Cox calibration coefficient. Ranger, CatBoost SurvivalAft,
#' XGBoost AFT, mboost AFT,
#' Bujar, SurvReg, randomForestSRC, partykit
#' Ctree/CForest, survival BART, Aalen, and proportional odds instead read
#' survival probability directly from
#' their fitted distributions, residual Kaplan-Meier curve, or coefficient
#' processes and do not invent a Cox baseline. Aalen uses its native
#' `exp(-H(t | x))`; an out-of-range native
#' estimate is reported as a probability failure rather than clipped or
#' recalibrated. Survival Rpart
#' uses an analysis-fold Breslow calibration matched to its relative-rate
#' score. The
#' corresponding IPCW Brier
#' score estimates each censoring distribution from the matching outer
#' analysis fold, orders tied events before censoring, and uses the
#' left-continuous censoring survival at event and prediction times. Censoring
#' exactly at the horizon is treated as known event-free at that instant. A
#' no-covariate event Kaplan-Meier reference is fitted in each outer analysis
#' fold and scored with the same assessment patients and IPCW weights. Model
#' and reference weighted losses are each divided by all assessment patients;
#' their aggregate loss ratio supplies a descriptive survival Brier skill
#' score.
#'
#' The same locked survival probabilities also produce a cumulative/dynamic
#' time-dependent AUC. The marker is `1 - S(t)`: observed events at or before
#' the horizon are cases, patients followed strictly beyond it are controls,
#' and patients censored at or before it are excluded. A case is weighted by
#' the inverse marginal censoring survival immediately before its event time;
#' exact marker ties receive half credit. Case-control pairs are formed only
#' inside the matching outer assessment fold, never across predictions from
#' different fitted models. The headline estimate divides summed successful
#' fold pair credit by the corresponding summed pair denominator and is
#' withheld unless every requested fold is evaluable. This uses the marginal
#' non-informative censoring assumption. Time-dependent AUC measures
#' discrimination, not calibration, transportability, clinical utility,
#' significance, or a decision threshold. These development-set OOF quantities
#' are not external validation.
#'
#' The saved risk scores separately produce Uno's truncated IPCW concordance at
#' the declared horizon. An index patient must have an observed event strictly
#' before the horizon; patients observed later and patients censored at the same
#' recorded time are comparators, while events tied at the same time are not.
#' Each eligible event receives 1 / G(T-)^2 from the matching outer-analysis
#' censoring curve, and exact risk ties receive half credit. Comparisons stay
#' within the matching outer assessment fold. The headline estimate is the
#' ratio of summed fold weighted pair credit to summed fold weighted pair
#' denominator, not a mean of fold estimates, and is withheld unless every fold
#' is evaluable. This statistic uses risk scores rather than survival
#' probabilities and does not compare scores from different outer models.
#'
#' When a horizon-enabled survival task also declares event-risk thresholds,
#' the same saved OOF probabilities enter a survival decision curve using
#' `1 - survival_probability >= threshold`. At each threshold, the selected
#' fraction is multiplied by the Kaplan-Meier event and event-free
#' probabilities among selected OOF patients to estimate true- and false-
#' positive rates. Treat-all uses the Kaplan-Meier survival of all OOF patients
#' and treat-none is zero. A nonempty group whose follow-up ends before the
#' horizon is unavailable unless its Kaplan-Meier survival has already reached
#' the absorbing value zero. This published estimator assumes independent
#' non-informative censoring in the relevant risk groups and censoring
#' independent of model predictors. It neither uses the IPCW Brier/time-AUC
#' weights nor refits a model. The curve is a point estimate at predeclared
#' thresholds; it does not choose a threshold or prove clinical effectiveness.
#'
#' Horizon-enabled survival results additionally fit one pooled moderate-
#' calibration curve to the saved OOF probabilities. A three-knot restricted
#' natural cubic spline of `log(-log(S(t)))` enters an Efron Cox model using the
#' complete right-censored follow-up; the calibrated survival probability is
#' then read at the locked horizon. ICI is the mean patient-level absolute
#' difference between predicted and calibrated survival probability, while E50
#' and E90 are its type-7 50th and 90th percentiles. The curve is descriptive,
#' assumes proportional hazards and independent non-informative censoring, and
#' has no universal acceptable cutoff. It is not applied as recalibration,
#' does not clip or change a prediction, and is distinct from mean or weak
#' calibration.
#'
#' The same horizon predictions separately report weak calibration. Mean
#' calibration is the pooled Kaplan-Meier observed event probability at the
#' horizon divided by mean predicted event probability (O:E); its target is 1.
#' The calibration slope is the coefficient from an Efron Cox model of the
#' endpoint administratively censored at the horizon on
#' `log(-log(survival_probability))`; its target is 1. Strictly later outcomes
#' cannot change that slope. The two components retain independent status:
#' exact zero/one or constant predictions can leave the slope unavailable
#' without erasing the mean ratio, and unsupported Kaplan-Meier follow-up can
#' leave the mean ratio unavailable without erasing an otherwise estimable
#' slope. No probability is clipped or recalibrated, and no survival calibration
#' intercept, p value, or universal cutoff is reported. The pooled OOF values
#' describe predictions from different outer fits; they are not external
#' validation or refit-uncertainty estimates. Any calibration failure remains
#' optional and does not erase otherwise valid OOF evidence.
#' Failed fits remain failed in the run, selection, prediction, fitted-state,
#' and metric audits; no replacement learner is used.
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
#' Ishwaran H, Kogalur UB, Blackstone EH, Lauer MS (2008). Random survival
#'   forests. *The Annals of Applied Statistics*, 2(3), 841-860.
#'   \doi{10.1214/08-AOAS169}.
#'
#' Wright MN, Ziegler A (2017). ranger: A fast implementation of random forests
#'   for high dimensional data in C++ and R. *Journal of Statistical Software*,
#'   77(1), 1-17. \doi{10.18637/jss.v077.i01}.
#'
#' @param task A task returned by [sb_surv_task()] or [sb_binary_task()]. A
#'   survival task may declare one study-defined `prediction_horizon` and
#'   clinically meaningful event-risk `decision_thresholds`; a binary task may
#'   likewise declare action thresholds. The package never chooses either from
#'   observed data or model performance.
#'   A separate `group_col` keeps repeated observations together in generated
#'   folds and glmnet learner tuning; supplied folds must respect these groups.
#'   Prediction rows and metric weights still refer to observations. Other
#'   engines' internal resampling is not configured by this declaration.
#' @param pipeline A matching pipeline returned by [sb_pipeline()]. `NULL`
#'   selects the minimum default for the supplied task.
#' @param folds Either one positive whole-number cross-validation fold count,
#'   for which each event or outcome stratum must contain at least that many
#'   patients, or a multi-element numeric vector named by the complete unique
#'   task patient IDs. Supplied values must be finite positive whole numbers
#'   with labels exactly `1, ..., K` for `K >= 2`. Vector order is ignored and
#'   the aligned assignment is preserved exactly. A length-one value, even if
#'   named, retains scalar-count behavior.
#' @param seed Positive integer used to generate scalar-count outer folds. It
#'   does not choose a supplied patient-ID assignment, whose effective outer
#'   seed is recorded as `NA`.
#' @param progress Show outer-fold start/completion messages; defaults to
#'   `interactive()`. Set TRUE for script logs.
#'
#' @return A `stormbreaker_cv_result` with:
#' * `predictions`: patient-level outcomes, OOF risk or class-1 probability,
#'   fold identity, and status; a survival task with a declared horizon also
#'   includes OOF survival probability and its status;
#' * `metrics`: fold-stratified Harrell concordance, or binary ROC AUC,
#'   non-interpolated average precision, and patient-level OOF Brier score and
#'   log loss with outer-analysis prevalence references and skill scores, plus
#'   descriptive OOF weak-calibration intercept and slope plus a separate pooled
#'   binary moderate-calibration curve and ICI/E50/E90, plus an optional
#'   patient-count decision curve with threshold operating characteristics at
#'   declared thresholds; per-fold audits and failures are retained, and no
#'   confidence interval is inferred from fold variability; a survival task
#'   with a declared horizon additionally reports
#'   analysis-fold IPCW Brier score, its censoring curves, a fold-local event
#'   Kaplan-Meier reference, the aggregate Brier skill score, and a
#'   within-fold cumulative/dynamic time-dependent AUC and truncated Uno IPCW
#'   concordance, plus a pooled
#'   three-knot Cox moderate-calibration curve and ICI/E50/E90, a pooled OOF
#'   observed/expected mean-calibration ratio and calibration slope, and, when
#'   declared, a pooled OOF survival decision curve at the locked horizon;
#' * `runs`: one audit row per fold;
#' * `selection_stability`: feature selection counts and frequencies over
#'   successful selector folds, plus all unique outer-analysis fold-pair
#'   Jaccard audits; pipelines without a learned selector are not applicable;
#' * `glmnet_active_set_stability`: for glmnet with `alpha > 0`, exact nonzero
#'   coefficient counts and frequencies plus all outer-fold active-set Jaccard
#'   audits; Ridge and non-glmnet learners are not applicable, and empty-union
#'   fold pairs remain not evaluable;
#' * `selections`: for supervised selectors, selected features and complete
#'   per-feature selection audit tables for every fold; a fixed selector retains its
#'   unchanged feature block without a score table;
#' * `tuning`: for glmnet learners, the selected `lambda`, predeclared lambda
#'   rule and standardization setting, lambda alignment, inner patient-to-fold
#'   assignment, error curve, compact fold-local
#'   selection audits when requested, and any tuning failure for every outer
#'   fold;
#' * `coefficients`: every fitted linear-survival coefficient, including
#'   zero-valued glmnet and SuperPC coefficients, every binary logistic
#'   coefficient including its intercept, or the cancerclass inner-OOF
#'   calibration intercept and score slope; Gaussian Naive Bayes and the fixed
#'   tree learners mark this audit as not applicable;
#' * `baseline_hazards`: when a survival horizon is declared, each Cox fold's
#'   `H0(t)`, `S0(t)`, engine-specific ties, source, failure, and any
#'   calibration coefficient multiplying the saved risk before applying
#'   `H0(t)`; ranger split-rule forest folds are transparently marked not
#'   applicable because their probabilities come directly from the forest;
#' * `distribution_parameters`: for Gaussian Naive Bayes, every outer-fold
#'   class prior, smoothing epsilon, mean, variance, and standard deviation;
#' * `tree_models`: for fixed LightGBM, fixed CatBoost (binary, Cox, or AFT), fixed
#'   GBM-Cox, fixed ranger binary or survival split-rule forests, fixed AdaBoost.M1, or
#'   fixed LogitBoost, every outer-fold serialized model, declared setting,
#'   fitted feature name, and all-feature importance audit; binary boosting also
#'   records requested and fitted iteration counts with split-count importance;
#' * `folds`: exact analysis and assessment patient IDs;
#' * `provenance`: the outer-fold source, effective count and seed, task-order
#'   patient-to-fold assignment, and R and dependency versions;
#' * `complete_oof`: whether every patient has one finite OOF risk or binary
#'   probability;
#' * `complete_survival_probability_oof`: when requested, whether every patient
#'   has one finite OOF survival probability at the declared horizon.
#' @export
#'
#' @examples
#' d <- sb_example_survival(n = 50, seed = 1)
#' task <- sb_surv_task(
#'   d, "patient_id", "time", "event",
#'   c("gene_a", "gene_b"), "month",
#'   prediction_horizon = 6
#' )
#' result <- sb_run_cv(task, folds = 5, seed = 2026)
#' result
#' outer_plan <- setNames(
#'   result$provenance$fold_assignment$fold,
#'   result$provenance$fold_assignment$patient_id
#' )
#' reused <- sb_run_cv(task, folds = rev(outer_plan), seed = 999)
#' identical(
#'   reused$provenance$fold_assignment$fold,
#'   result$provenance$fold_assignment$fold
#' )
#'
#' binary_data <- sb_example_binary(n = 50, seed = 1)
#' binary_task <- sb_binary_task(
#'   binary_data, "patient_id", "outcome",
#'   c("gene_a", "gene_b")
#' )
#' sb_run_cv(binary_task, folds = 5, seed = 2026)
#' sb_run_cv(
#'   binary_task,
#'   sb_pipeline(learner = "gaussian_naive_bayes"),
#'   folds = 5,
#'   seed = 2026
#' )
#' if (requireNamespace("lightgbm", quietly = TRUE)) {
#'   sb_run_cv(
#'     binary_task,
#'     sb_pipeline(
#'       learner = "lightgbm_binary",
#'       lightgbm_nrounds = 10
#'     ),
#'     folds = 5,
#'     seed = 2026
#'   )
#' }
#' if (requireNamespace("catboost", quietly = TRUE)) {
#'   sb_run_cv(
#'     binary_task,
#'     sb_pipeline(
#'       learner = "catboost_binary",
#'       catboost_iterations = 10
#'     ),
#'     folds = 5,
#'     seed = 2026
#'   )
#' }
sb_run_cv <- function(task, pipeline = NULL, folds = 5L,
                      seed = 1L, progress = interactive()) {
  sb_check_progress(progress)
  if (inherits(task, "stormbreaker_binary_task")) {
    return(sb_run_binary_cv(task, pipeline, folds, seed, progress))
  }
  if (is.null(pipeline)) {
    pipeline <- sb_pipeline()
  }
  if (!inherits(task, "stormbreaker_surv_task")) {
    sb_abort(
      "`task` must be created by `sb_surv_task()`.",
      "stormbreaker_bad_task"
    )
  }
  task <- sb_validate_task_contract(task)
  prediction_horizon <- if ("prediction_horizon" %in% names(task)) {
    task[["prediction_horizon"]]
  } else {
    NULL
  }
  pipeline <- sb_validate_pipeline_contract(pipeline, "survival")
  is_catboost <- pipeline$learner %in% c(
    "catboost_cox",
    "catboost_aft_survival"
  )
  is_gbm <- pipeline$learner %in% c("gbm_cox", "gbm_cox_cv")
  is_ranger <- pipeline$learner %in% sb_ranger_survival_learner_ids()
  runtime <- sb_resolve_pipeline_runtime(pipeline)

  split_object <- sb_make_outer_folds(task, folds, seed)
  fold_assignment <- data.frame(
    patient_id = as.character(task$data[[task$id_col]]),
    fold = split_object$assignment,
    stringsAsFactors = FALSE
  )
  fold_results <- lapply(
    split_object$splits,
    function(split) {
      sb_progress_message(progress, "Fold ", split$fold, "/", length(split_object$splits), " started")
      result <- sb_run_one_fold(
        task,
        pipeline,
        runtime,
        split,
        split_object$seed
      )
      sb_progress_message(progress, "Fold ", split$fold, "/", length(split_object$splits),
        " complete: ", result$run$status)
      result
    }
  )
  predictions <- do.call(
    rbind,
    lapply(fold_results, function(result) result$predictions)
  )
  runs <- do.call(rbind, lapply(fold_results, function(result) result$run))
  selections <- lapply(fold_results, function(result) result$selection)
  selection_stability <- sb_selection_stability(
    pipeline$selector,
    selections,
    task$predictors
  )
  tuning <- lapply(fold_results, function(result) result$tuning)
  coefficients <- lapply(
    fold_results,
    function(result) result$coefficients
  )
  glmnet_active_set_stability <- sb_glmnet_active_set_stability(
    pipeline,
    coefficients,
    task$predictors
  )
  baseline_hazards <- lapply(
    fold_results,
    function(result) result$baseline_hazard
  )
  tree_models <- lapply(
    fold_results,
    function(result) result$tree_model
  )
  rownames(predictions) <- NULL
  rownames(runs) <- NULL

  original_ids <- as.character(task$data[[task$id_col]])
  predictions <- predictions[
    match(original_ids, predictions$patient_id), ,
    drop = FALSE
  ]
  rownames(predictions) <- NULL
  exactly_once <- nrow(predictions) == task$n &&
    !anyNA(predictions$patient_id) &&
    anyDuplicated(predictions$patient_id) == 0L &&
    setequal(predictions$patient_id, original_ids)
  complete_oof <- exactly_once && all(predictions$status == "success") &&
    all(is.finite(predictions$risk))
  complete_survival_probability_oof <- if (is.null(prediction_horizon)) {
    TRUE
  } else {
    exactly_once &&
      all(predictions$survival_probability_status == "success") &&
      all(is.finite(predictions$survival_probability)) &&
      all(
        predictions$survival_probability >= 0 &
          predictions$survival_probability <= 1
      )
  }
  metrics <- list(
    harrell_c = sb_harrell_c_metric(predictions)
  )
  if (!is.null(prediction_horizon)) {
    analysis_inputs <- lapply(split_object$splits, function(split) {
      list(
        fold = split$fold,
        time = task$data[[task$time_col]][split$analysis_rows],
        event = as.integer(
          task$data[[task$event_col]][split$analysis_rows]
        )
      )
    })
    metrics$ipcw_brier <- sb_survival_ipcw_brier_metric(
      predictions,
      analysis_inputs,
      prediction_horizon
    )
    metrics$time_auc <- sb_survival_time_auc_metric(
      predictions,
      analysis_inputs,
      prediction_horizon,
      metrics$ipcw_brier$censoring_curves
    )
    metrics$uno_c <- sb_survival_uno_c_metric(
      predictions,
      analysis_inputs,
      prediction_horizon,
      metrics$ipcw_brier$censoring_curves
    )
    calibration_status <- ifelse(
      predictions$status == "success" &
        predictions$survival_probability_status == "success",
      "success",
      "failed"
    )
    metrics$calibration <- tryCatch(
      sb_survival_calibration_metric(
        patient_id = predictions$patient_id,
        time = predictions$time,
        event = predictions$event,
        survival_probability = predictions$survival_probability,
        prediction_status = calibration_status,
        prediction_horizon = prediction_horizon,
        time_unit = task$time_unit,
        scope = "pooled_patient_level_out_of_fold",
        fold = predictions$fold
      ),
      error = function(error) {
        sb_survival_calibration_error_metric(
          error = error,
          prediction_horizon = prediction_horizon,
          time_unit = task$time_unit,
          scope = "pooled_patient_level_out_of_fold",
          n = nrow(predictions),
          fold = predictions$fold
        )
      }
    )
    metrics$weak_calibration <- tryCatch(
      sb_survival_weak_calibration_metric(
        patient_id = predictions$patient_id,
        time = predictions$time,
        event = predictions$event,
        survival_probability = predictions$survival_probability,
        prediction_status = calibration_status,
        prediction_horizon = prediction_horizon,
        time_unit = task$time_unit,
        scope = "pooled_patient_level_out_of_fold",
        fold = predictions$fold
      ),
      error = function(error) {
        sb_survival_weak_calibration_error_metric(
          error = error,
          prediction_horizon = prediction_horizon,
          time_unit = task$time_unit,
          scope = "pooled_patient_level_out_of_fold",
          n = nrow(predictions),
          fold = predictions$fold
        )
      }
    )
    if ("decision_thresholds" %in% names(task)) {
      metrics$decision_curve <- sb_survival_decision_curve_metric(
        predictions,
        task$decision_thresholds,
        prediction_horizon,
        task$time_unit,
        scope = "pooled_patient_level_out_of_fold"
      )
    }
  }
  decision_curve_complete <-
    !("decision_thresholds" %in% names(task)) ||
      identical(metrics$decision_curve$status, "success")
  requested_metrics_complete <- is.null(prediction_horizon) ||
    (
      identical(metrics$ipcw_brier$status, "success") &&
        identical(metrics$time_auc$status, "success") &&
        identical(metrics$uno_c$status, "success") &&
        decision_curve_complete
    )

  result <- list(
    status = if (all(runs$status == "success") && complete_oof &&
      complete_survival_probability_oof && requested_metrics_complete) {
      "success"
    } else {
      "failed"
    },
    complete_oof = complete_oof,
    predictions = predictions,
    metrics = metrics,
    runs = runs,
    selection_stability = selection_stability,
    folds = split_object$splits,
    task = list(
      n = task$n,
      n_events = task$n_events,
      n_censored = task$n_censored,
      id_col = task$id_col,
      time_col = task$time_col,
      event_col = task$event_col,
      time_unit = task$time_unit,
      predictors = task$predictors,
      independent_unit = task$independent_unit
    ),
    pipeline = pipeline,
    provenance = list(
      seed = split_object$seed,
      folds = split_object$n_folds,
      fold_source = split_object$source,
      fold_assignment = fold_assignment,
      package = "stormbreaker",
      package_version = as.character(utils::packageVersion("stormbreaker")),
      r_version = R.version.string,
      survival_version = as.character(utils::packageVersion("survival"))
    )
  )
  if (!is.null(prediction_horizon)) {
    result <- append(
      result,
      list(
        complete_survival_probability_oof =
          complete_survival_probability_oof
      ),
      after = match("complete_oof", names(result))
    )
    result <- append(
      result,
      list(baseline_hazards = baseline_hazards),
      after = match("selection_stability", names(result))
    )
    result$task <- append(
      result$task,
      list(prediction_horizon = prediction_horizon),
      after = match("time_unit", names(result$task))
    )
    if ("decision_thresholds" %in% names(task)) {
      result$task <- append(
        result$task,
        list(decision_thresholds = task$decision_thresholds),
        after = match("prediction_horizon", names(result$task))
      )
    }
  }
  if (!is.null(task$group_col)) {
    result$task$group_col <- task$group_col
    result$provenance$group_ids <- stats::setNames(sb_task_group_ids(task), original_ids)
  }
  if (!identical(pipeline$selector, "none")) {
    result <- append(
      result,
      list(selections = selections),
      after = match("runs", names(result))
    )
  }
  if (pipeline$learner %in% c(
    sb_linear_survival_learner_ids(),
    sb_ahaz_survival_learner_ids()
  )) {
    result <- append(
      result,
      list(coefficients = coefficients),
      after = match("selection_stability", names(result))
    )
  }
  if (identical(pipeline$learner, "glmnet_cox")) {
    result$tuning <- tuning
    result$provenance$glmnet_version <- if (sb_glmnet_available()) {
      as.character(utils::packageVersion("glmnet"))
    } else {
      NA_character_
    }
  }
  classical_dependency <- if (
    pipeline$learner %in% sb_survivalsvm_learner_ids()) {
    c("survivalsvm_version", "survivalsvm")
  } else {
    switch(pipeline$learner,
      coxboost = c("coxboost_version", "CoxBoost"),
      plsrcox = c("plsrcox_version", "plsRcox"),
      superpc = c("superpc_version", "superpc"),
      NULL
    )
  }
  if (!is.null(classical_dependency)) {
    dependency_available <- if (
      pipeline$learner %in% sb_survivalsvm_learner_ids()) {
      sb_survivalsvm_available()
    } else {
      switch(pipeline$learner,
        coxboost = sb_coxboost_available(),
        plsrcox = sb_plsrcox_available(),
        superpc = sb_superpc_available()
      )
    }
    result$provenance[[classical_dependency[[1L]]]] <- if (
      dependency_available) {
      as.character(utils::packageVersion(classical_dependency[[2L]]))
    } else {
      NA_character_
    }
  }
  if (is_catboost || is_gbm || is_ranger) {
    result <- append(
      result,
      list(tree_models = tree_models),
      after = match("selection_stability", names(result))
    )
    if (is_catboost) {
      result$provenance$catboost_version <- if (sb_catboost_available()) {
        as.character(utils::packageVersion("catboost"))
      } else {
        NA_character_
      }
    } else if (is_gbm) {
      result$provenance$gbm_version <- if (sb_gbm_available()) {
        as.character(utils::packageVersion("gbm"))
      } else {
        NA_character_
      }
    } else {
      result$provenance$ranger_version <- if (sb_ranger_available()) {
        as.character(utils::packageVersion("ranger"))
      } else {
        NA_character_
      }
    }
  }
  if (pipeline$learner %in% sb_xgboost_survival_learner_ids() || pipeline$selector == "xgboost_cox_importance") {
    result$provenance$xgboost_version <- if (sb_xgboost_available()) {
      as.character(utils::packageVersion("xgboost"))
    } else {
      NA_character_
    }
  }
  if (pipeline$learner %in% sb_classical_binary_learner_ids()) {
    dependency <- sb_classical_binary_dependency(pipeline$learner)
    if (!is.na(dependency$package)) {
      field <- paste0(tolower(dependency$package), "_version")
      result$provenance[[field]] <- if (dependency$available) {
        as.character(utils::packageVersion(dependency$package))
      } else {
        NA_character_
      }
    }
  }
  if (pipeline$learner %in% sb_mboost_cox_learner_ids()) {
    dependency <- sb_mboost_cox_dependency(pipeline$learner)
    for (package in dependency$packages) {
      field <- paste0(tolower(package), "_version")
      result$provenance[[field]] <- if (
        requireNamespace(package, quietly = TRUE)) {
        as.character(utils::packageVersion(package))
      } else {
        NA_character_
      }
    }
  }
  if (pipeline$learner %in% sb_mboost_aft_learner_ids()) {
    dependency <- sb_mboost_aft_dependency(pipeline$learner)
    for (package in dependency$packages) {
      field <- paste0(tolower(package), "_version")
      result$provenance[[field]] <- if (
        requireNamespace(package, quietly = TRUE)) {
        as.character(utils::packageVersion(package))
      } else {
        NA_character_
      }
    }
  }
  if (pipeline$learner %in% sb_mboost_ranking_learner_ids()) {
    dependency <- sb_mboost_ranking_dependency(pipeline$learner)
    for (package in dependency$packages) {
      field <- paste0(tolower(package), "_version")
      result$provenance[[field]] <- if (
        requireNamespace(package, quietly = TRUE)) {
        as.character(utils::packageVersion(package))
      } else {
        NA_character_
      }
    }
  }
  if (pipeline$learner %in% sb_additional_survival_learner_ids()) {
    dependency <- sb_additional_survival_dependency(pipeline$learner)
    for (package in dependency$packages) {
      field <- paste0(tolower(package), "_version")
      if (!identical(field, "survival_version")) {
        result$provenance[[field]] <- if (
          requireNamespace(package, quietly = TRUE)) {
          as.character(utils::packageVersion(package))
        } else {
          NA_character_
        }
      }
    }
  }
  if (pipeline$learner %in% c(
    sb_flexsurv_parametric_learner_ids(),
    sb_flexsurvspline_learner_ids()
  )) {
    result$provenance$flexsurv_version <- if (requireNamespace(
      "flexsurv",
      quietly = TRUE
    )) {
      as.character(utils::packageVersion("flexsurv"))
    } else {
      NA_character_
    }
  }
  if (pipeline$learner %in% sb_timereg_survival_learner_ids()) {
    result$provenance$timereg_version <- if (sb_timereg_available()) {
      as.character(utils::packageVersion("timereg"))
    } else {
      NA_character_
    }
  }
  if (pipeline$learner %in% sb_abess_learner_ids() ||
    pipeline$selector %in% sb_abess_selector_ids()) {
    result$provenance$abess_version <- if (sb_abess_available()) {
      as.character(utils::packageVersion("abess"))
    } else {
      NA_character_
    }
  }
  if (pipeline$learner %in% sb_ncvreg_learner_ids()) {
    result$provenance$ncvreg_version <- if (sb_ncvreg_available()) {
      as.character(utils::packageVersion("ncvreg"))
    } else {
      NA_character_
    }
  }
  if (pipeline$learner %in% sb_rlt_survival_learner_ids()) {
    result$provenance$rlt_version <- if (sb_rlt_available()) {
      as.character(utils::packageVersion("RLT"))
    } else {
      NA_character_
    }
  }
  if (pipeline$learner %in% sb_coxphf_learner_ids()) {
    result$provenance$coxphf_version <- if (sb_coxphf_available()) {
      as.character(utils::packageVersion("coxphf"))
    } else {
      NA_character_
    }
  }
  if (pipeline$learner %in% sb_akritas_learner_ids()) {
    for (package in c("survivalmodels", "distr6"))
      result$provenance[[paste0(package, "_version")]] <- if (requireNamespace(package, quietly = TRUE))
        as.character(utils::packageVersion(package)) else NA_character_
  }
  if (pipeline$learner %in% sb_bagged_knn_survival_learner_ids()) {
    result$provenance$bnn_survival_version <- if (
      sb_bagged_knn_survival_available()
    ) {
      as.character(utils::packageVersion("bnnSurvival"))
    } else {
      NA_character_
    }
  }
  if (pipeline$learner %in% sb_ahaz_survival_learner_ids()) {
    result$provenance$ahaz_version <- if (sb_ahaz_available()) {
      as.character(utils::packageVersion("ahaz"))
    } else {
      NA_character_
    }
  }
  if (pipeline$learner %in% sb_pamr_survival_learner_ids()) {
    result$provenance$pamr_version <- if (sb_pamr_survival_available()) {
      as.character(utils::packageVersion("pamr"))
    } else {
      NA_character_
    }
  }
  if (pipeline$learner %in% sb_bujar_linear_learner_ids()) {
    result$provenance$bujar_version <- if (sb_bujar_available()) {
      as.character(utils::packageVersion("bujar"))
    } else {
      NA_character_
    }
  }
  if (pipeline$learner %in% sb_rule_aft_learner_ids()) {
    dependency <- sb_rule_aft_dependency(pipeline$learner)
    for (package in dependency$packages) {
      field <- paste0(tolower(package), "_version")
      result$provenance[[field]] <- if (
        requireNamespace(package, quietly = TRUE)) {
        as.character(utils::packageVersion(package))
      } else {
        NA_character_
      }
    }
  }
  if (pipeline$learner %in% sb_mgcv_gam_learner_ids()) {
    result$provenance$mgcv_version <- if (sb_mgcv_available()) {
      as.character(utils::packageVersion("mgcv"))
    } else {
      NA_character_
    }
  }
  if (pipeline$selector %in% sb_forest_selector_ids()) {
    selector_provenance <- sb_forest_selector_provenance(pipeline$selector)
    for (field in names(selector_provenance)) {
      result$provenance[[field]] <- selector_provenance[[field]]
    }
  }
  if (pipeline$selector %in% sb_pamr_selector_ids()) {
    result$provenance$pamr_version <- if (requireNamespace("pamr", quietly = TRUE))
      as.character(utils::packageVersion("pamr")) else NA_character_
  }
  high_dimensional_methods <- intersect(
    c(pipeline$selector, pipeline$learner),
    sb_high_dimensional_method_ids()
  )
  for (method in high_dimensional_methods) {
    method_provenance <- sb_high_dimensional_method_provenance(method)
    for (field in names(method_provenance)) {
      result$provenance[[field]] <- method_provenance[[field]]
    }
  }
  result$glmnet_active_set_stability <- glmnet_active_set_stability
  class(result) <- "stormbreaker_cv_result"
  result
}

#' @export
print.stormbreaker_cv_result <- function(x, ...) {
  cat("<stormbreaker_cv_result:", x$status, ">\n")
  cat("Pipeline:", x$pipeline$id, "| Folds:", nrow(x$runs), "\n")
  fold_label <- if (identical(
    x$provenance$fold_source,
    "preassigned_patient_id"
  )) {
    "preassigned patient-ID folds"
  } else {
    "generated stratified folds"
  }
  cat("Outer resampling:", fold_label, "\n")
  cat(
    "Patients:", nrow(x$predictions),
    "| Successful folds:", sum(x$runs$status == "success"), "\n"
  )
  cat("Complete OOF:", x$complete_oof, "\n")
  if (!is.null(x$complete_survival_probability_oof)) {
    cat(
      "Complete OOF survival probability at",
      format(x$task[["prediction_horizon"]]),
      paste0(x$task$time_unit, ":"),
      x$complete_survival_probability_oof, "\n"
    )
    ipcw_brier <- x$metrics$ipcw_brier
    if (identical(ipcw_brier$status, "success")) {
      cat(
        "IPCW Brier score:",
        format(round(ipcw_brier$estimate, 3), nsmall = 3), "\n"
      )
    } else {
      cat("IPCW Brier score: unavailable (incomplete)\n")
    }
    if (identical(ipcw_brier$reference$status, "success")) {
      cat(
        "IPCW Brier reference:",
        format(round(ipcw_brier$reference$estimate, 3), nsmall = 3), "\n"
      )
    } else {
      cat("IPCW Brier reference: unavailable (incomplete)\n")
    }
    if (identical(ipcw_brier$skill_score$status, "success")) {
      cat(
        "IPCW Brier skill score:",
        format(round(ipcw_brier$skill_score$estimate, 3), nsmall = 3), "\n"
      )
    } else {
      cat("IPCW Brier skill score: unavailable (incomplete)\n")
    }
    time_auc <- x$metrics$time_auc
    if (identical(time_auc$status, "success")) {
      cat(
        "Time-dependent AUC:",
        format(round(time_auc$estimate, 3), nsmall = 3),
        "| Weighted pair denominator:",
        format(time_auc$pair_denominator), "\n"
      )
    } else {
      cat("Time-dependent AUC: unavailable (incomplete)\n")
    }
    uno_c <- x$metrics$uno_c
    if (identical(uno_c$status, "success")) {
      cat(
        "Fold-stratified Uno C:",
        format(round(uno_c$estimate, 3), nsmall = 3),
        "| Weighted pair denominator:",
        format(uno_c$weighted_pair_denominator), "\n"
      )
    } else {
      cat("Fold-stratified Uno C: unavailable (incomplete)\n")
    }
    calibration <- x$metrics$calibration
    if (identical(calibration$status, "success")) {
      cat(
        "Survival calibration ICI:",
        format(round(calibration$ici, 3), nsmall = 3),
        "| E50:", format(round(calibration$e50, 3), nsmall = 3),
        "| E90:", format(round(calibration$e90, 3), nsmall = 3), "\n"
      )
    } else {
      cat("Survival calibration: unavailable (incomplete)\n")
    }
    weak_calibration <- x$metrics$weak_calibration
    if (identical(weak_calibration$mean_calibration$status, "success")) {
      cat(
        "Mean calibration O:E:",
        format(
          round(weak_calibration$mean_calibration_ratio, 3),
          nsmall = 3
        ), "\n"
      )
    } else {
      cat("Mean calibration O:E: unavailable (incomplete)\n")
    }
    if (identical(weak_calibration$slope_fit$status, "success")) {
      cat(
        "Survival calibration slope:",
        format(round(weak_calibration$calibration_slope, 3), nsmall = 3),
        "\n"
      )
    } else {
      cat("Survival calibration slope: unavailable (incomplete)\n")
    }
    if ("decision_curve" %in% names(x$metrics)) {
      decision_curve <- x$metrics$decision_curve
      if (identical(decision_curve$status, "success")) {
        cat(
          "OOF survival decision curve:",
          decision_curve$n_thresholds,
          "declared thresholds\n"
        )
      } else {
        cat("OOF survival decision curve: unavailable (incomplete)\n")
      }
    }
  }
  harrell_c <- x$metrics$harrell_c
  if (identical(harrell_c$status, "success")) {
    cat(
      "Fold-stratified Harrell C:",
      format(round(harrell_c$estimate, 3), nsmall = 3),
      "| Comparable pairs:", harrell_c$comparable_pairs, "\n"
    )
  } else {
    cat("Fold-stratified Harrell C: unavailable (incomplete)\n")
  }
  selection_stability <- x$selection_stability
  if (identical(selection_stability$status, "success")) {
    cat(
      "Selection stability (mean pairwise Jaccard):",
      format(
        round(selection_stability$mean_pairwise_jaccard, 3),
        nsmall = 3
      ), "\n"
    )
  } else if (identical(selection_stability$status, "incomplete")) {
    cat("Selection stability: unavailable (incomplete)\n")
  }
  active_set_stability <- x$glmnet_active_set_stability
  if (identical(active_set_stability$status, "success")) {
    cat(
      "Glmnet active-set stability (mean pairwise Jaccard):",
      format(
        round(active_set_stability$mean_pairwise_jaccard, 3),
        nsmall = 3
      ), "\n"
    )
  } else if (identical(active_set_stability$status, "incomplete")) {
    cat("Glmnet active-set stability: unavailable (incomplete)\n")
  } else if (identical(active_set_stability$status, "not_evaluable")) {
    cat("Glmnet active-set stability: not evaluable (empty union)\n")
  }
  invisible(x)
}
