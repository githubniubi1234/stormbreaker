test_that("cancerclass reproduces native z scores and fold-local calibration", {
  skip_if_not_installed("cancerclass", minimum_version = "1.52.0")
  skip_if_not_installed("Biobase", minimum_version = "2.68.0")

  d <- make_sb_binary_fixture(n = 60L, seed = 8901L)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  outcome <- d$outcome
  sample_names <- paste0("sample_", seq_len(nrow(x)))
  make_eset <- function(rows) {
    expression <- t(data.matrix(x[rows, , drop = FALSE]))
    colnames(expression) <- sample_names[rows]
    phenotype <- data.frame(
      class = factor(outcome[rows], levels = c(0L, 1L)),
      row.names = sample_names[rows]
    )
    Biobase::ExpressionSet(
      assayData = expression,
      phenoData = Biobase::AnnotatedDataFrame(phenotype)
    )
  }

  adapter <- sb_cancerclass_binary_learner_adapter()
  search_before_fit <- search()
  model <- adapter$fit(x, outcome)
  expect_identical(search(), search_before_fit)
  assessment_rows <- c(3L, 14L, 27L, 45L, 58L)
  withr::local_package("generics", quietly = TRUE)
  withr::local_package("BiocGenerics", quietly = TRUE)
  withr::local_package("Biobase", quietly = TRUE)
  native <- suppressMessages(cancerclass::fit(
    make_eset(seq_len(nrow(x))),
    class = "class",
    method = "welch.test",
    hparam = 0.75
  ))
  native_prediction <- suppressMessages(cancerclass::predict(
    native,
    make_eset(assessment_rows),
    positive = "1",
    class = "class",
    ngenes = ncol(x),
    dist = "cor"
  ))
  expected_score <- as.numeric(
    methods::slot(native_prediction, "prediction")[, "z"]
  )
  raw_score <- sb_cancerclass_raw_score(
    model,
    x[assessment_rows, , drop = FALSE]
  )
  probability <- adapter$predict(
    model,
    x[assessment_rows, , drop = FALSE]
  )

  expect_equal(model$predictor, methods::slot(native, "predictor"),
    tolerance = 0
  )
  expect_equal(raw_score, expected_score, tolerance = 1e-14)
  expect_equal(
    probability,
    stats::plogis(
      model$calibration_coefficients[[1L]] +
        model$calibration_coefficients[[2L]] * raw_score
    ),
    tolerance = 0
  )

  expected_fold_id <- sb_make_stratified_fold_id(
    outcome,
    folds = 5L,
    seed = 1L,
    stratum_label = "class"
  )
  expected_oof <- numeric(nrow(x))
  for (fold in seq_len(5L)) {
    analysis_rows <- which(expected_fold_id != fold)
    held_out_rows <- which(expected_fold_id == fold)
    inner_fit <- suppressMessages(cancerclass::fit(
      make_eset(analysis_rows),
      class = "class",
      method = "welch.test",
      hparam = 0.75
    ))
    predictor <- methods::slot(inner_fit, "predictor")
    held_out <- t(data.matrix(x[held_out_rows, , drop = FALSE]))
    held_out <- held_out[rownames(predictor), , drop = FALSE]
    distance_zero <- cancerclass::get.d(
      held_out,
      predictor[, "0"],
      method = "cor"
    )
    distance_one <- cancerclass::get.d(
      held_out,
      predictor[, "1"],
      method = "cor"
    )
    class_distance <- cancerclass::get.d(
      predictor[, "0"],
      predictor[, "1"],
      method = "cor"
    )
    expected_oof[held_out_rows] <-
      (distance_zero - distance_one) / class_distance
  }

  expect_identical(model$calibration_fold_id, expected_fold_id)
  expect_equal(model$calibration_oof_score, expected_oof, tolerance = 0)
  state <- unserialize(serialize(adapter$extract_state(model), NULL))
  replay <- adapter$predict_state(
    state,
    x[assessment_rows, , drop = FALSE]
  )
  expect_equal(replay$probability, probability, tolerance = 0)
  expect_identical(state$training_data_retention, "none")
  expect_identical(
    intersect(
      c("outcome", "calibration_fold_id", "calibration_oof_score"),
      names(state)
    ),
    character()
  )
})

test_that("cancerclass exposes only leakage-safe routes and runs binary CV", {
  skip_if_not_installed("cancerclass", minimum_version = "1.52.0")
  skip_if_not_installed("Biobase", minimum_version = "2.68.0")

  registry <- sb_pipeline_registry()
  rows <- registry[
    registry$learner == "cancerclass_binary",
    ,
    drop = FALSE
  ]
  expect_setequal(rows$selector, c("none", "fixed_features"))
  expect_identical(unique(rows$task_type), "binary")
  expect_identical(unique(rows$optional_dependency), "cancerclass + Biobase")
  unsafe <- capture_sb_error(sb_pipeline(
    selector = "univ_binary_auc_top_k",
    learner = "cancerclass_binary",
    selector_n_features = 2L
  ))
  expect_s3_class(unsafe, "stormbreaker_unadmitted_pipeline")

  d <- make_sb_binary_fixture(n = 90L, seed = 8911L)
  predictors <- c("gene_a", "gene_b", "gene_c")
  task <- sb_binary_task(d, "patient_id", "outcome", predictors)
  pipeline <- sb_pipeline(learner = "cancerclass_binary")
  result <- sb_run_cv(task, pipeline, folds = 3L, seed = 8912L)

  expect_identical(pipeline$prediction_type, "probability")
  expect_identical(pipeline$learner_params$feature_method, "welch.test")
  expect_identical(pipeline$learner_params$calibration, "oof_logistic")
  expect_identical(result$status, "success")
  expect_identical(result$complete_oof, TRUE)
  expect_true(all(is.finite(result$predictions$probability)))
  expect_identical(
    vapply(result$coefficients, `[[`, "", "status"),
    rep("success", 3L)
  )
  expect_identical(
    result$provenance$cancerclass_version,
    as.character(utils::packageVersion("cancerclass"))
  )
  expect_identical(
    result$provenance$biobase_version,
    as.character(utils::packageVersion("Biobase"))
  )

  fit <- sb_fit_final(task, pipeline)
  new_data <- d[c(7L, 28L, 74L),
    c("patient_id", "outcome", rev(predictors)),
    drop = FALSE
  ]
  altered <- new_data
  altered$outcome <- 1L - altered$outcome
  original_prediction <- sb_predict_final(fit, new_data, "patient_id")
  altered_prediction <- sb_predict_final(fit, altered, "patient_id")

  expect_identical(fit$status, "success")
  expect_identical(fit$learner_state$learner, "cancerclass_binary")
  expect_identical(original_prediction$status, "success")
  expect_equal(
    original_prediction$predictions$probability,
    altered_prediction$predictions$probability,
    tolerance = 0
  )
})

test_that("cancerclass fails explicitly without dependencies or feasible folds", {
  x <- data.frame(
    gene_a = seq_len(12L),
    gene_b = rev(seq_len(12L)),
    gene_c = rep(c(-1, 1), length.out = 12L)
  )
  adapter <- sb_cancerclass_binary_learner_adapter()

  local_mocked_bindings(sb_cancerclass_available = function() FALSE)
  missing <- capture_sb_error(adapter$fit(x, rep(c(0L, 1L), each = 6L)))
  expect_s3_class(missing, "stormbreaker_missing_cancerclass")

  local_mocked_bindings(sb_cancerclass_available = function() TRUE)
  single_class <- capture_sb_error(adapter$fit(x, rep(0L, nrow(x))))
  infeasible <- capture_sb_error(adapter$fit(
    x[seq_len(8L), , drop = FALSE],
    rep(c(0L, 1L), each = 4L)
  ))
  expect_s3_class(single_class, "stormbreaker_cancerclass_single_class")
  expect_s3_class(infeasible, "stormbreaker_infeasible_inner_folds")
})
