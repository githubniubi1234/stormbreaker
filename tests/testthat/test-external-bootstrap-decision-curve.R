make_sb_dca_bootstrap_evaluation <- function(
    outcome,
    probability,
    decision_thresholds = NULL) {
  ids <- sprintf("d%02d", seq_along(outcome))
  prediction <- list(
    status = "success",
    task_type = "binary",
    pipeline_id = "none__logistic",
    learner = "logistic",
    id_col = "id",
    predictions = data.frame(
      patient_id = ids,
      pipeline_id = rep("none__logistic", length(ids)),
      learner = rep("logistic", length(ids)),
      positive_class = rep.int(1L, length(ids)),
      probability = probability,
      status = rep("success", length(ids)),
      error_class = rep(NA_character_, length(ids)),
      error_message = rep(NA_character_, length(ids)),
      stringsAsFactors = FALSE
    )
  )
  if (!is.null(decision_thresholds)) {
    prediction$decision_thresholds <- decision_thresholds
  }
  class(prediction) <- "stormbreaker_final_prediction"
  endpoint <- sb_binary_endpoint(
    data.frame(id = ids, outcome = outcome),
    "id",
    "outcome"
  )
  sb_evaluate_external(prediction, endpoint, reverse = FALSE)
}

test_that("DCA bootstrap reuses paired patient draws", {
  outcome <- c(0L, 0L, 1L, 1L, 0L, 1L)
  probability <- c(0.10, 0.30, 0.55, 0.75, 0.45, 0.90)
  thresholds <- c(0.20, 0.50, 0.80)
  reps <- 7L
  seed <- 41L
  confidence_level <- 0.8
  plain_evaluation <- make_sb_dca_bootstrap_evaluation(
    outcome,
    probability
  )
  dca_evaluation <- make_sb_dca_bootstrap_evaluation(
    outcome,
    probability,
    thresholds
  )

  plain <- sb_bootstrap_external(
    plain_evaluation,
    reps = reps,
    confidence_level = confidence_level,
    seed = seed
  )
  result <- sb_bootstrap_external(
    dca_evaluation,
    reps = reps,
    confidence_level = confidence_level,
    seed = seed
  )
  indices <- withr::with_seed(seed, {
    lapply(seq_len(reps), function(unused) {
      sample.int(length(outcome), length(outcome), replace = TRUE)
    })
  })

  expect_false("decision_curve" %in% names(plain))
  expect_identical(result$point_estimates, plain$point_estimates)
  expect_identical(result$intervals, plain$intervals)
  expect_identical(result$replicates, plain$replicates)
  expect_identical(result$resample_audit, plain$resample_audit)
  expect_identical(result$decision_curve$status, "success")
  expect_identical(
    result$decision_curve$interval_scope,
    "pointwise_not_simultaneous"
  )
  expect_identical(
    result$decision_curve$comparison_targets,
    c("model_minus_treat_all", "model_minus_treat_none")
  )
  expect_identical(result$decision_curve$n_thresholds, 3L)
  expect_identical(
    result$decision_curve$point_curve$threshold_probability,
    thresholds
  )
  expect_identical(
    nrow(result$decision_curve$replicates),
    reps * length(thresholds)
  )
  expect_identical(
    result$decision_curve$replicates$replicate,
    rep(seq_len(reps), each = length(thresholds))
  )
  expect_false(any(c("patient_id", "indices") %in%
    names(result$decision_curve$replicates)))

  first_index <- indices[[1L]]
  sampled_outcome <- outcome[first_index]
  sampled_probability <- probability[first_index]
  selected <- sampled_probability >= 0.5
  true_positive <- as.integer(sum(selected & sampled_outcome == 1L))
  false_positive <- as.integer(sum(selected & sampled_outcome == 0L))
  true_negative <- as.integer(sum(!selected & sampled_outcome == 0L))
  false_negative <- as.integer(sum(!selected & sampled_outcome == 1L))
  threshold_odds <- 0.5 / (1 - 0.5)
  model_net_benefit <- true_positive / length(outcome) -
    false_positive / length(outcome) * threshold_odds
  treat_all_net_benefit <- sum(sampled_outcome == 1L) / length(outcome) -
    sum(sampled_outcome == 0L) / length(outcome) * threshold_odds
  model_minus_treat_all <-
    (true_negative * threshold_odds - false_negative) / length(outcome)
  row <- result$decision_curve$replicates[
    result$decision_curve$replicates$replicate == 1L &
      result$decision_curve$replicates$threshold_probability == 0.5,
    ,
    drop = FALSE
  ]

  expect_identical(row$true_positive, true_positive)
  expect_identical(row$false_positive, false_positive)
  expect_identical(row$true_negative, true_negative)
  expect_identical(row$false_negative, false_negative)
  expect_equal(row$model_net_benefit, model_net_benefit, tolerance = 0)
  expect_equal(
    row$treat_all_net_benefit,
    treat_all_net_benefit,
    tolerance = 0
  )
  expect_identical(row$treat_none_net_benefit, 0)
  expect_equal(
    row$model_minus_treat_all,
    model_minus_treat_all,
    tolerance = 0
  )
  expect_equal(
    row$model_minus_treat_none,
    model_net_benefit,
    tolerance = 0
  )

  probabilities <- c(
    (1 - confidence_level) / 2,
    1 - (1 - confidence_level) / 2
  )
  for (threshold in thresholds) {
    point <- result$decision_curve$point_comparisons[
      result$decision_curve$point_comparisons$threshold_probability ==
        threshold,
      ,
      drop = FALSE
    ]
    for (comparison in result$decision_curve$comparison_targets) {
      interval <- result$decision_curve$intervals[
        result$decision_curve$intervals$threshold_probability == threshold &
          result$decision_curve$intervals$comparison == comparison,
        ,
        drop = FALSE
      ]
      replicate_values <- result$decision_curve$replicates[
        result$decision_curve$replicates$threshold_probability == threshold,
        comparison
      ]
      expected <- unname(stats::quantile(
        replicate_values,
        probs = probabilities,
        type = 7
      ))
      expect_identical(interval$n_requested, reps)
      expect_equal(interval$point_estimate, point[[comparison]], tolerance = 0)
      expect_equal(
        c(interval$lower, interval$upper),
        expected,
        tolerance = 0
      )
    }
  }
  scalar_targets <- unique(c(
    result$point_estimates$metric,
    result$intervals$metric,
    result$replicates$metric
  ))
  expect_false(any(grepl(
    "decision|net_benefit|threshold",
    scalar_targets,
    ignore.case = TRUE
  )))
  expect_match(
    paste(capture.output(print(result)), collapse = "\n"),
    "DCA paired pointwise intervals: 3 declared thresholds"
  )
})

test_that("paired DCA contrast avoids near-one cancellation", {
  outcome <- c(rep.int(0L, 50L), rep.int(1L, 50L))
  probability <- c(rep.int(1, 50L), rep.int(1, 49L), 0)
  threshold <- 1 - .Machine$double.eps
  evaluation <- make_sb_dca_bootstrap_evaluation(
    outcome,
    probability,
    threshold
  )

  result <- sb_bootstrap_external(
    evaluation,
    reps = 2L,
    confidence_level = 0.8,
    seed = 7L
  )
  row <- result$decision_curve$point_curve
  naive_difference <-
    row$model_net_benefit - row$treat_all_net_benefit
  direct_difference <-
    (row$true_negative * row$threshold_odds - row$false_negative) / row$n

  expect_identical(row$true_positive, 49L)
  expect_identical(row$false_positive, 50L)
  expect_identical(row$true_negative, 0L)
  expect_identical(row$false_negative, 1L)
  expect_equal(naive_difference, 0, tolerance = 0)
  expect_equal(direct_difference, -0.01, tolerance = 0)
  expect_equal(row$model_minus_treat_all, -0.01, tolerance = 0)
  expect_equal(
    row$model_minus_treat_none,
    row$model_net_benefit,
    tolerance = 0
  )
})

test_that("DCA bootstrap recomputes curves and rejects invalid thresholds", {
  outcome <- c(0L, 0L, 1L, 1L, 0L, 1L)
  probability <- c(0.10, 0.30, 0.55, 0.75, 0.45, 0.90)
  thresholds <- c(0.2, 0.5)
  evaluation <- make_sb_dca_bootstrap_evaluation(
    outcome,
    probability,
    thresholds
  )
  mutated_values <- evaluation
  mutated_values$metrics$decision_curve$curve$model_net_benefit <- 999
  mutated_values$metrics$decision_curve$curve$treat_all_net_benefit <- 998
  mutated_values$metrics$decision_curve$curve$treat_none_net_benefit <- 997

  result <- sb_bootstrap_external(
    mutated_values,
    reps = 3L,
    confidence_level = 0.8,
    seed = 17L
  )
  expected <- sb_binary_decision_curve_metric(
    evaluation$scored_predictions,
    thresholds,
    scope = "test_oracle"
  )$curve

  expect_identical(
    result$decision_curve$point_curve$model_net_benefit,
    expected$model_net_benefit
  )
  expect_identical(
    result$decision_curve$point_curve$treat_all_net_benefit,
    expected$treat_all_net_benefit
  )
  expect_identical(
    result$decision_curve$point_curve$treat_none_net_benefit,
    expected$treat_none_net_benefit
  )
  expect_false(any(
    result$decision_curve$point_curve$model_net_benefit == 999
  ))

  invalid <- evaluation
  invalid$metrics$decision_curve$curve$threshold_probability <- c(0.2, 1)
  failed <- sb_bootstrap_external(invalid, reps = 2L, seed = 17L)

  expect_s3_class(failed, "stormbreaker_external_bootstrap")
  expect_identical(failed$status, "failed")
  expect_identical(failed$error_stage, "evaluation_contract")
  expect_identical(
    failed$error_class,
    "stormbreaker_bad_decision_thresholds"
  )
})

test_that(
  "one-class external DCA remains estimable inside incomplete bootstrap",
  {
  evaluation <- make_sb_dca_bootstrap_evaluation(
    rep.int(1L, 5L),
    c(0.55, 0.65, 0.75, 0.85, 0.95),
    c(0.6, 0.8)
  )
  expect_identical(evaluation$status, "incomplete")

  result <- sb_bootstrap_external(
    evaluation,
    reps = 5L,
    confidence_level = 0.8,
    seed = 23L
  )

  expect_identical(result$status, "incomplete")
  expect_identical(result$decision_curve$status, "success")
  expect_true(all(is.finite(
    result$decision_curve$replicates$model_minus_treat_all
  )))
  expect_true(all(is.finite(
    result$decision_curve$replicates$model_minus_treat_none
  )))
  expect_true(all(
    result$point_estimates$status[
      result$point_estimates$metric %in% c("roc_auc", "average_precision")
    ] != "success"
  ))
  }
)
