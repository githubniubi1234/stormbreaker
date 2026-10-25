make_binary_moderate_scored <- function(n = 24L) {
  probability <- seq(0.03, 0.97, length.out = n)
  data.frame(
    patient_id = paste0("patient_", seq_len(n)),
    fold = rep(seq_len(4L), length.out = n),
    outcome = rep(c(0L, 1L), length.out = n),
    probability = probability,
    status = rep("success", n),
    stringsAsFactors = FALSE
  )
}

manual_binary_rcs3_basis <- function(x, knots) {
  d <- function(value, knot) {
    (pmax(value - knot, 0)^3 - pmax(value - knots[[3L]], 0)^3) /
      (knots[[3L]] - knot)
  }
  cbind(
    rcs_linear = x,
    rcs_nonlinear = d(x, knots[[1L]]) - d(x, knots[[2L]])
  )
}

test_that("binary moderate calibration agrees with an independent RCS oracle", {
  scored <- make_binary_moderate_scored()
  result <- sb_binary_moderate_calibration_metric(
    scored,
    "pooled_patient_level_out_of_fold"
  )

  logit_probability <- stats::qlogis(scored$probability)
  knots <- unname(stats::quantile(
    logit_probability,
    c(0.10, 0.50, 0.90),
    names = FALSE,
    type = 7
  ))
  manual_basis <- manual_binary_rcs3_basis(logit_probability, knots)
  manual_design <- cbind("(Intercept)" = 1, manual_basis)
  manual_fit <- stats::glm.fit(
    x = manual_design,
    y = scored$outcome,
    family = stats::binomial()
  )
  manual_probability <- stats::plogis(drop(
    manual_design %*% manual_fit$coefficients
  ))
  manual_error <- abs(manual_probability - scored$probability)

  expect_identical(result$status, "success")
  expect_identical(result$input_status, "success")
  expect_identical(result$metric, "binary_calibration_curve")
  expect_identical(result$calibration_level, "moderate")
  expect_identical(
    result$method,
    "logistic_three_knot_restricted_natural_cubic_spline"
  )
  expect_identical(result$prediction_transform, "logit")
  expect_identical(result$probability_clipping, "none")
  expect_identical(result$spline$n_knots, 3L)
  expect_equal(result$knots$logit_probability, knots, tolerance = 1e-14)
  expect_equal(result$spline$boundary_knots, knots[c(1L, 3L)])
  expect_equal(result$spline$interior_knots, knots[[2L]])
  expect_identical(
    result$patient_calibration$patient_id,
    scored$patient_id
  )
  expect_identical(
    result$patient_calibration$row_instance,
    seq_len(nrow(scored))
  )
  expect_equal(
    result$patient_calibration$calibrated_probability,
    manual_probability,
    tolerance = 1e-10
  )
  expect_equal(result$ici, mean(manual_error), tolerance = 1e-10)
  expect_equal(
    result$e50,
    unname(stats::quantile(manual_error, 0.50, names = FALSE, type = 7)),
    tolerance = 1e-10
  )
  expect_equal(
    result$e90,
    unname(stats::quantile(manual_error, 0.90, names = FALSE, type = 7)),
    tolerance = 1e-10
  )
  expect_identical(result$fit$terms, result$fit$coefficients$term)
  expect_identical(result$fit$rank, 3L)
  expect_identical(result$fit$converged, TRUE)
  expect_identical(result$fit$boundary, FALSE)
  expect_identical(result$fit$warnings, character())
  expect_identical(result$fit$separation_status, "not_detected")
  expect_identical(result$fit$separation_detected, FALSE)
  expect_identical(result$per_fold$status, rep("success", 4L))
})

test_that("saved audit alone reconstructs the binary calibration evidence", {
  scored <- make_binary_moderate_scored()
  result <- sb_binary_moderate_calibration_metric(
    scored,
    "pooled_patient_level_out_of_fold"
  )
  saved_knots <- result$knots$logit_probability
  saved_basis <- splines::ns(
    result$patient_calibration$logit_probability,
    knots = result$spline$interior_knots,
    Boundary.knots = result$spline$boundary_knots,
    intercept = result$spline$intercept
  )
  colnames(saved_basis) <- c("calibration_ns1", "calibration_ns2")
  saved_design <- cbind("(Intercept)" = 1, saved_basis)
  saved_coefficients <- stats::setNames(
    result$fit$coefficients$estimate,
    result$fit$coefficients$term
  )[result$fit$terms]
  reconstructed <- stats::plogis(drop(saved_design %*% saved_coefficients))
  reconstructed_error <- abs(reconstructed - scored$probability)
  expected_curve <- data.frame(
    probability = sort(unique(scored$probability)),
    stringsAsFactors = FALSE
  )
  first <- match(expected_curve$probability, scored$probability)
  expected_curve$logit_probability <- stats::qlogis(
    expected_curve$probability
  )
  expected_curve$calibrated_probability <- reconstructed[first]
  expected_curve$n_row_instances <- as.integer(tabulate(
    match(scored$probability, expected_curve$probability),
    nbins = nrow(expected_curve)
  ))

  expect_equal(saved_knots, c(
    result$spline$boundary_knots[[1L]],
    result$spline$interior_knots[[1L]],
    result$spline$boundary_knots[[2L]]
  ))
  expect_equal(
    reconstructed,
    result$patient_calibration$calibrated_probability,
    tolerance = 1e-14
  )
  expect_equal(result$curve, expected_curve, tolerance = 1e-14)
  expect_equal(result$ici, mean(reconstructed_error), tolerance = 1e-14)
  expect_equal(
    result$e50,
    unname(stats::quantile(
      reconstructed_error,
      0.50,
      names = FALSE,
      type = 7
    )),
    tolerance = 1e-14
  )
  expect_equal(
    result$e90,
    unname(stats::quantile(
      reconstructed_error,
      0.90,
      names = FALSE,
      type = 7
    )),
    tolerance = 1e-14
  )
})

test_that("binary moderate calibration preserves ties rows and pooled meaning", {
  scored <- make_binary_moderate_scored()
  scored$probability[c(3L, 4L)] <- scored$probability[[2L]]
  scored$patient_id[c(7L, 8L)] <- "resampled_patient"
  result <- sb_binary_moderate_calibration_metric(
    scored,
    "patient_level_external_bootstrap_sample"
  )

  expect_identical(result$status, "success")
  expect_identical(nrow(result$patient_calibration), nrow(scored))
  expect_identical(
    sum(result$patient_calibration$patient_id == "resampled_patient"),
    2L
  )
  expect_identical(
    sum(result$curve$n_row_instances),
    as.integer(nrow(scored))
  )
  expect_identical(max(result$curve$n_row_instances), 3L)

  one_class_fold <- make_binary_moderate_scored()
  one_class_fold$outcome[one_class_fold$fold == 1L] <- 0L
  pooled <- sb_binary_moderate_calibration_metric(
    one_class_fold,
    "pooled_patient_level_out_of_fold"
  )
  expect_identical(pooled$status, "success")
  expect_identical(pooled$per_fold$status, rep("success", 4L))
  expect_identical(pooled$per_fold$n_positive[[1L]], 0L)
})

test_that("binary moderate calibration is invariant to patient row order", {
  scored <- make_binary_moderate_scored()
  original <- sb_binary_moderate_calibration_metric(
    scored,
    "patient_level_external_cohort"
  )
  order <- c(seq(2L, nrow(scored), by = 2L), seq(1L, nrow(scored), by = 2L))
  shuffled_scored <- scored[order, , drop = FALSE]
  shuffled <- sb_binary_moderate_calibration_metric(
    shuffled_scored,
    "patient_level_external_cohort"
  )

  expect_identical(shuffled$status, "success")
  expect_equal(shuffled[c("ici", "e50", "e90")], original[c(
    "ici", "e50", "e90"
  )], tolerance = 1e-13)
  expect_equal(shuffled$knots, original$knots, tolerance = 1e-13)
  expect_equal(shuffled$curve, original$curve, tolerance = 1e-13)
  aligned <- shuffled$patient_calibration[
    match(original$patient_calibration$patient_id,
      shuffled$patient_calibration$patient_id
    ),
    ,
    drop = FALSE
  ]
  expect_equal(
    aligned$calibrated_probability,
    original$patient_calibration$calibrated_probability,
    tolerance = 1e-13
  )
})

test_that("binary moderate calibration retains concrete unavailable states", {
  base <- make_binary_moderate_scored()

  failed <- base
  failed$status[[1L]] <- "prediction_failed"
  failed_result <- sb_binary_moderate_calibration_metric(failed, "test")
  expect_identical(failed_result$status, "incomplete")
  expect_identical(failed_result$input_status, "prediction_failed")

  one_class <- base
  one_class$outcome <- 0L
  class_result <- sb_binary_moderate_calibration_metric(one_class, "test")
  expect_identical(class_result$status, "incomplete")
  expect_identical(class_result$input_status, "not_evaluable")
  expect_match(class_result$error_message, "Both outcome classes", fixed = TRUE)

  boundary <- base
  boundary$probability[[1L]] <- 0
  boundary_result <- sb_binary_moderate_calibration_metric(boundary, "test")
  expect_identical(boundary_result$status, "incomplete")
  expect_identical(boundary_result$input_status, "not_evaluable")
  expect_identical(boundary_result$n_boundary_probability, 1L)
  expect_identical(boundary_result$probability_clipping, "none")

  tied <- base
  tied$probability <- 0.5
  tied_result <- sb_binary_moderate_calibration_metric(tied, "test")
  expect_identical(tied_result$status, "incomplete")
  expect_identical(tied_result$input_status, "not_evaluable")
  expect_match(tied_result$error_message, "strictly increasing", fixed = TRUE)
})

test_that("extreme finite probabilities are neither clipped nor rejected", {
  probability <- c(
    1e-300, 1e-20, 1e-8, 1e-4, 0.01, 0.05, 0.15, 0.30,
    0.50, 0.70, 0.85, 0.95, 0.99, 1 - 1e-4, 1 - 1e-8,
    1 - .Machine$double.eps
  )
  scored <- data.frame(
    patient_id = paste0("extreme_", seq_along(probability)),
    outcome = rep(c(0L, 1L), length.out = length(probability)),
    probability = probability,
    status = rep("success", length(probability)),
    stringsAsFactors = FALSE
  )
  result <- sb_binary_moderate_calibration_metric(scored, "test")

  expect_identical(result$status, "success")
  expect_identical(result$n_boundary_probability, 0L)
  expect_identical(
    result$patient_calibration$probability,
    probability
  )
  expect_equal(
    result$patient_calibration$logit_probability,
    stats::qlogis(probability),
    tolerance = 0
  )
  expect_identical(
    all(is.finite(result$patient_calibration$logit_probability)),
    TRUE
  )
  expect_identical(result$probability_clipping, "none")
  expect_identical(all(is.finite(c(result$ici, result$e50, result$e90))), TRUE)
})

test_that("probability warnings retain an estimable calibration fit", {
  scored <- make_binary_moderate_scored()
  expected <- sb_binary_moderate_calibration_metric(scored, "test")
  native_fit <- stats::glm.fit
  probability_warning <-
    "glm.fit: fitted probabilities numerically 0 or 1 occurred"
  local_mocked_bindings(
    glm.fit = function(...) {
      fit <- native_fit(...)
      warning(probability_warning, call. = FALSE)
      fit
    },
    .package = "stats"
  )

  result <- sb_binary_moderate_calibration_metric(scored, "test")

  expect_identical(result$status, "success")
  expect_identical(result$fit$status, "success")
  expect_identical(result$fit$warnings, probability_warning)
  expect_identical(result$fit$converged, TRUE)
  expect_identical(result$fit$rank, 3L)
  expect_identical(result$fit$separation_status, "not_detected")
  expect_identical(result$fit$separation_detected, FALSE)
  expect_equal(result$fit$coefficients, expected$fit$coefficients, tolerance = 0)
  expect_equal(result$patient_calibration, expected$patient_calibration,
    tolerance = 0
  )
  expect_equal(result[c("ici", "e50", "e90")],
    expected[c("ici", "e50", "e90")], tolerance = 0
  )
})

test_that("calibration warnings do not rescue invalid native fits", {
  scored <- make_binary_moderate_scored()
  native_fit <- stats::glm.fit
  probability_warning <-
    "glm.fit: fitted probabilities numerically 0 or 1 occurred"
  cases <- list(
    list(warning = "unrelated native warning", changes = list(),
      error = "glm_warning"),
    list(warning = "glm.fit: algorithm did not converge", changes = list(),
      error = "glm_warning"),
    list(warning = probability_warning, changes = list(converged = FALSE),
      error = "nonconvergence"),
    list(warning = probability_warning, changes = list(boundary = TRUE),
      error = "boundary_fit"),
    list(warning = probability_warning, changes = list(rank = 2L),
      error = "singular_fit")
  )
  for (case in cases) {
    result <- with_mocked_bindings(
      sb_binary_moderate_calibration_metric(scored, "test"),
      glm.fit = function(...) {
        fit <- native_fit(...)
        fit[names(case$changes)] <- case$changes
        warning(case$warning, call. = FALSE)
        fit
      },
      .package = "stats"
    )
    expect_identical(result$status, "fit_failed")
    expect_identical(result$fit$status, "failed")
    expect_identical(result$fit$warnings, case$warning)
    expect_identical(result$error_class, paste0(
      "stormbreaker_binary_moderate_calibration_", case$error
    ))
    expect_identical(is.na(result$ici), TRUE)
  }
})

test_that("a separated spline fit remains visible without fallback", {
  probability <- seq(0.03, 0.97, length.out = 40L)
  scored <- data.frame(
    patient_id = paste0("separated_", seq_along(probability)),
    outcome = as.integer(probability >= 0.5),
    probability = probability,
    status = rep("success", length(probability)),
    stringsAsFactors = FALSE
  )
  result <- sb_binary_moderate_calibration_metric(scored, "test")

  expect_identical(result$status, "fit_failed")
  expect_identical(result$input_status, "success")
  expect_identical(result$fit$status, "failed")
  expect_gt(length(result$fit$warnings), 0L)
  expect_true(is.na(result$ici))
  expect_identical(result$method, "logistic_three_knot_restricted_natural_cubic_spline")
})

test_that("silent finite IRLS output cannot hide nonlinear separation", {
  probability <- c(
    0.013529309201985598,
    0.232556497720070199,
    0.292985442220233372,
    0.376448372132144882,
    0.875341706429608180,
    0.887932648812420622
  )
  scored <- data.frame(
    patient_id = paste0("silent_separation_", seq_along(probability)),
    outcome = c(1L, 0L, 0L, 0L, 1L, 1L),
    probability = probability,
    status = rep("success", length(probability)),
    stringsAsFactors = FALSE
  )

  result <- sb_binary_moderate_calibration_metric(scored, "test")

  expect_identical(result$status, "fit_failed")
  expect_identical(result$input_status, "success")
  expect_identical(result$fit$status, "failed")
  expect_identical(result$fit$converged, TRUE)
  expect_identical(result$fit$boundary, FALSE)
  expect_identical(result$fit$warnings, character())
  expect_identical(result$fit$separation_status, "detected")
  expect_identical(result$fit$separation_detected, TRUE)
  expect_identical(result$fit$separation_solver_status, 0L)
  expect_gt(
    result$fit$separation_objective,
    result$fit$separation_tolerance
  )
  expect_gte(
    result$fit$separation_minimum_signed_margin,
    -result$fit$separation_tolerance
  )
  expect_gt(
    result$fit$separation_maximum_signed_margin,
    result$fit$separation_tolerance
  )
  expect_identical(
    result$error_class,
    "stormbreaker_binary_moderate_calibration_separation"
  )
  expect_identical(all(is.finite(result$fit$coefficients$estimate)), TRUE)
  expect_identical(is.na(result$ici), TRUE)
})

test_that("the signed-margin program detects quasi-complete separation", {
  probability <- c(0.05, 0.20, 0.50, 0.50, 0.80, 0.95)
  outcome <- c(0L, 0L, 0L, 1L, 1L, 1L)
  logit_probability <- stats::qlogis(probability)
  knots <- sb_binary_moderate_calibration_knots(
    logit_probability
  )$logit_probability
  basis <- sb_binary_moderate_calibration_basis(logit_probability, knots)
  design <- cbind("(Intercept)" = 1, basis)

  result <- sb_binary_moderate_calibration_separation(design, outcome)

  expect_identical(result$status, "detected")
  expect_identical(result$detected, TRUE)
  expect_identical(result$solver_status, 0L)
  expect_gt(result$objective, result$tolerance)
  expect_lte(abs(result$minimum_signed_margin), result$tolerance)
  expect_gt(result$maximum_signed_margin, result$tolerance)
})
