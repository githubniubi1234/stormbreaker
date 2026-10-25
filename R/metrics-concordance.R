sb_harrell_c_stratum <- function(time, event, risk) {
  data <- data.frame(
    time = time,
    event = event,
    risk = risk
  )
  result <- list(
    status = "not_evaluable",
    estimate = NA_real_,
    n = as.integer(length(time)),
    n_events = as.integer(sum(event)),
    concordant = 0,
    discordant = 0,
    tied_risk = 0,
    tied_time = 0,
    tied_both = 0,
    comparable_pairs = 0,
    error_class = NA_character_,
    error_message = "No comparable patient pairs."
  )
  concordance <- tryCatch(
    survival::concordance(
      survival::Surv(time, event) ~ risk,
      data = data,
      reverse = TRUE,
      timewt = "n"
    ),
    error = function(error) error
  )
  if (inherits(concordance, "error")) {
    result$status <- "metric_failed"
    result$error_class <- class(concordance)[1L]
    result$error_message <- conditionMessage(concordance)
    return(result)
  }

  counts <- concordance$count
  result$concordant <- unname(counts[["concordant"]])
  result$discordant <- unname(counts[["discordant"]])
  result$tied_risk <- unname(counts[["tied.x"]])
  result$tied_time <- unname(counts[["tied.y"]])
  result$tied_both <- unname(counts[["tied.xy"]])
  result$comparable_pairs <- unname(
    result$concordant + result$discordant + result$tied_risk
  )
  if (!is.finite(result$comparable_pairs) ||
      result$comparable_pairs <= 0) {
    return(result)
  }

  result$status <- "success"
  result$estimate <- unname(concordance$concordance)
  result$error_message <- NA_character_
  result
}

sb_harrell_c_metric <- function(predictions) {
  fold_ids <- sort(unique(as.integer(predictions$fold)))
  per_fold <- lapply(
    fold_ids,
    function(fold_id) {
      fold_predictions <- predictions[
        predictions$fold == fold_id,
        c("time", "event", "risk", "status"),
        drop = FALSE
      ]
      audit <- data.frame(
        fold = fold_id,
        status = "prediction_failed",
        n = nrow(fold_predictions),
        n_events = as.integer(sum(fold_predictions$event)),
        estimate = NA_real_,
        concordant = 0,
        discordant = 0,
        tied_risk = 0,
        tied_time = 0,
        tied_both = 0,
        comparable_pairs = 0,
        error_class = NA_character_,
        error_message = NA_character_,
        stringsAsFactors = FALSE
      )
      prediction_failed <- anyNA(fold_predictions$status) ||
        any(fold_predictions$status != "success") ||
        anyNA(fold_predictions$risk) ||
        any(!is.finite(fold_predictions$risk))
      if (prediction_failed) {
        audit$error_message <- paste0(
          "OOF predictions for this fold are failed or non-finite."
        )
        return(audit)
      }

      stratum <- sb_harrell_c_stratum(
        time = fold_predictions$time,
        event = fold_predictions$event,
        risk = fold_predictions$risk
      )
      audit$n <- stratum$n
      audit$n_events <- stratum$n_events
      audit$concordant <- stratum$concordant
      audit$discordant <- stratum$discordant
      audit$tied_risk <- stratum$tied_risk
      audit$tied_time <- stratum$tied_time
      audit$tied_both <- stratum$tied_both
      audit$comparable_pairs <- stratum$comparable_pairs
      if (identical(stratum$status, "metric_failed")) {
        audit$status <- "metric_failed"
        audit$error_class <- stratum$error_class
        audit$error_message <- stratum$error_message
        return(audit)
      }
      if (identical(stratum$status, "not_evaluable")) {
        audit$status <- "not_evaluable"
        audit$error_message <- "No comparable patient pairs in this fold."
        return(audit)
      }

      audit$status <- "success"
      audit$estimate <- stratum$estimate
      audit
    }
  )
  per_fold <- do.call(rbind, per_fold)
  rownames(per_fold) <- NULL
  evaluable <- per_fold$status == "success"
  count_columns <- c(
    "concordant", "discordant", "tied_risk", "tied_time",
    "tied_both", "comparable_pairs"
  )
  totals <- stats::setNames(numeric(length(count_columns)), count_columns)
  if (any(evaluable)) {
    totals <- colSums(per_fold[evaluable, count_columns, drop = FALSE])
  }
  complete <- length(fold_ids) > 0L && all(evaluable)
  estimate <- if (complete) {
    unname(
      (totals["concordant"] + 0.5 * totals["tied_risk"]) /
        totals["comparable_pairs"]
    )
  } else {
    NA_real_
  }

  list(
    metric = "harrell_c",
    scope = "within_outer_assessment_fold",
    aggregation = "summed_comparable_pairs",
    time_weight = "n",
    risk_direction = "higher_is_higher_hazard",
    status = if (complete) "success" else "incomplete",
    estimate = estimate,
    n_folds = as.integer(length(fold_ids)),
    n_evaluable_folds = as.integer(sum(evaluable)),
    concordant = unname(totals["concordant"]),
    discordant = unname(totals["discordant"]),
    tied_risk = unname(totals["tied_risk"]),
    tied_time = unname(totals["tied_time"]),
    tied_both = unname(totals["tied_both"]),
    comparable_pairs = unname(totals["comparable_pairs"]),
    per_fold = per_fold
  )
}
