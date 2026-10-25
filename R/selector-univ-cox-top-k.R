sb_score_univariable_cox_feature <- function(feature, x, survival_response,
                                              weights, control,
                                              original_order) {
  empty <- data.frame(
    feature = feature,
    original_order = as.integer(original_order),
    coefficient = NA_real_,
    standard_error = NA_real_,
    wald_z = NA_real_,
    p_value = NA_real_,
    status = "failed",
    error_class = NA_character_,
    error_message = NA_character_,
    stringsAsFactors = FALSE
  )
  tryCatch(
    {
      values <- x[[feature]]
      n_above_minimum <- sum(values > min(values))
      n_below_maximum <- sum(values < max(values))
      if (n_above_minimum < 2L || n_below_maximum < 2L) {
        sb_abort(
          paste0(
            "Univariable Cox scoring cannot use `", feature,
            "`: fewer than two patients lie away from one expression ",
            "extreme, so a single patient could define the feature score."
          ),
          "stormbreaker_selector_insufficient_feature_support"
        )
      }
      model <- withCallingHandlers(
        survival::coxph.fit(
          x = matrix(values, ncol = 1L),
          y = survival_response,
          strata = NULL,
          offset = NULL,
          init = NULL,
          control = control,
          weights = weights,
          method = "efron",
          rownames = NULL,
          resid = FALSE,
          nocenter = c(-1, 0, 1)
        ),
        warning = function(warning) {
          sb_abort(
            paste0(
              "Univariable Cox scoring raised a warning for `",
              feature, "`: ", conditionMessage(warning)
            ),
            "stormbreaker_selector_cox_warning",
            parent = warning
          )
        }
      )
      coefficient <- unname(model$coefficients[[1L]])
      standard_error <- unname(sqrt(model$var[[1L]]))
      if (!is.finite(coefficient) || !is.finite(standard_error) ||
          standard_error <= 0) {
        sb_abort(
          paste0(
            "Univariable Cox scoring produced invalid statistics for `",
            feature, "`."
          ),
          "stormbreaker_invalid_selector_score"
        )
      }
      wald_z <- coefficient / standard_error
      p_value <- 2 * stats::pnorm(abs(wald_z), lower.tail = FALSE)
      if (!is.finite(wald_z) || !is.finite(p_value)) {
        sb_abort(
          paste0(
            "Univariable Cox scoring produced a non-finite Wald score for `",
            feature, "`."
          ),
          "stormbreaker_invalid_selector_score"
        )
      }
      empty$coefficient <- coefficient
      empty$standard_error <- standard_error
      empty$wald_z <- wald_z
      empty$p_value <- p_value
      empty$status <- "success"
      empty
    },
    error = function(error) {
      empty$error_class <- class(error)[1]
      empty$error_message <- conditionMessage(error)
      empty
    }
  )
}

sb_fit_univ_cox_top_k <- function(x, time, event, n_features) {
  x <- sb_validate_predictor_frame(x, "x")
  time <- sb_validate_endpoint_vector(time, "time", nrow(x))
  event <- sb_validate_endpoint_vector(event, "event", nrow(x))
  if (any(time < 0)) {
    sb_abort(
      "`time` must contain non-negative values.",
      "stormbreaker_bad_adapter_time"
    )
  }
  if (!all(event %in% c(0, 1))) {
    sb_abort(
      "`event` must contain only 0/1 values.",
      "stormbreaker_bad_adapter_event"
    )
  }
  if (n_features > ncol(x)) {
    sb_abort(
      paste0(
        "The selector requested ", n_features,
        " features, but the task has only ", ncol(x),
        " candidate features."
      ),
      "stormbreaker_insufficient_candidate_features"
    )
  }

  survival_response <- survival::Surv(time, event)
  weights <- rep.int(1, nrow(x))
  control <- survival::coxph.control()

  scores <- do.call(
    rbind,
    lapply(seq_along(x), function(index) {
      sb_score_univariable_cox_feature(
        names(x)[index],
        x,
        survival_response,
        weights,
        control,
        index
      )
    })
  )
  rownames(scores) <- NULL
  scores$rank <- NA_integer_
  scores$selected <- FALSE
  valid <- which(scores$status == "success" & is.finite(scores$p_value))
  ordered_valid <- valid[order(
    scores$p_value[valid],
    scores$original_order[valid]
  )]
  if (length(ordered_valid) < n_features) {
    sb_abort(
      paste0(
        "The selector requested ", n_features,
        " features, but only ", length(ordered_valid),
        " features had valid univariable Cox scores."
      ),
      "stormbreaker_insufficient_valid_features"
    )
  }
  scores$rank[ordered_valid] <- seq_along(ordered_valid)
  selected_rows <- ordered_valid[seq_len(n_features)]
  scores$selected[selected_rows] <- TRUE
  selected_features <- scores$feature[selected_rows]

  structure(
    list(
      id = "univ_cox_top_k",
      n_features = n_features,
      source_features = names(x),
      selected_features = selected_features,
      scores = scores
    ),
    class = "stormbreaker_selector_fit"
  )
}

sb_apply_univ_cox_top_k <- function(fitted, x) {
  if (!inherits(fitted, "stormbreaker_selector_fit") ||
      !identical(fitted$id, "univ_cox_top_k") ||
      !is.character(fitted$source_features) ||
      !is.character(fitted$selected_features) ||
      length(fitted$selected_features) != fitted$n_features) {
    sb_abort(
      "The fitted selector is not a valid univariable Cox Top-K state.",
      "stormbreaker_bad_selector_fit"
    )
  }
  x <- sb_validate_predictor_frame(x, "x")
  if (!identical(names(x), fitted$source_features)) {
    sb_abort(
      paste0(
        "Predictor columns must exactly match the fitted selector schema: ",
        paste(fitted$source_features, collapse = ", "),
        "."
      ),
      "stormbreaker_selector_schema_mismatch"
    )
  }
  x[fitted$selected_features]
}

sb_univ_cox_top_k_selector_adapter <- function(n_features) {
  n_features <- sb_check_positive_integer(
    n_features,
    "n_features",
    minimum = 1L
  )
  list(
    id = "univ_cox_top_k",
    params = list(n_features = n_features),
    needs_outcome = TRUE,
    fit = function(x, time, event) {
      sb_fit_univ_cox_top_k(x, time, event, n_features)
    },
    apply = sb_apply_univ_cox_top_k
  )
}
