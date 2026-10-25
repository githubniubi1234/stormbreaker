sb_score_univariable_binary_auc_features <- function(x, outcome) {
  n_positive <- as.integer(sum(outcome == 1L))
  n_negative <- as.integer(sum(outcome == 0L))
  comparable_pairs <- as.numeric(n_positive) * as.numeric(n_negative)
  n_features <- ncol(x)
  auc <- rep(NA_real_, n_features)
  discrimination <- rep(NA_real_, n_features)
  direction <- rep(NA_character_, n_features)
  status <- rep("failed", n_features)
  error_class <- rep(NA_character_, n_features)
  error_message <- rep(NA_character_, n_features)

  for (index in seq_len(n_features)) {
    feature <- names(x)[index]
    values <- x[[index]]
    if (min(values) == max(values)) {
      error_class[index] <- "stormbreaker_selector_constant_feature"
      error_message[index] <- paste0(
        "Binary AUC scoring cannot rank constant predictor `",
        feature, "`."
      )
      next
    }

    ranks <- rank(values, ties.method = "average")
    pair_credit <- sum(ranks[outcome == 1L]) -
      n_positive * (n_positive + 1) / 2
    feature_auc <- unname(pair_credit / comparable_pairs)
    if (!is.finite(feature_auc) || feature_auc < 0 || feature_auc > 1) {
      error_class[index] <- "stormbreaker_invalid_selector_score"
      error_message[index] <- paste0(
        "Binary AUC scoring produced an invalid value for `",
        feature, "`."
      )
      next
    }

    auc[index] <- feature_auc
    discrimination[index] <- max(
      pair_credit,
      comparable_pairs - pair_credit
    ) / comparable_pairs
    direction[index] <- if (feature_auc > 0.5) {
      "higher_in_positive"
    } else if (feature_auc < 0.5) {
      "lower_in_positive"
    } else {
      "tied"
    }
    status[index] <- "success"
  }

  data.frame(
    feature = names(x),
    original_order = seq_len(n_features),
    auc = auc,
    discrimination = discrimination,
    direction = direction,
    n_positive = rep.int(n_positive, n_features),
    n_negative = rep.int(n_negative, n_features),
    comparable_pairs = rep.int(comparable_pairs, n_features),
    status = status,
    error_class = error_class,
    error_message = error_message,
    stringsAsFactors = FALSE
  )
}

sb_fit_univ_binary_auc_top_k <- function(x, outcome, n_features) {
  x <- sb_validate_predictor_frame(x, "x")
  outcome <- sb_validate_binary_outcome_vector(outcome, nrow(x))
  if (!all(c(0L, 1L) %in% outcome)) {
    sb_abort(
      "Binary AUC selection requires both outcome classes.",
      "stormbreaker_selector_single_class"
    )
  }
  n_features <- sb_check_positive_integer(
    n_features,
    "n_features",
    minimum = 1L
  )
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

  scores <- sb_score_univariable_binary_auc_features(x, outcome)
  scores$rank <- NA_integer_
  scores$selected <- FALSE
  valid <- which(
    scores$status == "success" & is.finite(scores$discrimination)
  )
  ordered_valid <- valid[order(
    -scores$discrimination[valid],
    scores$original_order[valid]
  )]
  if (length(ordered_valid) < n_features) {
    sb_abort(
      paste0(
        "The selector requested ", n_features,
        " features, but only ", length(ordered_valid),
        " features had valid binary AUC scores."
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
      id = "univ_binary_auc_top_k",
      n_features = n_features,
      source_features = names(x),
      selected_features = selected_features,
      scores = scores
    ),
    class = "stormbreaker_selector_fit"
  )
}

sb_apply_univ_binary_auc_top_k <- function(fitted, x) {
  valid_fit <- inherits(fitted, "stormbreaker_selector_fit") &&
    identical(fitted$id, "univ_binary_auc_top_k") &&
    is.character(fitted$source_features) &&
    is.character(fitted$selected_features) &&
    length(fitted$selected_features) == fitted$n_features
  if (!valid_fit) {
    sb_abort(
      "The fitted selector is not a valid binary AUC Top-K state.",
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

sb_univ_binary_auc_top_k_selector_adapter <- function(n_features) {
  n_features <- sb_check_positive_integer(
    n_features,
    "n_features",
    minimum = 1L
  )
  list(
    id = "univ_binary_auc_top_k",
    params = list(n_features = n_features),
    needs_outcome = TRUE,
    fit = function(x, outcome) {
      sb_fit_univ_binary_auc_top_k(x, outcome, n_features)
    },
    apply = sb_apply_univ_binary_auc_top_k
  )
}
