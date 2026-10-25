sb_empty_selection_pairwise <- function() {
  data.frame(
    fold_i = integer(),
    fold_j = integer(),
    status = character(),
    n_i = integer(),
    n_j = integer(),
    intersection = integer(),
    union = integer(),
    jaccard = numeric(),
    stringsAsFactors = FALSE
  )
}

sb_selection_stability <- function(selector, selections, candidate_features) {
  fold_ids <- vapply(
    selections,
    function(selection) as.integer(selection$fold[[1L]]),
    integer(1)
  )
  ordered <- order(fold_ids)
  fold_ids <- fold_ids[ordered]
  selections <- selections[ordered]
  n_folds <- as.integer(length(selections))
  not_applicable <- selector %in% c("none", "fixed_features")

  saved_success <- vapply(
    selections,
    function(selection) identical(selection$status, "success"),
    logical(1)
  )
  successful <- if (not_applicable) {
    rep(FALSE, n_folds)
  } else {
    saved_success
  }
  selected_sets <- lapply(
    selections,
    function(selection) unique(as.character(selection$selected_features))
  )

  error_field <- function(selection, field) {
    value <- selection[[field]]
    if (is.null(value) || length(value) == 0L || is.na(value[[1L]])) {
      return(NA_character_)
    }
    as.character(value[[1L]])
  }
  per_fold <- data.frame(
    fold = as.integer(fold_ids),
    status = if (not_applicable) {
      rep("not_applicable", n_folds)
    } else {
      vapply(selections, function(selection) selection$status, character(1))
    },
    n_selected = as.integer(vapply(seq_along(selections), function(index) {
      if (successful[[index]]) length(selected_sets[[index]]) else NA_integer_
    }, integer(1))),
    error_class = vapply(
      selections,
      error_field,
      character(1),
      field = "error_class"
    ),
    error_message = vapply(
      selections,
      error_field,
      character(1),
      field = "error_message"
    ),
    stringsAsFactors = FALSE
  )
  if (not_applicable) {
    per_fold$error_class <- NA_character_
    per_fold$error_message <- NA_character_
  }

  n_successful <- as.integer(sum(successful))
  selected_count <- vapply(candidate_features, function(feature) {
    as.integer(sum(vapply(
      selected_sets[successful],
      function(selected) feature %in% selected,
      logical(1)
    )))
  }, integer(1))
  feature_frequency <- data.frame(
    feature = as.character(candidate_features),
    selected_count = as.integer(selected_count),
    selection_frequency = if (n_successful > 0L) {
      unname(selected_count / n_successful)
    } else {
      rep(NA_real_, length(candidate_features))
    },
    stringsAsFactors = FALSE
  )

  pair_indices <- if (n_folds >= 2L) {
    t(utils::combn(seq_len(n_folds), 2L))
  } else {
    matrix(integer(), ncol = 2L)
  }
  pairwise <- if (nrow(pair_indices) == 0L) {
    sb_empty_selection_pairwise()
  } else {
    rows <- lapply(seq_len(nrow(pair_indices)), function(row) {
      index_i <- pair_indices[row, 1L]
      index_j <- pair_indices[row, 2L]
      status <- if (not_applicable) {
        "not_applicable"
      } else if (successful[[index_i]] && successful[[index_j]]) {
        "success"
      } else {
        "incomplete"
      }
      n_i <- per_fold$n_selected[[index_i]]
      n_j <- per_fold$n_selected[[index_j]]
      intersection_size <- NA_integer_
      union_size <- NA_integer_
      jaccard <- NA_real_
      if (identical(status, "success")) {
        intersection_size <- as.integer(length(intersect(
          selected_sets[[index_i]],
          selected_sets[[index_j]]
        )))
        union_size <- as.integer(length(union(
          selected_sets[[index_i]],
          selected_sets[[index_j]]
        )))
        if (union_size > 0L) {
          jaccard <- unname(intersection_size / union_size)
        } else {
          status <- "not_evaluable"
        }
      }
      data.frame(
        fold_i = as.integer(fold_ids[[index_i]]),
        fold_j = as.integer(fold_ids[[index_j]]),
        status = status,
        n_i = as.integer(n_i),
        n_j = as.integer(n_j),
        intersection = intersection_size,
        union = union_size,
        jaccard = jaccard,
        stringsAsFactors = FALSE
      )
    })
    combined <- do.call(rbind, rows)
    rownames(combined) <- NULL
    combined
  }

  n_pairs <- as.integer(nrow(pairwise))
  n_evaluable_pairs <- as.integer(sum(pairwise$status == "success"))
  complete <- !not_applicable && n_folds >= 2L && all(successful) &&
    n_evaluable_pairs == n_pairs
  status <- if (not_applicable) {
    "not_applicable"
  } else if (complete) {
    "success"
  } else {
    "incomplete"
  }
  error_message <- if (not_applicable || complete) {
    NA_character_
  } else if (n_folds < 2L) {
    "At least two outer-fold selections are required for stability."
  } else if (!all(successful)) {
    paste0(
      "Selection stability is incomplete because not every outer-fold ",
      "selection succeeded."
    )
  } else {
    "Selection stability is incomplete because a fold pair was not evaluable."
  }

  list(
    selector = selector,
    status = status,
    scope = "outer_analysis_fold_selected_feature_sets",
    aggregation = "mean_of_unique_fold_pair_jaccard",
    frequency_denominator = "successful_selection_folds",
    n_candidate_features = as.integer(length(candidate_features)),
    n_folds = n_folds,
    n_successful_folds = n_successful,
    n_pairs = n_pairs,
    n_evaluable_pairs = n_evaluable_pairs,
    mean_pairwise_jaccard = if (complete) {
      unname(mean(pairwise$jaccard))
    } else {
      NA_real_
    },
    feature_frequency = feature_frequency,
    per_fold = per_fold,
    pairwise = pairwise,
    error_message = error_message
  )
}
