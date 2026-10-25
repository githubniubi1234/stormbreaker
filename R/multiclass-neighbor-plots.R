sb_multiclass_neighbor_plots <- function(fit, new_data) {
  if (is.null(new_data)) stop("Supply new_data for prediction-neighborhood plots.", call. = FALSE)
  x <- sb_validate_predictor_frame(new_data[fit$features], "new_data")
  m <- fit$engine
  if (fit$method == "knn") {
    z <- sweep(sweep(as.matrix(x), 2, fit$center, "-"), 2, fit$scale, "/")
    rows <- lapply(seq_len(nrow(z)), function(i) {
      distance <- rowSums(sweep(m$training, 2, z[i, ], "-")^2)
      index <- order(distance, seq_along(distance))[seq_len(m$k)]
      data.frame(query = i, rank = seq_along(index), reference_row = index,
        distance = sqrt(distance[index]), class = as.character(m$outcome[index]),
        weight = 1 / m$k)
    })
    tab <- do.call(rbind, rows)
  } else {
    names(x) <- fit$feature_codes
    p <- fit$parameters
    native <- kknn::kknn(stats::as.formula(".outcome ~ .", env = baseenv()),
      train = m$training, test = x, k = p$k, distance = p$distance,
      kernel = p$kernel, scale = p$scale)
    tab <- data.frame(query = rep(seq_len(nrow(x)), each = ncol(native$C)),
      rank = rep(seq_len(ncol(native$C)), nrow(x)), reference_row = as.integer(t(native$C)),
      distance = as.numeric(t(native$D)), class = as.character(t(native$CL)),
      weight = as.numeric(t(native$W / rowSums(native$W))))
  }
  tab$class <- fit$classes[match(tab$class, fit$codes)]
  rank <- distance <- class <- weight <- NULL
  p <- ggplot2::ggplot(tab, ggplot2::aes(rank, distance, colour = class, size = weight)) +
    ggplot2::geom_point(alpha = .8) + ggplot2::facet_wrap(~query) +
    ggplot2::scale_x_continuous(breaks = function(z) unique(round(pretty(z)))) +
    sb_viz_theme() + ggplot2::labs(title = "Prediction neighborhoods", x = "Neighbor rank",
      y = "Distance", colour = "Class", size = "Voting weight",
      caption = "Original neighbor selection, distances and normalized class-voting weights for the supplied queries. Training queries retain self-neighbors; these are explanations, not held-out or leave-one-out performance. All queries and selected neighbors are retained across pages. Reference-row identities and classes may be sensitive.")
  list(plots = list(neighborhoods = p), tables = list(neighborhoods = tab))
}
