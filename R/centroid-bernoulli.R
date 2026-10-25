sb_predict_centroid_binary <- function(engine, x) {
  # Equivalent to (squared distance to class 0 - distance to class 1)/2,
  # avoiding cancellation of the common x^2 term.
  c0 <- engine$centroids[1L, ]
  c1 <- engine$centroids[2L, ]
  stats::plogis(as.numeric(sweep(x, 2L, (c0 + c1) / 2, "-") %*% (c1 - c0)))
}

sb_fit_bernoulli_binary <- function(x, outcome, params) {
  binary <- x > params$binarize
  counts <- c(sum(outcome == 0L), sum(outcome == 1L))
  probability <- rbind(
    (colSums(binary[outcome == 0L, , drop = FALSE]) + params$alpha) /
      (counts[1L] + 2 * params$alpha),
    (colSums(binary[outcome == 1L, , drop = FALSE]) + params$alpha) /
      (counts[2L] + 2 * params$alpha))
  list(feature_probability = probability, prior = counts / sum(counts))
}

sb_predict_bernoulli_binary <- function(engine, x, params) {
  binary <- 1 * (x > params$binarize)
  p <- engine$feature_probability
  joint <- binary %*% t(log(p)) + (1 - binary) %*% t(log1p(-p))
  joint <- sweep(joint, 2L, log(engine$prior), "+")
  stats::plogis(joint[, 2L] - joint[, 1L])
}
