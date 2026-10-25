sb_loclda_neighbors <- function(params, p, n, classes) {
  k <- params$k_max
  if (identical(params$k_rule, "dimension_aware")) {
    k <- max(k, 2L * (p + classes))
  }
  min(n, k)
}
