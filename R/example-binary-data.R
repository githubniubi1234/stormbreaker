#' Generate synthetic binary classification data
#'
#' This generator is only for examples and automated tests. It creates no
#' files and does not represent a clinical cohort. The outcome is balanced,
#' class `1` is positive, and all predictors are numeric.
#'
#' @param n Number of synthetic patients. Must be at least 10.
#' @param seed Positive integer controlling all random generation.
#'
#' @return A data frame containing `patient_id`, `outcome`, and three numeric
#'   predictor columns.
#' @export
#'
#' @examples
#' d <- sb_example_binary(n = 50, seed = 1)
#' head(d)
sb_example_binary <- function(n = 60L, seed = 20260818L) {
  n <- sb_check_positive_integer(n, "n", minimum = 10L)
  seed <- sb_check_positive_integer(seed, "seed", minimum = 1L)

  withr::with_seed(seed, {
    outcome <- sample(rep(c(0L, 1L), length.out = n))
    data.frame(
      patient_id = sprintf("B%03d", seq_len(n)),
      outcome = outcome,
      gene_a = stats::rnorm(n, mean = 0.5 * outcome),
      gene_b = stats::rnorm(n),
      gene_c = stats::rnorm(n),
      check.names = FALSE
    )
  })
}
