#' Generate synthetic right-censored survival data
#'
#' This generator is only for examples and automated tests. It creates no
#' files and does not represent a clinical cohort. `gene_a` influences the
#' simulated event hazard, while `gene_b` and `gene_c` are weaker noise
#' predictors. A common administrative censoring cutoff yields enough events
#' and censored patients for the package's five-fold example.
#'
#' @param n Number of synthetic patients. Must be at least 10.
#' @param seed Positive integer controlling all random generation.
#'
#' @return A data frame containing `patient_id`, `time`, `event`, and three
#'   numeric predictor columns. Time is expressed in arbitrary month-like
#'   units for demonstration.
#' @export
#'
#' @examples
#' d <- sb_example_survival(n = 50, seed = 1)
#' head(d)
sb_example_survival <- function(n = 60L, seed = 20260812L) {
  n <- sb_check_positive_integer(n, "n", minimum = 10L)
  seed <- sb_check_positive_integer(seed, "seed", minimum = 1L)

  withr::with_seed(seed, {
    gene_a <- stats::rnorm(n)
    gene_b <- stats::rnorm(n)
    gene_c <- stats::rnorm(n)
    latent_rate <- exp(0.55 * gene_a - 0.15 * gene_b)
    event_time <- stats::rexp(n, rate = 0.055 * latent_rate) + 0.1
    target_events <- floor(n / 2)
    ordered_time <- sort(event_time)
    censor_time <- mean(ordered_time[c(target_events, target_events + 1L)])
    event <- as.integer(event_time <= censor_time)
    observed_time <- pmin(event_time, censor_time)

    data.frame(
      patient_id = sprintf("P%03d", seq_len(n)),
      time = round(observed_time, 3),
      event = event,
      gene_a = gene_a,
      gene_b = gene_b,
      gene_c = gene_c,
      check.names = FALSE
    )
  })
}
