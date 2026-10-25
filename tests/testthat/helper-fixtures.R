make_sb_fixture <- function(n = 60L, seed = 20260812L) {
  withr::with_seed(seed, {
    event <- sample(rep(c(0L, 1L), length.out = n))
    data.frame(
      patient_id = sprintf("P%03d", seq_len(n)),
      time = round(stats::rexp(n, rate = 0.08) + 0.1, 3),
      event = event,
      gene_a = stats::rnorm(n),
      gene_b = stats::rnorm(n),
      gene_c = stats::rnorm(n),
      check.names = FALSE
    )
  })
}

make_sb_binary_fixture <- function(n = 60L, seed = 20260818L) {
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

make_sb_high_dim_binary_fixture <- function(n = 60L, p = 80L,
                                            seed = 20260819L) {
  withr::with_seed(seed, {
    outcome <- sample(rep(c(0L, 1L), length.out = n))
    expression <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
    expression[, 1L] <- expression[, 1L] + 0.8 * outcome
    expression[, 2L] <- expression[, 2L] - 0.5 * outcome
    colnames(expression) <- sprintf("gene_%03d", seq_len(p))
    data.frame(
      patient_id = sprintf("HD%03d", seq_len(n)),
      outcome = outcome,
      as.data.frame(expression, check.names = FALSE),
      check.names = FALSE
    )
  })
}

capture_sb_error <- function(code) {
  tryCatch(
    {
      force(code)
      NULL
    },
    stormbreaker_error = function(error) error
  )
}

population_variance <- function(x) {
  mean((x - mean(x))^2)
}

manual_gaussian_nb_probability <- function(x, priors, means, variances) {
  vapply(seq_len(nrow(x)), function(row) {
    log_joint <- vapply(c("0", "1"), function(class) {
      log(priors[[class]]) + sum(stats::dnorm(
        x[row, ],
        mean = means[class, ],
        sd = sqrt(variances[class, ]),
        log = TRUE
      ))
    }, numeric(1))
    unname(stats::plogis(log_joint[["1"]] - log_joint[["0"]]))
  }, numeric(1))
}
