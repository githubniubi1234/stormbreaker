sb_nomogram_aft_probability <- function(aft, score) {
  z <- (log(aft$prediction_horizon) + score) / aft$scale
  switch(aft$distribution,
    weibull = exp(-exp(z)), exponential = exp(-exp(z)),
    lognormal = stats::pnorm(z, lower.tail = FALSE),
    loglogistic = stats::plogis(z, lower.tail = FALSE))
}

sb_nomogram_aft_score <- function(aft, probability) {
  z <- switch(aft$distribution,
    weibull = log(-log(probability)), exponential = log(-log(probability)),
    lognormal = stats::qnorm(probability, lower.tail = FALSE),
    loglogistic = stats::qlogis(probability, lower.tail = FALSE))
  aft$scale * z - log(aft$prediction_horizon)
}
