sb_predict_fda_native <- function(engine, x, type, dimension = ncol(engine$means)) {
  n <- nrow(x)
  # mda::polybasis drops the row dimension for singleton interaction queries.
  # Its polynomial evaluation is row-local; duplicate only to retain dimensions.
  singleton <- n == 1L && inherits(engine$fit, "polyreg") &&
    engine$fit$degree > 1L && !isTRUE(engine$fit$monomial) && ncol(x) > 1L
  if (singleton) x <- x[c(1L, 1L), , drop = FALSE]
  native <- get("predict.fda", asNamespace("mda"))
  result <- native(engine, newdata = x, type = type, dimension = dimension)
  result[seq_len(n), , drop = FALSE]
}

sb_predict_fda_probability <- function(engine, x) {
  distances <- sb_predict_fda_native(engine, x, type = "distances")
  # The native posterior divides exponentiated class weights. Their log odds
  # give the same binary posterior without a zero/zero underflow.
  log_odds <- -0.5 * (distances[, "1"] - distances[, "0"]) +
    log(engine$prior[["1"]]) - log(engine$prior[["0"]])
  stats::plogis(log_odds)
}
