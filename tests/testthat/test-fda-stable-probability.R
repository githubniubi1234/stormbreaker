test_that("FDA probabilities match native posteriors on ordinary observations", {
  skip_if_not_installed("mda")
  d <- sb_example_binary(80, seed = 9301)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  for (id in c("fda_binary", "pda_binary")) {
    adapter <- sb_classical_binary_learner_adapter(id)
    fit <- adapter$fit(x, d$outcome)
    expected <- unname(predict(fit$engine, newdata = x, type = "posterior")[, "1"])
    expect_equal(adapter$predict(fit, x), expected, tolerance = 1e-12)
  }
})

test_that("FDA computes finite posteriors when native exponentials underflow", {
  skip_if_not_installed("mda")
  d <- sb_example_binary(80, seed = 9301)
  x <- d[c("gene_a", "gene_b", "gene_c")]
  extreme <- x[1:3, , drop = FALSE] * 1000
  for (id in c("fda_binary", "pda_binary")) {
    adapter <- sb_classical_binary_learner_adapter(id)
    fit <- adapter$fit(x, d$outcome)
    native <- predict(fit$engine, newdata = extreme, type = "posterior")
    distances <- predict(fit$engine, newdata = extreme, type = "distances")
    expect_true(any(!is.finite(native)))
    expected <- plogis(-0.5 * (distances[, "1"] - distances[, "0"]) +
      log(fit$engine$prior["1"] / fit$engine$prior["0"]))
    probability <- adapter$predict(fit, extreme)
    expect_equal(probability, unname(expected), tolerance = 1e-12)
    expect_true(all(is.finite(probability) & probability >= 0 & probability <= 1))
    restored <- unserialize(serialize(adapter$extract_state(fit), NULL))
    expect_equal(adapter$predict_state(restored, extreme)$probability, probability)
    expect_equal(adapter$predict(fit, extreme[3:1, , drop = FALSE]), rev(probability))
  }
})
