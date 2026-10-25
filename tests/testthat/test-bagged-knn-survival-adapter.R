test_that("bagged distance screening preserves exact selected neighbors", {
  withr::local_seed(491)
  for (p in c(1L, 5L, 30L)) {
    for (scenario in c("ordinary", "duplicates", "offset", "ill_conditioned",
                       "near_ties")) {
      x <- matrix(rnorm(100L * p), 100L, p)
      if (scenario == "duplicates") x[11:20, ] <- x[1:10, ]
      if (scenario == "offset") x <- x + 1e9
      if (scenario == "ill_conditioned" && p > 1L) {
        x[, p] <- x[, 1L] + 1e-6 * x[, p]
      }
      if (scenario == "near_ties") {
        x[2L, ] <- x[1L, ] + 1e-14
        x[3L, ] <- x[1L, ] - 1e-14
      }
      rownames(x) <- paste0("row", seq_len(nrow(x)))
      covariance <- cov(x)
      inverse <- solve(covariance, tol = 1e-25)
      distance <- sb_bagged_knn_distance_function(10L)
      queries <- rbind(x[1:4, , drop = FALSE], colMeans(x))
      for (i in seq_len(nrow(queries))) {
        expected <- sort(stats::mahalanobis(x, queries[i, ], inverse,
          inverted = TRUE), index.return = TRUE)
        actual <- sort(distance(x, queries[i, ], covariance, tol = 1e-25),
          index.return = TRUE)
        expect_identical(actual$ix[1:10], expected$ix[1:10])
        expect_identical(actual$x[1:10], expected$x[1:10])
      }
    }
  }
})

test_that("bagged full-neighborhood distances and singular errors stay native", {
  x <- matrix(seq_len(20), ncol = 1L)
  distance <- sb_bagged_knn_distance_function(nrow(x))
  expect_identical(distance(x, 2.5, cov(x), tol = 1e-25),
    stats::mahalanobis(x, 2.5, cov(x), tol = 1e-25))
  singular <- cbind(x, x)
  condition <- function(expr) tryCatch(expr, error = identity)
  native <- condition(stats::mahalanobis(singular, c(1, 1), cov(singular)))
  got <- condition(sb_bagged_knn_distance_function(2L)(
    singular, c(1, 1), cov(singular)))
  expect_identical(class(got), class(native))
  expect_identical(conditionMessage(got), conditionMessage(native))
})

test_that("bagged single-feature single-query native curves remain unchanged", {
  skip_if_not_installed("bnnSurvival", minimum_version = "0.1.5")
  withr::local_seed(731)
  x <- data.frame("gene one" = rnorm(60L), check.names = FALSE)
  time <- rexp(60L) + 0.01
  event <- rep(c(1L, 0L, 1L), 20L)
  adapter <- sb_bagged_knn_survival_learner_adapter("bagged_knn_survival")
  model <- adapter$fit(x, time, event)
  model <- unserialize(serialize(model, NULL))
  query <- x[4L, , drop = FALSE]
  prepared <- sb_bagged_knn_prediction_data(model, query)
  native <- get("predict", asNamespace("bnnSurvival"))
  before <- .Random.seed
  reference <- sb_with_bagged_knn_single_core(
    native(prepared$model$engine, prepared$internal))
  actual <- sb_predict_bagged_knn_survival_curves(model, query)
  expect_identical(actual$survival,
    as.matrix(bnnSurvival::predictions(reference)))
  expect_identical(.Random.seed, before)
})
