enet_gbm_task <- function() {
  withr::local_seed(42)
  d <- as.data.frame(matrix(rnorm(240 * 8), 240, 8))
  names(d) <- c("HLA-DRA", "gene with space", paste0("g", 3:8))
  d$id <- paste0("p", seq_len(nrow(d)))
  d$y <- rbinom(nrow(d), 1, plogis(2 * d[[1]] - 1.5 * d[[2]] + d[[3]]))
  sb_binary_task(d, "id", "y", names(d)[1:8])
}

test_that("Enet selection uses its alpha and preserves patient-level folds", {
  task <- enet_gbm_task()
  x <- task$data[task$predictors]
  groups <- rep(seq_len(120), each = 2)
  params <- sb_binary_model_selector_options("enet_binary_nonzero")
  state <- sb_fit_binary_model_selector(x, task$data$y, "enet_binary_nonzero", params, groups)
  work <- x
  names(work) <- paste0("feature_", seq_len(ncol(x)))
  native <- glmnet::glmnet(as.matrix(work), task$data$y, family = "binomial",
    alpha = .5, lambda = state$evidence$lambda, standardize = TRUE)
  expect_equal(state$scores$score, as.numeric(stats::coef(native))[-1])
  expect_identical(state$selected_features, names(x)[state$scores$score != 0])
  expect_identical(all(vapply(split(state$evidence$fold_id, groups),
    function(g) length(unique(g)) == 1L, logical(1))), TRUE)
  expect_identical(bi_selector("ElasticNet"), bi_selector("Enet"))
  expect_error(sb_binary_model_selector_options("enet_binary_nonzero", list(alpha = 0)),
    class = "stormbreaker_bad_selector_options")
})

test_that("GBM selection matches native relative influence without invented cutoffs", {
  task <- enet_gbm_task()
  x <- task$data[task$predictors]
  state <- sb_fit_binary_model_selector(x, task$data$y, "gbm_binary_importance",
    sb_binary_model_selector_options("gbm_binary_importance"))
  work <- x
  names(work) <- paste0("feature_", seq_len(ncol(x)))
  withr::local_seed(1)
  native <- gbm::gbm.fit(as.matrix(work), task$data$y, distribution = "bernoulli",
    n.trees = 100, interaction.depth = 1, n.minobsinnode = 10, shrinkage = .1,
    bag.fraction = 1, nTrain = nrow(x), keep.data = FALSE, verbose = FALSE,
    var.names = names(work))
  imp <- summary(native, n.trees = 100, plotit = FALSE)
  expect_equal(state$scores$score, imp$rel.inf[match(names(work), imp$var)])
  expect_identical(state$selected_features, names(x)[state$scores$score > 0])
  expect_equal(state$evidence$training_loss$loss, native$train.error)
})

test_that("new selectors run in ordered chains with saved genes and native plots", {
  task <- enet_gbm_task()
  for (selectors in list("Enet", "GBM", c("Enet", "GBM"), c("GBM", "Lasso"))) {
    p <- bi_pipeline("RF", bi_selector(selectors))
    fit <- sb_fit_final(task, p, capture_diagnostics = TRUE)
    expect_identical(fit$status, "success")
    expect_identical(sb_lint_pipeline(p)$status, "admitted")
    z <- sb_plots(fit)
    expect_contains(names(z$plots), "selector_1")
    if (selectors[1] == "GBM") {
      expect_contains(names(z$plots), "selector_1_training_loss")
      expect_identical(z$plots$selector_1$labels$x, "Importance")
    }
    if (length(selectors) == 2) {
      states <- fit$selection$state$stages
      expect_identical(states[[2]]$source_features, states[[1]]$selected_features)
    }
    cv <- sb_run_cv(task, p, folds = 3)
    expect_identical(cv$status, "success")
    expect_equal(sb_predict_final(fit, task$data, "id")$predictions,
      sb_predict_final(unserialize(serialize(fit, NULL)), task$data, "id")$predictions)
  }
})
