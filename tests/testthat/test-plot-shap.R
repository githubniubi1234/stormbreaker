make_sb_shap_plot_fixture <- function() {
  m <- rbind(c(a = 1, b = -0.3, c = 0.2), c(a = -0.5, b = 0.4, c = -0.1))
  rownames(m) <- c("patient-secret-1", "patient-secret-2")
  structure(list(status = "success", pipeline_id = "logistic", score_scale = "log_odds",
    reference_value = -0.4, contributions = m,
    scores = data.frame(patient_id = rownames(m), score = -0.4 + rowSums(m)),
    background = list(n = 10L)), class = "stormbreaker_additive_shap")
}

test_that("SHAP importance is mean absolute saved contribution on the score scale", {
  skip_if_not_installed("ggplot2")
  x <- make_sb_shap_plot_fixture()
  before <- x
  p <- sb_plot_shap_importance(x)
  expect_equal(p$data$mean_absolute_shap, c(0.75, 0.35, 0.15))
  expect_identical(p$data$feature, c("a", "b", "c"))
  expect_identical(p$labels$x, "Mean |SHAP| (log odds)")
  expect_identical(x, before)
})

test_that("groupby reuses contributions and aligns named groups by ID", {
  skip_if_not_installed("ggplot2")
  x <- make_sb_shap_plot_fixture()
  before <- x
  p <- sb_plot_shap_importance(x, groupby = c("Normal", "Disease"), top_n = 2)
  expect_equal(p$data$mean_absolute_shap, c(1, .3, .5, .4))
  expect_identical(p$data$feature, rep(c("a", "b"), 2))
  expect_equal(p$data$n_observations, rep(1, 4))
  groups <- setNames(c("Disease", "Normal"), rev(rownames(x$contributions)))
  expect_identical(sb_plot_shap_importance(x, top_n = 2, groupby = groups)$data, p$data)
  q <- sb_plot_shap_importance(x, groupby = c("All", "All"))
  expect_equal(q$data$mean_absolute_shap, sb_plot_shap_importance(x)$data$mean_absolute_shap)
  one <- sb_plot_shap_importance(x, top_n = 1, groupby = factor(c("Normal", "Disease"),
    levels = c("Disease", "Normal", "Unused")))
  expect_equal(one$data$mean_absolute_shap, c(.5, 1))
  expect_identical(levels(one$data$group), c("Disease", "Normal"))
  built <- ggplot2::ggplot_build(p)
  expect_equal(length(unique(built$data[[1]]$y)), 4)
  expect_equal(built$data[[1]]$y, built$data[[2]]$y)
  expect_null(p$labels$subtitle)
  expect_identical(x, before)
  expect_error(sb_plot_shap_importance(x, groupby = "Normal"), class = "stormbreaker_bad_visualization_argument")
  expect_error(sb_plot_shap_importance(x, groupby = c(NA, "Normal")), class = "stormbreaker_bad_visualization_argument")
  expect_error(sb_plot_shap_importance(x, groupby = c("", "Normal")), class = "stormbreaker_bad_visualization_argument")
  expect_error(sb_plot_shap_importance(x, groupby = c(wrong = "Normal", other = "Disease")),
    class = "stormbreaker_bad_visualization_argument")
})

test_that("waterfall retains baseline and remainder and reconstructs the full score", {
  skip_if_not_installed("ggplot2")
  x <- make_sb_shap_plot_fixture()
  p <- sb_plot_shap_waterfall(x, row = 1, top_n = 1)
  expect_identical(p$data$kind, c("reference", "feature", "remainder", "prediction"))
  expect_equal(p$data$contribution, c(-0.4, 1, -0.1, 0.5))
  expect_equal(p$data$end, c(-0.4, 0.6, 0.5, 0.5))
  expect_equal(tail(p$data$end, 1), x$scores$score[1])
  expect_match(p$labels$subtitle, "Observation 1", fixed = TRUE)
  expect_identical(any(grepl("patient-secret", unlist(p$data))), FALSE)
  expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
})

test_that("beeswarm is deterministic and feature values align by explicit IDs", {
  skip_if_not_installed("ggplot2")
  x <- make_sb_shap_plot_fixture()
  set.seed(29)
  before <- .Random.seed
  p <- sb_plot_shap_beeswarm(x)
  q <- sb_plot_shap_beeswarm(x)
  expect_identical(p$data, q$data)
  expect_identical(.Random.seed, before)
  expect_equal(nrow(p$data), 6L)
  feature_values <- data.frame(id = rev(rownames(x$contributions)), a = c(20, 10), b = c(40, 30), c = c(60, 50))
  coloured <- sb_plot_shap_beeswarm(x, feature_values = feature_values, id_col = "id")
  expect_equal(coloured$data$feature_value, c(10, 20, 30, 40, 50, 60))
  expect_identical(any(grepl("patient-secret", unlist(coloured$data))), FALSE)
  expect_error(sb_plot_shap_beeswarm(x, feature_values = feature_values),
    class = "stormbreaker_bad_visualization_argument")
})

test_that("bad attribution states are not silently relabelled as valid explanations", {
  skip_if_not_installed("ggplot2")
  x <- make_sb_shap_plot_fixture()
  expect_error(sb_plot_shap_waterfall(x, row = 3), class = "stormbreaker_bad_visualization_argument")
  x$scores$score[1] <- 100
  expect_error(sb_plot_shap_waterfall(x), class = "stormbreaker_bad_shap_visualization")
})
