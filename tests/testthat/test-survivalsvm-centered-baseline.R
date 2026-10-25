test_that("kernel probabilities are invariant to large risk offsets", {
  x <- data.frame(a = c(-1, 0, 1, 2), b = c(1, -1, 0, 2))
  time <- c(1, 2, 3, 4)
  event <- c(1, 1, 0, 1)
  for (learner in c("survivalsvm_rbf", "survivalsvm_additive")) {
    params <- sb_survivalsvm_params(1, learner)
    model <- sb_new_survivalsvm_kernel_model(
      learner, names(x), as.matrix(x), c(.1, -.1, .2, -.2), 0,
      1, if (params$kernel == "rbf_kernel") 1 else NULL,
      params, list(center = c(0, 0), scale = c(1, 1)), NULL, "test"
    )
    raw <- sb_predict_survivalsvm_kernel(model, x)
    expected <- exp(-vapply(raw, function(r) {
      sum(vapply(c(1, 2), function(t) {
        sum(event == 1 & time == t) / sum(exp(raw[time >= t] - r))
      }, numeric(1)))
    }, numeric(1)))
    context <- list(x = x, time = time, event = event,
                    prediction_horizon = 2, time_unit = "day")
    for (offset in c(0, -2000, 2000)) {
      shifted <- model
      shifted$intercept <- offset
      direct <- sb_predict_survivalsvm_kernel_survival(shifted, x, time, event, x, 2)
      expect_equal(direct$probability, expected, tolerance = 1e-10)
      state <- sb_extract_survivalsvm_kernel_state(shifted, context)
      saved <- sb_predict_survivalsvm_kernel_state(unserialize(serialize(state, NULL)), x)
      expect_equal(saved$risk, raw - offset)
      expect_equal(saved$survival_probability, expected, tolerance = 1e-10)
      expect_equal(state$baseline_hazard$risk_reference, mean(raw - offset))
      for (bad in list(NA_real_, Inf, "zero", c(0, 1))) {
        broken <- state
        broken$baseline_hazard$risk_reference <- bad
        expect_error(sb_predict_survivalsvm_kernel_state(broken, x),
                     class = "stormbreaker_bad_learner_state")
      }
      broken <- state
      broken$baseline_hazard$prediction_horizon <- -1
      expect_error(sb_predict_survivalsvm_kernel_state(broken, x))
      broken <- state
      broken$baseline_hazard$cumulative_baseline_hazard <- -1
      expect_error(sb_predict_survivalsvm_kernel_state(broken, x))
    }
    legacy <- sb_extract_survivalsvm_kernel_state(model, context)
    legacy$baseline_hazard$risk_reference <- NULL
    legacy$baseline_hazard$cumulative_baseline_hazard <-
      sb_breslow_baseline_hazard(time, event, raw, 2)
    expect_equal(sb_predict_survivalsvm_kernel_state(legacy, x)$survival_probability,
                 expected, tolerance = 1e-10)
  }
})
