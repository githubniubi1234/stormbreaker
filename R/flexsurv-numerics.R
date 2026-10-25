# Numerical coordinates change neither the declared model nor the feature set.
# Compact prediction parameters are mapped back to original input units.
sb_flexsurv_numerical_settings <- function(x) {
  center <- colMeans(x)
  scale <- vapply(x, stats::sd, numeric(1))
  if (any(!is.finite(scale) | scale <= 0)) {
    sb_abort("Flexible survival fitting requires positive finite training predictor scales.",
      "stormbreaker_flexsurv_invalid_scale")
  }
  list(center = center, scale = scale, maxit = 1000L,
    standardization = "training_center_sd", initialization = "native")
}

sb_flexsurv_numerical_predictors <- function(x, settings) {
  as.data.frame(sweep(sweep(as.matrix(x), 2L, settings$center, "-"),
    2L, settings$scale, "/"), check.names = FALSE)
}

sb_flexsurvspline_numerical_fit <- function(formula, data, k, scale,
                                           timescale, spline, hessian, maxit) {
  baseline <- flexsurv::flexsurvspline(stats::update(formula, . ~ 1),
    data = data, k = k, scale = scale, timescale = timescale, spline = spline,
    hessian = TRUE, control = list(maxit = maxit))
  if (!identical(baseline$opt$convergence, 0L) || any(!is.finite(stats::coef(baseline)))) {
    sb_abort("The training-fold intercept-only spline initialization did not converge.",
      "stormbreaker_flexsurvspline_nonconvergence")
  }
  mm <- stats::model.matrix(formula, stats::model.frame(formula, data = data))
  init <- c(unname(stats::coef(baseline)), rep(0, ncol(mm) - 1L))
  fit <- flexsurv::flexsurvspline(formula, data = data, k = k, scale = scale,
    timescale = timescale, spline = spline, hessian = hessian, inits = init,
    control = list(maxit = maxit))
  if (!isTRUE(all.equal(baseline$knots, fit$knots))) {
    sb_abort("Spline initialization changed the training-fold knots.",
      "stormbreaker_flexsurvspline_fit_failure")
  }
  fit
}
