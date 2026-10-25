#' Inspect original HDA and HDDA discriminant subspaces
#' @param fit Successful HDA or HDDA final fit.
#' @param type HDA: loadings or densities. HDDA: loadings, spectrum or centers.
#'   NULL chooses loadings.
#' @return Original retained subspace loadings, HDA Gaussian coordinate
#'   densities, or HDDA eigenvalues, retained dimensions and original-unit
#'   centers. HDDA loadings and eigenvalues use the fitted standardized space.
#'   Spectra are not a newly computed PCA or explained-variance analysis.
#'   Axis signs are basis conventions, not feature importance or effects.
#' @export
sb_high_dimensional_diagnostics <- function(fit, type = NULL) {
  fit <- sb_validate_final_fit_for_prediction(fit)
  s <- fit$learner_state
  if (!s$learner %in% c("hda_binary", "hdda_binary_cattell")) sb_abort(
    "These diagnostics require HDA or HDDA.", "stormbreaker_unsupported_diagnostic")
  hda <- s$learner == "hda_binary"
  choices <- if (hda) c("loadings", "densities") else c("loadings", "spectrum", "centers")
  if (is.null(type)) type <- "loadings"
  if (!is.character(type) || length(type) != 1L || is.na(type) || !type %in% choices) sb_abort(
    paste("Choose", paste(choices, collapse = ", ")), "stormbreaker_unsupported_diagnostic")
  e <- s$engine
  p <- length(s$features)
  parameters <- NULL
  if (hda) {
    if (!inherits(e, "hda") || !identical(rownames(e$hda.loadings), s$features)) sb_abort(
      "Original HDA loadings are invalid.", "stormbreaker_invalid_diagnostic")
    k <- e$reduced.dimension
    if (type == "loadings") {
      m <- e$hda.loadings[, seq_len(k), drop = FALSE]
      tab <- data.frame(feature = rep(s$features, k), axis = rep(seq_len(k), each = p), loading = as.numeric(m), class = "Shared")
    } else {
      parameters <- do.call(rbind, lapply(seq_len(k), function(j) {
        m <- e$naivebayes$tables[[j]]
        data.frame(axis = j, class = rownames(m), mean = m[, 1L], sd = m[, 2L])
      }))
      tab <- do.call(rbind, lapply(seq_len(k), function(j) {
        z <- parameters[parameters$axis == j, ]
        grid <- seq(min(z$mean - 4 * z$sd), max(z$mean + 4 * z$sd), length.out = 201L)
        do.call(rbind, lapply(seq_len(nrow(z)), function(i) data.frame(axis = j, class = z$class[i],
          coordinate = grid, density = stats::dnorm(grid, z$mean[i], z$sd[i]))))
      }))
    }
  } else {
    if (!inherits(e, "hdc") || e$model != "AKJBKQKDK" || !identical(e$kname, c("0", "1"))) sb_abort(
      "Original HDDA model is invalid.", "stormbreaker_invalid_diagnostic")
    tab <- do.call(rbind, lapply(seq_along(e$kname), function(i) {
      k <- as.integer(e$d[i])
      if (type == "loadings") return(data.frame(class = e$kname[i],
        axis = rep(seq_len(k), each = p), feature = rep(s$features, k), loading = as.numeric(e$Q[[i]])))
      if (type == "spectrum") return(data.frame(class = e$kname[i], axis = seq_len(p),
        eigenvalue = as.numeric(e$ev[i, ]), retained = seq_len(p) <= k,
        model_variance = c(as.numeric(e$a[i, seq_len(k)]), rep(as.numeric(e$b[i]), p - k))))
      mu <- as.numeric(e$mu[i, ])
      if (length(e$scaling)) mu <- mu * e$scaling$sd + e$scaling$mu
      data.frame(class = e$kname[i], feature = s$features, mean = as.numeric(mu))
    }))
  }
  structure(list(type = paste0("high_dimensional_", type), table = tab, parameters = parameters,
    learner = s$learner, pipeline_id = fit$pipeline$id,
    source = "Original HDA loadings/naive-Bayes parameters or HDDA Q/ev/d/a/b/mu",
    interpretation = "Only retained native loadings are displayed. HDA densities use reduced coordinates; HDDA axes use standardized features and centers use original units. Signs are basis conventions, not effects or importance. No PCA recomputation or refitting."),
    class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_high_dimensional_diagnostic <- function(x) {
  loading <- feature <- coordinate <- density <- class <- axis <- eigenvalue <- retained <- mean <- NULL
  if (x$type == "high_dimensional_loadings") {
    p <- ggplot2::ggplot(x$table, ggplot2::aes(loading, feature)) +
      ggplot2::geom_vline(xintercept = 0, colour = "grey65", linewidth = .3) +
      ggplot2::geom_point(colour = "#2878A0", size = 2.5) +
      ggplot2::facet_wrap(~class + axis, labeller = ggplot2::label_both) +
      ggplot2::labs(title = "Retained discriminant loadings", x = "Loading", y = NULL)
  } else if (x$type == "high_dimensional_densities") {
    p <- ggplot2::ggplot(x$table, ggplot2::aes(coordinate, density, colour = class)) +
      ggplot2::geom_line(linewidth = .7) + ggplot2::facet_wrap(~axis, labeller = ggplot2::label_both) +
      ggplot2::scale_colour_manual(values = c("0" = "#2878A0", "1" = "#B86236")) +
      ggplot2::labs(title = "HDA reduced-space class densities", x = "Reduced coordinate", y = "Density", colour = "Class")
  } else if (x$type == "high_dimensional_spectrum") {
    p <- ggplot2::ggplot(x$table, ggplot2::aes(axis, eigenvalue)) +
      ggplot2::geom_line(colour = "grey60", linewidth = .4) + ggplot2::geom_point(ggplot2::aes(colour = retained), size = 2.5) +
      ggplot2::scale_colour_manual(values = c("FALSE" = "grey65", "TRUE" = "#2878A0"), labels = c("FALSE" = "Not retained", "TRUE" = "Retained")) +
      ggplot2::scale_x_continuous(breaks = unique(x$table$axis)) +
      ggplot2::facet_wrap(~class, labeller = ggplot2::label_both) +
      ggplot2::labs(title = "HDDA class spectra", x = "Axis", y = "Eigenvalue", colour = NULL)
  } else {
    p <- ggplot2::ggplot(x$table, ggplot2::aes(mean, feature)) + ggplot2::geom_point(colour = "#2878A0", size = 2.5) +
      ggplot2::facet_wrap(~class, labeller = ggplot2::label_both) +
      ggplot2::labs(title = "HDDA class centers", x = "Mean in original units", y = NULL)
  }
  p <- p + sb_viz_theme()
  attr(p, "stormbreaker_suggested_size") <- c(width = 11, height = 5)
  p
}
