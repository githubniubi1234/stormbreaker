sb_nomogram_overlays <- function(p, nomogram, reference_data, patient) {
  axes <- attr(p, "stormbreaker_nomogram_axes")
  features <- nomogram$features$feature
  shown <- features[features %in% axes$name]
  raised <- p$data$label_y > p$data$y
  p$data$label_y[raised] <- p$data$y[raised] - .45
  density_rows <- list()
  if (!is.null(reference_data)) {
    values <- sb_nomogram_values(nomogram, reference_data)
    if (any(values$outside_display_range)) stop("Reference observations extend beyond the nomogram display ranges.")
    add_density <- function(v, y) {
      if (length(v) < 2L || length(unique(v)) < 2L) return(NULL)
      d <- stats::density(v, from = min(v), to = max(v), n = 128)
      data.frame(x = d$x, y = y + .55 * d$y / max(d$y), axis = y)
    }
    for (f in shown) density_rows[[length(density_rows) + 1L]] <- add_density(values$points[, f], axes$y[match(f, axes$name)])
    if (nomogram$total_points > 0) density_rows[[length(density_rows) + 1L]] <- add_density(values$total_points / nomogram$total_points * 100, 2)
  }
  x <- y <- axis <- label <- NULL
  if (length(density_rows)) {
    d <- do.call(rbind, density_rows)
    if (!is.null(d)) p <- p + ggplot2::geom_line(data = d, ggplot2::aes(x, y, group = axis),
      inherit.aes = FALSE, colour = "#4477cc", linewidth = .45)
    attr(p, "stormbreaker_nomogram_density") <- d
  }
  if (!is.null(patient)) {
    if (!is.data.frame(patient) || nrow(patient) != 1L) stop("Supply exactly one patient row.")
    value <- sb_nomogram_values(nomogram, patient)
    if (value$outside_display_range) stop("Patient is outside the displayed ranges; recreate the nomogram with suitable ranges.")
    markers <- data.frame(x = as.numeric(value$points[1, shown]), y = axes$y[match(shown, axes$name)])
    p <- p + ggplot2::geom_segment(data = markers, ggplot2::aes(x = x, xend = x, y = y, yend = max(axes$y)),
      inherit.aes = FALSE, colour = "#d73027", linetype = 3, linewidth = .35)
    total_x <- if (nomogram$total_points > 0) value$total_points / nomogram$total_points * 100 else 0
    markers <- rbind(markers, data.frame(x = total_x, y = c(2, 1)))
    labels <- data.frame(x = total_x, y = c(2, 1) + .20,
      label = c(sprintf("%.1f", value$total_points), sprintf("%.3f", if (is.null(value$probability)) value$score else value$probability)))
    p <- p + ggplot2::geom_point(data = markers, ggplot2::aes(x, y), inherit.aes = FALSE, colour = "#d73027", size = 2) +
      ggplot2::geom_text(data = labels, ggplot2::aes(x, y, label = label), inherit.aes = FALSE, colour = "#d73027", size = 3)
    attr(p, "stormbreaker_nomogram_patient") <- value
  }
  p + ggplot2::labs(title = "Nomogram", subtitle = NULL, caption = NULL)
}
