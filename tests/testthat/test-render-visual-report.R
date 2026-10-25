test_that("visual gallery is offline, escaped and refuses existing outputs", {
  skip_if_not_installed("ggplot2")
  skip_if_not(capabilities("cairo"))
  p <- ggplot2::ggplot(data.frame(x = 1:3, y = 3:1), ggplot2::aes(x, y)) +
    ggplot2::geom_point()
  file <- tempfile(fileext = ".html")
  device <- grDevices::dev.cur()
  result <- sb_render_visual_report(list("A < B" = p, Other = p), file,
    title = "<script>not executable</script>", width = c(4, 5), height = 4)
  html <- paste(readLines(file), collapse = "\n")
  expect_identical(result, normalizePath(file, winslash = "/"))
  expect_identical(grDevices::dev.cur(), device)
  expect_match(html, "&lt;script&gt;not executable&lt;/script&gt;", fixed = TRUE)
  expect_match(html, "A &lt; B", fixed = TRUE)
  expect_equal(length(gregexpr("data:image/svg+xml", html, fixed = TRUE)[[1]]), 2L)
  expect_identical(grepl("<script|https?://|<table", html), FALSE)
  expect_error(sb_render_visual_report(list(A = p), file), class = "stormbreaker_report_exists")
  expect_identical(paste(readLines(file), collapse = "\n"), html)
  expect_error(sb_render_visual_report(list(A = p), tempfile(fileext = ".html"), width = NA_real_),
    class = "stormbreaker_bad_report_input")
  expect_error(sb_render_visual_report(list(A = 1), tempfile(fileext = ".html")),
    class = "stormbreaker_bad_report_input")
})

test_that("gallery restores graphics device after a drawing error", {
  skip_if_not_installed("ggplot2")
  skip_if_not(capabilities("cairo"))
  p <- ggplot2::ggplot(data.frame(x = 1), ggplot2::aes(x, missing_y)) +
    ggplot2::geom_point()
  device <- grDevices::dev.cur()
  file <- tempfile(fileext = ".html")
  expect_error(sb_render_visual_report(list(Broken = p), file))
  expect_identical(grDevices::dev.cur(), device)
  expect_identical(file.exists(file), FALSE)
})

test_that("grouped reports retain every figure and move explanations outside images", {
  skip_if_not_installed("ggplot2")
  skip_if_not(capabilities("cairo"))
  p <- ggplot2::ggplot(data.frame(x = 1:3, y = 3:1), ggplot2::aes(x, y)) +
    ggplot2::geom_point() + ggplot2::labs(title = "Native result", subtitle = "Saved fit",
      caption = "Original scale < not probability")
  original <- list(data = p$data, labels = p$labels, built = ggplot2::ggplot_build(p)$data)
  file <- tempfile(fileext = ".html")
  sb_render_visual_report(list(A = p, B = p, C = p), file,
    groups = c("Clinical", "Algorithm < details", "Clinical"),
    captions = c("Additional legend", "", ""), width = 4, height = 3)
  html <- paste(readLines(file), collapse = "\n")
  expect_match(html, "<details id='topic-1' open>", fixed = TRUE)
  expect_match(html, "<details id='topic-2'>", fixed = TRUE)
  expect_match(html, "Algorithm &lt; details", fixed = TRUE)
  expect_match(html, "Original scale &lt; not probability", fixed = TRUE)
  expect_match(html, "Additional legend", fixed = TRUE)
  expect_equal(length(gregexpr("data:image/svg+xml", html, fixed = TRUE)[[1]]), 3L)
  expect_identical(p$data, original$data)
  expect_identical(p$labels, original$labels)
  expect_identical(ggplot2::ggplot_build(p)$data, original$built)
  file_all <- tempfile(fileext = ".html")
  sb_render_visual_report(list(A = p, B = p), file_all, groups = c("A", "B"), expanded = c("A", "B"))
  all_html <- paste(readLines(file_all), collapse = "\n")
  expect_match(all_html, "<details id='topic-2' open>", fixed = TRUE)
  expect_error(sb_render_visual_report(list(A = p), tempfile(fileext = ".html"), groups = "A", expanded = "missing"),
    class = "stormbreaker_bad_report_input")
  expect_error(sb_render_visual_report(list(A = p), tempfile(fileext = ".html"), captions = c("a", "b")),
    class = "stormbreaker_bad_report_input")
})
