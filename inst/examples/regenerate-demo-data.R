# Synthetic inputs only. Source this file, then call regenerate_demo_data("new-directory").
regenerate_demo_data <- function(directory) {
  if (file.exists(directory)) stop("Choose a new directory; existing files are not overwritten.")
  specs <- data.frame(
    file = c("demo-binary-development.csv", "demo-binary-external-a.csv",
      "demo-binary-external-b.csv", "demo-survival-development.csv",
      "demo-survival-external.csv"),
    task = c("binary", "binary", "binary", "survival", "survival"),
    n = c(120, 80, 70, 160, 100), seed = c(42, 17, 23, 42, 17),
    prefix = c("DEV_", "EXTA_", "EXTB_", "SDEV_", "SEXT_")
  )
  if (!dir.create(directory, recursive = TRUE)) stop("Could not create the requested directory.")
  for (i in seq_len(nrow(specs))) {
    d <- if (specs$task[i] == "binary") {
      stormbreaker::bi_example(specs$n[i], specs$seed[i])
    } else {
      stormbreaker::surv_example(specs$n[i], specs$seed[i])
    }
    d$patient_id <- paste0(specs$prefix[i], d$patient_id)
    utils::write.csv(d, file.path(directory, specs$file[i]),
      row.names = FALSE, fileEncoding = "UTF-8")
  }
  invisible(specs)
}
