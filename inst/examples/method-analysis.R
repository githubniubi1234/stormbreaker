# Runnable synthetic example; no files are exported automatically.
library(stormbreaker)
set.seed(881)
n <- 120L
d <- data.frame(patient_id = paste0("s", seq_len(n)))
genes <- paste0("gene_", letters[1:8])
for (gene in genes) d[[gene]] <- rnorm(n)
d$outcome <- rbinom(n, 1, plogis(2 * d$gene_a - 1.3 * d$gene_b))
task <- bi_task(d, "patient_id", "outcome", genes)

# Each requested method is fitted independently, once.
analysis <- analyze_methods(task, c("LASSO", "Boruta", "random_forest"))
analysis$status
selected_genes(analysis, "LASSO")
sb_method_features(analysis, "random_forest")
analysis$results$lasso$plots$coefficient_path

# Optional: combine saved plots. No method is fitted again.
if (requireNamespace("patchwork", quietly = TRUE)) {
  panel <- patchwork::wrap_plots(list(
    analysis$results$random_forest$plots$oob_error,
    analysis$results$random_forest$plots$importance,
    analysis$results$boruta$plots$boruta_summary,
    analysis$results$boruta$plots$boruta_history,
    analysis$results$lasso$plots$coefficient_path,
    analysis$results$lasso$plots$cv_curve
  ), ncol = 2) + patchwork::plot_annotation(tag_levels = "A")
  print(panel)
}
