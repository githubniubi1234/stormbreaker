library(stormbreaker)

make_binary_task <- function(prefix, n, seed) {
  d <- bi_example(n = n, seed = seed)
  d$patient_id <- paste0(prefix, d$patient_id)
  bi_task(d, "patient_id", "outcome",
    predictors = c("gene_a", "gene_b", "gene_c"))
}

multi_cohort_result <- sb_compare_cohorts(
  development = make_binary_task("dev_", 100, 42),
  cohorts = list(
    Validation1 = make_binary_task("val1_", 60, 43),
    Validation2 = make_binary_task("val2_", 60, 44),
    Test = make_binary_task("test_", 80, 45)
  ),
  roles = c(Validation1 = "validation", Validation2 = "validation", Test = "test"),
  models = c("logistic", "gaussian_naive_bayes"),
  development_name = "Training", folds = 5, seed = 42
)

print(multi_cohort_result)
multi_cohort_result$cohorts$Test$comparison

# For a persistent report, choose a new filename:
# sb_view_results(multi_cohort_result, path = "cohort-results.html", open = FALSE)
