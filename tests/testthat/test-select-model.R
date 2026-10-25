selection_example <- function() {
  make_task <- function(prefix, seed) {
    d <- sb_example_binary(60, seed)
    d$patient_id <- paste0(prefix, d$patient_id)
    sb_binary_task(d, "patient_id", "outcome", c("gene_a", "gene_b"))
  }
  sb_compare_cohorts(make_task("train", 1),
    list(Validation = make_task("val", 2), Test = make_task("test", 3)),
    roles = c(Validation = "validation", Test = "test"),
    models = c("logistic", "gaussian_naive_bayes"), folds = 3, view = FALSE
  )
}

test_that("validation selection returns a directly predictable existing fit", {
  compared <- selection_example()
  table <- compared$cohorts$Validation$comparison
  selected <- sb_select_model(compared)
  id <- table$pipeline_id[which.max(table$roc_auc)]
  index <- match(id, vapply(compared$fits, function(x) x$pipeline$id, character(1)))
  expect_s3_class(selected, "stormbreaker_final_fit")
  expect_identical(selected$pipeline$id, id)
  expect_identical(selected$learner_state, compared$fits[[index]]$learner_state)
  d <- sb_example_binary(20, 4)
  expect_identical(
    sb_predict_final(selected, d, "patient_id")$predictions,
    sb_predict_final(unserialize(serialize(selected, NULL)), d, "patient_id")$predictions
  )
  expect_identical(selected$model_selection$cohorts, "Validation")
})

test_that("validation scores are aligned by model ID and equally weighted", {
  compared <- selection_example()
  compared$cohorts$Second <- compared$cohorts$Validation
  compared$cohort_roles <- c(compared$cohort_roles, Second = "validation")
  compared$cohorts$Validation$comparison$roc_auc <- c(0.9, 0.6)
  compared$cohorts$Second$comparison$roc_auc <- c(0.2, 0.8)
  compared$cohorts$Second$comparison <- compared$cohorts$Second$comparison[2:1, ]
  selected <- sb_select_model(compared)
  expect_identical(selected$pipeline$id, compared$model_pipelines[[2]]$id)
  expect_equal(unname(selected$model_selection$ranking$validation_mean), c(0.7, 0.55))
  expect_identical(selected$model_selection$cohorts, c("Validation", "Second"))
})

test_that("test scores never change selection and ties keep declaration order", {
  compared <- selection_example()
  before <- sb_select_model(compared)$pipeline$id
  compared$cohorts$Test$comparison$roc_auc <- c(999, -999)
  expect_identical(sb_select_model(compared)$pipeline$id, before)
  compared$cohorts$Validation$comparison$roc_auc <- c(0.7, 0.7)
  expect_identical(sb_select_model(compared)$pipeline$id, compared$model_pipelines[[1]]$id)
  compared$fits[[1]]$status <- "failed"
  expect_identical(sb_select_model(compared)$pipeline$id, compared$model_pipelines[[2]]$id)
})

test_that("selection requires validation and complete eligible candidates", {
  compared <- selection_example()
  compared$cohort_roles["Validation"] <- "external"
  expect_error(sb_select_model(compared), "validation", class = "stormbreaker_no_validation_cohort")
  compared$cohort_roles["Validation"] <- "validation"
  compared$cohorts$Validation$comparison$roc_auc <- NA_real_
  expect_error(sb_select_model(compared), "eligible", class = "stormbreaker_no_eligible_model")
})

test_that("survival selection uses Harrell C and preserves the fitted state", {
  make_task <- function(prefix, seed) {
    d <- sb_example_survival(55, seed)
    d$patient_id <- paste0(prefix, d$patient_id)
    sb_surv_task(d, "patient_id", "time", "event", c("gene_a", "gene_b"), "month")
  }
  compared <- sb_compare_cohorts(make_task("train", 10),
    list(Validation = make_task("val", 11)), roles = c(Validation = "validation"),
    pipelines = list(sb_pipeline(), sb_pipeline(selector = "univ_cox_top_k", selector_n_features = 1)),
    folds = 3, view = FALSE, progress = FALSE)
  selected <- sb_select_model(compared)
  tab <- compared$cohorts$Validation$comparison
  expected <- tab$pipeline_id[which.max(tab$harrell_c)]
  expect_identical(selected$pipeline$id, expected)
  expect_identical(selected$model_selection$metric, "harrell_c")
})
