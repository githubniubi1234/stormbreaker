test_that("new defaults keep historical pipeline settings readable", {
  nn <- sb_pipeline(learner = "neuralnet_multilayer_binary")
  expect_identical(nn$learner_params$hidden, c(3L, 2L))
  expect_identical(nn$learner_params$error_function, "sse")
  expect_identical(nn$learner_params$stepmax, 500000L)
  old <- nn
  old$learner_params$hidden <- c(5L, 3L)
  old$learner_params$error_function <- "ce"
  old$learner_params$stepmax <- 100000L
  expect_identical(sb_validate_pipeline_contract(old, "binary"), old)
  fda <- sb_pipeline(learner = "fda_binary")
  expect_identical(fda$learner_params$degree, 1L)
  fda$learner_params$degree <- 2L
  expect_identical(sb_validate_pipeline_contract(fda, "binary"), fda)
  lda <- sb_pipeline(learner = "loclda_binary")
  expect_identical(lda$learner_params$k_rule, "dimension_aware")
  lda$learner_params$k_rule <- NULL
  expect_identical(sb_validate_pipeline_contract(lda, "binary"), lda)
  svm <- sb_pipeline(learner = "survivalsvm_rbf")
  expect_identical(svm$learner_params$kernel_parameter_rule, "analysis_median_squared_distance")
  svm$learner_params$kernel_parameter_rule <- "native_inverse_analysis_feature_count"
  expect_identical(sb_validate_pipeline_contract(svm, "survival"), svm)
  svm$learner_params$opt_meth <- "unregistered"
  expect_s3_class(tryCatch(sb_validate_pipeline_contract(svm, "survival"), error = identity),
    "stormbreaker_unresolved_pipeline")
})

test_that("local LDA neighbourhood accounts for dimension and class count", {
  p <- sb_classical_binary_params("loclda_binary")
  expect_equal(sb_loclda_neighbors(p, 3L, 100L, 2L), 30L)
  expect_equal(sb_loclda_neighbors(p, 60L, 400L, 2L), 124L)
  expect_equal(sb_loclda_neighbors(p, 60L, 400L, 3L), 126L)
  expect_equal(sb_loclda_neighbors(p, 60L, 80L, 2L), 80L)
  p$k_rule <- NULL
  expect_equal(sb_loclda_neighbors(p, 60L, 400L, 2L), 30L)
})
