test_that("the full learner registry uses clean canonical identifiers", {
  registry <- sb_pipeline_registry()
  expect_length(unique(registry$learner), 176L)
  expect_length(grep("_fixed$", registry$learner), 0L)
  for (old in names(sb_legacy_learner_names)) {
    clean <- sb_legacy_learner_names[[old]]
    expect_identical(sb_pipeline(learner = old), sb_pipeline(learner = clean))
  }
})

test_that("name compatibility preserves feature names and parameter values", {
  pipeline <- sb_pipeline(learner = "ranger_binary", ranger_num_trees = 23L)
  legacy <- pipeline
  legacy$learner <- "ranger_binary_fixed"
  legacy$id <- sub("ranger_binary", "ranger_binary_fixed", legacy$id, fixed = TRUE)
  expect_identical(sb_upgrade_legacy_names(legacy), pipeline)
  expect_identical(sb_ranger_binary_fixed_params(), sb_ranger_binary_params())
  genes <- data.frame(ranger_binary_fixed = c(1, 2), check.names = FALSE)
  expect_identical(sb_upgrade_legacy_names(genes), genes)
  expect_identical(sb_canonical_learner_id("invented_fixed"), "invented_fixed")
  expect_identical(sb_canonical_learner_id("fixed_features"), "fixed_features")
})
