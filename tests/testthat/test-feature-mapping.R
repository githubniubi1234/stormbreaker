test_that("explicit feature mapping is outcome-blind and ordered", {
  new_data <- data.frame(
    patient_id = c("p1", "p2", "p3"),
    outcome = c(0L, 1L, 0L),
    ENSG_B = c(2.0, 2.5, 3.0),
    unrelated = c(10, 20, 30),
    ENSG_A = c(-1.0, 0.0, 1.0),
    check.names = FALSE
  )
  feature_map <- data.frame(
    source_feature = c("ENSG_B", "ENSG_A"),
    target_feature = c("gene_b", "gene_a"),
    annotation_release = c("declared_release", "declared_release"),
    stringsAsFactors = FALSE
  )

  observed <- sb_apply_feature_map(
    new_data,
    "patient_id",
    feature_map,
    target_features = c("gene_a", "gene_b")
  )
  changed <- new_data
  changed$outcome <- 1L - changed$outcome
  changed$unrelated <- -changed$unrelated
  repeated <- sb_apply_feature_map(
    changed,
    "patient_id",
    feature_map,
    target_features = c("gene_a", "gene_b")
  )

  expect_s3_class(observed, "stormbreaker_mapped_features")
  expect_identical(observed$status, "success")
  expect_identical(observed$id_col, "patient_id")
  expect_identical(observed$n, 3L)
  expect_identical(
    observed$mapping$target_feature,
    c("gene_a", "gene_b")
  )
  expect_identical(
    observed$mapping$source_feature,
    c("ENSG_A", "ENSG_B")
  )
  expect_identical(
    observed$mapping$annotation_release,
    rep("declared_release", 2L)
  )
  expect_identical(
    names(observed$data),
    c("patient_id", "gene_a", "gene_b")
  )
  expect_identical(observed$data$patient_id, new_data$patient_id)
  expect_identical(observed$data$gene_a, new_data$ENSG_A)
  expect_identical(observed$data$gene_b, new_data$ENSG_B)
  expect_identical(repeated, observed)

  signature <- sb_linear_signature(
    "mapped_signature",
    "binary",
    coefficients = c(gene_a = 0.5, gene_b = -0.25),
    intercept = -0.1
  )
  prediction <- sb_predict_signature(
    signature,
    observed$data,
    "patient_id"
  )
  expected_score <- -0.1 + 0.5 * new_data$ENSG_A -
    0.25 * new_data$ENSG_B
  expect_identical(prediction$status, "success")
  expect_equal(
    prediction$predictions$probability,
    stats::plogis(expected_score),
    tolerance = 0
  )
})

test_that("mapped features feed an unchanged final prediction", {
  development <- make_sb_binary_fixture(n = 60L, seed = 9951L)
  task <- sb_binary_task(
    development,
    "patient_id",
    "outcome",
    c("gene_a", "gene_b", "gene_c")
  )
  fit <- sb_fit_final(task, sb_pipeline(learner = "logistic"))
  direct_data <- make_sb_binary_fixture(n = 6L, seed = 9952L)
  direct_data$patient_id <- paste0("new_", direct_data$patient_id)
  source_data <- direct_data
  source_names <- paste0("source_", seq_along(fit$required_features))
  names(source_data)[match(fit$required_features, names(source_data))] <-
    source_names
  feature_map <- data.frame(
    source_feature = rev(source_names),
    target_feature = rev(fit$required_features),
    stringsAsFactors = FALSE
  )

  mapped <- sb_apply_feature_map(
    source_data,
    "patient_id",
    feature_map,
    fit$required_features
  )
  direct <- sb_predict_final(fit, direct_data, "patient_id")
  through_map <- sb_predict_final(fit, mapped$data, "patient_id")

  expect_identical(through_map, direct)
})

test_that("feature mapping refuses ambiguity and missing sources", {
  data <- data.frame(
    patient_id = c("p1", "p2"),
    source_a = c(1, 2),
    source_b = c(3, 4),
    check.names = FALSE
  )
  one_to_many <- data.frame(
    source_feature = c("source_a", "source_a"),
    target_feature = c("gene_a", "gene_b")
  )
  many_to_one <- data.frame(
    source_feature = c("source_a", "source_b"),
    target_feature = c("gene_a", "gene_a")
  )
  incomplete <- data.frame(
    source_feature = "source_a",
    target_feature = "gene_a"
  )
  missing_source <- data.frame(
    source_feature = c("source_a", "source_c"),
    target_feature = c("gene_a", "gene_b")
  )

  expect_s3_class(
    capture_sb_error(sb_apply_feature_map(
      data, "patient_id", one_to_many, c("gene_a", "gene_b")
    )),
    "stormbreaker_non_bijective_feature_map"
  )
  expect_s3_class(
    capture_sb_error(sb_apply_feature_map(
      data, "patient_id", many_to_one, "gene_a"
    )),
    "stormbreaker_non_bijective_feature_map"
  )
  expect_s3_class(
    capture_sb_error(sb_apply_feature_map(
      data, "patient_id", incomplete, c("gene_a", "gene_b")
    )),
    "stormbreaker_feature_map_target_mismatch"
  )
  expect_s3_class(
    capture_sb_error(sb_apply_feature_map(
      data, "patient_id", missing_source, c("gene_a", "gene_b")
    )),
    "stormbreaker_missing_feature_map_source"
  )
})

test_that("feature mapping never aggregates or repairs predictor values", {
  nonnumeric <- data.frame(
    patient_id = c("p1", "p2"),
    source_a = c("low", "high")
  )
  missing <- data.frame(
    patient_id = c("p1", "p2"),
    source_a = c(1, NA_real_)
  )
  mapping <- data.frame(
    source_feature = "source_a",
    target_feature = "gene_a"
  )

  expect_s3_class(
    capture_sb_error(sb_apply_feature_map(
      nonnumeric, "patient_id", mapping, "gene_a"
    )),
    "stormbreaker_bad_adapter_predictor_type"
  )
  expect_s3_class(
    capture_sb_error(sb_apply_feature_map(
      missing, "patient_id", mapping, "gene_a"
    )),
    "stormbreaker_bad_adapter_predictor_values"
  )
})
