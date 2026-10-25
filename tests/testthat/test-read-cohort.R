csv_fixture <- function(lines) {
  path <- tempfile(fileext = ".csv")
  writeLines(lines, path)
  path
}

test_that("CSV cohort joins by preserved IDs, not row position", {
  expression <- csv_fixture(c("patient_id,EGFR,TP53", "001,1.2,4", "002,2.3,5"))
  clinical <- csv_fixture(c("patient_id,outcome", "002,0", "001,1"))
  d <- sb_read_cohort(expression, clinical)
  expect_identical(d$patient_id, c("001", "002"))
  expect_equal(d$outcome, c(1, 0))
  expect_equal(d$EGFR, c(1.2, 2.3))
  expect_named(sb_read_cohort(expression), c("patient_id", "EGFR", "TP53"))
  expect_s3_class(sb_binary_task(d, "patient_id", "outcome", c("EGFR", "TP53")),
    "stormbreaker_binary_task")
})

test_that("unmatched samples identify both missing and extra clinical IDs", {
  expression <- csv_fixture(c("patient_id,EGFR", "A,1", "B,2"))
  clinical <- csv_fixture(c("patient_id,outcome", "B,0", "C,1"))
  expect_error(sb_read_cohort(expression, clinical), "A.*C",
    class = "stormbreaker_unmatched_cohort_ids")
})

test_that("CSV input rejects duplicate IDs, columns and overlapping values", {
  expression <- csv_fixture(c("patient_id,EGFR", "A,1", "A,2"))
  expect_error(sb_read_cohort(expression), "A", class = "stormbreaker_duplicate_id")
  duplicate <- csv_fixture(c("patient_id,EGFR,EGFR", "A,1,2"))
  expect_error(sb_read_cohort(duplicate), "EGFR", class = "stormbreaker_duplicate_column_names")
  expression <- csv_fixture(c("patient_id,EGFR", "A,1", "B,2"))
  clinical <- csv_fixture(c("patient_id,EGFR", "A,0", "B,1"))
  expect_error(sb_read_cohort(expression, clinical), "EGFR",
    class = "stormbreaker_overlapping_cohort_columns")
  no_id <- csv_fixture(c("sample,EGFR", "A,1", "B,2"))
  expect_error(sb_read_cohort(no_id), "patient_id", class = "stormbreaker_missing_endpoint_columns")
})

test_that("missing measurements and character outcomes are not silently repaired", {
  expression <- csv_fixture(c("patient_id,EGFR", "A,NA", "B,2"))
  clinical <- csv_fixture(c("patient_id,outcome", "A,control", "B,cancer"))
  d <- sb_read_cohort(expression, clinical)
  expect_identical(is.na(d$EGFR), c(TRUE, FALSE))
  expect_identical(d$outcome, c("control", "cancer"))
  expect_error(sb_binary_task(d, "patient_id", "outcome", "EGFR"),
    "0/1", class = "stormbreaker_bad_binary_outcome")
})
