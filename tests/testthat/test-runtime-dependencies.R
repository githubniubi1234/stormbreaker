test_that("dependency descriptions separate R packages and Python modules", {
  expect_identical(
    sb_dependency_components("mboost + partykit"),
    list(
      r_packages = c("mboost", "partykit"),
      python_modules = character()
    )
  )
  expect_identical(
    sb_dependency_components(
      "reticulate + Python numpy + pandas + torch + torchtuples + pycox"
    ),
    list(
      r_packages = "reticulate",
      python_modules = c(
        "numpy", "pandas", "torch", "torchtuples", "pycox"
      )
    )
  )
})

test_that("dependency status distinguishes missing R and unchecked Python", {
  local_mocked_bindings(
    sb_r_package_installed = function(package) {
      package %in% c("mboost", "reticulate")
    }
  )

  missing <- sb_dependency_status_row("mboost + partykit")
  python <- sb_dependency_status_row(
    "reticulate + Python numpy + pandas + torch + torchtuples + pycox"
  )

  expect_identical(missing$missing_r_packages, "partykit")
  expect_identical(missing$dependency_status, "missing_r_packages")
  expect_identical(python$missing_r_packages, NA_character_)
  expect_identical(python$dependency_status, "python_not_checked")
})

test_that("installed checks require one logical value", {
  error <- tryCatch(
    sb_pipeline_registry(check_installed = 1),
    error = function(error) error
  )

  expect_s3_class(error, "stormbreaker_bad_registry_check_installed")
})
