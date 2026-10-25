# Contributing

Bug reports should include the package and R versions, the relevant engine
version, expected behavior and a small reproducible example using synthetic data.
Use the repository issue tracker. Do not include patient records or credentials.

## Development checks

Run these commands from the package root:

```r
testthat::test_local(
  filter = "^(shap-prediction-reuse|logrank-status|short-scenarios|read-cohort|example-data)$"
)
testthat::test_local()
```

The first command runs the core regressions and tutorial assertions. The second
runs the complete suite and requires the optional engines used by those tests.
Report failures, skips and unavailable dependencies separately.

When changing documentation comments, regenerate the help files:

```r
roxygen2::roxygenise()
```

Build and check the package:

```text
R CMD build --no-manual .
R CMD check --no-manual stormbreaker_0.3.47.tar.gz
```

Tutorial building needs knitr, rmarkdown, markdown and Pandoc. The six runnable
scenarios use ggplot2, patchwork and glmnet. Run tools/check-scenarios.R after
loading the source package or an isolated installation. On a shared workstation,
use one R test process and OMP_NUM_THREADS=1.

The core GitHub workflow runs regression tests, tutorial assertions, the scenario
scripts and package/documentation checks. A separate Windows job exercises
selected optional engines. A successful subset does not replace a complete local
suite or establish support on an untested platform.

## Changes and documentation

Keep changes focused. Model changes should be checked against native predictions
and saved-fit replay. New behavior needs tests. Short names and original names
must remain consistent. Write user documentation in English and add a NEWS entry
for changes users will notice.

The optional package website is separate from these checks. Before deploying a
website, configure its canonical url in _pkgdown.yml and run
pkgdown::check_pkgdown(). The reference index can be maintained without claiming
that a website has already been deployed.

Retain LICENSE, LICENSE.md and inst/COPYRIGHTS. Do not commit session histories,
.RData, .Rproj.user, credentials, generated logs, local libraries or patient-level
outputs. Built installation archives belong in Release attachments rather than
the source tree.
