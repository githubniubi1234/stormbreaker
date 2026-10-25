# Compare cohorts and view results

sb_compare_cohorts() accepts a development task and named cohort tasks with
compatible endpoints and predictors. Survival tasks must also agree on units and
horizons. Patient identifiers must be unique across these cohorts; do not rename
repeated records to manufacture independent samples.

Each route is cross-validated in development, fitted on development data and
applied to the supplied cohorts. Development CV, validation and test roles remain
separate. The result viewer displays summary metrics and configurations, not
patient-level records. Existing evaluation results can be combined without fitting.

sb_view_results() writes a standalone HTML page. Use a new file path and
open=FALSE for noninteractive workflows. Model names are readable while original
IDs and settings are retained. A historical _fixed suffix identifies a predefined
configuration, not a new algorithm.

Sorting and importing are display operations, not model fitting. Validation means
use only the displayed validation cohorts with compatible metrics and horizons;
they exclude development and test cohorts. Missing or unavailable metrics remain
unavailable. CSV export retains all filtered routes; SVG/print follow the displayed
page. Imported files are read locally in the browser.

## Code examples

The following fragments use the input and object names described above. Supply your prepared data and the required engines before running them.

### Example 1

```r
source(system.file("examples", "multi-cohort.R", package = "stormbreaker"))
```

### Example 2

```r
sb_view_results(multi_cohort_result)
sb_view_results(multi_cohort_result, path = "my-results.html", open = FALSE)
```

### Example 3

```r

sb_view_results(
  list(Training = dev_result, Validation = validation_result, Test = test_result),
  roles = c(Training = "development_cv", Validation = "validation", Test = "test")
)
```
