# STORMBREAKER 0.3.47

## Bug fixes

- `model_shap()` and `tree_shap()` reuse validated model state across permutation
  predictions, avoiding repeated full-model validation.
- `surv_risk_plots()` returns log-rank computation status and diagnostic messages
  alongside P values.
- `finish_run()` includes log-rank status and calculation errors in survival
  risk-group summaries and exported issue tables.

## Documentation

- Updated installation instructions and binary and survival quick-start examples.
- Updated the six-scenario illustrated tutorial and user guides in English.
- Added executable help examples for the main workflows and common views.
- Extended core automated checks to cover tutorial behavior.
- Updated the data-reading example to use a synthetic cohort.
- Removed obsolete example CSVs and pre-generated result pages from the
  distribution.
