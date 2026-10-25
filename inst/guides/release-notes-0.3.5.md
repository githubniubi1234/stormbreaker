# stormbreaker 0.3.5

stormbreaker connects candidate-gene modeling with independent-cohort evaluation
for binary and survival outcomes.

Choose models from a task-specific menu, compare them on shared cross-validation
folds, and inspect the results in an offline leaderboard. The same workflow
supports final fitting, new-patient predictions, calibration, decision curves
and cohort summaries. Identical fold-local Boruta and varSelRF selections can be
reused across learners to reduce repeated computation.

## Changes in 0.3.5

- Stable FDA/PDA probabilities when extreme discriminant distances cause numerical underflow.
- XGBoost state loading through both older and current R interfaces, with JSON storage for newly saved states.
- English and Chinese quick starts, runnable examples and installation guidance for optional engines.

## Installation

R 4.3 or newer is required.

```r
install.packages("remotes")
remotes::install_github("githubniubi1234/stormbreaker", dependencies = NA)
```

Use the attached `.tar.gz` for R source installation or the `.zip` for Windows
with R 4.5. Install optional model engines as needed; see the repository's
installation guide.
