library(stormbreaker)

d <- bi_example(n = 160, seed = 42)
genes <- c("gene_a", "gene_b", "gene_c")
task <- bi_task(d, "patient_id", "outcome", genes,
  decision_thresholds = seq(.05, .8, by = .05))
pipeline <- sb_pipeline(learner = "logistic")
fit <- fit_final(task, pipeline)
pred <- predict_final(fit, d[c("patient_id", genes)], "patient_id")
training_roc <- plot_pred_roc(d$outcome, pred$predictions$probability,
  kind = "apparent", cohort = "Training", ci_reps = 1000)

external <- bi_example(n = 120, seed = 43)
external$patient_id <- paste0("external_", external$patient_id)
ep <- predict_final(fit, external[c("patient_id", genes)], "patient_id")
evaluation <- eval_external(ep, bi_endpoint(external, "patient_id", "outcome"))
external_roc <- plot_roc(evaluation, cohort_label = "Validation", ci_reps = 1000)

cv <- run_cv(task, pipeline, folds = 3, seed = 42)
oof_roc <- plot_roc(cv)
nomogram <- plot_nomogram(nomogram(fit, d), reference_data = d,
  patient = d[1, , drop = FALSE])
calibration <- sb_bootstrap_calibration(task, pipeline, reps = 1000, seed = 42)
calibration_plot <- plot_boot_calib(calibration)
dca <- plot_dca(cv, standardized = TRUE, cost_benefit_axis = TRUE)


performance <- data.frame(model = rep(c("Model A", "Model B", "Model C"), each = 3),
  cohort = rep(c("Development", "Validation1", "Validation2"), 3),
  auc = c(.85, .80, .82, .90, .78, .76, .82, .79, .80))
ranking <- plot_rank(performance, c("Validation1", "Validation2"))


training_roc
ranking$plot
nomogram
calibration_plot
dca
