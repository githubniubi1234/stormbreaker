test_that("rotation forest diagnostics preserve original bases and native predictions", {
  skip_if_not_installed("rotationForest")
  d <- sb_example_binary(160L,seed=108L)
  features <- c("gene_c","gene_a","gene_b")
  fit <- sb_fit_final(sb_binary_task(d,"patient_id","outcome",features),sb_pipeline(learner="rotation_forest_binary"))
  expect_identical(fit$status,"success")
  before <- serialize(fit,NULL,version=2L)
  e <- unserialize(fit$learner_state$engine_raw)
  loadings <- sb_rotation_forest_diagnostics(fit)
  sizes <- sb_rotation_forest_diagnostics(fit,"tree_size")
  query <- d[c(1L,7L), ]; query$gene_c <- query$gene_c+.17
  predictions <- sb_rotation_forest_diagnostics(fit,"tree_predictions",query)
  for(i in seq_along(e$models)) {
    tab <- loadings$table[loadings$table$tree==i,]
    m <- matrix(tab$loading,nrow=length(features))
    expect_equal(m,unname(as.matrix(e$loadings[[i]])))
    expect_identical(tab$feature,rep(features,length(features)))
    ids <- as.integer(rownames(e$models[[i]]$frame))
    expect_equal(sizes$table$terminal_nodes[i],sum(!((2L*ids) %in% ids)))
    rotated <- as.data.frame(as.matrix(query[features]) %*% as.matrix(e$loadings[[i]]))
    native <- predict(e$models[[i]],rotated,type="prob")[,"1"]
    expect_equal(predictions$table$probability[predictions$table$tree==i],unname(native))
  }
  expect_equal(as.numeric(tapply(predictions$table$probability,predictions$table$query,mean)),
    sb_predict_final(fit,query,"patient_id")$predictions$probability,tolerance=1e-12)
  for(z in list(loadings,sizes,predictions)) expect_s3_class(ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(z)),"ggplot_built")
  expect_identical(sb_native_diagnostics(fit,complete=TRUE)$status$status,c("available","available","unavailable"))
  expect_identical(sb_native_diagnostics(fit,complete=TRUE,newdata=query)$status$status,rep("available",3))
  expect_identical(serialize(fit,NULL,version=2L),before)
})
