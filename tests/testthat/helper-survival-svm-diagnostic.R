make_survival_svm_diagnostic_task <- function() {
  set.seed(11032L)
  x <- as.data.frame(matrix(rnorm(72L * 4L), nrow = 72L))
  names(x) <- c("gene A", "gene-B", "gene_c", "gene_d")
  lp <- .7 * x[[1]] - .5 * x[[2]] + .35 * x[[3]]^2 - .25 * x[[1]] * x[[4]]
  et <- rexp(72L, rate = .07 * exp(lp))
  ct <- rexp(72L, rate = .03)
  d <- data.frame(id = seq_len(72L), time = pmax(pmin(et, ct), .01), event = as.integer(et <= ct), x, check.names = FALSE)
  sb_surv_task(d, "id", "time", "event", names(x), "month")
}
