make_component_survival_task <- function() {
  set.seed(8901L)
  x <- as.data.frame(matrix(stats::rnorm(90L * 12L), nrow = 90L))
  names(x) <- sprintf("gene_%02d", seq_len(12L))
  lp <- x[[1]] - x[[2]] + .6 * x[[3]]
  et <- stats::rexp(90L, rate = exp(lp) / 15)
  ct <- stats::rexp(90L, rate = 1 / 25)
  d <- data.frame(id = seq_len(90L), time = pmax(pmin(et, ct), .05), event = as.integer(et <= ct), x)
  sb_surv_task(d, "id", "time", "event", names(x), "month")
}
