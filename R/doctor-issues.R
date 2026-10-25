#' Explain issues from a doctor-facing model run
#' @param x Result of bi_run or surv_run.
#' @return Data frame with model, stage, category, Chinese explanation and
#'   original reason. Plot unavailability is separate from fitting failures.
#'   Unknown failures remain unclassified; raw evidence is never replaced.
#' @export
sb_run_issues <- function(x) {
  if (!is.list(x) || !inherits(x$comparison, "stormbreaker_cv_comparison")) stop("Supply a bi_run or surv_run result.")
  rules <- list(
    dependency = c("Package .*required|Install |not installed|there is no package", "\u7f3a\u5c11\u7b97\u6cd5\u4f9d\u8d56\uff0c\u8bf7\u5b89\u88c5\u63d0\u793a\u7684\u8f6f\u4ef6\u5305\u540e\u518d\u8fd0\u884c\u3002"),
    features = c("no features|no predictors|no genes|empty.*select|retain 1 or less|No features met|nonempty set", "\u7b5b\u9009\u540e\u53d8\u91cf\u4e0d\u8db3\uff1b\u53ef\u51cf\u5c11\u7b5b\u9009\u6b65\u9aa4\uff0c\u4e0d\u80fd\u968f\u610f\u8865\u57fa\u56e0\u51d1\u6570\u3002"),
    evidence = c("not retained|not captured|not saved|capture_diagnostics|disabled importance|OOB.*disabled", "\u6a21\u578b\u6ca1\u6709\u4fdd\u5b58\u8fd9\u9879\u7ed8\u56fe\u6570\u636e\uff1b\u4e0d\u5f71\u54cd\u5df2\u6210\u529f\u751f\u6210\u7684\u5176\u4ed6\u7ed3\u679c\u3002"),
    data = c("single.class|both.*classes|0/1|non-finite|finite|missing values|schema|columns.*match|same.*columns", "\u68c0\u67e5\u5206\u7ec4\u3001\u7f3a\u5931\u503c\u53ca\u8bad\u7ec3\u4e0e\u9a8c\u8bc1\u6570\u636e\u7684\u5217\u540d\uff1b\u4e5f\u53ef\u80fd\u662f\u8be5\u7b97\u6cd5\u5728\u8fd9\u4efd\u6570\u636e\u4e0a\u65e0\u6cd5\u4f30\u8ba1\u3002"),
    unsupported = c("not admitted|not supported|not implemented|No native|nonlinear kernel|cannot.*group|No specialized|no single exact|not a variable.selection", "\u8fd9\u4e2a\u7ec4\u5408\u6216\u56fe\u5f62\u4e0d\u9002\u7528\uff1b\u4e0d\u4f1a\u7528\u53e6\u4e00\u79cd\u7b97\u6cd5\u6216\u865a\u6784\u56fe\u5f62\u4ee3\u66ff\u3002")
  )
  rows <- list()
  add <- function(model, stage, reason) {
    if (is.null(reason) || is.na(reason) || !nzchar(reason)) reason <- "No detailed reason was saved."
    category <- "other"
    explanation <- "\u9700\u8981\u67e5\u770b\u539f\u59cb\u539f\u56e0\uff0c\u6682\u4e0d\u80fd\u81ea\u52a8\u5224\u65ad\u3002"
    for (k in names(rules)) {
      if (grepl(rules[[k]][1], reason, ignore.case = TRUE)) {
        category <- k
        explanation <- rules[[k]][2]
        break
      }
    }
    rows[[length(rows) + 1L]] <<- data.frame(
      model = model, stage = stage,
      category = category, explanation = explanation, reason = reason
    )
  }
  failures <- sb_failure_table(x$comparison)
  labels <- stats::setNames(as.character(names(x$pipelines)), vapply(x$pipelines, function(p) p$id, character(1)))
  for (i in seq_len(nrow(failures))) {
    f <- failures[i, ]
    if (f$issue_type == "pipeline" && any(failures$pipeline_id == f$pipeline_id & failures$issue_type == "run")) next
    model <- unname(labels[f$pipeline_id])
    if (is.na(model)) model <- f$pipeline_id
    reason <- f$error_message
    if (is.na(reason) || !nzchar(reason)) reason <- paste(f$metric, f$status)
    add(model, if (is.na(f$fold)) f$stage else paste0(f$stage, " ", f$fold), reason)
  }
  for (n in names(x$unavailable)) add(n, "combination", x$unavailable[[n]])
  for (n in names(x$plot_errors)) add("", "overview_plot", x$plot_errors[[n]])
  for (n in names(x$details)) {
    d <- x$details[[n]]
    if (d$fit$status != "success") add(n, "final_fit", d$fit$error_message)
    for (v in names(d$figures$unavailable)) {
      issue <- d$figures$unavailable[[v]]
      if (is.data.frame(issue)) {
        for (i in which(issue$status != "available")) {
          add(n, paste("plot", issue$view[i]), issue$reason[i])
        }
      } else {
        add(n, paste("plot", v), issue)
      }
    }
    native <- d$figures$native$status
    if (is.null(d$figures$unavailable$native) && is.data.frame(native)) for (i in which(native$status != "available")) add(n, paste("plot", native$view[i]), native$reason[i])
  }
  if (!length(rows)) {
    return(data.frame(
      model = character(), stage = character(),
      category = character(), explanation = character(), reason = character()
    ))
  }
  do.call(rbind, rows)
}
