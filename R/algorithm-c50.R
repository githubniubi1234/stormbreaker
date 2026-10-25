#' Inspect captured C5.0 tree or rule ensembles
#' @param fit Successful C5.0 final fit with original diagnostic capture.
#' @param type usage, splits or size.
#' @return Native C5imp percentages or original trial sizes (leaves for trees,
#'   rules for rule ensembles). Native unparseable importance is unavailable,
#'   not replaced with fabricated zeros. Compact prediction state deliberately
#'   omits these summaries, so original capture is required. No refitting occurs.
#' @export
sb_c50_diagnostics <- function(fit,type="usage") {
  fit <- sb_validate_final_fit_for_prediction(fit)
  s <- fit$learner_state
  if (!s$learner %in% c("c50_binary","c50_rules_binary") ||
      !is.character(type) || length(type)!=1L || is.na(type) || !type %in% c("usage","splits","size")) sb_abort(
    "Choose usage, splits or size for a C5.0 ensemble.","stormbreaker_unsupported_diagnostic")
  e <- fit$diagnostic_artifacts$learner$engine
  if(is.null(e)) sb_abort("Original C5.0 summaries were not captured.","stormbreaker_missing_diagnostic_artifacts")
  if(!inherits(e,"C5.0")) sb_abort("Captured C5.0 engine is invalid.","stormbreaker_invalid_diagnostic")
  if(type=="size") {
    if(length(e$size)!=e$trials[["Actual"]] || any(!is.finite(e$size)) || any(e$size<0)) sb_abort(
      "Native C5.0 sizes are invalid.","stormbreaker_invalid_diagnostic")
    tab <- data.frame(trial=seq_along(e$size),size=e$size)
  } else {
    native <- tryCatch(C50::C5imp(e,metric=type),error=identity)
    if(inherits(native,"error")) sb_abort(paste("Native C5.0 importance unavailable:",conditionMessage(native)),
      "stormbreaker_unavailable_diagnostic")
    if(is.data.frame(native) && "Overall" %in% names(native) && all(is.nan(native$Overall))) sb_abort(
      "Native C5.0 split percentages are undefined because no split importance was reported.","stormbreaker_unavailable_diagnostic")
    if(!is.data.frame(native) || !"Overall" %in% names(native) ||
        !setequal(rownames(native),s$features) || any(!is.finite(native$Overall))) sb_abort(
      "Native C5.0 feature summary is invalid.","stormbreaker_invalid_diagnostic")
    tab <- data.frame(feature=s$features,value=native[s$features,"Overall"])
  }
  structure(list(type=paste0("c50_",type),table=tab,learner=s$learner,pipeline_id=fit$pipeline$id,
    rules=isTRUE(e$rbm),source="Original captured C5.0 summaries and native C5imp",
    interpretation="Native usage and split percentages are different summaries, not permutation importance or causal effects. Trial sizes are not performance. No refitting or missing-value reconstruction."),
    class="stormbreaker_algorithm_diagnostic")
}

sb_plot_c50_diagnostic <- function(x) {
  feature <- value <- trial <- size <- NULL
  if(x$type=="c50_size") {
    p <- ggplot2::ggplot(x$table,ggplot2::aes(trial,size)) + ggplot2::geom_point(colour="#2878A0") +
      ggplot2::scale_x_continuous(breaks=x$table$trial) +
      ggplot2::scale_y_continuous(breaks=function(z) unique(round(pretty(z)))) +
      ggplot2::labs(title="C5.0 trial sizes",x="Trial",y=if(x$rules) "Number of rules" else "Number of leaves")
  } else {
    d <- x$table; d$feature <- factor(d$feature,levels=d$feature[order(d$value)])
    p <- ggplot2::ggplot(d,ggplot2::aes(value,feature)) + ggplot2::geom_col(fill="#2878A0") +
      ggplot2::labs(title="C5.0 native feature summary",x=if(x$type=="c50_usage") "Native usage (%)" else "Native splits (%)",y=NULL)
  }
  p + sb_viz_theme()
}
