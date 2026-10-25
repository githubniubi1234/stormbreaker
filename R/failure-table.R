sb_failure_empty <- function() {
  data.frame(cohort_label=character(),validation=character(),task_type=character(),
    pipeline_id=character(),issue_type=character(),stage=character(),fold=integer(),
    metric=character(),status=character(),run_status=character(),
    source_error_stage=character(),error_class=character(),error_message=character(),
    estimate=numeric(),source_row=integer(),prediction_horizon=numeric(),time_unit=character(),
    stringsAsFactors=FALSE)
}

sb_failure_results <- function(x) {
  external <- inherits(x,c("stormbreaker_external_evaluation","stormbreaker_external_comparison"))
  is_cv <- function(value) inherits(value,c("stormbreaker_cv_result","stormbreaker_binary_cv_result"))
  results <- if (inherits(x,"stormbreaker_cv_comparison")) x$results else if (
    inherits(x,"stormbreaker_external_comparison")) x$evaluations else if (
      external || is_cv(x)) list(x) else NULL
  predicate <- if (external) function(value) inherits(value,"stormbreaker_external_evaluation") else is_cv
  if (!is.list(results) || !length(results) || !all(vapply(results,predicate,logical(1)))) {
    sb_abort("Supply supported CV or external results/comparisons.","stormbreaker_bad_failure_view")
  }
  list(results=results,external=external)
}

sb_failure_rows <- function(x,cohort_label) {
  if (sb_is_composite_view(x) || !is.null(attr(x, "stormbreaker_view_context"))) {
    return(sb_composite_failures(x, cohort_label))
  }
  parts <- sb_failure_results(x)
  rows <- list()
  value_or <- function(value,default) if(is.null(value) || length(value)!=1L) default else value
  problem <- function(status) !is.null(status) && length(status)==1L && !is.na(status) &&
    !status %in% c("success","not_requested","not_applicable")
  for (result in parts$results) {
    descriptor <- if(parts$external) result$endpoint else result$task
    type <- if(parts$external) result$task_type else if(inherits(result,"stormbreaker_binary_cv_result")) "binary" else "survival"
    pipeline <- if(parts$external) result$pipeline_id else result$pipeline$id
    add <- function(issue,stage,status,source,fold=NA_integer_,metric=NA_character_,
                    row=NA_integer_,estimate=NA_real_) {
      rows[[length(rows)+1L]] <<- data.frame(
        cohort_label=cohort_label,validation=if(parts$external) "external_evaluation" else "development_cv",
        task_type=type,pipeline_id=pipeline,issue_type=issue,stage=stage,fold=as.integer(fold),
        metric=metric,status=value_or(status,NA_character_),run_status=result$status,
        source_error_stage=value_or(source$error_stage,NA_character_),
        error_class=value_or(source$error_class,NA_character_),
        error_message=value_or(source$error_message,NA_character_),estimate=estimate,
        source_row=as.integer(row),
        prediction_horizon=value_or(source$prediction_horizon,value_or(descriptor$prediction_horizon,NA_real_)),
        time_unit=value_or(source$time_unit,value_or(descriptor$time_unit,NA_character_)),
        stringsAsFactors=FALSE)
    }
    if (problem(result$status)) {
      add("pipeline",if(parts$external) "external_evaluation" else "cv",result$status,result)
    }
    if (!parts$external && is.data.frame(result$runs)) {
      for(i in which(vapply(as.list(result$runs$status),problem,logical(1)))) {
        source <- as.list(result$runs[i,,drop=FALSE])
        add("run","cv_fold",source$status,source,fold=value_or(source$fold,NA_integer_))
      }
    }
    visit <- function(value,path) {
      if (!is.list(value) || is.data.frame(value)) return(invisible(NULL))
      if(problem(value$status)) {
        add("metric","metric",value$status,value,metric=path,
            estimate=value_or(value$estimate,NA_real_))
      } else if(identical(value$status,"success") && is.numeric(value$estimate) &&
                length(value$estimate)==1L && !is.finite(value$estimate)) {
        add("metric_nonfinite","metric",value$status,value,metric=path,estimate=value$estimate)
      }
      for(name in names(value)) {
        child <- value[[name]]
        if(is.data.frame(child) && "status" %in% names(child)) {
          # Only saved aggregate curve/fold diagnostic rows, never predictions.
          if(!name %in% c("curve","per_fold")) next
          for(i in which(vapply(as.list(child$status),problem,logical(1)))) {
            source <- as.list(child[i,,drop=FALSE])
            add(if(name=="curve") "curve_row" else "metric_fold",name,
                source$status,source,fold=value_or(source$fold,NA_integer_),metric=path,row=i)
          }
        } else if(is.list(child) && !is.data.frame(child)) {
          visit(child,paste(path,name,sep="."))
        }
      }
      invisible(NULL)
    }
    for(metric in names(result$metrics)) visit(result$metrics[[metric]],metric)
  }
  if(!length(rows)) return(sb_failure_empty())
  result <- do.call(rbind,rows); rownames(result) <- NULL
  result
}

#' Extract saved pipeline and metric failures
#'
#' Read failure records without executing, modifying or dropping model results.
#' @param x A CV result/comparison, external evaluation/comparison, whole
#'   repeated CV, IECV or stacking CV result, or a uniquely
#'   named list of these objects representing distinct report sections.
#' @param cohort_label Optional non-empty cohort display label for a single
#'   result or comparison. For a named list, leave NULL and use its names.
#' @return A data frame with one row per saved pipeline/run/metric issue.
#'   Pipeline summaries and detailed fold or metric issues are distinct rows;
#'   the number of rows is not the number of failed pipelines. Includes cohort,
#'   validation boundary, pipeline, issue type, structural stage, fold, metric,
#'   source status, complete-run status, source error stage/class/message, scalar
#'   estimate and saved horizon/unit. No patient prediction rows are included.
#' @details Explicitly unrequested/not-applicable components are not failures.
#'   Composite results add repetition, held-out cohort, member/stack role,
#'   validation scheme and structural `view_id` columns. IECV is development
#'   validation. Mixed named lists keep these columns with missing values for
#'   ordinary results. Saved failures are neither averaged nor deduplicated.
#'   A saved success status with a nonfinite scalar estimate is labeled through
#'   issue_type = "metric_nonfinite" without changing its source status.
#'   Missing source stages and reasons remain NA. Structural stage "cv_fold"
#'   identifies a fold record, not an inferred fitting or prediction error stage.
#'   Empty results return the same columns with zero rows. Raw error messages
#'   can contain sensitive information: review before sharing or opening in
#'   software that may interpret cell contents as formulas.
#' @seealso [sb_export_failures()]
#' @export
sb_failure_table <- function(x,cohort_label=NULL) {
  label <- sb_result_cohort_label(cohort_label)
  supported <- inherits(x,c("stormbreaker_cv_result","stormbreaker_binary_cv_result",
    "stormbreaker_cv_comparison","stormbreaker_external_evaluation","stormbreaker_external_comparison"))
  if(supported || sb_is_composite_view(x)) return(sb_failure_rows(x,label))
  if(!is.null(cohort_label) || !is.list(x) || !length(x) || is.null(names(x)) ||
     anyNA(names(x)) || any(!nzchar(trimws(names(x)))) || anyDuplicated(names(x))) {
    sb_abort("Supply a supported result or a uniquely named non-empty result list.",
             "stormbreaker_bad_failure_view")
  }
  rows <- lapply(seq_along(x),function(i) sb_failure_rows(x[[i]],names(x)[i]))
  output <- sb_view_bind(rows); rownames(output) <- NULL
  output
}

sb_new_output_path <- function(file,extension) {
  if(!is.character(file) || length(file)!=1L || is.na(file) || !nzchar(trimws(file)) ||
     !endsWith(tolower(file),paste0(".",extension)) || !dir.exists(dirname(file))) {
    sb_abort(paste0("Supply a new .",extension," file in an existing directory."),
             "stormbreaker_bad_output_path")
  }
  if(file.exists(file)) sb_abort("Output already exists; choose a new filename.","stormbreaker_output_exists")
  file.path(normalizePath(dirname(file),winslash="/",mustWork=TRUE),basename(file))
}

#' Export saved failures to a separate CSV
#' @param x A result or named result list accepted by [sb_failure_table()].
#' @param file A new CSV path in an existing directory. Existing files are refused.
#' @param cohort_label Optional label, forwarded to [sb_failure_table()].
#' @return Invisibly, the normalized absolute CSV path. A header-only CSV is
#'   written when there are no issues. Uses base R and does not require ggplot2.
#' @details The CSV preserves saved error text, including quoted commas and
#'   line breaks, plus repetition, held-out cohort and model-role columns for
#'   composite objects. It excludes patient-level prediction rows but raw error text
#'   still requires privacy review. Results are not retrained or recomputed.
#' @export
sb_export_failures <- function(x,file,cohort_label=NULL) {
  file <- sb_new_output_path(file,"csv")
  tab <- sb_failure_table(x,cohort_label)
  file <- sb_new_output_path(file,"csv")
  utils::write.csv(tab,file,row.names=FALSE,na="NA",fileEncoding="UTF-8")
  invisible(normalizePath(file,winslash="/",mustWork=TRUE))
}
