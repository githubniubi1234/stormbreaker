#' Inspect the original tree-package classifier
#' @param fit Successful tree_binary final fit.
#' @param type structure or leaf_probabilities.
#' @return All native nodes and split text, or terminal node class probabilities
#'   and original sample counts. Probabilities are training-leaf proportions,
#'   not held-out calibration. Native split strings retain original precision;
#'   the complete frame is returned, without surrogate fitting or pruning.
#' @export
sb_classification_tree_diagnostics <- function(fit, type = "structure") {
  fit <- sb_validate_final_fit_for_prediction(fit)
  s <- fit$learner_state
  if (s$learner != "tree_binary" || !is.character(type) || length(type) != 1L ||
      is.na(type) || !type %in% c("structure", "leaf_probabilities")) sb_abort(
    "Choose structure or leaf_probabilities for the tree-package classifier.", "stormbreaker_unsupported_diagnostic")
  e <- unserialize(s$engine_raw)
  f <- e$frame
  ids <- as.integer(rownames(f))
  if (!inherits(e, "tree") || !identical(attr(e,"ylevels"), c("0","1")) ||
      !is.matrix(f$yprob) || !identical(colnames(f$yprob), c("0","1"))) sb_abort(
    "Native classification tree or class probabilities are invalid.", "stormbreaker_invalid_diagnostic")
  internal <- sb_tree_binary_internal_names(length(s$features))
  vars <- as.character(f$var)
  terminal <- vars == "<leaf>"
  if (any(!vars[!terminal] %in% internal)) sb_abort("A native split feature is unknown.", "stormbreaker_invalid_diagnostic")
  feature <- s$features[match(vars,internal)]
  parent <- ifelse(ids == 1L, NA_integer_, ids %/% 2L)
  tab <- data.frame(node=ids,parent=parent,depth=floor(log2(ids)),terminal=terminal,
    incoming_rule="Root",split_feature=feature,n=f$n,probability=f$yprob[,"1"],
    x=NA_real_,y=-floor(log2(ids)))
  for (i in which(ids != 1L)) {
    j <- match(parent[i],ids)
    tab$incoming_rule[i] <- paste(feature[j], f$splits[j,if(ids[i] %% 2L == 0L) 1L else 2L])
  }
  position <- 0L
  layout <- function(id) {
    i <- match(id,ids)
    if (terminal[i]) { position <<- position+1L; tab$x[i] <<- position } else {
      layout(id*2L); layout(id*2L+1L)
      tab$x[i] <<- mean(range(tab$x[match(c(id*2L,id*2L+1L),ids)]))
    }
  }
  layout(1L)
  if (type == "leaf_probabilities") tab <- tab[tab$terminal, , drop=FALSE]
  structure(list(type=if(type=="structure") "ctree_structure" else "classification_tree_leaf_probabilities",
    table=tab,native_frame=f,learner=s$learner,pipeline_id=fit$pipeline$id,
    source="Original tree frame, native split strings and class probabilities",
    interpretation="Original training node counts and leaf probabilities; not held-out performance. No pruning or refitting."),
    class="stormbreaker_algorithm_diagnostic")
}

sb_plot_classification_tree_probabilities <- function(x) {
  node <- probability <- n <- NULL
  ggplot2::ggplot(x$table,ggplot2::aes(factor(node),probability,size=n)) +
    ggplot2::geom_point(colour="#2878A0") + ggplot2::scale_y_continuous(limits=c(0,1)) +
    ggplot2::labs(title="Terminal-node probabilities",x="Node",y="Class-1 probability",size="Training rows") + sb_viz_theme()
}
