#' Inspect original rotation-forest geometry and predictions
#' @param fit Successful rotation-forest final fit.
#' @param type loadings, tree_size or tree_predictions.
#' @param newdata Explicit query rows required for tree_predictions.
#' @return All original rotation matrices, tree sizes or per-tree probabilities.
#'   Rotated axes are not the original features despite native reused names.
#'   Loading signs are basis conventions, not effect directions or importance.
#'   Per-tree spread is not a confidence interval. No new PCA, fitting or pruning.
#' @export
sb_rotation_forest_diagnostics <- function(fit,type="loadings",newdata=NULL) {
  fit <- sb_validate_final_fit_for_prediction(fit)
  s <- fit$learner_state
  if(s$learner!="rotation_forest_binary" || !is.character(type) || length(type)!=1L ||
      is.na(type) || !type %in% c("loadings","tree_size","tree_predictions")) sb_abort(
    "Choose loadings, tree_size or tree_predictions for rotation forest.","stormbreaker_unsupported_diagnostic")
  e <- unserialize(s$engine_raw)
  if(!inherits(e,"rotationForest") || !identical(e$columnnames,s$features) || length(e$models)!=length(e$loadings)) sb_abort(
    "Original rotation forest is invalid.","stormbreaker_invalid_diagnostic")
  if(type=="loadings") {
    tab <- do.call(rbind,lapply(seq_along(e$loadings),function(i) {
      m <- as.matrix(e$loadings[[i]])
      if(!identical(dim(m),rep(length(s$features),2)) || any(!is.finite(m))) sb_abort(
        "An original rotation matrix is invalid.","stormbreaker_invalid_diagnostic")
      data.frame(tree=i,feature=rep(s$features,ncol(m)),axis=rep(seq_len(ncol(m)),each=nrow(m)),
        native_axis=rep(colnames(m),each=nrow(m)),loading=as.numeric(m))
    }))
  } else if(type=="tree_size") {
    tab <- do.call(rbind,lapply(seq_along(e$models),function(i) {
      f <- e$models[[i]]$frame
      data.frame(tree=i,nodes=nrow(f),terminal_nodes=sum(f$var=="<leaf>"))
    }))
  } else {
    if(is.null(newdata)) sb_abort("Supply explicit queries for rotation-forest tree predictions.","stormbreaker_unavailable_diagnostic")
    x <- sb_select_final_feature_frame(newdata,s$features,"Query data","query fitted features")
    if(!nrow(x)) sb_abort("Queries must contain at least one row.","stormbreaker_invalid_diagnostic")
    native_predict <- get("predict.rotationForest",asNamespace("rotationForest"))
    m <- native_predict(e,x,all=TRUE)
    if(!is.matrix(m) || nrow(m)!=nrow(x) || ncol(m)!=length(e$models) || any(!is.finite(m)) || any(m<0 | m>1)) sb_abort(
      "Native per-tree probabilities are invalid.","stormbreaker_invalid_diagnostic")
    tab <- data.frame(query=rep(seq_len(nrow(m)),ncol(m)),tree=rep(seq_len(ncol(m)),each=nrow(m)),probability=as.numeric(m))
  }
  structure(list(type=paste0("rotation_forest_",type),table=tab,learner=s$learner,pipeline_id=fit$pipeline$id,
    source="Original saved rotation matrices, rpart frames and native per-tree probabilities",
    interpretation="Rotated axes are not original features. Loading signs are basis conventions, not feature effects or importance. Tree variation is not a confidence interval. No PCA refit, pruning or validation."),
    class="stormbreaker_algorithm_diagnostic")
}

sb_plot_rotation_forest_diagnostic <- function(x) {
  axis <- feature <- loading <- tree <- terminal_nodes <- probability <- query <- NULL
  if(x$type=="rotation_forest_loadings") {
    p <- ggplot2::ggplot(x$table,ggplot2::aes(factor(axis),feature,fill=loading)) +
      ggplot2::geom_tile(colour="white",linewidth=.2) + ggplot2::facet_wrap(~tree,ncol=5,labeller=ggplot2::label_both) +
      ggplot2::scale_fill_gradient2(low="#2878A0",mid="white",high="#B86236",midpoint=0,
        guide=ggplot2::guide_colourbar(display="rectangles")) +
      ggplot2::labs(title="Original rotation matrices",x="Rotated axis",y=NULL,fill="Loading")
    size <- c(width=14,height=7)
  } else if(x$type=="rotation_forest_tree_size") {
    p <- ggplot2::ggplot(x$table,ggplot2::aes(tree,terminal_nodes)) + ggplot2::geom_point(colour="#2878A0") +
      ggplot2::scale_x_continuous(breaks=x$table$tree) +
      ggplot2::labs(title="Rotation-forest tree sizes",x="Tree",y="Terminal nodes")
    size <- c(width=11,height=5)
  } else {
    p <- ggplot2::ggplot(x$table,ggplot2::aes(tree,probability)) + ggplot2::geom_point(colour="#2878A0") +
      ggplot2::facet_wrap(~query,labeller=ggplot2::label_both) + ggplot2::scale_y_continuous(limits=c(0,1)) +
      ggplot2::scale_x_continuous(breaks=unique(x$table$tree)) +
      ggplot2::labs(title="Native per-tree predictions",x="Tree",y="Class-1 probability")
    size <- c(width=11,height=5)
  }
  p <- p + sb_viz_theme(); attr(p,"stormbreaker_suggested_size") <- size; p
}
