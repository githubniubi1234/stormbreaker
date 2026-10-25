sb_validate_show_unavailable <- function(value) {
  if(!is.logical(value) || length(value)!=1L || is.na(value)) {
    sb_abort("`show_unavailable` must be TRUE or FALSE.","stormbreaker_bad_display_option")
  }
  invisible(value)
}

sb_empty_display_plot <- function(title,subtitle=NULL,caption=NULL) {
  ggplot2::ggplot()+
    ggplot2::annotate("text",x=0,y=0,label="No usable results for this display",size=4)+
    ggplot2::labs(title=title,subtitle=subtitle,caption=caption)+
    ggplot2::theme_void(base_size=11)+
    ggplot2::theme(plot.caption=ggplot2::element_text(hjust=0))
}

sb_performance_display_rows <- function(tab) {
  !is.na(tab$status) & !is.na(tab$run_status) & tab$status=="success" &
    tab$run_status=="success" & is.finite(tab$estimate)
}

sb_curve_display_view <- function(view,kind,show_unavailable=FALSE) {
  sb_validate_show_unavailable(show_unavailable)
  if(show_unavailable) return(view)
  usable <- if(kind=="calibration") {
    unique(sb_curve_view_lines(view$data)$pipeline_id)
  } else {
    unique(view$data$pipeline_id[view$data$draw & view$data$strategy=="model"])
  }
  view$status <- view$status[view$status$pipeline_id %in% usable,,drop=FALSE]
  view$data <- view$data[view$data$pipeline_id %in% usable,,drop=FALSE]
  view$data$panel <- factor(view$data$panel,levels=view$status$panel)
  view
}
