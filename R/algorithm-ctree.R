#' Inspect an original conditional inference tree
#' @param fit A successful binary or survival ctree final fit.
#' @return An algorithm diagnostic containing every original node, branch rule,
#'   training-node count, exact native split and available native node information.
#'   Counts describe training membership, not validation performance. No new
#'   hypothesis tests, probabilities or survival estimates are computed. The
#'   compact saved engine suffices; no refitting or pruning is performed.
#' @export
sb_ctree_structure <- function(fit) {
  fit <- sb_validate_final_fit_for_prediction(fit)
  if (!fit$pipeline$learner %in% c("ctree_binary", "ctree_survival")) sb_abort(
    "This diagnostic requires a single conditional inference tree.",
    "stormbreaker_unsupported_diagnostic")
  sb_party_tree_structure(fit)
}

#' Inspect the original evolutionary classification tree
#' @param fit Successful evolutionary tree final fit.
#' @return All original nodes, split rules and retained training memberships.
#'   The compact predictor schema maps native split indices to feature names;
#'   fitted membership columns are never used as predictor labels. No new
#'   evolutionary search, pruning, hypothesis test or surrogate tree is fitted.
#' @export
sb_evtree_structure <- function(fit) {
  fit <- sb_validate_final_fit_for_prediction(fit)
  if (fit$pipeline$learner != "evtree_binary") sb_abort(
    "This diagnostic requires an evolutionary classification tree.", "stormbreaker_unsupported_diagnostic")
  sb_party_tree_structure(fit)
}

sb_party_tree_structure <- function(fit) {
  e <- unserialize(fit$learner_state$engine_raw)
  sb_party_native_structure(e, fit$pipeline$learner, fit$pipeline$id)
}

sb_party_native_structure <- function(e, learner, pipeline_id = NULL) {
  evolutionary <- learner %in% c("evtree_binary", "evtree_multiclass")
  if (!inherits(e, "constparty")) sb_abort("The original conditional tree is missing.",
    "stormbreaker_bad_learner_state")
  frame <- partykit::data_party(e)
  label_frame <- if (evolutionary) e$data else frame
  membership <- frame[["(fitted)"]]
  weights <- frame[["(weights)"]]
  if (is.null(membership) || is.null(weights) || length(membership) != length(weights)) sb_abort(
    "Original tree membership or weights are missing.", "stormbreaker_bad_learner_state")
  rows <- splits <- information <- list()
  leaf_position <- 0L
  visit <- function(node, parent = NA_integer_, depth = 0L, rule = "Root") {
    id <- partykit::id_node(node)
    key <- as.character(id)
    kids <- partykit::kids_node(node)
    split <- partykit::split_node(node)
    terminal <- !length(kids)
    rules <- if (terminal) NULL else partykit::character_split(split, data = label_frame, digits = 7)
    if (!terminal && length(rules$levels) != length(kids)) sb_abort(
      "Native split labels do not match child branches.", "stormbreaker_bad_learner_state")
    rows[[key]] <<- data.frame(node = id, parent = parent, depth = depth,
      terminal = terminal, incoming_rule = rule,
      split_feature = if (terminal) NA_character_ else rules$name,
      n = NA_integer_, weight = NA_real_, x = NA_real_, y = -depth)
    splits[key] <<- list(split)
    information[key] <<- list(partykit::info_node(node))
    if (terminal) {
      leaf_position <<- leaf_position + 1L
      descendants <- id
      position <- leaf_position
    } else {
      descendants <- unlist(lapply(seq_along(kids), function(j) {
        visit(kids[[j]], id, depth + 1L, paste(rules$name, rules$levels[j]))
      }), use.names = FALSE)
      positions <- vapply(kids, function(kid) rows[[as.character(partykit::id_node(kid))]]$x, 0)
      position <- mean(range(positions))
    }
    included <- membership %in% descendants
    rows[[key]]$n <<- sum(included)
    rows[[key]]$weight <<- sum(weights[included])
    rows[[key]]$x <<- position
    descendants
  }
  leaves <- visit(partykit::node_party(e))
  if (!all(membership %in% leaves)) sb_abort("Native membership contains unknown terminal nodes.",
    "stormbreaker_bad_learner_state")
  tab <- do.call(rbind, unname(rows))
  rownames(tab) <- NULL
  structure(list(type = "ctree_structure", table = tab, native_splits = splits,
    native_node_information = information, learner = learner,
    pipeline_id = pipeline_id,
    source = if (evolutionary) "original evolutionary tree and retained training memberships" else
      "all original conditional inference tree nodes and saved training memberships",
    interpretation = paste("Counts are original training-node memberships, not validation performance.",
      "Exact native split objects are retained separately from rounded display rules.",
      "Only originally retained node statistics are available; missing tests are not recomputed.")),
    class = "stormbreaker_algorithm_diagnostic")
}

sb_plot_ctree_structure <- function(x) {
  d <- x$table
  parent <- match(d$parent, d$node)
  edges <- data.frame(x = d$x[parent], y = d$y[parent] - .22,
    xend = d$x, yend = d$y + .22)
  edges <- edges[!is.na(parent), , drop = FALSE]
  d$label <- vapply(seq_len(nrow(d)), function(i) paste0("Node ", d$node[i],
    " | n = ", d$n[i], "\n",
    paste(strwrap(d$incoming_rule[i], width = 28), collapse = "\n")), "")
  y <- xend <- yend <- label <- terminal <- NULL
  p <- ggplot2::ggplot(d, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_segment(data = edges, ggplot2::aes(xend = xend, yend = yend), colour = "grey60") +
    ggplot2::geom_label(ggplot2::aes(label = label, fill = terminal), size = 3.1,
      linewidth = .2, show.legend = FALSE) +
    ggplot2::scale_fill_manual(values = c("FALSE" = "#EDF3F7", "TRUE" = "#D7EBE0")) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(add = .8)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(add = .5)) +
    ggplot2::labs(title = if (x$learner == "evtree_binary") "Evolutionary classification tree" else
      if (x$learner == "tree_binary") "Classification tree" else
      "Conditional inference tree", x = NULL, y = NULL) +
    sb_viz_theme() + ggplot2::theme(axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(), panel.grid = ggplot2::element_blank())
  attr(p, "stormbreaker_suggested_size") <- c(width = max(9, 2.3 * sum(d$terminal)),
    height = max(5, 1.6 * (max(d$depth) + 1)))
  p
}
