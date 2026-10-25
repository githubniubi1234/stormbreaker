sb_multiclass_tree_structure <- function(fit) {
  engine <- fit$engine
  frame <- engine$frame
  node <- as.numeric(rownames(frame))
  leaf <- frame$var == "<leaf>"
  rules <- if (inherits(engine, "tree")) labels(engine, pretty = FALSE) else
    labels(engine, digits = 5, minlength = 0)
  for (i in order(nchar(fit$feature_codes), decreasing = TRUE))
    rules <- gsub(fit$feature_codes[i], fit$features[i], rules, fixed = TRUE)
  codes <- attr(engine, "ylevels")
  predicted <- fit$classes[match(codes[frame$yval], fit$codes)]
  x <- rep(NA_real_, length(node))
  x[leaf] <- seq_len(sum(leaf))
  place <- function(n) {
    i <- match(n, node)
    if (is.na(x[i])) x[i] <<- mean(c(place(2 * n), place(2 * n + 1)))
    x[i]
  }
  place(1)
  tab <- data.frame(node = node, parent = ifelse(node == 1, NA_real_, floor(node / 2)),
    x = x, y = -floor(log2(node)), terminal = leaf, n = frame$n,
    incoming_rule = rules, prediction = predicted,
    display_label = paste0("n = ", frame$n, "\n", rules, "\nClass: ", predicted))
  list(table = tab, prediction_scale = "Training-node class",
    learner = if (inherits(engine, "tree")) "multiclass_tree" else "multiclass_rpart",
    suggested_width = max(11, 1.6 * sum(leaf)),
    suggested_height = max(7, 1.2 * (max(-tab$y) + 1)))
}
