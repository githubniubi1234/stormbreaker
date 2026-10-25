sb_multiclass_tree_catalog <- function() {
  data.frame(method = c("ranger", "ranger_extratrees", "c50", "c50_rules", "ctree",
    "cforest", "bagging", "rrf", "tree", "evtree"),
    label = c("Ranger forest", "Extremely randomized forest", "C5.0 trees", "C5.0 rules",
      "Conditional inference tree", "Conditional inference forest", "Bagged trees",
      "Regularized random forest", "Classification tree", "Evolutionary tree"),
    supported = TRUE, engine = c("ranger", "ranger", "C50", "C50", "partykit",
      "partykit", "ipred", "RRF", "tree", "evtree"),
    feature_result = c("ranking", "ranking", rep("model_features", 8)),
    selection_rule = "No standalone gene-selection rule")
}

sb_fit_multiclass_tree <- function(fit, x, y) {
  catalog <- sb_multiclass_tree_catalog()
  method <- fit$method
  package <- catalog$engine[match(method, catalog$method)]
  if (!requireNamespace(package, quietly = TRUE)) sb_abort(
    paste("Install", package, "for this method."), "stormbreaker_missing_dependency")
  fit$package <- package
  fit$feature_codes <- paste0("feature", seq_along(fit$features))
  names(x) <- fit$feature_codes
  d <- x
  d$.outcome <- y
  formula <- stats::as.formula(".outcome ~ .", env = baseenv())
  params <- if (method %in% c("ranger", "ranger_extratrees"))
    sb_ranger_binary_params(learner = paste0(method, "_binary")) else
    sb_classical_binary_params(paste0(method, "_binary"))
  params$seed <- fit$seed
  record_warning <- function(w) fit$warnings <<- c(fit$warnings, conditionMessage(w))
  fit$engine <- withCallingHandlers(withr::with_seed(fit$seed, {
    if (method %in% c("ranger", "ranger_extratrees")) {
      ranger::ranger(x = x, y = y, probability = TRUE, num.trees = params$num_trees,
        mtry = max(1L, floor(sqrt(ncol(x)))), min.node.size = params$min_node_size,
        splitrule = params$splitrule, importance = "impurity", replace = params$replace,
        sample.fraction = params$sample_fraction, num.threads = 1L,
        seed = fit$seed, oob.error = FALSE)
    } else if (method %in% c("c50", "c50_rules")) {
      C50::C5.0(x, y, trials = params$trials, rules = params$rules,
        control = C50::C5.0Control(winnow = FALSE, earlyStopping = FALSE, seed = fit$seed))
    } else if (method %in% c("ctree", "cforest")) {
      control <- partykit::ctree_control(teststat = params$teststat, testtype = params$testtype,
        mincriterion = params$mincriterion, minsplit = params$minsplit,
        minbucket = params$minbucket, maxdepth = params$maxdepth, saveinfo = params$saveinfo)
      if (method == "ctree") partykit::ctree(formula, data = d, control = control) else
        partykit::cforest(formula, data = d, ntree = params$num_trees,
          mtry = max(1L, ceiling(sqrt(ncol(x)))),
          perturb = list(replace = params$replace, fraction = params$sample_fraction),
          control = control, cores = 1L, trace = FALSE)
    } else if (method == "bagging") {
      ipred::bagging(formula, data = d, nbagg = params$num_trees, coob = FALSE,
        keepX = params$keep_x, control = rpart::rpart.control(cp = params$cp,
          minsplit = params$minsplit, minbucket = params$minbucket,
          maxdepth = params$maxdepth, xval = 0L))
    } else if (method == "rrf") {
      RRF::RRF(x, y, ntree = params$num_trees, mtry = max(1L, floor(sqrt(ncol(x)))),
        nodesize = params$min_node_size, replace = params$replace,
        coefReg = params$regularization_coefficient, flagReg = params$regularization_mode,
        importance = FALSE, keep.forest = TRUE)
    } else if (method == "tree") {
      tree::tree(formula, data = d, split = params$split,
        control = tree::tree.control(nrow(d), mincut = params$mincut,
          minsize = params$minsize, mindev = params$mindev),
        model = FALSE, x = FALSE, y = FALSE, wts = FALSE)
    } else {
      evtree::evtree(formula, data = d, control = evtree::evtree.control(
        minbucket = params$minbucket, minsplit = params$minsplit, maxdepth = params$maxdepth,
        niterations = params$iterations, ntrees = params$population_size,
        alpha = params$complexity_weight, seed = fit$seed))
    }
  }), warning = record_warning)
  fit$parameters <- params
  structure(fit, class = "stormbreaker_multiclass_fit")
}

sb_predict_multiclass_tree <- function(fit, x) {
  if (!requireNamespace(fit$package, quietly = TRUE)) sb_abort(
    paste("Install", fit$package, "to predict with this saved fit."), "stormbreaker_missing_dependency")
  names(x) <- fit$feature_codes
  method <- fit$method
  probability <- if (method %in% c("ranger", "ranger_extratrees")) {
    stats::predict(fit$engine, data = x, num.threads = 1L)$predictions
  } else if (method == "cforest") {
    stats::predict(fit$engine, newdata = x, type = "prob", simplify = TRUE, OOB = FALSE)
  } else if (method == "tree") {
    stats::predict(fit$engine, newdata = x, type = "vector")
  } else stats::predict(fit$engine, newdata = x, type = "prob")
  if (is.list(probability)) probability <- do.call(rbind, probability)
  probability[, fit$codes, drop = FALSE]
}

sb_multiclass_tree_evidence <- function(fit) {
  tab <- data.frame(feature = fit$features, selected = NA)
  if (fit$method %in% c("ranger", "ranger_extratrees")) {
    tab$importance <- unname(fit$engine$variable.importance[fit$feature_codes])
    tab$rank <- rank(-tab$importance, ties.method = "min")
  }
  list(features = tab, by_class = NULL)
}

sb_multiclass_tree_plots <- function(fit) {
  if (fit$method == "tree") {
    structure <- sb_multiclass_tree_structure(fit)
    return(list(plots = list(tree_structure = sb_plot_rpart_tree(structure)),
      tables = list(tree_structure = structure$table)))
  }
  if (fit$method %in% c("ctree", "evtree")) {
    structure <- sb_party_native_structure(fit$engine, paste0(fit$method, "_multiclass"))
    for (i in order(nchar(fit$feature_codes), decreasing = TRUE)) {
      structure$table$incoming_rule <- gsub(fit$feature_codes[i], fit$features[i],
        structure$table$incoming_rule, fixed = TRUE)
      structure$table$split_feature <- gsub(fit$feature_codes[i], fit$features[i],
        structure$table$split_feature, fixed = TRUE)
    }
    p <- sb_plot_ctree_structure(structure) + ggplot2::labs(
      title = if (fit$method == "evtree") "Evolutionary tree" else "Conditional inference tree",
      caption = structure$interpretation)
    return(list(plots = list(tree_structure = p), tables = list(tree_structure = structure$table)))
  }
  if (fit$method %in% c("c50", "c50_rules", "rrf")) return(sb_multiclass_tree_native_plots(fit))
  if (!fit$method %in% c("ranger", "ranger_extratrees")) return(list(plots = list(), tables = list()))
  tab <- sb_multiclass_tree_evidence(fit)$features
  tab$feature <- factor(tab$feature, levels = tab$feature[order(tab$importance)])
  feature <- importance <- NULL
  plot <- ggplot2::ggplot(tab, ggplot2::aes(importance, feature)) +
    ggplot2::geom_col(fill = "#2878A0") + sb_viz_theme() +
    ggplot2::labs(title = "Forest: Importance", x = "Impurity importance", y = NULL)
  list(plots = list(importance = plot), tables = list(importance = tab))
}
