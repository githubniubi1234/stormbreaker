test_that("conditional trees retain native branches, membership and compact state", {
  skip_if_not_installed("partykit")
  for (id in c("ctree_binary", "ctree_survival")) {
    binary <- id == "ctree_binary"
    d <- if (binary) sb_example_binary(160L, seed = 108L) else sb_example_survival(160L, seed = 8612L)
    if (binary) d$gene_a <- d$gene_a + 1.5 * d$outcome
    features <- c("gene_c", "gene_a", "gene_b")
    task <- if (binary) sb_binary_task(d, "patient_id", "outcome", features) else
      sb_surv_task(d, "patient_id", "time", "event", features, "month")
    fit <- sb_fit_final(task, sb_pipeline(learner = id), capture_diagnostics = TRUE)
    expect_identical(fit$status, "success")
    before <- serialize(fit, NULL, version = 2L)
    e <- fit$diagnostic_artifacts$learner$engine
    view <- sb_ctree_structure(fit)
    tab <- view$table
    expect_equal(tab$node, partykit::nodeids(e))
    expect_equal(tab$node[tab$terminal], partykit::nodeids(e, terminal = TRUE))
    expect_gt(nrow(tab), 1L)
    membership <- partykit::fitted_node(partykit::node_party(e), partykit::data_party(e))
    expect_equal(tab$n[tab$terminal], as.integer(table(factor(membership, levels = tab$node[tab$terminal]))))
    expect_equal(tab$n[1], nrow(d))
    expect_equal(sum(tab$n[tab$terminal]), nrow(d))
    inspect <- function(node) {
      key <- as.character(partykit::id_node(node))
      split <- partykit::split_node(node)
      expect_identical(view$native_splits[[key]], split)
      expect_identical(view$native_node_information[[key]], partykit::info_node(node))
      kids <- partykit::kids_node(node)
      if (length(kids)) {
        labels <- partykit::character_split(split, partykit::data_party(e), digits = 7)
        children <- match(vapply(kids, partykit::id_node, 1L), tab$node)
        expect_equal(tab$parent[children], rep(as.integer(key), length(kids)))
        expect_identical(tab$incoming_rule[children], paste(labels$name, labels$levels))
        expect_equal(sum(tab$n[children]), tab$n[match(as.integer(key), tab$node)])
        invisible(lapply(kids, inspect))
      }
    }
    inspect(partykit::node_party(e))
    compact <- fit
    compact$diagnostic_artifacts <- NULL
    expect_identical(sb_ctree_structure(compact), view)
    expect_identical(sb_native_diagnostics(compact)$status$status, "available")
    p <- sb_plot_algorithm_diagnostic(view)
    expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
    expect_null(p$labels$subtitle)
    expect_null(p$labels$caption)
    expect_identical(serialize(fit, NULL, version = 2L), before)
  }
})

test_that("an unsplit conditional tree remains a valid one-node view", {
  skip_if_not_installed("partykit")
  d <- sb_example_binary(160L, seed = 108L)
  task <- sb_binary_task(d, "patient_id", "outcome", c("gene_c", "gene_a", "gene_b"))
  fit <- sb_fit_final(task, sb_pipeline(learner = "ctree_binary"))
  view <- sb_ctree_structure(fit)
  expect_equal(nrow(view$table), 1L)
  expect_identical(view$table$terminal, TRUE)
  expect_identical(view$table$incoming_rule, "Root")
  expect_equal(view$table$n, 160L)
  expect_null(view$native_splits[[1]])
  expect_s3_class(ggplot2::ggplot_build(sb_plot_algorithm_diagnostic(view)), "ggplot_built")
})
