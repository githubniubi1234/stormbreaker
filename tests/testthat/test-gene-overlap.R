test_that("gene intersections preserve exact exclusive members", {
  x <- sb_gene_overlap(list(A = c("a", "b", "b"), B = c("b", "c"), C = character()))
  expect_identical(x$common, character())
  expect_equal(sum(x$intersections$n), 3)
  expect_identical(x$intersections$genes[[match("A & B", x$intersections$intersection)]], "b")
  expect_identical(colnames(x$membership), c("A", "B", "C"))
  expect_s3_class(ggplot2::ggplotGrob(sb_plot_gene_overlap(x, "bar")), "gtable")
  empty <- sb_gene_overlap(list(A = character(), B = character()))
  expect_equal(nrow(empty$intersections), 0)
  expect_equal(dim(empty$membership), c(0L, 2L))
})

test_that("Venn diagrams render two to four gene sets", {
  skip_if_not_installed("ggvenn")
  sets <- list(A = c("a", "b"), B = c("b", "c"), C = c("b", "d"), D = "b")
  for (k in 2:4) {
    x <- sb_gene_overlap(sets[seq_len(k)])
    p <- sb_plot_gene_overlap(x)
    expect_identical(attr(p, "gene_overlap")$common, "b")
    expect_s3_class(ggplot2::ggplotGrob(p), "gtable")
  }
})
