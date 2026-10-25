make_gene_rank_result <- function(id = 'p1', selector = 'screen', counts = c(a=2L,b=1L,c=0L), n = 2L, learner = 'logistic') {
  summary <- list(status = if(n == 2L) 'success' else 'incomplete', n_successful_folds = n,
    feature_frequency = data.frame(feature=names(counts), selected_count=unname(counts), selection_frequency=unname(counts)/n),
    per_fold=data.frame(fold=1:2,status=c('success',if(n==2L)'success' else 'failed'),n_selected=c(2L,1L)))
  structure(list(status='success', task=list(predictors=names(counts)),
    pipeline=list(id=id, selector=selector, selector_params=list(k=2L), learner=learner, learner_params=list(alpha=1)),
    selection_stability=summary, glmnet_active_set_stability=summary,
    provenance=list(fold_assignment=data.frame(patient_id=c('secret-patient-1','secret-patient-2'),fold=1:2))),
    class='stormbreaker_binary_cv_result')
}
make_gene_rank_comparison <- function(...) structure(list(results=list(...)),class='stormbreaker_cv_comparison')

test_that('shared selector copies count once and retain zero and unobserved genes', {
  a <- make_gene_rank_result()
  b <- a; b$pipeline$id <- 'p2'; b$pipeline$learner <- 'another_learner'
  c <- make_gene_rank_result('p3','other',c(a=1L,b=2L,c=0L))
  r <- sb_gene_selection_ranking(make_gene_rank_comparison(a,b,c), candidate_features=c('a','b','c','missing'))
  expect_identical(r$ranking$feature,c('a','b','c','missing'))
  expect_equal(r$ranking$selected_count,c(3,3,0,NA))
  expect_equal(r$ranking$eligible_count,c(4,4,4,0))
  expect_equal(r$ranking$rank,c(1,1,3,NA))
  expect_equal(r$ranking$selection_frequency,c(.75,.75,0,NA))
  expect_equal(r$groups$pipeline_copies,c(2L,1L))
  expect_equal(r$status,'success')
  expect_identical(any(grepl('secret-patient',unlist(r))),FALSE)
  expect_s3_class(sb_plot_gene_selection_ranking(r), 'ggplot')
})

test_that('failed folds do not enter the eligible denominator', {
  a <- make_gene_rank_result(counts=c(a=1L,b=0L,c=0L),n=1L)
  a$status <- 'failed'
  r <- sb_gene_selection_ranking(a)
  expect_equal(r$ranking$selected_count,c(1,0,0))
  expect_equal(r$ranking$eligible_count,c(1,1,1))
  expect_equal(r$ranking$selection_frequency,c(1,0,0))
  expect_identical(r$status,'incomplete')
})

test_that('none fixed panels and all-failed summaries do not invent zero evidence', {
  for(selector in c('none','fixed_features')) {
    r <- sb_gene_selection_ranking(make_gene_rank_result(selector=selector))
    expect_identical(r$status,'not_available')
    expect_identical(all(is.na(r$ranking$selected_count)),TRUE)
  }
  a <- make_gene_rank_result(counts=c(a=0L,b=0L,c=0L),n=0L)
  r <- sb_gene_selection_ranking(a)
  expect_identical(all(is.na(r$ranking$selected_count)),TRUE)
})

test_that('conflicting saved counts are visible and not arbitrarily selected', {
  a <- make_gene_rank_result()
  b <- make_gene_rank_result('p2',counts=c(a=1L,b=1L,c=0L))
  c <- make_gene_rank_result('p3','other')
  r <- sb_gene_selection_ranking(make_gene_rank_comparison(a,b,c))
  expect_identical(r$groups$status,c('conflicting_counts','success'))
  expect_equal(r$ranking$selected_count,c(2,1,0))
  expect_equal(r$ranking$eligible_count,c(2,2,2))
  expect_identical(r$status,'incomplete')
})

test_that('different split assignments and CV repetitions remain separate', {
  a <- make_gene_rank_result()
  b <- a; b$pipeline$id <- 'p2'; b$provenance$fold_assignment$fold <- 2:1
  r <- sb_gene_selection_ranking(make_gene_rank_comparison(a,b))
  expect_equal(r$ranking$selected_count,c(4,2,0))
  repeated <- structure(list(repetitions=list(make_gene_rank_comparison(a),make_gene_rank_comparison(a))),class='stormbreaker_repeated_cv_comparison')
  q <- sb_gene_selection_ranking(repeated)
  expect_equal(q$ranking$selected_count,c(4,2,0))
  expect_equal(nrow(q$groups),2L)
})

test_that('glmnet nonzero counts exclude ridge and stay separate from selectors', {
  a <- make_gene_rank_result(selector='none',learner='glmnet_logistic')
  b <- a; b$pipeline$id <- 'ridge'; b$pipeline$learner_params$alpha <- 0
  r <- sb_gene_selection_ranking(make_gene_rank_comparison(a,b),source='glmnet_active_set')
  expect_equal(r$ranking$selected_count,c(2,1,0))
  expect_equal(r$ranking$eligible_count,c(2,2,2))
  expect_identical(sb_gene_selection_ranking(a)$status,'not_available')
})

test_that('malformed counts and external evaluation inputs are rejected', {
  a <- make_gene_rank_result(counts=c(a=3L,b=1L,c=0L))
  expect_error(sb_gene_selection_ranking(a),class='stormbreaker_bad_gene_ranking')
  expect_error(sb_gene_selection_ranking(structure(list(),class='stormbreaker_external_evaluation')),class='stormbreaker_bad_gene_ranking')
})
