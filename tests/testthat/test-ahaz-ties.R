test_that("ahaz fit-only offsets preserve endpoints, order and caller RNG", {
  time <- c(1, 2, 2, 2, 3, 3, 4)
  event <- c(1L, 1L, 1L, 0L, 0L, 0L, 1L)
  original <- time
  withr::local_seed(62)
  before <- .Random.seed
  got <- sb_ahaz_fit_times(time, event)
  expect_identical(time, original)
  expect_identical(.Random.seed, before)
  expect_identical(got, sb_ahaz_fit_times(time, event))
  expect_identical(got$time[event == 0L], time[event == 0L])
  expect_equal(anyDuplicated(got$time[event == 1L]), 0L)
  expect_equal(all(got$time <= time), TRUE)
  expect_equal(all(got$time[time == 2] > 1), TRUE)
  expect_equal(got$audit$tied_event_groups, 1L)
  expect_equal(got$audit$tied_events, 2L)
  expect_equal(got$audit$changed_events, 1L)
  expect_equal(got$audit$max_abs_shift > 0, TRUE)
})

test_that("ahaz leaves supported times unchanged and breaks all-event ties", {
  time <- c(1, 1, 2, 2)
  event <- c(1L, 0L, 1L, 0L)
  expect_identical(sb_ahaz_fit_times(time, event)$time, time)
  same <- sb_ahaz_fit_times(rep(10, 20), rep(1L, 20))
  expect_equal(anyDuplicated(same$time), 0L)
  expect_equal(all(same$time > 0 & same$time <= 10), TRUE)
  expect_equal(max(10 - same$time) <= 1e-7, TRUE)
  alternative <- sb_ahaz_fit_times(rep(10, 20), rep(1L, 20), seed = 10102L)
  expect_equal(identical(same$time, alternative$time), FALSE)
  for (scale in c(1e-10, 1, 1e10)) {
    x <- c(1, 2, 2, 2.0001) * scale
    z <- sb_ahaz_fit_times(x, rep(1L, 4))$time
    expect_equal(anyDuplicated(z), 0L)
    expect_equal(z[1] < min(z[2:3]) && max(z[2:3]) < z[4], TRUE)
  }
})

test_that("unrepresentable ahaz offsets fail rather than reorder distinct times", {
  time <- c(1, 1 + .Machine$double.eps, 1 + .Machine$double.eps)
  error <- tryCatch(sb_ahaz_fit_times(time, rep(1L, 3)), error = identity)
  expect_s3_class(error, "stormbreaker_ahaz_time_precision")
})
