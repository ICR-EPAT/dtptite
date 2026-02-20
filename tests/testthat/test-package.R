test_that("dtptite package loads", {
  expect_true(is.character(utils::packageDescription("dtptite")$Package))
})
