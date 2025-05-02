library(testthat)
library(here)

test_that("function returns ggplot object", {
  source(here("R", "plot_plate_zfactors.R"))

  test_df <- data.frame(
    experiment_id = c("0"),
    experiment_type = c("0"),
    assay_type = c("0"),
    start_date = c(0),
    plate_zprime = c(0)
  )
  expect(is_ggplot(plot_plate_zfactors(test_df)),
         "function does not return a ggplot object")
})
