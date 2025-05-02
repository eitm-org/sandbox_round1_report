library(here)

test_that("function returns ggplot object", {
  source(here("R", "plot_zfactor_outliers.R"))

  test_df <- data.frame(
    plate_zprime = c("0"),
    assay_type = c("0"),
    experiment_id = c("0"),
    data.plate_id = c("0"),
    ctrl = c("0"),
    data.result = c(0),
    coi = c("0")
  )
  expect(is_ggplot(plot_zfactor_outliers(test_df)),
         "function does not return a ggplot object")
})
