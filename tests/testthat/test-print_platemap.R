library(here)
source(here("R", "print_platemap.R"))

test_that("returns nothing", {
  test_df <- data.frame(
    data.plate_id = c("0"),
    assay_type = c("0"),
    col = c("0"),
    row = c("0"),
    data.result = c(0),
    nice_coi = c("0"),
    coi = c("0"),
    experiment_id = c("0")
  )
  expect(is.null(print_platemap("0", test_df)), "why did this return something")
})
