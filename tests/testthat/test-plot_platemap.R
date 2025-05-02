library(here)

test_that("function returns ggplot object", {
  source(here("R", "plot_platemap.R"))

  test_df <- data.frame(
    data.plate_id = c("0"),
    assay_type = c("0"),
    col = c("0"),
    row = c("0"),
    data.result = c(0),
    nice_coi = c("0"),
    coi = c("0")
  )
  expect(is_ggplot(plot_platemap("0", "test", platemap_df = test_df)),
         "function does not return a ggplot object")
})
