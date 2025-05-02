library(here)
library(ggplot2)
library(tidyverse)
source(here("R", "plot_platemap.R"))

#' print_platemap
#'
#' @param plate character string designating plate (value for data.plate_id) to plot
#' @param platemap_df dataframe from data_clean with the columns for well, col, row, and nice_coi added. must contain rows for designated plate, or the plot won't work.
#'
#' @returns nothing-- prints with cat() [so this function will require the results = "asis" argument in the markdown chunk you run it in]
#' @export
#'
#' @examples
#' lapply(plates, print_platemap, platemap_df = platemap_df)
print_platemap <- function(plate, platemap_df) {
  #find experiment number for this plate
  exp <- platemap_df %>%
    filter(data.plate_id == plate) %>%
    distinct(experiment_id)
  if (nrow(exp) > 1) {
    warning(paste("more than one experiment designated for plate", exp))
  }
  exp <- exp[[1]]
  cat("#### Platemap for Plate", plate, ", Experiment ", exp)
  cat("\n")
  print(plot_platemap(plate, exp, platemap_df))
  cat("\\newpage")
}
