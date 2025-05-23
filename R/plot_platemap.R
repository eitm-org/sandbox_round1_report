library(ggplot2)
library(tidyverse)

#' plot_platemap
#'
#' @param plate string that designates which plate you're creating a platemap for
#' @param platemap_df dataframe from data_clean with the columns for well, col, row, and nice_coi added
#'
#' @returns a ggplot object. a platemap for the plate you plug in
#' @export
#'
#' @examples
#' print(plot_platemap(plate, platemap_df))
plot_platemap <- function(plate, exp, platemap_df) {
  #filter platemap_df to include this plate only
  forplot_df <- platemap_df %>%
    filter(data.plate_id == plate)
  #designate y axis label
  plate_assay_type <- unique(platemap_df$assay_type)
  if (length(plate_assay_type) > 1) {
    warning(paste("more than one assay type designated for plate", plate))
  }
  #if it's a polar binding assay, the signal is in polarization (mP)
  ylabel <- case_when(plate_assay_type[[1]] == "Polar Binding" ~ "Polarization (mP)",
                      #if it's a screen or dose response curve assay, the signal is in rlu
                      TRUE ~ "Raw RLU")
  #create platemap ggplot
  plate_plot <- ggplot(data = forplot_df, aes(x = col, y = row, fill = data.result)) +
    geom_tile(color = "black", linewidth = 1) +
    scale_fill_viridis_c(option = "magma") +
    geom_text(color = "white",
              size = 2,
              aes(label = nice_coi)) +
    coord_fixed() +
    labs(
      title = paste("Experiment", exp, ", Plate", plate),
      x = "Column",
      y = "Row",
      fill = ylabel
    ) +
    theme(text = element_text(family = "serif"))
  return(plate_plot)
}
