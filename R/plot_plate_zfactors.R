library(ggplot2)
library(tidyverse)
library(ggbeeswarm)
#' plot_plate_zfactors
#'
#' @param df dataframe of all experiment results, produced by data_clean.R
#'
#' @returns ggplot object. boxplot of zfactors by plate
#' @export
#'
#' @examples
#' print(plot_plate_zfactors(alldf))
plot_plate_zfactors <- function(df) {
  #designate yaxis limits
    #[0,1], unless your input dataframe has values outside that range
  minz <- min(df$plate_zprime, na.rm = TRUE)
  maxz <- max(df$plate_zprime, na.rm = TRUE)
  ymin <- min(minz, 0)
  ymax <- max(maxz, 1)
  #data cleaning for plot
  qcplot <- df %>%
    select(experiment_id,
           experiment_type,
           assay_type,
           start_date,
           plate_zprime) %>%
    arrange(assay_type, experiment_id) %>%
    distinct() %>%
    mutate(
      qc_flag = case_when(plate_zprime < 0 ~ "Z'Factor < 0", TRUE ~ "Passed QC"),
      qc_flag = factor(qc_flag, levels = c("Passed QC", "Z'Factor < 0"))
    ) %>%
    #time to plot!
    ggplot(aes(x = assay_type, y = plate_zprime)) +
    geom_boxplot(outliers = FALSE) +
    geom_quasirandom(size = 3,
                     alpha = .6,
                     aes(color = experiment_id, shape = qc_flag)) +
    scale_color_viridis_d(end = .9) +
    theme_bw() +
    #add labels
    labs(
      x = "Assay Type",
      y = "Plate Z'Factor",
      color = "Experiment ID",
      shape = "QC Flag",
      title = "Z'Factors by Plate"
    ) +
    #remove gridlines (they're distracting)
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      text = element_text(family = "serif")
    ) +
    ylim(c(ymin, ymax))
  return(qcplot)
}
