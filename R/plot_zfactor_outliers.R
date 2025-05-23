library(tidyverse)
library(ggbeeswarm)
#' plot_zfactor_outliers
#'
#' @param dfdataframe of all experiment results, produced by data_clean.R
#'
#' @returns ggplot object. boxplot of raw signal by control for plates with outlier z'factors
#' @export
#'
#' @examples
#' print(plot_zfactor_outliers(alldf))
#'
plot_zfactor_outliers <- function(df) {
  zplot <- df %>%
    #select outlier z'factors
    filter(plate_zprime < 0 | plate_zprime > 20) %>%
    #add columns that will help make your plot nicer
    mutate(
      nice_label = paste(assay_type, "Exp", experiment_id, "Plate", data.plate_id),
      nice_ctrl = case_when(
        ctrl == "100% control" ~ "100%\ncontrol",
        ctrl == "0% control" ~ "0%\ncontrol",
        ctrl == "control compound" ~ "control\ncompound",
        TRUE ~ ctrl
      )
    ) %>%
    #start plotting!
    ggplot(aes(x = nice_ctrl, y = data.result)) +
    facet_grid(cols = vars(nice_label), scales = "fixed", margins = FALSE) +
    # facet_wrap( ~ nice_label, scales = "free") +
    geom_boxplot(outliers = FALSE) +
    geom_quasirandom(size = 3, alpha = .5, aes(color = coi)) +
    scale_color_viridis_d(guide = "none", end = .8) +
    theme_bw() +
    labs(
      title = "Raw Signal for Z'Factor Outliers",
      y = "Raw RLU",
      x = "Control",
      caption = "Points colored by compound.",
      color = "Z Factor"
    ) +
    theme(text = element_text(family = "serif"))
  return(zplot)
}
