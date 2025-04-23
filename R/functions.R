get_coi <- function(data, ctrl_list, other_reagents = c()) {
  drug_list <- unique(data$data.drug_name)
  drug_list <- drug_list[!(drug_list %in% other_reagents)]
  ndrug <- length(drug_list)
  noctrls_list <- drug_list[!(drug_list %in% ctrl_list)]
  if (length(noctrls_list) > 1) {
    return("what the ")
  } else if (length(noctrls_list) == 1) {
    return(unlist(noctrls_list))
  } else if (length(noctrls_list) == 0) {
    if (ctrl100 %in% drug_list) {
      return(ctrl100)
    } else if (ctrl0 %in% drug_list) {
      return(ctrl0)
    } else {
      return("blank well")
    }
  }
}

get_coi_concentration <- function(data, coi) {
  coi_conc <- data[data$data.drug_name == coi, "data.concentration"][[1]]
  return(coi_conc)
}

get_raw_drc_plot <- function(mcule, drc_exps) {
  cat("## ", mcule)
  cat("\n")
  #raw plot
  #get the experiment that contains the curve for this coi
  exps <- drc_exps %>%
    filter(coi == mcule) %>%
    distinct(experiment_id) %>%
    as.vector() %>%
    unname() %>%
    unlist() %>%
    sort()
  exps_df <- drc_exps %>%
    filter(experiment_id %in% exps) %>%
    filter(coi == mcule | ctrl != "sample")
  #make a list of shapes to graph with
  #(to make sure you have enough shapes for all the compounds in these experiments)
  ncois <- length(unique(exps_df$coi))
  if (ncois == 4) {
    shapes <- c(15, 17, 18, 19)
  } else if (ncois == 3) {
    shapes <- c(15, 17, 19)
  } else if (ncois == 2) {
    shapes <- c(15, 19)
  } else if (ncois == 5) {
    shapes <- c(15, 17, 6, 1, 19)
  } else if (ncois == 6) {
    shapes <- c(15, 17, 6, 1, 0, 19)
  } else {
    shapes <- c(15, 17, 6, 1, 0, 2, 19)
  }
  raw_plot <- ggplot(data = exps_df, aes(x = coi_dose, y = data.result, shape = coi_lab, color = experiment_id)) +
    geom_point(size = 3, alpha = .5, position = "jitter") +
    geom_smooth(se = FALSE) +
    theme_bw() +
    scale_x_log10() +
    scale_color_viridis_d(end = .8) +
    labs(
      title = paste(mcule, "Raw Dose Response"),
      x = "Dose",
      y = "Raw RLU",
      color = "Experiment\nID",
      shape = "Compound"
    ) +
    scale_shape_manual(values = shapes) +
    ylim(c(0, max(exps_df$data.result, na.rm = TRUE)))
  print(raw_plot)
  cat("\n")
}

