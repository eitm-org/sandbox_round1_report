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

#TODO: write this function
# overlay_drcs <- function() {
#
# }

get_normed_plot <- function(df, min_y, max_y) {
  #TODO: should i factor the compounds somewhere?
  color_df <- data.frame(cmpnds = unique(df$coi_lab)) %>%
    arrange(cmpnds)
  color_df["color"] <- viridis(nrow(color_df), end = .8)
  outplot <- ggplot(data = df, aes(x = coi_dose, y = normed_results, color = coi_lab)) +
    geom_point(size = 2, alpha = .5, aes(shape = data.plate_id)) +
    theme_bw() +
    scale_x_log10() +
    scale_color_manual(breaks = color_df$cmpnds, values = color_df$color) +
    ylim(c(min_y, max_y)) +
    theme(legend.position = "bottom",
          legend.box = "vertical",
          text = element_text(family = "serif")) +
    labs(x = "Dose",
         y = "Percent Control (%)",
         shape = "Plate",
         color = "Compound") +
    guides(col = guide_legend(ncol = 2))
  models <- stats %>%
    #filter for relevant compounds
    filter(coi %in% df$coi & data.plate_id %in% df$data.plate_id & ctrl %in% df$ctrl) %>%
    #filter out 0 and 100% controls-- we don't need to see the drcs for those!
    filter(!(ctrl %in% c("100% control", "0% control")))
  they_model <- models %>%
    filter(does_it_drc == TRUE)
  # map2(they_model$coi, they_model$experiment_id, overlay_drcs)
  #TODO: one day do it like this
  for (cmp in they_model$coi) {
    cmp_lab <- str_remove(cmp, "MCULE-")
    color1 <- color_df[color_df$cmpnds == cmp_lab, "color"]
    #TODO: make sure the experiment u loop on is also plotted
      #or u plot with experiment instead of date
    for (exp in they_model$experiment_id) {
      #TODO: add test to make sure you only get one model here
      #we're grabbing out the model for this compound (of interest..... coi) and experiment
      #and models were generated for each unique combination of coi and experiment
      #so there SHOULD only EVER be one model popping up here
      #but add a test to make sure!
      #TODO: filter this instead
      model_df <- they_model[they_model$coi == cmp & they_model$experiment_id == exp,]
      #TODO: loop better to avoid having to do this:
      if (nrow(model_df) == 0) {
        next
      }
      model <- model_df$drm_col[[1]]
      #get max and min doses this drug was tested at in this experiment
      #TODO: do this more elegantly
      min_dose <- min(model_df$data[[1]]$coi_dose, na.rm = TRUE)
      #TODO: add a test that checks that none of these doses are NA
      max_dose <- max(model_df$data[[1]]$coi_dose, na.rm = TRUE)
      #get points to draw the curve
      curve.data.drc <- PR(model, xVec = pracma::logseq(min_dose, max_dose, n = 100))
      curve.data <- data.frame(y = as.vector(unlist(curve.data.drc)), x = as.numeric(names(curve.data.drc)))
      #add the drc curve to the plot
      outplot <- outplot +
        geom_line(data = curve.data, aes(x = x, y = y), color = color1)
      #           , linetype = linetype) +
      # annotate("label", x = label_x_plotter, hjust = 1, y = label_y_plotter, label = this_date, color = this_color, size = 4,
      #          family = "serif")
      #add vertical line for absolute ec50
      abs_ec50 <- model_df$abs_ec50
      if (!is.nan(abs_ec50 & !is.na(abs_ec50))) {
        outplot <- outplot +
          geom_vline(xintercept = abs_ec50, alpha = .5, linetype = 5, linewidth = .5)
      }
    }
  }
  return(outplot)
}

print_drcs <- function(mcule, drc_exps, raw_plot_ylab = "input raw plot ylab") {
  cat("## ", mcule)
  cat("\n")
  cat(mcule, "Raw Plot")
  #raw plot
  #get the experiment that contains the curve for this coi
  plates <- drc_exps %>%
    filter(coi == mcule) %>%
    distinct(data.plate_id) %>%
    as.vector() %>%
    unname() %>%
    unlist() %>%
    sort()
  exps_df <- drc_exps %>%
    filter(data.plate_id %in% plates) %>%
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
  raw_plot <- ggplot(data = exps_df,
                     aes(
                       x = coi_dose,
                       y = data.result,
                       shape = coi_lab,
                       color = data.plate_id
                     )) +
    geom_point(size = 3,
               alpha = .5,
               position = "jitter") +
    geom_smooth(se = FALSE, aes(linetype = ctrl), linewidth = 1) +
    theme_bw() +
    scale_x_log10() +
    scale_color_viridis_d(option = "plasma", end = .9) +
    labs(
      title = paste(mcule, "Raw Dose Response"),
      x = "Dose",
      y = raw_plot_ylab,
      color = "Plate\nID",
      shape = "Compound"
    ) +
    scale_shape_manual(values = shapes) +
    ylim(c(0, max(exps_df$data.result, na.rm = TRUE))) +
    theme(text = element_text(family = "serif"))
  print(raw_plot)
  cat("\n")
  #make table with qc stats
  qc_tab <- exps_df %>%
    select(
      experiment_id,
      start_date,
      data.plate_id,
      data.plate_cv,
      data.plate_z_prime,
      plate_mean_ctrl100,
      plate_sd_ctrl100,
      plate_mean_ctrl0,
      plate_sd_ctrl0,
      cv_cutoff
    ) %>%
    distinct()
  names(qc_tab) <- c(
    "Experiment ID",
    "Start Date",
    "Plate ID",
    "Plate %CV",
    "Plate Z'Factor",
    paste("Mean", ctrl100, "RLU"),
    paste("SD", ctrl100, "RLU"),
    paste("Mean", ctrl0, "RLU"),
    paste("SD", ctrl0, "RLU"),
    "Plate %CV Cutoff"
  )
  print(kable(qc_tab, digits = 3))
  cat("\n")
  #get normalized plots
  ctrl_cmpnds <- exps_df %>%
    filter(ctrl != "sample" & n_dose_per_cmpnd_plate > 1)
  samples <- exps_df %>%
    filter(ctrl == "sample" | n_dose_per_cmpnd_plate < 2)

  #get y axis bounds
  min_y <- min(exps_df$normed_results, na.rm = TRUE)
  max_y <- max(exps_df$normed_results, na.rm = TRUE)

  cat("\n")
  cat("\\newpage")
  cat("\n")

  if (nrow(ctrl_cmpnds) > 0) {
    ctrls_plot <- get_normed_plot(ctrl_cmpnds, min_y, max_y)
    coi_plot <- get_normed_plot(samples, min_y, max_y)
    printer <- ctrls_plot + coi_plot
  } else {
    coi_plot <- get_normed_plot(samples, min_y, max_y)
    printer <- coi_plot
  }

  cat(mcule, "Normalized Plot")
  cat("\n")
  print(printer)
  cat("\n")

  #print stats tab
  mcule_stats <- stats %>%
    filter(data.plate_id %in% plates) %>%
    filter(coi == mcule) %>%
    mutate(abs_ec50_ci = paste0("(", round(abs_ec50_lower, 2), ", ", round(abs_ec50_upper, 2), ")"),
           rel_ec50_ci = paste0("(", round(rel_ec50_lower, 2), ", ", round(rel_ec50_upper, 2), ")")
    ) %>%
    select(experiment_id, data.plate_id, abs_ec50, abs_ec50_ci,
           rel_ec50, rel_ec50_ci, rse, rse_cutoff, extrapolated_rel_ec50, extrapolated_abs_ec50)
  names(mcule_stats) <- c("Experiment", "Plate", "Absolute EC50", "Absolute EC50 CI",
                          "Relative EC50", "Relative EC50 CI", "RSE", "RSE Cutoff", "Extrapolated Relative C50", "Extrapolated Absolute EC50")
  cat("\n")

  print(kable(mcule_stats, digits = 3))

  cat("\n")
  cat("\\newpage")
  cat("\n")

}

does_it_model <- function(df) {
  tryCatch({
    suppressMessages(drm(
      normed_results ~ coi_dose,
      data = df,
      fct = LL.4(names = c(
        "hill", "min_value", "max_value", "ec_50"
      ))
    ))
    return(TRUE)
  }, error = function(cond) {
    return(FALSE)
  })
}

model_drcs <- function(df) {
  drm(
    normed_results ~ coi_dose,
    data = df,
    fct = LL.4(names = c(
      "hill", "min_value", "max_value", "ec_50"
    ))
  )
}

get_rse <- function(drm_col1) {
  mod_summary <- summary(drm_col1)
  rse <- mod_summary$rseMat[[1]]
  return(rse)
}


