#TODO: throw this into a script in the munge directory, functionalize, add tests

library(drc)
library(tidyverse)
library(conflicted)
library(here)

source(here("R", "functions.R"))

#resolve conflicts
conflicts_prefer(dplyr::filter(), stats::gaussian(), dplyr::select())

ctrl0 <<- "DMSO"
ctrl100 <<- "R1881"

screen <- readRDS(here("data", "screen_data.rds")) %>%
  unnest(cols = data, names_sep = ".") %>%
  select(
    study_name,
    study_rid,
    experiment_name,
    experiment_rid,
    experiment_id,
    experiment_type,
    assay_type,
    container,
    start_date,
    instrument,
    sop,
    data.plate_id,
    data.replicate_group,
    data.drug_name,
    data.drug_id,
    data.concentration,
    data.units,
    sample_id,
    data.well_id,
    data.result,
    data.plate_cv,
    data.plate_z_prime,
    data.norm_result
  ) %>%
  distinct()

drc <- readRDS(here("data", "drc_data.rds")) %>%
  unnest(cols = data, names_sep = ".") %>%
  select(
    study_name,
    study_rid,
    experiment_name,
    experiment_rid,
    experiment_id,
    experiment_type,
    assay_type,
    container,
    start_date,
    instrument,
    sop,
    data.plate_id,
    data.replicate_group,
    data.drug_name,
    data.drug_id,
    data.concentration,
    data.units,
    sample_id,
    data.well_id,
    data.result,
    data.plate_cv,
    data.plate_z_prime,
    data.norm_result
  ) %>%
  distinct()

binding <- readRDS(here("data", "binding_data.rds")) %>%
  unnest(cols = data, names_sep = ".") %>%
  select(
    study_name,
    study_rid,
    experiment_name,
    experiment_rid,
    experiment_id,
    experiment_type,
    assay_type,
    container,
    start_date,
    instrument,
    sop,
    data.plate_id,
    data.replicate_group,
    data.drug_name,
    data.drug_id,
    data.concentration,
    data.units,
    sample_id,
    data.well_id,
    data.result,
    data.plate_cv,
    data.plate_z_prime,
    data.norm_result
  ) %>%
  distinct() %>%
  #ignore the first column for the polar binding assays
  filter(!grepl("01$", data.well_id))

#add dmso to binding
meta <- binding %>%
  select(
    study_name,
    study_rid,
    experiment_name,
    experiment_rid,
    experiment_id,
    experiment_type,
    assay_type,
    container,
    start_date,
    instrument,
    sop,
    sample_id,
    data.plate_id,
    data.replicate_group,
    data.well_id,
    data.result,
    data.plate_cv,
    data.plate_z_prime,
    data.norm_result
  ) %>%
  distinct() %>%
  mutate(
    data.drug_name = "DMSO",
    data.concentration = .1,
    data.units = "%"
  )

ctrl_list <- c(ctrl0, ctrl100)
other_reagents <- c("FL AL Green", "AR  Screeining Buffer", "AR-LBD (GST)")

alldf <- bind_rows(screen, drc, binding, meta) %>%
  #first, designate controls and designate "compound of interest"
  #to do that, we have to nest by well (and keep all other columns outside of the nested part)
  group_by(
    data.well_id,
    study_name,
    study_rid,
    experiment_name,
    experiment_rid,
    experiment_id,
    experiment_type,
    assay_type,
    container,
    start_date,
    instrument,
    sop,
    data.plate_id,
    data.replicate_group,
    data.result,
    data.plate_cv,
    data.plate_z_prime,
    data.norm_result
  ) %>%
  nest() %>%
  mutate(
    coi = map(data, get_coi, ctrl_list, other_reagents),
    coi_dose = map2(data, coi, get_coi_concentration)
  ) %>%
  unnest(cols = c(coi, coi_dose)) %>%
  mutate(
    ctrl = case_when(
      #this will get fixed when you get the ctrl designations straight away
      coi == ctrl0 & assay_type != "Polar Binding" ~ "0% control",
      coi == ctrl100 & assay_type != "Polar Binding" ~ "100% control",
      coi == ctrl0 & assay_type == "Polar Binding" ~ "100% control",
      #controls are in the second column of hte polar binding assayss
      coi == ctrl100 &
        assay_type == "Polar Binding" &
        grepl("02$", data.well_id) ~ "0% control",
      #R1881 is in this list because it is also used as a control compound for the polar binding assays
      #(and the 20um dose is used for the 0%)
      grepl("VPC|VCP|Enzalutamide|Darolutamide|R1881", coi) ~ "control compound",
      coi == "what the " ~ "what the ",
      TRUE ~ "sample"
    )
  ) %>%
  unnest(cols = ctrl) %>%
  ungroup() %>%
  #plate 1498, 1499, 1500. rows f and g were mislabelled "r1881" when they should be blank wells
  mutate(
    coi = case_when(
      data.plate_id %in% c(1498, 1499, 1500) &
        grepl("F|G", data.well_id) ~ "blank well",
      TRUE ~ coi
    ),
    ctrl = case_when(
      data.plate_id %in% c(1498, 1499, 1500) &
        grepl("F|G", data.well_id) ~ "blank well",
      TRUE ~ ctrl
    )
  ) %>%
  #next, calculate mean and sd result for each replicate group
  group_by(experiment_id, data.plate_id, coi, coi_dose, ctrl) %>%
  mutate(
    plate_mean_res = mean(data.result, na.rm = TRUE),
    plate_sd_res = sd(data.result, na.rm = TRUE),
    plate_cv_result = (plate_sd_res / plate_mean_res) * 100
  ) %>%
  ungroup() %>%
  group_by(experiment_id, coi, coi_dose, ctrl) %>%
  mutate(
    exp_mean_res = mean(data.result, na.rm = TRUE),
    exp_sd_res = sd(data.result, na.rm = TRUE),
    exp_cv_result = (exp_sd_res / exp_mean_res) * 100
  ) %>%
  ungroup() %>%
  #use those to calculate plate %CV and zfactor
  group_by(data.plate_id) %>%
  mutate(
    plate_cv = mean(plate_cv_result, na.rm = TRUE),
    plate_mean_ctrl100 = mean(data.result[ctrl == "100% control"], na.rm = TRUE),
    plate_sd_ctrl100 = sd(data.result[ctrl == "100% control"], na.rm = TRUE),
    plate_mean_ctrl0 = mean(data.result[ctrl == "0% control"], na.rm = TRUE),
    plate_sd_ctrl0 = sd(data.result[ctrl == "0% control"], na.rm = TRUE),
    normed_results = (data.result - plate_mean_ctrl0) / (plate_mean_ctrl100 - plate_mean_ctrl0) *
      100
  ) %>%
  ungroup() %>%
  group_by(experiment_id) %>%
  mutate(
    exp_cv = mean(exp_cv_result, na.rm = TRUE),
    exp_mean_ctrl100 = mean(data.result[ctrl == "100% control"], na.rm = TRUE),
    exp_sd_ctrl100 = sd(data.result[ctrl == "100% control"], na.rm = TRUE),
    exp_mean_ctrl0 = mean(data.result[ctrl == "0% control"], na.rm = TRUE),
    exp_sd_ctrl0 = sd(data.result[ctrl == "0% control"], na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    plate_zfactor = 1 - (
      3 * (plate_sd_res + plate_sd_ctrl0) / (plate_mean_res - plate_mean_ctrl0)
    ),
    plate_zprime = 1 - (
      3 * (plate_sd_ctrl100 + plate_sd_ctrl0) / (plate_mean_ctrl100 - plate_mean_ctrl0)
    ),
    exp_zfactor = 1 - (3 * (exp_sd_res + exp_sd_ctrl0) / (exp_mean_res - exp_mean_ctrl0)),
    exp_zprime = 1 - (
      3 * (exp_sd_ctrl100 + exp_sd_ctrl0) / (exp_mean_ctrl100 - exp_mean_ctrl0)
    ),
    summary_zfactor = case_when(
      assay_type == "Screen" ~ exp_zfactor,
      assay_type == "Dose Response Curve" |
        assay_type == "Polar Binding" ~ plate_zfactor
    ),
    summary_zprime = case_when(
      assay_type == "Screen" ~ exp_zprime,
      assay_type == "Dose Response Curve" |
        assay_type == "Polar Binding" ~ plate_zprime
    ),
    summary_cv = case_when(
      assay_type == "Screen" ~ exp_cv,
      assay_type == "Dose Response Curve" |
        assay_type == "Polar Binding" ~ plate_cv
    ),
    assay_type = factor(
      assay_type,
      levels = c("Screen", "Dose Response Curve", "Polar Binding")
    )
  ) %>%
  #calculate cv cutoffs
  group_by(experiment_id, summary_cv, assay_type) %>%
  nest() %>%
  group_by(assay_type) %>%
  mutate(cv_cutoff = quantile(summary_cv, prob = .75) + 1.5 * IQR(summary_cv)) %>%
  unnest(cols = c(data)) %>%
  ungroup() %>%
  #designate hits
  #group by replicate group (over experiments)
  group_by(experiment_id, coi, coi_dose, ctrl) %>%
  mutate(exp_mean_pc = mean(normed_results, na.rm = TRUE)) %>%
  group_by(experiment_id) %>%
  mutate(
    exp_mean_ctrl100_pc = mean(normed_results[ctrl == "100% control"], na.rm = TRUE),
    exp_sd_ctrl100_pc = sd(normed_results[ctrl == "100% control"], na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    scr_hit = case_when(
      assay_type == "Screen" &
        exp_mean_pc < (exp_mean_ctrl100_pc - 9 * exp_sd_ctrl100_pc) &
        ctrl == "sample" ~ "9SD hit",
      assay_type == "Screen" &
        exp_mean_pc < (exp_mean_ctrl100_pc - 8 * exp_sd_ctrl100_pc) &
        ctrl == "sample" ~ "8SD hit",
      assay_type == "Screen" &
        exp_mean_pc < (exp_mean_ctrl100_pc - 7 * exp_sd_ctrl100_pc) &
        ctrl == "sample" ~ "7SD hit",
      assay_type == "Screen" &
        exp_mean_pc < (exp_mean_ctrl100_pc - 6 * exp_sd_ctrl100_pc) &
        ctrl == "sample" ~ "6SD hit",
      assay_type == "Screen" &
        exp_mean_pc < (exp_mean_ctrl100_pc - 5 * exp_sd_ctrl100_pc) &
        ctrl == "sample" ~ "5SD hit",
      assay_type == "Screen" &
        exp_mean_pc < (exp_mean_ctrl100_pc - 4 * exp_sd_ctrl100_pc) &
        ctrl == "sample" ~ "4SD hit",
      assay_type == "Screen" &
        exp_mean_pc < (exp_mean_ctrl100_pc - 3 * exp_sd_ctrl100_pc) &
        ctrl == "sample" ~ "3SD hit",
      assay_type == "Screen" ~ "not a hit",
      assay_type != "Screen" ~ "not a screen"
    ),
    scr_hit = factor(
      scr_hit,
      levels = c(
        "not a hit",
        "9SD hit",
        "8SD hit",
        "7SD hit",
        "6SD hit",
        "5SD hit",
        "4SD hit",
        "3SD hit",
        "not a screen"
      )
    ),
    #calculate fold change from 0% control
    fc_ctrl0 = data.result / plate_mean_ctrl0,
    #calculate separation from 100% control
    sep_ctrl100 = case_when(
      assay_type == "Screen" ~ 1 - (
        3 * (exp_sd_ctrl100 + exp_sd_res) / (exp_mean_ctrl100 - exp_mean_res)
      ),
      assay_type %in% c("Dose Response Curve", "Polar Binding") ~ 1 - (
        3 * (plate_sd_ctrl100 + plate_sd_res) / (plate_mean_ctrl100 - plate_mean_res)
      )
    )
  ) %>%
  #see if it models with drc
  #group by coi and plate
  group_by(
    study_name,
    study_rid,
    experiment_name,
    experiment_id,
    experiment_rid,
    experiment_type,
    container,
    start_date,
    instrument,
    sop,
    data.plate_id,
    assay_type,
    data.plate_cv,
    data.plate_z_prime,
    plate_cv,
    plate_mean_ctrl100,
    plate_sd_ctrl100,
    plate_mean_ctrl0,
    plate_sd_ctrl0,
    exp_cv,
    exp_mean_ctrl100,
    exp_sd_ctrl100,
    exp_mean_ctrl0,
    exp_sd_ctrl0,
    plate_zprime,
    exp_zprime,
    summary_zprime,
    cv_cutoff,
    coi,
    ctrl,
    summary_cv
  ) %>%
  nest() %>%
  mutate(does_it_drc = map_lgl(data, does_it_model),
         n_dose_per_cmpnd_plate = map_dbl(data, function(data) length(unique(data$coi_dose))))
# %>%
#   unnest(cols = c(data, does_it_drc)) %>%
#   ungroup()


#ok now make stats table
# alldf2 <- alldf %>%
#   group_by(data.plate_id, coi, assay_type, does_it_drc) %>%
#   nest() %>%
#   mutate(drm_col = map_if(data, does_it_drc, model_drcs, .else = NULL),
#          drm_coefs = map_if(drm_col, is.list, function(drm_col) drm_col$coefficients, .else = NULL),
#          hill = map_if(drm_coefs, function(coefs) !is.null(coefs), function(drm_coefs) drm_coefs[[1]], .else = 0),
#          min_val = map_if(drm_coefs, function(coefs) !is.null(coefs), function(drm_coefs) drm_coefs[[2]], .else = 0),
#          max_val = map_if(drm_coefs, function(coefs) !is.null(coefs), function(drm_coefs) drm_coefs[[3]], .else = 0)
#          ) %>%
#   unnest(cols = c(hill, min_val, max_val))


#TODO: read this guy https://www.tidyverse.org/blog/2019/02/purrr-0-3-0/
#TODO ADD TEST: there shouldn't be any drugs that model if the assay_type == "Screening"

this_models <- alldf %>%
  filter(does_it_drc == TRUE) %>%
  mutate(
    drm_col = map(data, model_drcs),
    drm_coefs = map(drm_col, function(drm_col)
      drm_col$coefficients),
    hill = map_dbl(drm_coefs, function(coefs)
      coefs[[1]]),
    min_val = map_dbl(drm_coefs, function(coefs)
      coefs[[2]]),
    max_val = map_dbl(drm_coefs, function(coefs)
      coefs[[3]]),
    abs_ec50_obj = map(drm_col, function(drm_col)
      ED(
        drm_col,
        respLev = 50,
        type = "absolute",
        interval = "delta"
      )),
    abs_ec50 = map_dbl(abs_ec50_obj, function(abs_ec50_obj)
      abs_ec50_obj[[1]]),
    abs_ec50_lower = map_dbl(abs_ec50_obj, function(abs_ec50_obj)
      abs_ec50_obj[[3]]),
    abs_ec50_upper = map_dbl(abs_ec50_obj, function(abs_ec50_obj)
      abs_ec50_obj[[4]]),
    rel_ec50_obj = map(drm_col, function(drm_col)
      ED(
        drm_col,
        respLev = 50,
        type = "relative",
        interval = "delta"
      )),
    rel_ec50 = map_dbl(rel_ec50_obj, function(rel_ec50_obj)
      rel_ec50_obj[[1]]),
    rel_ec50_lower = map_dbl(rel_ec50_obj, function(rel_ec50_obj)
      rel_ec50_obj[[3]]),
    rel_ec50_upper = map_dbl(rel_ec50_obj, function(rel_ec50_obj)
      rel_ec50_obj[[4]]),
    rse = map_dbl(drm_col, get_rse),
    max_dose = map_dbl(data, function(data)
      max(data$coi_dose)),
    extrapolated_rel_ec50 = case_when(rel_ec50 > max_dose ~ TRUE, TRUE ~ FALSE),
    extrapolated_abs_ec50 = case_when(abs_ec50 > max_dose ~ TRUE, TRUE ~ FALSE)
  ) %>%
  group_by(assay_type) %>%
  #find outlier rses
  mutate(rse_cutoff = quantile(rse, prob = .75) + 1.5 * IQR(rse)) %>%
  ungroup()

this_doesnt_model <- alldf %>%
  filter(does_it_drc == FALSE)

na_cols <- setdiff(names(this_models), names(this_doesnt_model))
this_doesnt_model[na_cols] <- NA

data <- alldf %>%
  unnest(cols = data) %>%
  ungroup()

stats <- bind_rows(this_models, this_doesnt_model)


saveRDS(data, here("data", "well_level_data.RDS"))
saveRDS(stats, here("data", "plate_level_stats.RDS"))
