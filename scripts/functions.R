library(tidyverse)
library(broom)
library(dplyr)
library(janitor)
library(survey)
library(microsynth)


# Notin function
`%notin%` <- Negate(`%in%`)


# Function to clean finger tips data
clean_fingertips_data <- function(df) {
  df %>%
    janitor::clean_names() %>%
    filter(area_type == "GPs") %>%
    mutate(year = as.numeric(substr(timeperiod_sortable, 1, 4))) %>%
    filter(year >= 2019) %>%
    select(area_code, value, year) %>%
    rename(gpcode = area_code)
}

adjust_indicator_name <- function(df) {
  df %>% 
    rename("{df$indicator_name[1]}" := value) %>% 
    select(-c(indicator_name)) %>% 
    clean_names()
}



##################################################### Synthetic control functions ##################################################### 



# Function to extract data from synth and return it as a datafrom -----------
extract_syn_results <- function(...) {
  
  # capture all model objects passed into the function
  models <- list(...)
  
  # helper to clean and format one model's Results
  format_results <- function(model) {
    model$Results %>%
      as.data.frame() %>%
      rownames_to_column("Outcome") %>%
      rename_with(~ str_remove(., "^X\\d+\\.")) %>%
      mutate(
        Pct.Chng     = if_else(is.na(Pct.Chng), "NA", glue("{round(Pct.Chng * 100, 1)}%")),
        # Linear.Lower = glue("{round(Linear.Lower * 100, 1)}%"),
        # Linear.Upper = glue("{round(Linear.Upper * 100, 1)}%"),
        Perm.Lower   = glue("{round(Perm.Lower * 100, 1)}%"),
        Perm.Upper   = glue("{round(Perm.Upper * 100, 1)}%"), 
        Outcome = case_when(Outcome == "i_mmr1_24m"   ~ "MMR1 (24 months)",
                            Outcome ==  "i_mmr1_5y"    ~ "MMR1 (5 years)",
                            Outcome ==  "i_mmr2_5y"    ~ "MMR2 (5 years)",
                            Outcome ==  "i_rota"       ~ "Rotavirus (12 months)",
                            Outcome ==  "i_pcv_boost"  ~ "PCV booster (24 months)",
                            Outcome ==  "i_six_in_one" ~ "6-in-1 (12 months)")
      ) %>%
      mutate(across(where(is.numeric), round, 4)) %>% 
      select(-c(Linear.pVal, Linear.Lower, Linear.Upper))
  }
  
  # apply to each model provided
  results <- models %>%
    purrr::map(format_results) %>%
    dplyr::bind_rows()
  
  return(results)
}  



# Svy ratios -------
extract_svy_ratios <- function(model, outcome_denom_map, model_name = NULL, synth.df = NULL) {
  if (is.null(model_name)) {
    model_name <- deparse(substitute(model))
  }
  
  if (is.null(synth.df)) {
    synth.df <- gp 
  }
  
  valid_times <- c(model$info$start.pre:model$info$end.post)
  
  # Extract weights
  weights_df <- as.data.frame(model$w$Weight) %>%
    tibble::rownames_to_column(var = "gpcode")  # only weights needed
  
  # Merge weights into main dataset
  gpweights <- synth.df %>%
    dplyr::filter(time %in% valid_times) %>% 
    dplyr::left_join(weights_df %>%
                       select(gpcode, Main), by = "gpcode")
  
  # Survey design
  svyfull <- svydesign(ids = ~0, weights = ~Main, data = gpweights)
  
  # Loop over outcome/denominator pairs
  results <- purrr::map_dfr(names(outcome_denom_map), function(outcome_var) {
    denom_var <- outcome_denom_map[[outcome_var]]
    
    svy_ratio <- svyby(
      formula = as.formula(paste0("~", outcome_var)),
      denominator = as.formula(paste0("~", denom_var)),
      by = ~interv + time,
      design = svyfull,
      FUN = svyratio,
      na.rm = TRUE
    )
    
    svy_ci <- confint(svyby( formula = as.formula(paste0("~", outcome_var)),
                             denominator = as.formula(paste0("~", denom_var)),
                             by=~interv+time ,design=svyfull,FUN=svyratio))
    
    # Combine and label
    result_df <- cbind(svy_ratio, svy_ci)
    colnames(result_df) <- c("interv", "time", "uptake", "se", "lcl", "ucl")
    
    result_df$Outcome <- outcome_var
    result_df$Model <- model_name
    
    result_df
  })
  
  return(results)
}



## svy ratios plot -------
plot_svy_ratios <- function(data, endpre, ncol = 2) {
  
  quarter_labels <- c(
    "Q1 (Apr–Jun 2019)",
    "Q2 (Jul–Sep 2019)",
    "Q3 (Oct–Dec 2019)",
    "Q4 (Jan–Mar 2020)",
    "Q5 (Apr–Jun 2020)",
    "Q6 (Jul–Sep 2020)",
    "Q7 (Oct–Dec 2020)",
    "Q8 (Jan–Mar 2021)",
    "Q9 (Apr–Jun 2021)",
    "Q10 (Jul–Sep 2021)",
    "Q11 (Oct–Dec 2021)",
    "Q12 (Jan–Mar 2022)",
    "Q13 (Apr–Jun 2022)",
    "Q14 (Jul–Sep 2022)",
    "Q15 (Oct–Dec 2022)",
    "Q16 (Jan–Mar 2023)",
    "Q17 (Apr–Jun 2023)",
    "Q18 (Jul–Sep 2023)",
    "Q19 (Oct–Dec 2023)",
    "Q20 (Jan–Mar 2024)",
    "Q21 (Apr–Jun 2024)",
    "Q22 (Jul–Sep 2024)",
    "Q23 (Oct–Dec 2024)",
    "Q24 (Jan–Mar 2025)"
  )
  
  outcome_labels <- c(
    "i_mmr1_24m"   = "MMR1 (24 months)",
    "i_mmr1_5y"    = "MMR1 (5 years)",
    "i_mmr2_5y"    = "MMR2 (5 years)",
    "i_rota"       = "Rotavirus (12 months)",
    "i_pcv_boost"  = "PCV booster (24 months)",
    "i_six_in_one" = "6-in-1 (12 months)"
  )
  
  y_min <- min(data$uptake*100, na.rm = TRUE)
  y_max <- max(data$uptake*100, na.rm = TRUE)
  
  ggplot(data, aes(x = time, y = uptake * 100, color = as.factor(interv))) +
    geom_line(aes(linetype = as.factor(interv)), size = 1) +
    geom_ribbon(aes(ymin = lcl * 100, ymax = ucl * 100, fill = as.factor(interv)), alpha = 0.3) +
    geom_vline(xintercept = endpre, linetype = "dashed", color = "black") +
    facet_wrap(~Outcome, scales = "free_y", ncol = ncol,  
               labeller = labeller(Outcome = outcome_labels)) +
    scale_color_manual(
      values = c("0" = "#440154FF", "1" = "#FDE725FF"),
      labels = c("Synthetic Control", "Intervention")
    ) +
    scale_fill_manual(
      values = c("0" = "#440154FF", "1" = "#FDE725FF"),
      labels = c("Synthetic Control", "Intervention")
    ) +
    scale_linetype_manual(
      values = c("0" = "solid", "1" = "dashed"),
      labels = c("Synthetic Control", "Intervention")
    ) +
    
    scale_x_continuous(
      breaks = 1:24,
      labels = quarter_labels,
      expand = c(0,0)
    ) +
    
    labs(
      x = "Quarter",
      y = "Quarterly vaccination rate %",
      color = "Group",
      fill = "Group",
      linetype = "Group"
    ) +
    theme_classic() +
    theme(
      text = element_text(size=16),
      legend.position="bottom",
      legend.box = "horizontal",
      legend.text = element_text(size=24),
      legend.title = element_blank(),
      legend.key.size = unit(1, "cm"),
      strip.background = element_blank(),
      axis.text.x = element_text(angle = 70, hjust = 1)
    ) +
    theme(strip.background = element_blank(), legend.position = "bottom")
}








# Forest plotting  --------------------------------------------------------


forest_plot_results <- function(all_models, title, show_model = TRUE) {
  group_order <- names(all_models)
  
  results_df <- purrr::imap_dfr(
    all_models,
    function(models, model_label) {
      
      if (!is.list(models) || inherits(models, "microsynth")) models <- list(models)
      
      purrr::imap_dfr(models, function(model, sub_label) {
        if (is.null(model$Results)) return(NULL)
        
        model$Results %>%
          as.data.frame() %>%
          rename_with(~ stringr::str_remove(., "^X\\d+\\.")) %>%
          tibble::rownames_to_column("Outcome") %>%
          dplyr::mutate(
            Pct.Chng = dplyr::if_else(is.na(Pct.Chng), NA_real_, round(Pct.Chng * 100, 2)),
            Linear.Lower = Linear.Lower * 100,
            Linear.Upper = Linear.Upper * 100,
            Perm.Lower = Perm.Lower * 100,
            Perm.Upper = Perm.Upper * 100
          ) %>%
          dplyr::mutate(dplyr::across(where(is.numeric), ~ round(., 4))) %>%
          dplyr::mutate(
            Outcome = dplyr::case_when(
              Outcome %in% c("i_mmr1_24m", "mmr1_24m") ~ "MMR1 at 24 m",
              Outcome %in% c("i_mmr1_5y", "mmr1_5y") ~ "MMR1 at 5 y",
              Outcome %in% c("i_mmr2_5y", "mmr2_5y") ~ "MMR2 at 5 y",
              Outcome %in% c("i_rota", "rota") ~ "Rota at 12 m",
              Outcome %in% c("i_six_in_one", "six_in_one") ~ "6-in-1 at 12 m",
              Outcome %in% c("i_pcv_boost", "pcv_boost") ~ "PCV Booster at 24 m",
              TRUE ~ Outcome
            ),
            Model = factor(model_label, levels = group_order)
          )
      })
    }
  )
  
  results_df$` ` <- paste(rep(" ", 20), collapse = " ")
  
  if (show_model) {
    results_df <- results_df %>%
      dplyr::group_by(Model) %>%
      dplyr::group_split() %>%
      purrr::map_df(function(df) {
        header_row <- df[1, ]
        header_row <- header_row %>%
          dplyr::mutate(
            dplyr::across(where(is.numeric), ~ NA_real_),
            dplyr::across(where(is.character), ~ "")
          )
        header_row$Model <- unique(df$Model)
        dplyr::bind_rows(header_row, df %>% dplyr::mutate(Model = ""))
      }) %>%
      dplyr::ungroup()
  } else {
    results_df <- results_df %>%
      dplyr::mutate(Model = "")
  }
  
  results_df$`Pct.Chng (95% CI)` <- ifelse(
    is.na(results_df$Pct.Chng), "",
    sprintf("%.2f (%.2f to %.2f)",
            results_df$Pct.Chng, results_df$Perm.Lower, results_df$Perm.Upper)
  )
  
  library(forestploter)
  
  plot_table <- results_df %>%
    dplyr::select(
      dplyr::any_of(if (show_model) "Model" else NULL),
      Outcome,
      `Intervention` = Trt,
      `Synth Con` = Con,
      ` `,
      `Pct.Chng (95% CI)`,
      Perm.pVal
    ) %>%
    dplyr::mutate(
      Perm.pVal = dplyr::case_when(
        is.na(Perm.pVal) ~ "",
        Perm.pVal <= 0.001 ~ "<0.001*",
        Perm.pVal <= 0.05 ~ paste0(formatC(Perm.pVal, format = "f", digits = 3), "*"),
        TRUE ~ formatC(Perm.pVal, format = "f", digits = 3)
      )
    ) %>%
    dplyr::mutate(dplyr::across(where(is.numeric), ~ ifelse(is.na(.), "", .)))
  
  ci_col <- if (show_model) 5 else 4
  
  f_plot <- forest(
    plot_table,
    est = results_df$Pct.Chng,
    lower = results_df$Perm.Lower,
    upper = results_df$Perm.Upper,
    ci_column = ci_col,
    ref_line = 0,
    arrow_lab = c("Uptake Reduced", "Uptake Increased"),
    xlab = "Pct.Change",
    title = title
  )
  
  f_plot <- add_border(f_plot, part = "header", row = 1, gp = gpar(lwd = 1))
  f_plot <- wrap_elements(full = f_plot)
  
  return(f_plot)
}



# Matching quality -----
# Extract ε (epsilon) because microsynth reports this matching-quality metric in a warning
# message (not in the returned model object). We capture it so we can quantify and compare
# pre-intervention match quality across model specifications (smaller ε = better match).

run_microsynth <- function(...) {
  eps_value <- NA
  iter_value <- NA
  
  w <- withCallingHandlers(
    microsynth(...),
    warning = function(wrn) {
      msg <- conditionMessage(wrn)
      if (grepl("Failed to converge: eps=", msg)) {
        eps_value <<- as.numeric(sub(".*eps=([0-9.]+).*", "\\1", msg))
        iter_value <<- as.numeric(sub(".*in ([0-9]+) iterations.*", "\\1", msg))
      }
    }
  )
  
  w$info$eps <- eps_value
  w$info$iter <- iter_value
  
  return(w)
}


# Extracting effect for results section -----------------------------------
extract_effect <- function(model, outcome) {
  model$Results %>%
    as.data.frame() %>%
    rownames_to_column("Outcome") %>%
    rename_with(~ str_remove(., "^X\\d+\\.")) %>%
    mutate(
      Pct.Chng     = if_else(is.na(Pct.Chng), NA, Pct.Chng *100),
      lci   = Perm.Lower * 100,
      uci   = Perm.Upper * 100) %>%
    mutate(across(
      c(Trt, Con, Pct.Chng, Perm.Lower, Perm.Upper),~ as.numeric(.x))) %>% 
    filter(Outcome == outcome) %>% 
    summarise(
      pct  = round(Pct.Chng, 2),
      lci  = round(lci, 2),
      uci  = round(uci, 2),
      p_val = round(Perm.pVal, 3),
      treat_n = round(Trt),
      synth_n = round(Con),
      abs_diff = (treat_n - synth_n),
      abs_diff_lcl = round(synth_n * (lci / 100)),
      abs_diff_ucl = round(synth_n * (uci / 100))
    )
}
