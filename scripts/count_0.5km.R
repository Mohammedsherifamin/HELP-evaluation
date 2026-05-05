gp <- readRDS("data/analysis_data.rds")

interv_0.5km <- readRDS("data/interv_gp/0.5km_interv.rds")


gp <- gp %>%  
  mutate(
    interv = case_when(gpcode %in% interv_0.5km  ~ 1, 
                       TRUE ~ 0)
  )

## MMR -----
set.seed(12913)
gp_all <- readRDS("data/GP_MMR_IMD_ICB_POP_FING_ETHN_dist.rds")


c_n_gpcodes <- gp_all %>% filter(pcn_name %in% c("NORTH LIVERPOOL PCN", "CENTRAL LIVERPOOL PCN") ) %>%
  select(gpcode) %>% distinct() %>% pull()

c_n_gpcodes <- setdiff(c_n_gpcodes, interv_0.5km)

# Time-invariant covariates
cov.var <- c(
  'pct_boys_0_4', 'appoint_per_patient',
  'qof', 'satisfied', "caring", 'imd_2025',
  'pct_black', 'pct_asian', 'pct_mixed', 'distance')

# # Time-varying co-variates 
match.out <- c("i_mmr1_24m", "i_mmr1_5y", "i_mmr2_5y", "den_24m", "den_5y")    

# Result variable 
result.var <- c("i_mmr1_24m", "i_mmr1_5y", "i_mmr2_5y")    


count_0.5km_exgp <- run_microsynth(gp %>% 
                               filter(time >= 13 & gpcode %notin% c_n_gpcodes), 
                             idvar = "gpcode", timevar = "time", 
                             intvar = "interv", start.pre = 13, end.pre = 19, 
                             end.post = 24, match.out = match.out, match.covar = cov.var,
                             result.var = result.var, test = "twosided", n.cores = 1, use.backup = T,
                             omnibus.var = NULL, perm = 250, confidence = 0.95)


## Other vaccines -------
set.seed(12913)


# Time-invariant covariates
cov.var <- c(
  'pct_boys_0_4', 'appoint_per_patient',
  'qof', 'satisfied', "caring", 'imd_2025',
  'pct_black', 'pct_asian', 'pct_mixed', 'distance')

# Time-varinat covariates (including outcomes)
match.out <- c("i_rota", "i_six_in_one", "i_pcv_boost", "den_12m", "den_24m")    

# Result variable 
result.var <- c("i_rota", "i_six_in_one", "i_pcv_boost")     


count_0.5km_exgp_other <- run_microsynth(gp %>% 
                                     filter(time >= 13 & gpcode %notin% c_n_gpcodes), 
                                   idvar = "gpcode", timevar = "time", 
                                   intvar = "interv", start.pre = 13, end.pre = 19, 
                                   end.post = 24, match.out = match.out, match.covar = cov.var,
                                   result.var = result.var, test = "twosided", n.cores = 1, use.backup = T, 
                                   omnibus.var = NULL, perm = 250, confidence = 0.95)




### table -----
count_0.5km_exgp_tbl <- extract_syn_results(count_0.5km_exgp, count_0.5km_exgp_other) %>% 
  gt() %>% 
  tab_header(
    title = "Estimated Intervention Effects on Vaccination Counts",
    subtitle = "Model excludes practices from Central & North Liverpool PCNs"
  ) %>%
  tab_footnote(
    footnote = glue(
      "Model parameters: start.pre = {count_0.5km_exgp$info$start.pre}, ",
      "end.pre = {count_0.5km_exgp$info$end.pre}, ",
      "end.post = {count_0.5km_exgp$info$end.post}"),
    locations = cells_title(groups = "subtitle"))



### Plots -------
out_den <- c(
  "i_mmr1_24m" = "den_24m",
  "i_mmr1_5y"  = "den_5y",
  "i_mmr2_5y"  = "den_5y")

p_count_0.5km_exgp  <- extract_svy_ratios(
  count_0.5km_exgp, outcome_denom_map = out_den,  synth.df = gp %>% 
    filter(time >= 13 & gpcode %notin% c_n_gpcodes))

count_0.5km_exgp_plot <- plot_svy_ratios(p_count_0.5km_exgp, ncol = 1, endpre = count_0.5km_exgp$info$end.pre)



out_den <- c(
  "i_rota"     =  "den_12m", 
  "i_six_in_one"  =  "den_12m",
  "i_pcv_boost" = "den_24m"
)

p_count_0.5km_exgp_other  <- extract_svy_ratios(
  count_0.5km_exgp_other, outcome_denom_map = out_den,  synth.df = gp %>% 
    filter(time >= 13 & gpcode %notin% c_n_gpcodes))

count_0.5km_exgp_plot_other <- plot_svy_ratios(p_count_0.5km_exgp_other, ncol = 1, 
                                             endpre = count_0.5km_exgp_other$info$end.pre)

(count_0.5km_exgp_plot + count_0.5km_exgp_plot_other) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")
