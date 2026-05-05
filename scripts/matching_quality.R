gp <- readRDS("data/analysis_data.rds")

# Parameters
target_crs <- 27700  # British national grid (meters)

##  Load lookups and filter to Liverpool LSOAs ----
## LSOA 2011 to LSOA 2021 lookup. 
lsoa_lookup <- read_csv(
  "data/interv_gp/LSOA_(2011)_to_LSOA_(2021)_to_Local_Authority_District_(2022)_Best_Fit_Lookup_for_EW_(V2) (1).csv") %>%
  filter(LAD22NM == "Liverpool") %>%     # filtering to keep only Liverpool
  rename(lsoa11 = LSOA11CD, lsoa21 = LSOA21CD)

## loading LSOA 2021 shapfile 
lsoa_sf <- read_sf(
  "data/interv_gp/Lower_layer_Super_Output_Areas_(December_2021)_Boundaries_EW_BFC_(V10)/Lower_layer_Super_Output_Areas_(December_2021)_Boundaries_EW_BFC_(V10).shp"
) %>%
  rename(lsoa21 = LSOA21CD) %>%
  filter(lsoa21 %in% lsoa_lookup$lsoa21) %>%
  st_transform(target_crs) %>% # project to meters
  left_join(lsoa_lookup, by = "lsoa21")

# derive the centroid for each LSOA. 
lsoa_pts <- lsoa_sf %>%
  st_centroid() %>%            
  select(lsoa21, lsoa11)                        

# Extract the coordinates 
lsoa_xy <- st_coordinates(lsoa_pts)    
lsoa_pts <- lsoa_pts %>%
  mutate(X = lsoa_xy[, "X"], Y = lsoa_xy[, "Y"])

# Filtering practices that are physically located in a Liverpool LSOA11 to avoid having non-Liverpool practices contributing to the weighted centroids. 
# To do so, we use 2 look-up tables: 1) A GP to PC table and 2) a PC to LSOA one

## loading GP to postcode look-up table
pc_gp <- read_csv("data/interv_gp/epraccur.csv", col_names = F) %>%
  rename(gpcode = X1, gpname = X2, pc = X10)

## loading postcode to LSOA lookup table
pc_lsoa <- read_csv("data/interv_gp/PCD_OA_LSOA_MSOA_LAD_NOV22_UK_LU.csv") %>%
  select(pc = pcds, lsoa11 = lsoa11cd) %>% 
  left_join(lsoa_lookup, by = "lsoa11")

## filtering to keep only practices in Liverpool
gp_in_liverpool <- pc_gp %>%
  left_join(pc_lsoa, by = "pc") %>%
  filter(lsoa11 %in% lsoa_lookup$lsoa11) %>%    # practice location in Liverpool
  distinct(gpcode)

# Calculated the patients weighted centroids for the GPs based on their respective registered population

## Loading the GP-registered population for July 2022 (intervention starting period)

gp_pop <- read_csv(
  "data/interv_gp/gp-reg-pat-prac-lsoa-male-female-July-2022/gp-reg-pat-prac-lsoa-all.csv",
  show_col_types = FALSE) %>%
  rename(
    gpcode = PRACTICE_CODE,
    gpname = PRACTICE_NAME,
    lsoa21 = LSOA_CODE,
    pop    = `NUMBER_OF_PATIENTS`) %>% filter(gpcode %in% gp_in_liverpool$gpcode) %>% 
  left_join(lsoa_lookup, by = "lsoa21") %>%
  mutate(lsoa21 = coalesce(lsoa21, lsoa11)) %>%
  select(gpcode, gpname, lsoa21, pop)



## joining lsoa points (x, y) to GP-LSOA weights
gp_weights_xy <- gp_pop %>%
  left_join(lsoa_pts %>% st_drop_geometry() %>% select(lsoa21, X, Y), by = "lsoa21") %>%
  filter(!is.na(X), !is.na(Y))

# Weighted x and y per GP:
gp_weighted_xy <- gp_weights_xy %>%
  group_by(gpcode, gpname) %>%
  summarise(
    wsum = sum(pop, na.rm = TRUE),
    xbar = sum(pop * X, na.rm = TRUE) / wsum,
    ybar = sum(pop * Y, na.rm = TRUE) / wsum,
    .groups = "drop"
  ) %>%
  # create sf point for each GP weighted centroid
  st_as_sf(coords = c("xbar", "ybar"), crs = target_crs)

# Defining the intervention sites using the national statistics postcode lookup table
nspl <- read_sf("data/interv_gp//NSPL_Online_Latest_Centroids/NSPL_Online_Latest_Centroids.shp") %>%
  filter(CTY == "E11000002") %>%     # Liverpool
  st_transform(target_crs)

intervention_site_postcodes <- c("L6 1AE", "L8 2UU", "L5 2PY", "L8 1TH",
                                 "L8 0TP", "L7 8TQ", "L7 6HD", "L15 0EE", "L11 2RY")

intervention_sites <- nspl %>%
  filter(PCDS %in% intervention_site_postcodes)


# Setting a buffer zone around the intervention site, and assigning intervention status
radius_km <- seq(0.5, 3, by = 0.1)


# A function that returns the intervention GPs for a given radius
get_intervention_gps <- function(r_km) {
  sites_buff <- st_buffer(intervention_sites, dist = r_km * 1000) %>%
    mutate(radius_km = r_km)
  
  gp_flag <- gp_weighted_xy %>%
    mutate(in_buffer = lengths(st_intersects(., sites_buff)) > 0) %>%
    transmute(radius_km = r_km, gpcode, gpname, in_buffer)
  
  list(
    gp_flag = gp_flag,
    sites_buff = sites_buff
  )
}

# apply over all radii and bind results
all_results  <- map(radius_km, get_intervention_gps)

intervention_gps <- map_dfr(all_results, "gp_flag")
sites_buff <- map_dfr(all_results, "sites_buff")


buff_distances <- intervention_gps %>% filter(in_buffer == T) %>% st_drop_geometry() %>%
  count(radius_km, name = "n_gps_in_buffer") 





gtsave(buff_distances %>% rename(`Buffer Radius` = radius_km, `Number of Intervention GPs` =n_gps_in_buffer) %>% 
         gt() %>% 
         tab_caption("Number of intervention GPs within each buffer radius"), filename = "outputs/buff_distances.docx")


buffer_distances <- sort(unique(intervention_gps$radius_km))
buffer_distances


buffer_list <- lapply(buffer_distances, function(d){
  intervention_gps %>%
    filter(radius_km == d, in_buffer == TRUE) %>%
    pull(gpcode) %>%
    unique()
})

names(buffer_list) <- as.character(buffer_distances)



set.seed(12913)
gp_all <- readRDS("data/GP_MMR_IMD_ICB_POP_FING_ETHN_dist.rds")


cov.var <- c(
  'pct_boys_0_4', 'appoint_per_patient',
  'qof', 'satisfied', "caring", 'imd_2025',
  'pct_black', 'pct_asian', 'pct_mixed', 'distance')

# # Time-varying co-variates 
match.out <- c("i_mmr1_24m", "i_mmr1_5y", "i_mmr2_5y", "den_24m", "den_5y")    

# Result variable 
result.var <- c("i_mmr1_24m", "i_mmr1_5y", "i_mmr2_5y")    


eps_results <- list()

for(d in names(buffer_list)){
  
  cat("Running buffer distance:", d, "\n")
  
  interv_set <- buffer_list[[d]]
  
  c_n_gpcodes <- gp_all %>% filter(pcn_name %in% c("NORTH LIVERPOOL PCN", "CENTRAL LIVERPOOL PCN") ) %>%
    select(gpcode) %>% distinct() %>% pull()

  c_n_gpcodes <- setdiff(c_n_gpcodes, interv_set)

  # # Create intervention variable & drop excluded PCNs
  gp_sub <- gp %>%
    mutate(interv = ifelse(gpcode %in% interv_set, 1, 0)) %>%
    # filter(time >= 13)
    filter(time >= 13 & !(gpcode %in% c_n_gpcodes))
  
  
  # Run microsynth model
  model <- tryCatch(
    run_microsynth(
      gp_sub,
      idvar       = "gpcode",
      timevar     = "time",
      intvar      = "interv",
      start.pre   = 13,
      end.pre     = 19,   # lag applied via shifted pre-period
      end.post    = 24,
      match.out   = match.out,
      match.covar = cov.var,
      result.var  = result.var,
      test        = "twosided",
      n.cores     = 1,
      use.backup  = TRUE,
      perm        = 250,
      confidence  = 0.95
    ),
    error = function(e) NULL
  )
  
  # Extract epsilon
  if (is.null(model)) {
    eps_results[[d]] <- NA
  } else {
    eps_val <- model$info$eps
    if (is.na(eps_val)) eps_val <- 0   # perfect match
    eps_results[[d]] <- eps_val
  }
}



eps_df <- tibble(
  distance = as.numeric(names(eps_results)),
  epsilon  = unlist(eps_results)
)

ggplot(eps_df, aes(x = distance, y = epsilon)) +
  geom_line(size = 1, color = "#440154FF") +
  geom_point(size = 1, color = "#FDE725FF") +
  labs(
    x = "Buffer distance (km)",
    y = "Epsilon (matching imbalance)",
    title = "Matching quality across buffer distances"
  ) +
  scale_x_continuous(limits = c(0.5, 3),
                     breaks = seq(0.5, 3, by = 0.1), 
                     expand = c(0, 0)
  ) +
  theme_classic() + 
  theme(text = element_text(size = 16), 
        axis.text.x = element_text(angle = 45, hjust = 1))

eps_df %>% gt() %>% gtsave(filename = "outputs/eps.docx")

ggsave(plot = last_plot(), filename = "outputs/matching_quality_postcovid.pdf", width = 10)
