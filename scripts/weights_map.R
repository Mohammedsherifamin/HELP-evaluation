rm(list = ls())
library(microsynth)
library(janitor)
library(tidyverse)
library(readODS)
library(sf)
library(readxl)
library(patchwork)
library(ggspatial)
library(cowplot)
source("scripts/functions.R")

# Loading the analysis data and running the model 
gp <- readRDS("data/analysis_data.rds")
interv_1km <- readRDS("data/interv_gp/1km_interv.rds")


gp <- gp %>%  
  mutate(
    interv = case_when(gpcode %in% interv_1km ~ 1, 
                       TRUE ~ 0))


set.seed(12913)
gp_all <- readRDS("data/GP_MMR_IMD_ICB_POP_FING_ETHN_dist.rds")


c_n_gpcodes <- gp_all %>% filter(pcn_name %in% c("NORTH LIVERPOOL PCN", "CENTRAL LIVERPOOL PCN") ) %>%
  select(gpcode) %>% distinct() %>% pull()

c_n_gpcodes <- setdiff(c_n_gpcodes, interv_1km)

# Time-invariant covariates
cov.var <- c(
  'pct_boys_0_4', 'appoint_per_patient',
  'qof', 'satisfied', "caring", 'imd_2025',
  'pct_black', 'pct_asian', 'pct_mixed', 'distance')

# # Time-varying co-variates 
match.out <- c("i_mmr1_24m", "i_mmr1_5y", "i_mmr2_5y", "den_24m", "den_5y")    

# Result variable 
result.var <- c("i_mmr1_24m", "i_mmr1_5y", "i_mmr2_5y")    


count_1km_exgp <- microsynth(gp %>% 
                               filter(time >= 13 & gpcode %notin% c_n_gpcodes), 
                             idvar = "gpcode", timevar = "time", 
                             intvar = "interv", start.pre = 13, end.pre = 19, 
                             end.post = 24, match.out = match.out, match.covar = cov.var,
                             result.var = result.var, test = "twosided", n.cores = 1, use.backup = T,
                             omnibus.var = result.var, perm = 250, confidence = 0.95)


plot_microsynth(count_1km_exgp, end.pre = 19, legend.spot = "topleft", plot.var = "i_mmr2_5y")



# Extracting the weights
weights <- as.data.frame(count_1km_exgp$w$Weight) %>%
  tibble::rownames_to_column(var = "gpcode") %>% select(gpcode, Main)


gp_w <- gp %>% filter(gpcode %notin% c_n_gpcodes) %>% 
  left_join(weights, by = "gpcode")


# Parameters
target_crs <- 27700  # British national grid (meters)


# Loading LSOA lookup and Shape files
lsoa_lookup <- read_csv(
  "data/interv_gp/LSOA_(2011)_to_LSOA_(2021)_to_Local_Authority_District_(2022)_Best_Fit_Lookup_for_EW_(V2) (1).csv") %>%
  rename(lsoa11 = LSOA11CD, lsoa21 = LSOA21CD)


lsoa_sf <- read_sf(
  "data/interv_gp/Lower_layer_Super_Output_Areas_(December_2021)_Boundaries_EW_BFC_(V10)/Lower_layer_Super_Output_Areas_(December_2021)_Boundaries_EW_BFC_(V10).shp"
) %>%
  rename(lsoa21 = LSOA21CD) %>%
  filter(lsoa21 %in% lsoa_lookup$lsoa21) %>%
  st_transform(target_crs)


# loading GP to postcode look-up table
pc_gp <- read_csv("data/interv_gp/epraccur.csv", col_names = F) %>%
  rename(gpcode = X1, gpname = X2, pc = X10)


# NSPL to link the GPs to Physcial points on the a map 
nspl <- read_sf("data/interv_gp//NSPL_Online_Latest_Centroids/NSPL_Online_Latest_Centroids.shp") %>%
  st_transform(target_crs)




# helper to standardise postcodes for joining
pc_clean <- function(x) {
  x %>%
    str_trim() %>%
    str_to_upper() %>%
    str_replace_all("\\s+", "")   
}

# NSPL: pick the postcode column
nspl_pts <- nspl %>%
  mutate(pc = pc_clean(PCDS)) %>%   
  select(pc, geometry)


# GP postcodes 
pc_gp2 <- pc_gp %>%
  mutate(pc = pc_clean(pc)) %>%
  select(gpcode, pc)


# GP weights + geometry 
gp_pts <- gp_w %>%
  # keep only practices with a non-zero weight
  mutate(weight = as.numeric(Main), 
         is_intervention = gpcode %in% interv_1km) %>%
  filter(!is.na(weight), weight > 0) %>%
  left_join(pc_gp2, by = "gpcode") %>%
  filter(!is.na(pc)) %>%
  left_join(nspl_pts, by = "pc") %>%
  filter(!st_is_empty(geometry), !is.na(st_coordinates(geometry)[,1])) %>%
  st_as_sf(crs = target_crs)


# Loading IMD 2025 and joining to shp files 
imd <- read_csv("data/IMD/File_7_IoD2025_All_Ranks_Scores_Deciles_Population_Denominators.csv") %>% 
  select(lsoa21 =  "LSOA code (2021)" , imd_score = "Index of Multiple Deprivation (IMD) Score", 
         imd_decile = "Index of Multiple Deprivation (IMD) Decile (where 1 is most deprived 10% of LSOAs)")


lsoa_sf <- lsoa_sf %>%
  left_join(imd, by = "lsoa21") %>%
  mutate(
    imd_quintile = 6 - ntile(imd_score, 5),
    imd_quintile = factor(
      imd_quintile,
      levels = 1:5,
      labels = c("Q1 Most deprived", "Q2", "Q3", "Q4", "Q5 Least deprived")))



# England Map -------------------------------------------------------------
england <- ggplot() +
  geom_sf(
    data = lsoa_sf %>% filter(str_starts(lsoa21, "E")),
    fill = "gray99",
    linewidth = 0.01,
    alpha = 0.7) +
  
  # scale_fill_manual(
  #   name = "IMD quintile (2025)",
  #   values = c(
  #     "Q1 Most deprived"   = "#0B3C5D",
  #     "Q2"                 = "#2E6F95",
  #     "Q3"                 = "#5FA8D3",
  #     "Q4"                 = "#9BCDEB",
  #     "Q5 Least deprived"  = "#E6F2FA" )) +
  
  
  geom_sf(
    data = gp_pts %>% filter(!is_intervention),
    aes(color = weight),
    size = 0.55,
    alpha = 0.8) +
  
  
  scale_color_gradientn(
    name = "Synthetic control weight",
    colours = c(
      "yellow",  
      "#FD8D3C",  
      "#E31A1C", 
      "#800026" ) ) +
  
  # Intervention practices 
  geom_sf(
    data = gp_pts %>% filter(is_intervention),
    aes(shape = "Intervention practices"),
    color = "cyan",
    size = 0.55,
    alpha = 0.8) +
  
  scale_shape_manual(
    name = "",
    values = c("Intervention practices" = 16)) +
  
  theme_void() +
  
  guides(
    fill  = guide_legend(order = 1),
    color = guide_colorbar(order = 2),
    shape = guide_legend(
      order = 3,
      override.aes = list(color = "cyan", size = 3))) +
  
  coord_sf(crs = target_crs, expand = FALSE) +
  
  
  labs(title = "England") +
  theme(legend.position = "right", plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) + 
  
  annotation_scale(
    location = "br",  
    style = "bar") +
  
  annotation_north_arrow(
    location = "tl", 
    which_north = "true", 
    style = north_arrow_fancy_orienteering)




# Liverpool Map -----------------------------------------------------------
#Liverpool-specific data

liverpool_lsoa <- lsoa_sf %>% filter(str_detect(LSOA21NM, "Liverpool")) 
# liverpool_lsoa <- lsoa_sf %>% filter(str_detect(LSOA21NM, "Halton|Knowsley|Liverpool|St Helens|Sefton|Wirral")) 
liverpool_pts  <- gp_pts %>% st_intersection(st_union(liverpool_lsoa))

weight_limits <- range(gp_pts$weight, na.rm = TRUE)


# the Liverpool only 
liverpool <- ggplot() +
  geom_sf(
    data = liverpool_lsoa, fill = "gray99", linewidth = 0.02,  alpha = 0.7) +
  
  # scale_fill_manual(
  #   name = "IMD quintile (2025)",
  #   values = c(
  #     "Q1 Most deprived"   = "#0B3C5D",
  #     "Q2"                 = "#2E6F95",
  #     "Q3"                 = "#5FA8D3",
  #     "Q4"                 = "#9BCDEB",
  #     "Q5 Least deprived"  = "#E6F2FA" )) +
  
  
  # Non-intervention contributors 
  geom_sf(data = liverpool_pts %>% filter(!is_intervention), aes(color = weight), size = 2, alpha = 0.8) +
  
  scale_color_gradientn( name = "Synthetic control weight",
                         colours = c(
                           "yellow",  
                           "#FD8D3C",  
                           "#E31A1C", 
                           "#800026" ) ) +
  
  # Intervention practices 
  geom_sf(data = gp_pts %>% filter(gpcode %in% interv_1km), aes(shape = "Intervention practices"),
          color = "cyan",
          size = 2,
          alpha = 0.8) +
  
  scale_shape_manual(
    name = "",
    values = c("Intervention practices" = 16)) +
  
  guides(
    fill  = guide_legend(order = 1),
    color = guide_colorbar(order = 2),
    shape = guide_legend(
      order = 3,
      override.aes = list(color = "cyan", size = 3))) +
  
  coord_sf(crs = target_crs, expand = FALSE) + 
  labs(title = "Liverpool") +
  theme_void() +
  theme( plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
         panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
         legend.position = "none") 



liverpool

final_map <- england + 
  inset_element(
    liverpool, 
    left = 0.02,   
    bottom = 0.50, 
    right = 0.35,  
    top = 0.88,  
    align_to = "panel"
  )


ggsave(plot = last_plot(), filename = "outputs/weights_map.jpeg", dpi = 1200, width= 7, height = 7)
ggsave(plot = last_plot(), filename = "outputs/weights_map.pdf", width = 7, height = 7)
