remotes::install_github("pfmc-assessments/indexwc")
library(nwfscSurvey)
library(sdmTMB)
library(dplyr)
bio_dat <- pull_bio_cache(common_name = "lingcod")
names(bio_dat) <- tolower(names(bio_dat))
bio_dat$yday <- lubridate::yday(bio_dat$date)

#grabbing from the indexwc output
load("output/coastwide/2019/lingcod/wcgbts/lat_32.1-48.5_years_2003-2019/delta_gamma/diagnostics/predictions.rdata")
head(predictions)

# filter out only females
bio_dat <- dplyr::filter(bio_dat, sex == "F", !is.na(length_cm))
# add y_lin -- 152 in table 1 Andrea paper; 103.645 females for assessment
bio_dat$y_lin <- log(103.645 - bio_dat$length_cm)

haul_dat <- pull_haul_cache()
names(haul_dat) <- tolower(names(haul_dat))

bio_dat <- sdmTMB::add_utm_columns(
  bio_dat,
  ll_names = c("longitude_dd", "latitude_dd")
) |>
  dplyr::filter(!is.na(age))

bio_dat <-bio_dat |>
  dplyr::filter(!is.na(age))

mesh <- sdmTMB::make_mesh(bio_dat, xy_cols = c("X", "Y"), cutoff = 20)
mesh$mesh$n

# center and scale data
year_pars <- c(mean(bio_dat$year), sd(bio_dat$year))
bio_dat$year_scaled <- (bio_dat$year - year_pars[1]) / year_pars[2]
age_pars <- c(mean(bio_dat$age), sd(bio_dat$age))
bio_dat$age_scaled <- (bio_dat$age - age_pars[1]) / age_pars[2]
bio_dat$fyear <- as.factor(bio_dat$year)
doy_pars <- c(mean(bio_dat$yday), sd(bio_dat$yday))
bio_dat$doy_scaled <- (bio_dat$yday - doy_pars[1]) / doy_pars[2]
bio_dat$doy_scaled2 <- bio_dat$doy_scaled^2
bio_dat<-bio_dat%>%mutate(N_S=ifelse(latitude_dd > 40.166667, "N", "S"))
names(bio_dat)
map_list <- list(
  ln_tau_Z = factor(rep(1, times = length(unique(bio_dat$fyear))))
)
fit <- sdmTMB(
  y_lin ~ age + fyear + doy_scaled + doy_scaled2, # reminder fyear is what you would pull out for comparison with a cov
  spatial_varying = ~ age:fyear, #if fyear is dropped in main you want to keep it here!! (could put cov here too but less meaningful)
  mesh = mesh,
  spatiotemporal = "off", # off because random fields in age interaction
  data = bio_dat,
  time = "year",
  control = sdmTMBcontrol(newton_loops = 1, map = map_list)
)

# use indexwc grid
# grid <- indexwc::california_current_grid |>
#   dplyr::filter(area_km2_WCGBTS == 4)
pred_grid <- predictions
grid <- dplyr::mutate(pred_grid, age = 8, doy_scaled = 0, doy_scaled2 = 0)
unique(grid$year)
grid$fyear <- as.factor(grid$year)
#grid <- dplyr::rename(grid, X = x, Y = y)

# here is where we'd fit the biomass model, and use the predicted biomass as the RHS of the vector argument
indxN <- sdmTMB::get_weighted_average(
  fit,
  newdata = grid,
  vector =exp(grid$est), #you need to replace this with index WC
  bias_correct = TRUE
)

L_inf <- 103.645
pred_age <- 8

indx_k <- indxN |>
  dplyr::mutate(
    # Back-transform y_lin to (L_inf - L_at_age8)
    diff_est = exp(est),
    # Calculate k (assuming t0 = 0)
    k_est = -log(diff_est / L_inf) / pred_age,
    # Approximate SE via Delta method
    k_se = se / pred_age,
    k_lwr = -log(exp(upr) / L_inf) / pred_age,
    k_upr = -log(exp(lwr) / L_inf) / pred_age
  )



ggplot(indx_k, aes(x = year, y = est)) +
  geom_line(color = "#2c7fb8", linewidth = 1) +
  geom_point(color = "#2c7fb8", size = 2) +
  labs(
    title = "Time Series Plot",
    x = "Date",
    y = "Value"
  ) +
  theme_minimal()

write.csv(indx, "output/index_coastwide2019.csv", row.names = FALSE)

library(ggplot2)
library(dplyr)
library(tidyr)
library(sf)
library(rnaturalearth)

# 1. Base map projected to WGS84 / UTM Zone 10N (meters) to match sdmTMB's X*1000, Y*1000
map_data <- rnaturalearth::ne_states(
  returnclass = "sf", 
  country = "united states of america"
)

usa <- suppressWarnings(suppressMessages(
  sf::st_crop(map_data, c(xmin = -145, ymin = 25, xmax = -117, ymax = 48.5))
))
coast <- sf::st_transform(usa, crs = 32610) # Matches UTM Zone 10N in meters

# 2. Extract main 'age' fixed effect directly
age_coef <- sdmTMB::tidy(fit, effects = "fixed") |>
  dplyr::filter(term == "age") |>
  dplyr::pull(estimate)

# 3. Ungroup predictions (CRITICAL: removes rowwise grouping) & filter to ref_year
ref_year <- 2007

pred_grid_sub <- predictions |>
  dplyr::ungroup() |>
  dplyr::filter(year == ref_year) |>
  dplyr::mutate(
    age = 8,
    doy_scaled = 0,
    doy_scaled2 = 0,
    fyear = as.factor(year)
  ) |>
  # Add UTM columns using the exact column names ('longitude', 'latitude') and zone 10
  sdmTMB::add_utm_columns(
    ll_names = c("longitude", "latitude"),
    utm_zone = 10
  )

# 4. Predict from the fitted model
preds <- predict(fit, newdata = pred_grid_sub)

# 5. Reshape and calculate k_f, spatial mean, and annual anomalies
k_df <- preds |>
  tidyr::pivot_longer(
    cols = dplyr::starts_with("zeta_s_age:fyear"),
    names_to = "zeta_s_fyear",
    names_prefix = "zeta_s_age:fyear",
    values_to = "zeta_s_vals"
  ) |>
  dplyr::mutate(
    zeta_vals_age_adjusted = zeta_s_vals + age_coef,
    k_f = -zeta_vals_age_adjusted
  ) |>
  dplyr::group_by(X, Y) |>
  dplyr::mutate(
    mean_k_f = mean(k_f, na.rm = TRUE),
    sd_k_f = sd(k_f, na.rm = TRUE),
    anom_k_f = k_f - mean_k_f
  ) |>
  dplyr::ungroup()

# 6. Plot 1: Spatial Average k_f Map
p1 <- ggplot(coast) +
  geom_sf() +
  geom_tile(
    data = dplyr::distinct(k_df, X, Y, mean_k_f),
    aes(x = X * 1000, y = Y * 1000, fill = mean_k_f),
    width = 4000,
    height = 4000
  ) +
  scale_fill_viridis_c(name = expression(k[f])) +
  coord_sf(
    xlim = range(k_df$X * 1000),
    ylim = range(k_df$Y * 1000)
  ) +
  theme_bw() +
  xlab("") + ylab("") +
  labs(subtitle = "Lingcod (Female) - Mean Growth Rate") +
  theme(
    strip.background = element_rect(fill = "white"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# 7. Plot 2: Annual k_f Anomaly Faceted Map
p2 <- ggplot(coast) +
  geom_sf() +
  geom_tile(
    data = k_df,
    aes(x = X * 1000, y = Y * 1000, fill = anom_k_f),
    width = 4000,
    height = 4000
  ) +
  scale_fill_gradient2(name = expression(k[f] ~ "anom")) +
  coord_sf(
    xlim = range(k_df$X * 1000),
    ylim = range(k_df$Y * 1000)
  ) +
  theme_bw() +
  facet_wrap(~ zeta_s_fyear) +
  xlab("") + ylab("") +
  labs(subtitle = "Lingcod (Female) - Annual Growth Anomaly") +
  theme(
    strip.background = element_rect(fill = "white"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Save plots
dir.create("plots/lingcod", recursive = TRUE, showWarnings = FALSE)
ggsave(p1, file = "plots/lingcod/k_st_avg_length_lingcod.png", width = 5, height = 7)
ggsave(p2, file = "plots/lingcod/k_st_annual_length_lingcod.png", width = 9, height = 9)
