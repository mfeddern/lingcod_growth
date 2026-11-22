remotes::install_github("pfmc-assessments/indexwc")
library(nwfscSurvey)
library(sdmTMB)
library(dplyr)
library(dplyr)
library(tidyr)
library(sf)
library(rnaturalearth)
library(sdmTMB)
library(ggplot2)
bio_dat_raw <- pull_bio(common_name = "lingcod")
names(bio_dat_raw) <- tolower(names(bio_dat_raw))
bio_dat_raw$yday <- lubridate::yday(bio_dat_raw$date)

#grabbing from the indexwc output
load("output/coastwide/2019/lingcod/wcgbts/lat_32.1-48.5_years_2003-2019/delta_gamma/diagnostics/predictions.rdata")
#load("output/south/2019/lingcod/wcgbts/lat_32.1-40.2_years_2003-2019/delta_gamma/diagnostics/predictions.rdata")
head(predictions)
exp(predictions$est)
sex_selection='M'
location_selection="Coastwide"
l_inf<- 103.645

run_growth<-function(sex_selection, location_selection, l_inf){
  
  message("Step 1/6: Filtering data and preparing variables...")
  # filter out only females
  bio_dat <- dplyr::filter(bio_dat_raw, sex == sex_selection)
  # add y_lin -- 152 in table 1 Andrea paper
  lin<- l_inf
  bio_dat$y_lin <- log(lin - bio_dat$length_cm)
  
  bio_dat <- sdmTMB::add_utm_columns(
    bio_dat,
    ll_names = c("longitude_dd", "latitude_dd")
  ) |>
    dplyr::filter(!is.na(age))
  
  #pulling in haul and bio dat
  haul_dat <- pull_haul()
  names(haul_dat) <- tolower(names(haul_dat))
  
  message("Step 2/6: Generating SPDE mesh...")
  #generating the mesh from biodat
  mesh <- sdmTMB::make_mesh(bio_dat, xy_cols = c("X", "Y"), cutoff = 20)
  mesh$mesh$n
  
  message("Step 3/6: Centering and scaling covariates...")
  # center and scale data
  year_pars <- c(mean(bio_dat$year), sd(bio_dat$year))
  bio_dat$year_scaled <- (bio_dat$year - year_pars[1]) / year_pars[2]
  age_pars <- c(mean(bio_dat$age), sd(bio_dat$age))
  bio_dat$age_scaled <- (bio_dat$age - age_pars[1]) / age_pars[2]
  bio_dat$fyear <- as.factor(bio_dat$year)
  doy_pars <- c(mean(bio_dat$yday), sd(bio_dat$yday))
  bio_dat$doy_scaled <- (bio_dat$yday - doy_pars[1]) / doy_pars[2]
  bio_dat$doy_scaled2 <- bio_dat$doy_scaled^2
  
  map_list <- list(
    ln_tau_Z = factor(rep(1, times = length(unique(bio_dat$fyear))))
  )
  
  message("Step 4/6: Fitting sdmTMB model...")
  # fitting the sdm
  fit <- sdmTMB(
    y_lin ~ age + fyear + doy_scaled + doy_scaled2, # reminder fyear is what you would pull out for comparison with a cov
    spatial_varying = ~ age:fyear, #if fyear is dropped in main you want to keep it here!! (could put cov here too but less meaningful)
    mesh = mesh,
    spatiotemporal = "off", # off because random fields in age interaction
    data = bio_dat,
    time = "year",
    control = sdmTMBcontrol(newton_loops = 1, map = map_list)
  )
  
  message("Step 5/6: Filtering grid and calculating biomass weights...")
  # Filter predictions grid by location BEFORE calculating biomass weights
  # Generating the grid for the coastwide index
  
  pred_grid <- predictions
  if (location == "Coastwide") {
    # No filtering: keep coastwide predictions
    pred_grid <- pred_grid
  } else if (location == "N") {
    pred_grid <- dplyr::filter(pred_grid, split_mendocino == "N")
  } else if (location == "S") {
    pred_grid <- dplyr::filter(pred_grid, split_mendocino == "S")
  } else {
    stop("Invalid location: must be 'N', 'S', or NULL.")
  }
  grid <- dplyr::mutate(pred_grid, age = 8, doy_scaled = 0, doy_scaled2 = 0)
  grid <- replicate_df(
    grid,
    time_name = "year",
    time_values = unique(bio_dat$year)
  )
  grid$fyear <- as.factor(grid$year)
  
  # 1. Clean grid with annual normalized biomass weights (sums to 1 per year)
  grid_age8 <- predictions |>
    sf::st_drop_geometry() |>
    dplyr::ungroup() |>
    dplyr::filter(year %in% unique(bio_dat$year)) |>
    dplyr::group_by(year) |>
    dplyr::mutate(
      age = 8,
      doy_scaled = 0,
      doy_scaled2 = 0,
      fyear = as.factor(year),
      biomass_wt = exp(est) / sum(exp(est), na.rm = TRUE)
    ) |>
    dplyr::ungroup()
  
  grid_age9 <- dplyr::mutate(grid_age8, age = 9)
  
  message("Step 6/6: Running get_index() and compiling k_index...")
  # 2. Run sdmTMB::get_index() with area = biomass_wt on your existing `fit`
  idx_8FC <- sdmTMB::get_index(
    fit,
    newdata = grid_age8,
    area = grid_age8$biomass_wt,
    bias_correct = FALSE
  )
  
  idx_9FC <- sdmTMB::get_index(
    fit,
    newdata = grid_age9,
    area = grid_age9$biomass_wt,
    bias_correct = FALSE
  )
  
  # Final biomass-weighted k index using sdmTMB::get_index()
  k_index <- data.frame(
    year = idx_8FC$year,
    k_est = idx_8FC$est - idx_9FC$est
  )%>%
    mutate(sex=sex_selection, location=location_selection)
  
  years <- sort(unique(grid_age8$year))
  n_yrs <- length(years)
  start_time <- Sys.time()
  
  k_summary <- lapply(seq_along(years), function(i) {
    yr <- years[i]
    yr_start <- Sys.time()
    
    # Print progress before starting the year
    message(sprintf("[%d/%d] Processing year %s...", i, n_yrs, yr))
    
    # Subset grid to a single year
    g8_yr <- dplyr::filter(grid_age8, year == yr)
    g9_yr <- dplyr::filter(grid_age9, year == yr)
    
    # Normalize biomass weights for this year
    # (Replace 'YOUR_BIOMASS_COL' with your actual column, and remove exp() if already natural scale)
    w <- exp(g8_yr$est)
    w <- w / sum(w)
    
    # Draw simulations for just this single year
    set.seed(42)
    s8 <- predict(fit, newdata = g8_yr, nsim = 200)
    
    set.seed(42)
    s9 <- predict(fit, newdata = g9_yr, nsim = 200)
    
    # Compute weighted average of k = (y_8 - y_9) across space for each draw
    yr_draws <- colSums((s8 - s9) * w)
    
    # Free up memory before the next year
    rm(s8, s9)
    gc(verbose = FALSE)
    
    # Print completion time for this year
    yr_elapsed <- round(difftime(Sys.time(), yr_start, units = "secs"), 1)
    message(sprintf("   -> Done (%s sec) | k_est = %.4f", yr_elapsed, mean(yr_draws)))
    
    data.frame(
      year  = yr,
      k_est = mean(yr_draws),
      k_se  = sd(yr_draws),
      k_lwr = quantile(yr_draws, 0.025),
      k_upr = quantile(yr_draws, 0.975)
    )
  }) |> dplyr::bind_rows()
  return(list(
    k = k_summary,
    fit_model = fit
  ))
}

k_index_FC<- run_growth("F", "Coastwide",l_inf)


#k_summary_north <- k_summary%>%mutate(region="N")
#k_summary_south <- k_summary_south%>%mutate(region="S",sex="F")
 
 

# 1. Base map projected to WGS84 / UTM Zone 10N (meters)
coast <- ne_states(returnclass = "sf", country = "united states of america") |>
  st_crop(c(xmin = -145, ymin = 25, xmax = -117, ymax = 48.5)) |>
  st_transform(crs = 32610)

# 2. Extract main 'age' fixed effect
age_coef <- tidy(k_index_FC[["fit_model"]], effects = "fixed") |>
  filter(term == "age") |>
  pull(estimate)

# 3. Prepare prediction grid for reference year and add UTM coordinates
ref_year <- 2003

pred_grid_sub <- predictions |>
  ungroup() |>
  filter(year == ref_year) |>
  mutate(
    age = 8,
    doy_scaled = 0,
    doy_scaled2 = 0,
    fyear = factor(year)
  )

# 4. Predict from the fitted model
preds <- predict(fit, newdata = pred_grid_sub)

# 5. Reshape and calculate k_f, spatial mean, and annual anomalies
k_df <- preds |>
  pivot_longer(
    cols = starts_with("zeta_s_age:fyear"),
    names_to = "zeta_s_fyear",
    names_prefix = "zeta_s_age:fyear",
    values_to = "zeta_s_vals"
  ) |>
  mutate(
    year_anom = as.integer(zeta_s_fyear),
    zeta_vals_age_adjusted = zeta_s_vals + age_coef,
    k_f = -zeta_vals_age_adjusted
  ) |>
  group_by(X, Y) |>
  mutate(
    mean_k_f = mean(k_f, na.rm = TRUE),
    sd_k_f = sd(k_f, na.rm = TRUE),
    anom_k_f = k_f - mean_k_f
  ) |>
  ungroup()

# Max absolute anomaly to keep diverging color scale centered at 0
max_anom <- max(abs(k_df$anom_k_f), na.rm = TRUE)

# 6. Plot spatial anomalies by year
ggplot(coast) +
  geom_tile(
    data = k_df,
    aes(x = X * 1000, y = Y * 1000, fill = anom_k_f),
    width = 4000,
    height = 4000
  ) +
  geom_sf(fill = "grey85", color = "grey40", linewidth = 0.25) +
  scale_fill_distiller(
    palette = "RdBu",
    direction = 1,
    limits = c(-max_anom, max_anom),
    name = expression(Delta * k[f])
  ) +
  facet_wrap(~ year_anom) +
  coord_sf(
    xlim = range(k_df$X * 1000, na.rm = TRUE),
    ylim = range(k_df$Y * 1000, na.rm = TRUE),
    expand = FALSE
  ) +
  theme_bw() +
  xlab("") +
  ylab("") +
  labs(
    title = "Lingcod (Female) - Annual Spatial Growth Anomalies",
    subtitle = expression("Anomaly relative to site-level mean (" * Delta * k[f] == k[f] - bar(k)[f] * ")")
  ) +
  theme(
    strip.background = element_rect(fill = "white"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )


# 6. Plot spatial anomalies by year

ggplot(coast) +
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
