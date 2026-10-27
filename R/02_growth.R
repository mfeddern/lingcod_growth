remotes::install_github("pfmc-assessments/indexwc")
library(nwfscSurvey)
library(sdmTMB)
bio_dat <- pull_bio(common_name = "lingcod")
names(bio_dat) <- tolower(names(bio_dat))
bio_dat$yday <- lubridate::yday(bio_dat$date)

#grabbing from the indexwc output
load("output/coastwide/lingcod/wcgbts/lat_32.1-40.2_years_2003-2025/delta_gamma/diagnostics/predictions.rdata")
head(predictions)

# filter out only females
bio_dat <- dplyr::filter(bio_dat, sex == "F", !is.na(length_cm))
# add y_lin -- 152 in table 1 Andrea paper
bio_dat$y_lin <- log(152 - bio_dat$length_cm)

haul_dat <- pull_haul()
names(haul_dat) <- tolower(names(haul_dat))

bio_dat <- sdmTMB::add_utm_columns(
  bio_dat,
  ll_names = c("longitude_dd", "latitude_dd")
) |>
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
grid <- indexwc::california_current_grid |>
  dplyr::filter(area_km2_WCGBTS == 4)
grid <- dplyr::mutate(grid, age = 8, doy_scaled = 0, doy_scaled2 = 0)
grid <- replicate_df(
  grid,
  time_name = "year",
  time_values = unique(bio_dat$year)
)

grid$fyear <- as.factor(grid$year)
grid <- dplyr::rename(grid, X = x, Y = y)

# here is where we'd fit the biomass model, and use the predicted biomass as the RHS of the vector argument
indx <- sdmTMB::get_weighted_average(
  fit,
  newdata = grid,
  vector =exp(predictions$est), #you need to replace this with index WC
  bias_correct = TRUE
)

