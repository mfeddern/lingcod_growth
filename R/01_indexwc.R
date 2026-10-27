library(indexwc)
library(dplyr)
library(purrr)
library(sdmTMB)
library(tibble)


lingcod<- pull_catch(common_name = "lingcod")
names(lingcod) <- tolower(names(lingcod))
# Load from the rda file
configuration_to_run <- configuration |>
  dplyr::filter(species == "lingcod", source == "NWFSC.Combo", family == "sdmTMB::delta_gamma()")
configuration_coastwide <-configuration_to_run[1,] #coastwide configuration

configuration_NS <-configuration_to_run[3,]


#### Coastwide ####

#pulling in the configuration
my_data <- pull_and_format_data(
  configuration_to_run = configuration_coastwide,  common_name = "lingcod",
  survey = "NWFSC.Combo" 
)

names(my_data)
my_data$family
my_data$formula
my_data$knots
my_data$share_range
my_data$anisotropy
list(my_data$spatiotemporal1, my_data$spatiotemporal2)

## Estimate the sdmTMB model

fit <- run_sdmtmb(
  data = my_data$data_filtered[[1]],
  family = my_data$family,
  formula = my_data$formula,
  n_knots = my_data$knots,
  share_range = my_data$share_range,
  anisotropy = my_data$anisotropy,
  spatiotemporal = list(my_data$spatiotemporal1, my_data$spatiotemporal2)
)
fit

diag <- diagnose(
  fit = fit
)
args(diagnose)

### generate the index
index <- calc_index_areas(
  data = fit$data,
  fit = fit,
  boundaries = "Coastwide"
)

save_index_outputs(
  fit = fit,
  diagnostics = diag,
  indices = index,
  overwrite = TRUE,
  dir = "output/coastwide"
)

#### N/S Medncino ####

#pulling in the configuration
my_data <- pull_and_format_data(
  configuration_to_run = configuration_coastwide 
)

my_data[["data_filtered"]][[1]] <- my_data[["data_filtered"]][[1]] %>%
  mutate(split_mendocino = ifelse(latitude > 40.1666667, "N", "S"))

names(my_data)
my_data$family
my_data$formula
my_data$knots
my_data$share_range
my_data$anisotropy
list(my_data$spatiotemporal1, my_data$spatiotemporal2)

## Estimate the sdmTMB model

fit <- run_sdmtmb(
  data = my_data$data_filtered[[1]],
  family = my_data$family,
  formula = my_data$formula,
  n_knots = my_data$knots,
  share_range = my_data$share_range,
  anisotropy = my_data$anisotropy,
  spatiotemporal = list(my_data$spatiotemporal1, my_data$spatiotemporal2)
)
fit

diag <- diagnose(
  fit = fit
)
args(diagnose)

### generate the index
index <- calc_index_areas(
  data = fit$data,
  fit = fit,
  boundaries = c("North of Cape Mendocino",
                 "South of Cape Mendocino")
)

save_index_outputs(
  fit = fit,
  diagnostics = diag,
  indices = index,
  dir = "output/NS"
)
