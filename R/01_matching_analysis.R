# Research: Evaluating the impact of private land conservation
# with synthetic control design
# Author: Roshan Sharma
# Published: Conservation Biology
# Date: 04 July 2023

# This script contains the code for the matching analysis

# load packages
packages <- c("MatchIt", "data.table", "tidyverse", "furrr", "future")

lapply(packages, require, character.only = TRUE)

# get the path with the current script and move one level up to
# get to the project directory with the all files and folders
# if it does not work, set the path manually
pd <- dirname(rstudioapi::getActiveDocumentContext()$path) |> dirname()

# set global options
options(scipen = 999)
debuggingState(on = FALSE)

# set data and results path
data_path <- file.path(pd, "data")
results_path <- file.path(pd, "results")

# load data
df_raw <- data.table::fread(file.path(data_path, "df_impact2.csv"))

# scale logit transformation of the outcome data
lower_limit <- 0
upper_limit <- 100.01  #add small value to avoid Inf

# refer for the formula for transformation
# https://robjhyndman.com/hyndsight/forecasting-within-limits/

# function to transform
scaled_log_transform <- function(x) {
  (log((x - lower_limit) / (upper_limit - x)))
}

# function to back-transform
back_transform <- function(x) {
  (upper_limit - lower_limit) * exp(x) / (1 + exp(x)) + lower_limit
}

# log transform outcome - proportion of woody vegetation
# variables starting with w_ are the outcome variables
df <- df_raw |>
  mutate_at(vars(starts_with("w_")), scaled_log_transform)

# set parameters for matching analysis
caliper_size <- c(0.6, 0.8, 1.0, 1.2) # caliper size
maxcontrols <- 3000 # no of controls to match

# important note: set a high number of controls, this ensures selection of as
# many controls as possible within this threshold - ideally, we wanted max no
# of controls. However, there the control pool is large
# in this data (n ~ 45,000) due to which there is software implementation issue
# so I tested different max values to check feasibility
# by setting maximum controls to 3000. No covenants had
# max controls greater than 2900, so setting 3000 will still
# result in the max number of controls that could be matched

# covariates for matching
# biophysical covariates
biophysical_covs <- data.frame(
  cov = c("slope", "elev", "soc", "awc", "der",
          "riv", "vegden", "road", "area_ha"),
  year = 0)

# woody vegetation
# note: need the year data to filter data for different years
# as each covenants have different pre-intervention data

# select names starting with "w_"
woody_covs <- colnames(df) |>
  str_subset("w_") |>
  data.frame() |>
  set_names("cov") |>
  mutate(year = str_extract(cov, "\\d+")) |>
  arrange(year)

# join the biophysical and woody covariates
covs <- rbind(biophysical_covs, woody_covs)

# pull the ID of covenants for looping
covenant_id <- stringr::str_subset(df$id, "^C")

# create a df to get combinations of covenants id and
# and caliper size for looping
df_loop <- crossing(covenant_id, caliper_size) |>
  mutate(sn = parse_number(covenant_id)) |>
  arrange(sn) |>
  mutate(list_name = paste0(covenant_id, "_", caliper_size))

head(df_loop)

#' a function to run matching analysis for each covenant
#' @covenants_id a string of covenant ID
#' @caliper_size is a numeric value of caliper size
#' @m_output a dataframe of matched treatment and control groups
do_matching <- function(covenant_id, caliper_size) {
  id <- treat <- year <- NULL # binding variables/columns
  # prepare data for matching
  df_prematch <- df |> filter(id == covenant_id | treat == 0)
  # get covenanted year
  yr_protect <- df_prematch |> filter(id == covenant_id) |>
    dplyr::pull(yr_protect)
  # select the covariates but filter only preintervention woody names
  covs_select <- covs |> filter(year < yr_protect) |> dplyr::pull(cov)
  # set calipers on the selected covariates
  set_calipers <- rep(caliper_size, length(covs_select)) |>
    purrr::set_names(covs_select)
  # create a matching formula
  match_formula <- as.formula(
    paste("treat~",
          paste(paste(covs_select, collapse = "+"))))
  # run matching
  match_out <- MatchIt::matchit(
    data = df_prematch,
    formula = match_formula,
    method = "nearest", distance = "mahalanobis",
    caliper = set_calipers,
    ratio = maxcontrols)
  # get matched dataset from matched object
  matched_dataset <- MatchIt::match.data(match_out)
  # assign a match_id and caliper size
  matched_dataset <- matched_dataset |>
    dplyr::mutate(
      match_id = covenant_id,
      caliper_size = caliper_size)
  }

# use parallelization and all available cores
future::plan(multisession, workers = availableCores())

# apply the function to the list of covenants and caliper size and
# assign the covenant_id as the list name
matched_list <- list(df_loop$covenant_id, df_loop$caliper_size) |>
  furrr::future_pmap(possibly(do_matching,
    print("No matches found"),
    quiet = TRUE),
    .options = furrr_options(seed = TRUE)) |>
  set_names(df_loop$list_name)

names(matched_list)

# save to rds - incase the process crashes
saveRDS(matched_list, file.path(results_path, "matched_list.rds"))
