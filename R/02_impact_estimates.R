# Research: Evaluating the impact of private land conservation
# with synthetic control design
# Author: Roshan Sharma
# Published: Conservation Biology
# Date: 04 July 2023

# This script contains the code for the impact analysis

# load packages
packages <- c("CausalImpact", "tidyverse", "furrr", "future")
# packages <- c("CausalImpact", "brms", "bayestestR",
#               "data.table", "tidyverse", "purrr", "cowplot", "scico",
#               "bayesplot", "gridExtra", "future")

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

# function to back-transform

# scale logit transformation of the outcome data
lower_limit <- 0
upper_limit <- 100.01  #add small value to avoid Inf

back_transform <- function(x) {
  (upper_limit - lower_limit) * exp(x) / (1 + exp(x)) + lower_limit
}

# load the matched results
matched_list <- readRDS(file.path(results_path, "matched_list.rds"))

# remove the elements with no matched control these are characters
# so remove them from the list
matched_list <- matched_list |> purrr::keep(~ !is.character(.x))

# the results in matched_list are in wide format convert to long format
# for plotting the time series

# check the first element of a list it should be a dataframe
# with matched treatment and control groups
first_element <- matched_list[2]
head(first_element)

# create a function to convert wide to long format
convert_to_long <- function(df) {
  id <- prop_woody <- yr_protect <- year <- NULL
  df |>
    tidyr::pivot_longer(
      cols = dplyr::contains("w_", ignore.case = FALSE),
      names_to = "year",
      values_to = "prop_woody") |>
    dplyr::mutate(year = readr::parse_number(year),
                  yr_protect = max(yr_protect)) |>
    dplyr::ungroup() |>
    dplyr::select(id, prop_woody, year, yr_protect, caliper_size)
    }

# apply the function to the lists to get panel dataframes
matched_panel <- matched_list |> map(convert_to_long)

# save the matched panel object as an RDS file
saveRDS(matched_panel, file.path(results_path, "matched_panel.rds"))

# create a function to get a dataframe for CausalImpact
convert_for_analysis <- function(df) {
  id <- prop_woody <- NULL
   df |>
    tidyr::pivot_wider(
      names_from = id,
      values_from = prop_woody) |>
    # make first column response variable and the rest controls
    dplyr::relocate(dplyr::starts_with(c("C", "P"), ignore.case = FALSE))|>
    as.data.frame()
}

# apply the function to the list
analysis_list <- matched_panel |> map(convert_for_analysis)

df <- analysis_list[[10]]

# create a function to get impact estimates
# and back transform the estimates
get_impact <- function(df) {
  # run causal impact analysis
  impact <- CausalImpact::CausalImpact(
    df |>
      dplyr::select(dplyr::starts_with(c("C", "P"), ignore.case = FALSE)),
    model.args = list(prior.level.sd = 0.01),
    pre.period = c(1, which(df$year == df$yr_protect) - 1),
    post.period = c(which(df$year == df$yr_protect), nrow(df)))
  # get series from the impact also apply back transformation
  series <- impact$series |>
    as.data.frame() |>
    # back transform the columns in the series
    dplyr::mutate(across(everything(), .fns = back_transform)) |>
    dplyr::transmute(
      id = colnames(df) |> stringr::str_subset("^C"),
      year = df$year,
      yr_protect = df$yr_protect,
      age = max(df$year) - max(df$yr_protect),
      matches = colnames(df) |> stringr::str_subset("^P") |> length(),
      caliper = df$caliper_size,
      observed = response,
      cf = point.pred, # cf means counterfactual
      cf_lower = point.pred.lower,
      cf_upper = point.pred.upper,
      effect = response - point.pred,
      effect_lower = response - point.pred.upper, # cf upper is effect lower
      effect_upper = response - point.pred.lower)
  # get average effect per year (aepy)
  series <- series |>
    dplyr::mutate(
      aepy = series$effect[series$year == max(series$year)] / series$age,
      aepy_lower =
        series$effect_lower[series$year == max(series$year)] / series$age,
      aepy_upper =
        series$effect_upper[series$year == max(series$year)] / series$age,
      aepy_sd = (aepy_upper - aepy_lower) / 3.92) |>
    dplyr::mutate(across(where(is.numeric), round, 2))
  # return the series
  series
  }

# set up parallel processing
future::plan(multisession, workers = availableCores())

# apply the function
impact_list <- analysis_list |>
  furrr::future_map(possibly(get_impact,
    print("no_series_to_estimate_impact"),
    quiet = TRUE),
    .options = furrr_options(seed = TRUE))

# save the impact object as an RDS file
saveRDS(impact_list, file.path(results_path, "impact_list.rds"))
