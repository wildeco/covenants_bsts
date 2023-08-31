# Research: Evaluating the impact of private land conservation
# with synthetic control design
# Author: Roshan Sharma and Ascelin Gordon
# Published: Conservation Biology
# Date: 04 July 2023

# Objective: Run meta-analysis of the effects from individual covenants

# brms had issues with this CausalImpact package, unload this if it is loaded
# from previous run
devtools::unload("CausalImpact")

# load packages
packages <- c("tidyverse", "brms", "parallel", "bayestestR",
  "ggtext", "bayesplot")
lapply(packages, require, character.only = TRUE)

# project directory
pd <- dirname(rstudioapi::getActiveDocumentContext()$path) |> dirname()
results_path <- file.path(pd, "results")

# load the impact results
impact_list <- readRDS(file.path(results_path, "impact_list.rds"))

# filter the results with the caliper values
df_summary <- impact_list |>
  discard(~ is.character(.x)) |>
  map(~ .x |> drop_na(effect)) |> # remove effects with NA
  purrr::list_rbind() |> # flatten the list
  dplyr::group_by(id, caliper) |>
  slice(1) |> # get the first row for each id and each caliper
  ungroup()

# note: add a small value to SD as in meta-analysis the main
# effect is normalised with SD, if SD is 0 then the infinity
# result will crash the run
df_summary <- df_summary |>
  mutate(aepy_sd = case_when(aepy <= 0.001 ~ 0.01,
                     TRUE ~ aepy_sd))

# models with different caliper sizes
models <- c(0.6, 0.8, 1.0, 1.2)

# sensitivity for the effect priors (effect and variance)
# N(0,1),N(0,2),N(0,4),N(0,10)

# for convenience, we can loop over the variance
variance <- c(1, 2, 4, 10)

# looping list
model_list <- crossing(models, variance)

# note: to pass a variable in a function in stan requires defining
# a variable by stanvars

# function to run meta analysis
do_meta <- function(x, y) {
  stanvars <- stanvar(y, name = "y")
  brms::brm(aepy | se(aepy_sd) ~ 1 + (1 | id),
      data = df_summary |>
        filter(caliper == x),
      prior = c(prior(normal(0, y), class = Intercept),
                prior(cauchy(0, 0.5), class = sd)),
      stanvars = stanvars,
      iter = 5000,
      cores = parallel::detectCores(),
      seed = 2021)
}

# apply the function to the model list
meta_results <- pmap(list(
  x = model_list$models, y = model_list$variance),
  .f = do_meta)

# save the impact object as an RDS file
saveRDS(meta_results, file.path(results_path, "output_meta.RDS")) # no lint