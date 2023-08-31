# Research: Evaluating the impact of private land conservation
# with synthetic control design
# Author: Roshan Sharma
# Published: Conservation Biology
# Date: 04 July 2023

# This script contains the code for plotting the results

# load packages
packages <- c("tidyverse", "cowplot", "bayestestR")
lapply(packages, require, character.only = TRUE)

# project directory
pd <- dirname(rstudioapi::getActiveDocumentContext()$path) |> dirname()
results_path <- file.path(pd, "results")

# scale logit transformation of the outcome data
lower_limit <- 0
upper_limit <- 100.01  #add small value to avoid Inf

back_transform <- function(x) {
  (upper_limit - lower_limit) * exp(x) / (1 + exp(x)) + lower_limit
}

# load the panel data for plotting
matched_panel <- readRDS(file.path(results_path, "matched_panel.rds"))

# load the impact results
impact_list <- readRDS(file.path(results_path, "impact_list.rds"))

# create a custom theme
mytheme <- function() {
  theme_bw(base_size = 30) %+replace%
    theme(
      legend.position = "none",
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
}

# generate the plot for each covenant (e.g., figure 3 in manuscript)
# the plot is generated in three panels

make_panel_plots <- function(x) {
 panel <- matched_panel[[x]]
 df <- impact_list[[x]]

# vertical line to indicate the year of protection
 vline <- ggplot2::geom_vline(
    xintercept = df$yr_protect, colour = "darkgrey",
    size = 2, linetype = "dashed")

# panel a: plot the time series of the treatment and control groups
panel_a <- panel |>
    mutate(across(prop_woody, back_transform),
           type = as.factor(case_when(str_starts(id, "C") ~ 1, TRUE ~ 0))) |>
    mutate(id = fct_reorder(id, as.numeric(type))) |>
    ggplot() +
    geom_line(
        aes(x = year, y = prop_woody, col = type,
            group = id), linewidth = 2) +
    #scale_size_manual(values = c(1, 2)) +
    scale_color_manual(
        values = alpha(c("#97dbf2", "#000000"), c(0.8, 1))) +
    vline + ylim(0, 100) + labs(x = "", y = "") +
    mytheme() +
    theme(plot.margin = margin(2, 2, 0, 0, "cm"))

# panel b: plot the time series of the treatment and counterfactual
panel_b <- ggplot(df) +
    geom_ribbon(aes(
      x = year, ymin = cf_lower,
      ymax = cf_upper), fill = "#e8e4a6") +
    geom_line(aes(
        x = year, y = observed), linewidth = 2,
        colour = "#000000", na.rm = TRUE) +
    geom_line(aes(
        x = year, y = cf), size = 2, colour = "#0864a4",
        linetype = "dashed", na.rm = TRUE) + vline +
    ylim(0, 100) + labs(x = "", y = "") +
    mytheme() +
    theme(plot.margin = margin(2, 2, 0, 0, "cm"))

panel_c <- ggplot(df) +
    geom_hline(aes(
        yintercept = 0), colour = "darkgrey",
        linewidth = 2, linetype = "solid") +
    geom_ribbon(aes(
        x = year, ymin = effect_lower,
        ymax = effect_upper), fill = "#e8e4a6") +
    geom_line(aes(
        x = year, y = effect), linewidth = 2,
        colour = "#7f0f0f",
        linetype = "dashed", na.rm = TRUE) +
    vline + labs(x = "", y = "") + theme_bw(base_size = 30) +
    coord_cartesian(ylim = c(-100, 100))

panel_d <- cowplot::plot_grid(
    panel_a, panel_b, panel_c,
    nrow = 3, labels = "AUTO",
    label_size = 35, align = "v",
    rel_heights = c(1, 1, 1.3)
  )

# add the title
title <- cowplot::ggdraw() +
    cowplot::draw_label(
      str_c(
      "ID: ", str_extract(x, ".*(?=_)"), " ",
      "caliper:", str_extract(x, "(?<=_).*")),
      fontface = "bold", size = 35, x = 0, hjust = 0) +
    theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(0, 0, 0, 7))

plot_e <- cowplot::plot_grid(
    title, panel_d,
    ncol = 1,
    rel_heights = c(0.1, 1)
    )
}

# generate the plots for the selected covenants for the figure for the
# (manuscript figure 3)
selected_covs <- c("C79_0.6", "C80_0.6", "C99_0.6",
  "C47_0.6", "C65_0.6", "C138_0.6")

selected_covs_plots <- selected_covs |> map(make_panel_plots)

# arrange the plots in a grid
selected_covs_plots_grid <-  cowplot::plot_grid(
    plotlist = selected_covs_plots,
    ncol = 2, nrow = 3, align = "hv",
    rel_widths = c(1, 1), rel_heights = c(1, 1)) +
  theme(plot.margin = margin(3, 3, 3, 3, "cm"))

# save the plot to add the titles manually
ggsave(
  file.path(results_path, "plots", "fig_3_panel_plots.pdf"),
  selected_covs_plots_grid, width = 20, height = 26, dpi = 300,
  bg = "white"
)

# panel plots for all covenants (figures in supplementary information)
all_covs_plots <- map(names(impact_list), make_panel_plots)

# plot all plots in a single pdf
all_covs_plots_pdf <- marrangeGrob(
  all_covs_plots,
  layout_matrix = matrix(1:36, 4, 6, TRUE),
  as.table = FALSE)

ggsave(
  file.path(results_path, "plots", "si_panel_plots_all_covs.png"),
  selected_covs_plots_grid, width = 30, height = 40, dpi = 50,
  units = "in", bg = "white")

# plot 2 temporal distribution of the impact
# (manuscript figure 4)

# discard the list elements that are characters
impact_list_filter <- impact_list |>
    discard(~ is.character(.x)) |>
    map(~ .x |> drop_na(effect)) |>
    map(~ .x |> filter(year > yr_protect))

# flatten the list and assign the id as a new column using purrr
df_impact <- impact_list_filter |>
    map_dfr(~ .x, .id = "name") |>
    mutate(
      id = str_extract(name, ".*(?=_)"),
      caliper = as.numeric(str_extract(name, "(?<=_).*"))) |>
    group_by(id, caliper) |>
    mutate(seq = seq_along(year))

# plot the temporal distribution
n  <- df_impact |>
  group_by(caliper) |>
  drop_na(effect) |>
  summarise(
    n = n_distinct(id)) |>
  mutate(
    names = str_c("caliper: ",
     sprintf("%.1f", caliper), ", *n* = ", n)
  )

caliper_names <- n$names |>
  set_names(n$caliper)

# create temporal data
temporal_data <- df_impact |>
  filter(seq <= 15) |>
  group_modify(~ add_row(.)) |>
  mutate_all(~replace(., is.na(.), 0)) |>
  ungroup()

temporal_plot <- temporal_data |>
  ggplot() +
  geom_ribbon(
    aes(x = seq, group = id,
        ymin = effect_lower, ymax = effect_upper),
        fill = "gray", color = NA, alpha = 0.1) +
  geom_line(aes(x = seq, y = effect, group = id),
            alpha = 0.4, linewidth = 1, color = "#7E1900") +
  geom_line(data = temporal_data |>
              group_by(caliper, seq) |>
              summarise(
                estimate = mean(effect, na.rm = TRUE)),
            aes(x = seq, y = estimate), color = "gold", linewidth = 1) +
  geom_hline(yintercept = 0,
             linetype = "dashed",
             color = "black",
             linewidth = 1) +
  labs(x = "Years since intervention",
       y = "Impact (% woody vegetation)") +
  facet_wrap(~ caliper,
    labeller = as_labeller(caliper_names), nrow = 2) +
  scale_y_continuous(breaks = seq(-100, 100, by = 50)) +
  theme_minimal() +
  theme(
    text = element_text(size = 30),
    strip.text = element_markdown()) # ggtext for italic text)

# save the plot
ggsave(
  file.path(results_path, "plots", "fig_4_temporal_impact_plot.png"),
  temporal_plot, width = 16, height = 12, dpi = 300,
  bg = "white"
)

######## ************* plots for point estimates (figure 5) ********** ######
#############################################################################

# summary of the effects
df_summary <- impact_list |>
  discard(~ is.character(.x)) |>
  map(~ .x |> drop_na(effect)) |> # remove effects with NA
  purrr::list_rbind() |> # flatten the list
  dplyr::group_by(id, caliper) |>
  slice(1) |> # get the first row for each id and each caliper
  ungroup()

summary(df_summary$effect

# count the number of positive and negative effects (only statistically
# significant that is upper and lower bounds do not cross 0)
prop_positive <- df_summary |>
  group_by(caliper) |>
  summarise(
    n = n(),
    positive = sum(aepy_lower > 0), # lower to be above 0
    prop = round(positive / n * 100, 0)
  )

prop_negative <- df_summary |>
  group_by(caliper) |>
  summarise(
    n = n(),
    negative = sum(aepy_upper < 0), # upper to be below 0
    prop = round(negative / n * 100, 0)
  )

# set title names
stat_sig <- str_c("(*n** = ",
  prop_positive$positive, " +ve", " & ",
  prop_negative$negative, " -ve", ")") |>
  str_pad(40, "left", " ")

# title names
title_names <- n |>
  mutate(
    names = str_c(names, stat_sig)
  )

title_names <- title_names$names |> set_names(n$caliper)

# sort the effects based on the last panel with caliper 1.2

# get id levels from caliper 1.2
levs <- df_summary |>
  dplyr::filter(caliper == 1.2) |>
  dplyr::arrange(desc(aepy)) |>
  dplyr::mutate(id = factor(id, levels = id))

levs <- levels(levs$id)

# change the level of the summary dataframe
df_summary <- df_summary |>
  dplyr::arrange(desc(aepy)) |>
  dplyr::mutate(id = factor(id, levels = levs))

# plot the average effect per year (aepy)

effect_point_plot <- df_summary |>
  ggplot() +
  geom_errorbar(aes(
    x = id, ymin = aepy_lower, ymax = aepy_upper,
    col = case_when(
      aepy > 0 ~ "#009e73",
      TRUE ~ "#d55e00")),
    width = .1, size = 1, alpha = 0.3) +
  geom_point(aes(
    x = id, y = aepy,
    col = case_when(
      aepy > 0 ~ "#009e73",
      TRUE ~ "#d55e00")),
    size = 3.5) +
  scale_colour_identity() +
  facet_wrap(~ caliper,
    nrow = 4,
    labeller = as_labeller(title_names)) +
  geom_hline(
    yintercept = 0, linetype = "dashed", color = "black", size = 1) +
  labs(x = "", y = "Impact (% woody vegetation) per year") +
  ylim(-10, 10) +
  theme_minimal() +
  theme(
    text = element_text(size = 30),
    strip.text = element_markdown(), # ggtext for italic text
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

effect_point_plot

# density of the average effect per year

effect_density_plot <- df_summary |>
  ggplot(aes(x = aepy)) +
  geom_histogram(
    aes(y = after_stat(density)),
    binwidth = 1, col = "black", fill = "white") +
  geom_density(alpha = .5, fill = "gold", color = NA) +
  labs(x = " ", y = " ") +
  geom_vline(
    xintercept = 0, linetype = "dashed",
    color = "black", linewidth = 2) +
  labs(x = " ", y = " ") +
  facet_wrap(~ caliper, nrow = 4) +
  coord_flip() +
  theme_minimal() +
  theme(
    text = element_text(size = 20),
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_blank())

# merge the point plot and the density plot
merged_plot <- plot_grid(
  effect_point_plot, effect_density_plot, nrow = 1, rel_widths = c(6, 1))

# save the plot
ggsave(
  file.path(results_path, "plots", "fig_5_point_effect_plot.pdf"),
  merged_plot, width = 18, height = 12, dpi = 300,
  bg = "white"
)

######## ************* plots for meta-analysis results ************* ########
#############################################################################

# load data
meta_results <- readRDS(file.path(results_path, "output_meta.RDS"))

# summary of meta as lists
meta_summary <- meta_results |> map(~ summary(.x))

n <- length(meta_summary)
models <- c(0.6, 0.8, 1.0, 1.2) # models with different caliper sizes
variance <- c(1, 2, 4, 10)
model_list <- crossing(models, variance)

# extract the summary statistics
df_meta_summary <- model_list |>
  mutate(
    estimate = map_dbl(1:n, ~ meta_summary[[.x]] |>
      pluck("fixed") |>
      pluck("Estimate")),
    dist = map(1:n, ~ meta_results[[.x]] |>
                 as.data.frame() |>
                 pluck("b_Intercept") |>
                 as.list()),
    ci_low = map_dbl(1:n,
      ~ ci(unlist(dist[[.x]]), ci = 0.95, method = "HDI") |>
      pluck("CI_low")),
    ci_high = map_dbl(1:n,
      ~ ci(unlist(dist[[.x]]), ci = 0.95, method = "HDI") |>
      pluck("CI_high"))) |>
  mutate(across(where(is.numeric), round, 2))

# plots for meta-analysis

# color palette for density plot
den_color <- scico::scico(4, palette = "roma")

# label names
labels <- c("μ ~ N(0,1)", "μ ~ N(0,2)", "μ ~ N(0,4)", "μ ~ N(0,10)") |>
       set_names(variance)

# plot for the distribution of the effects across difference
# priors
effect_dist_sensitivity <- df_meta_summary |>
  unnest(cols = dist) |>
  mutate(dist = as.numeric(dist)) |>
  ggplot() +
  geom_density(
    aes(x = dist, fill = as.factor(models)), color = NA, alpha = .6) +
  scale_fill_manual(values = den_color,
    labels = c("0.6", "0.8", "1.0", "1.2")) +
  geom_vline(xintercept = 0, col = "red", linewidth = 1) +
  labs(title = "", x = "Impact per year", y = "Density",
       fill = "Models with calipers") +
  xlim(-1, 3) +
  theme_minimal(base_size = 20) +
  theme(legend.position = "top") +
  facet_wrap(~ variance, labeller = as_labeller(labels), nrow = 4)

# save the plot
ggsave(
  file.path(results_path, "plots", "sensitivity_density.png"),
  effect_dist_sensitivity, width = 8, height = 11, dpi = 300,
  units = "in", bg = "white"
)

# trace plot

meta_results_0_4 <- model_list |>
  mutate(
    meta_results = meta_results) |>
  dplyr::filter(variance == 4) |>
  pull(meta_results)

head(meta_results_0_4)

# get posterior in a list
posterior <- map(meta_results_0_4, as.array)

# plot the intervals
bayesplot::color_scheme_set("mix-blue-red")

get_traceplot <- function(x) {
  bayesplot::mcmc_trace(x, pars = c("b_Intercept"),
             facet_args = list(ncol = 1, strip.position = "left")) +
    ggplot2::labs(y = "intercept")
}

trace_plot <- map(posterior, get_traceplot)

# arrange in a grid
trace_plot_grid <- cowplot::plot_grid(
                 plotlist = trace_plot, nrow = 4,
                 labels = str_extract(n$names, "[^,]+"))

ggsave(
  file.path(results_path, "plots", "traceplot.png"),
  trace_plot_grid, width = 10, height = 10, dpi = 300,
  units = "in", bg = "white"
)

# plot the posterior predictive checks
ppc_plot <- map(meta_results_0_4, bayesplot::pp_check)

# arrange in grid
ppc_grid <- cowplot::plot_grid(
  plotlist = den_plots,
  labels = str_extract(n$names, "[^,]+"))

ggsave(
  file.path(results_path, "plots", "posterior_predictive_check.png"),
  ppc_grid, width = 14, height = 12, dpi = 300,
  units = "in", bg = "white"
)