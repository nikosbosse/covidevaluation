# Packages ----------------------------------------------------------------

require(EpiSoon)
require(EpiNow)
require(dplyr)
require(bsts)
require(fable)
require(fabletools)
require(feasts)
require(furrr)
require(future)
require(future.apply)
require(urca)
require(dplyr)
require(ggplot2)
require(patchwork)
require(cowplot)
require(EpiSoon)

# Get timeseries ----------------------------------------------------------

data_samples <- 100

dir <- "covid-global/national" # gl
dir <- "../covid-regional/united-states/regional" # regional in US

## Extract the Rt and case timeseries as produced by EpiNow
timeseries <- EpiNow::get_timeseries(dir)


## Extract rt timeseries and format
rt_timeseries <- timeseries$rt %>% 
  dplyr::filter(rt_type %in% "nowcast", 
                type %in% "nowcast") %>% 
  # dplyr::filter(region %in% "New York") %>%
  dplyr::mutate(sample = as.numeric(sample)) %>% 
  dplyr::group_by(region, date, sample) %>% 
  dplyr::mutate(rt_sample = 1:dplyr::n()) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(timeseries = region, rt = R, date, rt_sample, sample) %>% 
  dplyr::group_by(timeseries, date, sample) %>% 
  ## Each cases dataset has multiple Rt samples linked to it
  ## Take mean here to prevent duplicates
  dplyr::summarise(rt = mean(rt, na.rm = TRUE)) %>% 
  dplyr::ungroup() 

## Extract cases timeseries and format
case_timeseries <- timeseries$incidence %>% 
  dplyr::filter(import_status %in% "local") %>% 
#  dplyr::filter(region %in% "New York") %>%
  dplyr::mutate(sample = as.numeric(sample)) %>% 
  dplyr::select(timeseries = region, cases, date, sample)


## Sample timeseries
samples <- sample(1:max(rt_timeseries$sample), data_samples)

rt_timeseries <- rt_timeseries %>% 
  dplyr::filter(sample %in% samples)

case_timeseries <- case_timeseries %>% 
  dplyr::filter(sample %in% samples)

## Save samples
saveRDS(rt_timeseries, "forecast/model-choice/data/rt_timeseries.rds")

saveRDS(case_timeseries, "forecast/model-choice/data/case_timeseries.rds")

# Define a serial interval ------------------------------------------------

## In EpiNow we use a sampled serial interval - this is not yet supported 
## by EpiSoon and so instead here we define the mean serial interval
covid_serial_interval <- rowMeans(EpiNow::covid_serial_intervals)

# Define models for evaluation --------------------------------------------

## Inverse variance weighted
baseline_ensemble <- function(...) {
  EpiSoon::fable_model(model = fabletools::combination_model(fable::ARIMA(y), fable::ETS(y), fable::NAIVE(y),
                                                             fable::RW(y ~ drift()), cmbn_args = list(weights = "inv_var")), ...)
}

## Mean ensemble
baseline_mean_ensemble <- function(...) {
  EpiSoon::fable_model(model = fabletools::combination_model(fable::ARIMA(y), fable::ETS(y), fable::NAIVE(y),
                                                             fable::RW(y ~ drift()),
                                                             cmbn_args = list(weights = "equal")), ...)
}

## Last 14 days of data
current_ensemble <- function(y, ...) {
  y <- y[max(1, length(y) - 14):length(y)]
  EpiSoon::fable_model(y = y, model = fabletools::combination_model(fable::ARIMA(y), fable::ETS(y), 
                                                                    fable::NAIVE(y), fable::RW(y ~ drift()),
                                                                    cmbn_args = list(weights = "inv_var")), ...)
}

## Last 7 days of data
current_7_ensemble <- function(y, ...) {
  y <- y[max(1, length(y) - 7):length(y)]
  EpiSoon::fable_model(y = y, model = fabletools::combination_model(fable::ARIMA(y), fable::ETS(y), 
                                                                    fable::NAIVE(y), fable::RW(y ~ drift()),
                                                                    cmbn_args = list(weights = "inv_var")), ...)
}

## Inverse variance weighted
no_arima_ensemble <- function(...) {
  EpiSoon::fable_model(model = fabletools::combination_model(fable::ETS(y), fable::NAIVE(y),
                                                             fable::RW(y ~ drift()), 
                                                             cmbn_args = list(weights = "inv_var")), ...)
}

saveRDS(no_arima_ensemble, here::here("forecast/model-choice/data/ensemble.rds"))


models <- list( 
  "Naive" = function(...){EpiSoon::fable_model(model = fable::NAIVE(y), ...)},
  "ARIMA" = function(...){EpiSoon::fable_model(model = fable::ARIMA(y), ...)},
  "ETS" = function(...){EpiSoon::fable_model(model = fable::ETS(y), ...)},
  "Drift" = function(...){EpiSoon::fable_model(model = fable::RW(y ~ drift()), ...)},
  "Ensemble (mean)" = baseline_mean_ensemble,
  "Ensemble (weighted)" =  baseline_ensemble,
  "Ensemble (no ARIMA)" = no_arima_ensemble,
  "Ensemble (previous 14 days)" = current_ensemble,
  "Ensemble (previous 7 days)" = current_7_ensemble
)

# Set up parallelisation --------------------------------------------------
if (!interactive()) {
  ## If running as a script enable this
  options(future.fork.enable = TRUE)
}

future::plan("multiprocess")

# Evaluate models ---------------------------------------------------------


## Run iterative evaluation
forecast_eval <- EpiSoon::compare_timeseries(obs_rts = rt_timeseries,
                                             obs_cases = case_timeseries,
                                             models = models,
                                             horizon = 21,
                                             samples = 10,
                                             serial_interval = covid_serial_interval)


# Load in evaluation ------------------------------------------------------

## Forecasts
forecast_rts <- forecast_eval$forecast_rts  
saveRDS(forecast_rts, here::here("forecast/model-choice/data/forecast_rts.rds")) 

forecast_cases <- forecast_eval$forecast_cases  
saveRDS(forecast_cases, here::here("forecast/model-choice/data/forecast_cases.rds")) 

## Scores
rt_scores <- forecast_eval$rt_scores
case_scores  <- forecast_eval$case_scores 

## Summarised scores

scores_0_7 <- list(rt = rt_scores, case = case_scores) %>% 
  purrr::map_dfr( 
    ~ dplyr::filter(., horizon <= 7) %>% 
      summarise_scores("timeseries"), .id = "type")

scores_8_plus <- list(rt = rt_scores, case = case_scores) %>% 
  purrr::map_dfr( 
    ~ dplyr::filter(., horizon > 7) %>% 
      summarise_scores("timeseries"), .id = "type")

## Functions for plotting/summary
adjust_score <- function(df, group_var) {
  
  df %>% 
    dplyr::group_by(score, .dots = group_var) %>% 
    dplyr::mutate(upper_min = 10 * min(upper)) %>% 
    dplyr::ungroup() %>% 
    dplyr::filter(upper <= upper_min | score %in% c("bias", "calibration")) %>%
    dplyr::ungroup() %>% 
    dplyr::filter(!score %in% c("logs", "dss")) %>% 
    dplyr::mutate(score = score %>% 
                    factor(levels = c("crps", "calibration",
                                      "sharpness", "bias",
                                      "median", "iqr", "ci")) %>% 
                    dplyr::recode_factor(
                      "crps" = "CRPS",
                      "calibration" = "Calibration",
                      "sharpness" = "Sharpness",
                      "bias" = "Bias",
                      "median" = "Median",
                      "iqr" = "IQR",
                      "ci" = "CI"
                    ))
}



## Summarise scores
summarise_scores_by_horizon <- function(scores) {
  
  
  score_7 <- scores %>% 
    dplyr::filter(horizon <= 7) %>%
    EpiSoon::summarise_scores() %>% 
    dplyr::mutate(horizon = "0 -- 7")
  
  score_14 <- scores %>% 
    dplyr::filter(horizon <= 14, horizon > 7) %>%
    EpiSoon::summarise_scores() %>% 
    dplyr::mutate(horizon = "8 -- 14")
  
  # score_14_plus <- scores %>% 
  #   dplyr::filter(horizon > 14) %>%
  #   EpiSoon::summarise_scores() %>%
  #   dplyr::mutate(horizon = "14+")
  # 
  
  scores <- score_7 %>% 
    dplyr::bind_rows(score_14) %>% 
    # dplyr::bind_rows(score_14_plus) %>% 
    dplyr::mutate(horizon = horizon %>% 
                    factor(levels = c("0 -- 7", "8 -- 14")))
  
  return(scores)
}

## Summarise by horizon for all and recent data only
summarise_scores_recent_all <- function(scores, recent = NULL) {
  all_scores <- scores %>% 
    summarise_scores_by_horizon()
  
  recent_scores <- scores %>% 
    dplyr::filter(date >= max(date) - lubridate::days(recent)) %>% 
    summarise_scores_by_horizon()
  
  
  sum_scores <- all_scores %>% 
    dplyr::mutate(data = "All") %>% 
    dplyr::bind_rows(recent_scores %>% 
                       dplyr::mutate(data = paste0("Most recent ", recent, " days")))
  
  return(sum_scores)
}

sum_scores <- summarise_scores_recent_all(rt_scores, recent = 7) %>% 
  dplyr::mutate(type = "Rt") %>% 
  dplyr::bind_rows(summarise_scores_recent_all(case_scores, recent = 7) %>% 
                     dplyr::mutate(type = "Cases")) %>% 
  adjust_score(group_var = c("horizon", "type", "data")) %>% 
  dplyr::mutate(type = type %>% 
                  factor(levels = c("Rt", "Cases")))


saveRDS(sum_scores, here::here("forecast/model-choice/data/summarised_scores.rds"))



# Generate summary plots --------------------------------------------------


summary_plot <- function(scores, target_score) {
  scores %>% 
    dplyr::filter(score %in% target_score) %>% 
    ggplot2::ggplot(
      ggplot2::aes(x = model, y = mean, col = data)) +
    ggplot2::geom_point(size = 2, position = ggplot2::position_dodge(width = 1)) +
    ggplot2::geom_point(ggplot2::aes(y = median), shape = 2, size = 2,
                        position = ggplot2::position_dodge(width = 1)) +
    ggplot2::geom_linerange(ggplot2::aes(ymin = lower, ymax = upper),
                            alpha = 0.4, size = 1.5,
                            position = ggplot2::position_dodge(width = 1)) +
    ggplot2::scale_fill_viridis_d(option = "cividis", begin = 0.2, end = 0.6) +
    cowplot::theme_cowplot() +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::coord_flip() +
    facet_grid(horizon ~ type, scales = "free_x") +
    cowplot::panel_border()  +
    labs(x = "Model", y = target_score, col = "Data")
  
  
}

bias_summary_plot <- summary_plot(sum_scores, "Bias")

ggsave(paste0("forecast/model-choice/figures/bias.png"), bias_summary_plot, 
       dpi = 330, width = 8, height = 6)


crps_summary_plot <- summary_plot(sum_scores, "CRPS")

ggsave(paste0("forecast/model-choice/figures/crps.png"), crps_summary_plot, 
       dpi = 330, width = 8, height = 6)

sharpness_summary_plot <- summary_plot(sum_scores, "Sharpness") 

ggsave(paste0("forecast/model-choice/figures/sharpness.png"), sharpness_summary_plot, 
       dpi = 330, width = 8, height = 6)

calibration_summary_plot <- summary_plot(sum_scores, "Calibration") 

ggsave(paste0("forecast/model-choice/figures/calibration.png"), calibration_summary_plot, 
       dpi = 330, width = 8, height = 6)


median_summary_plot <- summary_plot(sum_scores, "Median") 

ggsave(paste0("forecast/model-choice/figures/median.png"), median_summary_plot, 
       dpi = 330, width = 8, height = 6)

iqr_summary_plot <- summary_plot(sum_scores, "IQR") 

ggsave(paste0("forecast/model-choice/figures/iqr.png"), iqr_summary_plot, 
       dpi = 330, width = 8, height = 6)

ci_summary_plot <- summary_plot(sum_scores, "CI") 

ggsave(paste0("forecast/model-choice/figures/ci.png"), ci_summary_plot, 
       dpi = 330, width = 8, height = 6)

# Plot evaluation ---------------------------------------------------------

plot_case_rt_comparison <- function(horizon = NULL) {
  
  ## Plot look ahead forecasts of rts
  rt_plot <- EpiSoon::plot_forecast_evaluation(forecast_rts, rt_timeseries,
                                               horizon_to_plot = horizon) +
    ggplot2::facet_grid(model ~ timeseries) +
    cowplot::panel_border()
  
  
  ggsave(paste0("forecast/model-choice/figures/rt_plot_", horizon, ".png"), rt_plot, 
         dpi = 330, width = 36, height = 24)
  
  
  ## Plot look ahead forecasts of cases
  case_plot <- EpiSoon::plot_forecast_evaluation(forecast_cases, case_timeseries,
                                                 horizon_to_plot = horizon) +
    ggplot2::facet_grid(model ~ timeseries, scales = "free_y") +
    cowplot::panel_border()
  
  ggsave(paste0("forecast/model-choice/figures/case_plot_", horizon, ".png"), case_plot, 
         dpi = 330, width = 36, height = 24)
  
  return(invisible(NULL))
}

plot_case_rt_comparison(horizon = 7)

plot_case_rt_comparison(horizon = 14)

# Summarise forecast scores -----------------------------------------------

make_plot_horizon_score <- function(df, label = NULL) {
  df <- df %>% 
    dplyr::filter(horizon %% 3 == 0, horizon <= 15)
  
  plot_internal <- function(df, label = NULL) {
    plot <- df %>% 
      ggplot2::ggplot(ggplot2::aes(x = horizon, y = mean, col = model,
                                   group = model)) +
      ggplot2::geom_line(size = 1.2, alpha = 0.6) +
      ggplot2::geom_point(size = 2) + 
      ggplot2::geom_point(ggplot2::aes(y = median), shape = 2, size = 2) +
      ggplot2::geom_linerange(ggplot2::aes(ymin = lower, ymax = upper),
                              alpha = 0.4, size = 1.5, position = position_dodge(width = 3)) +
      ggplot2::scale_color_viridis_d(option = "cividis") +
      ggplot2::scale_x_continuous(breaks = unique(df$horizon)) +
      cowplot::theme_cowplot() +
      ggplot2::theme(legend.position = "bottom") +
      labs(col = "Model", y = "", x = "", tag = label) +
      ggplot2::facet_wrap(~ score, scales = "free_y")
    
    return(plot)
  }
  
  score_plot <- df %>% 
    dplyr::filter(!score %in% "Bias") %>% 
    plot_internal(label = label) +
    ggplot2::scale_y_log10()
  
  bias_plot <- df %>% 
    dplyr::filter(score %in% "Bias") %>% 
    plot_internal(label = "") +
    ggplot2::theme(legend.position = "none")
  
  plot <- score_plot +
    bias_plot + 
    patchwork::plot_layout(widths = c(3, 1), nrow = 1)
  
  return(plot)
}

## Summarise scores across time horizons
plot_rt_horizon_score <- 
  summarise_scores(rt_scores, "horizon") %>% 
  adjust_score("horizon") %>% 
  make_plot_horizon_score(label = "Rt") &
  ggplot2::theme(legend.position = "none")

plot_cases_horizon_score <- 
  summarise_scores(case_scores, "horizon") %>% 
  adjust_score("horizon") %>%  
  make_plot_horizon_score(label = "Cases")

plot_horizon_score <- patchwork::wrap_plots(plot_rt_horizon_score,
                                            plot_cases_horizon_score) + 
  patchwork::plot_layout(ncol = 1)

ggsave("forecast/model-choice/figures/horizon_score.png", plot_horizon_score, 
       dpi = 330, width = 18, height = 18)

## Summarise scores across regions

make_region_score_plot <- function(df, title = NULL) {
  
  plot_region_score <- function(scores, label = NULL) {
    scores %>% 
      adjust_score("timeseries") %>%     
      ggplot2::ggplot(ggplot2::aes(x = timeseries, y = median, col = model)) +
      ggplot2::geom_point(size = 1.1, position = position_dodge(width = 1)) + 
      ggplot2::geom_point(ggplot2::aes(y = median), shape = 2, size = 1.1,
                          position = position_dodge(width = 1)) +
      ggplot2::geom_linerange(ggplot2::aes(ymin = lower, ymax = upper),
                              alpha = 0.4, size = 1.1, position = position_dodge(width = 1)) +
      ggplot2::scale_color_viridis_d(option = "cividis") +
      cowplot::theme_cowplot() +
      ggplot2::theme(legend.position = "bottom") +
      labs(col = "Model", y = "Score value", x = "Region", tag = label) +
      ggplot2::facet_wrap(~ score, scales = "free_x") +
      ggplot2::coord_flip()
    
  }
  
  plot_region_score_rt <- df %>% 
    dplyr::filter(type %in% "rt") %>% 
    plot_region_score(label = "Rt") +
    ggplot2::labs(title = title, x = "", y = "") 
  
  
  plot_region_score_cases <-  df %>%
    dplyr::filter(type %in% "case") %>% 
    plot_region_score(label = "Cases") +
    ggplot2::theme(legend.position = "none")
  
  plot_region_score <- plot_region_score_rt + 
    plot_region_score_cases +
    patchwork::plot_layout(nrow = 1)
  
  return(plot_region_score)
}

plot_regions_7 <- scores_0_7 %>%
  dplyr::filter(!score %in% "ci") %>% 
  make_region_score_plot(title = "Horizon: 0 to 7") 

plot_regions_7_up <- scores_8_plus %>%
  dplyr::filter(!score %in% "ci") %>% 
  make_region_score_plot(title = "Horizon: 7+") &
  ggplot2::labs(x = "", y = "")


ggsave("forecast/model-choice/figures/region_score_7.png", plot_regions_7, 
       dpi = 330, width = 18, height = 12)


ggsave("forecast/model-choice/figures/region_score_8_plus.png", plot_regions_7_up, 
       dpi = 330, width = 18, height = 12)