Forecast Evaluation
================

# Introduction

The purpose of this repo is to provide a quick overview over the
utilities used to generate predictions on epiforecasts.io/covid and to
provide tools to improve and evaluate these models. Forecast generation
and evaluation are somewhat separate for computational reasons:
Forecasts two weeks into the future are run daily for the website. Model
evaluation is much more costly as it entails iteratively fitting the
model(s) multiple times, each time adding one more data point. From this
we obtain a set of predicted values and true observations that can be
used for scoring.

### Infrastructure

The framework relies on the two packages EpiNow and EpiSoon. EpiNow
provides tools to do nowcasts of cases by date of infection or date of
symptom onset and for estimating the current effective Reproduction
Number. EpiSoon provides functionality to forecast the trajectory of
\(R_t\) values into the future and to generate case forecasts based on
the \(R_t\) forecasts using the renewal equation to link \(R_t\) back to
cases. EpiSoon provides wrappers for different classes of timeseries
models through R packages. Currently, bsts, fable and soon brms are
supported.

## Fitting

EpiNow:

  - input: reported cases
  - fit a delay distribution between symptom onset and reporting using
    exp or gamma distribution, evaluated using LOOIC.
  - use delay distribution to transform reported cases to cases on the
    day of symptom onset
  - number of cases on the day of symptom onset is upscaled to account
    for the fact that some might have symptoms but haven’t yet been
    reported. Upscaling is done using a negative binomial distribution.
  - use EpiEstim to estimate time-varying reproduction number on sympom
    onset cases. To account for uncertainty, the used mean and sd of the
    serial interval are sampled 1000 times, so that the Reproduction
    Number is estimated 1000 times.
  - This is repeated for different sliding moving average windows in the
    range of 1 to 7 days.
  - For every of the 1000 combinations of mean and sd of the serial
    interval, the optimal time window is chosen to yield the best one
    day ahead predicitions, implicitly resulting in a mixture of time
    windows.

EpiSoon:

  - forecast \(R_t\)-values using time series models provided by bsts,
    fable, and soon brms
  - convert \(R_t\)-forecasts back to case forecasts using a renewal
    equation
  - currently, every region is forecasted completely separately

# Forecasting

Forecasts and estimates can be obtained in two ways: cloning them from
the existing epiforecasts repositories or generating forecasts anew
using an R script.

### Existing Repos

The results displayed on epiforecasts.io/covid/ are saved in two repos
on github.com/epiforecasts. To access the global data, clone the
epiforecasts/covid-global repository. For access to regional breakdowns
in Germany, Italy, The UK or the US, clone the
epiforecasts/covid-regional repository (both repos are fairly large).

git clone <https://github.com/epiforecasts/covid-global>

git clone <https://github.com/epiforecasts/covid-regional>

### Creating New Forecasts

To create new forecasts with settings different from the ones we
currently use, one can ran the update\_nowcasts.R script. In the
covid-global repo it lies in the top directory,
covid-global/update\_nowcasts.R. In the regional repo, there is one
script per region in the regional folders, e.g. in
covid-regional/united-states/update\_nowcasts.R.

# Model Evaluation

Model evaluation can be performed with the evaluate\_models.R script in
this repo. As the forecasting is currently done on \(R_t\) values rather
than cases (and then transforming back to cases) the script needs the
\(R_t\) estimates and nowcasts created by EpiNow. The easiest way to use
it is to change `dir` in the script to the directory with the latest
\(R_t\) estimates. Then, a list of models can be provided. Those will be
iteratively fit to data and the forecasts can then be compared against
the observed values.
