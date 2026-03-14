# US Recession Forecasting with AR(1) Probit Model

This project implements an autoregressive probit model to forecast US recessions using yield spread and the National Financial Conditions Index (NFCI) as predictors. The model estimates parameters via maximum composite likelihood (MCL) and generates probabilistic recession forecasts at various time horizons.

The methodology is based on the paper by Poulin-Moore, Antoine, and Kerem Tuzcuoglu. "Forecasting recessions in Canada: an autoregressive probit model approach." Staff Working Paper No. 2024-10. Bank of Canada, 2024.

## Structure

- `analysis.qmd`: Main analysis document with data preparation and model estimation
- `R/functions.R`: Core MCL estimation functions (gradient, Hessian, forecasting)
- `data/`: Quarterly macroeconomic data for the US
