# =============================================================================
# Retail Sales Forecasting: OLS (HAC) vs. SARIMA vs. Random Forest
# Data: U.S. Census Bureau MRTS via FRED (Federal Reserve Bank of St. Louis)
# Author: Kristina Zogovic, University of Belgrade
# =============================================================================
# FRED series:
#   MRTSSM44000USS  Retail Sales: Retail Trade (NSA, monthly, millions USD)
#   CPIAUCSL        CPI All Urban Consumers (SA, 1982-84 = 100)
#   UNRATE          Unemployment Rate (SA, %)
#   DCOILWTICO      WTI Crude Oil Price (USD/barrel)
# =============================================================================

# -- 0. Packages ---------------------------------------------------------------
required_packages <- c(
  "fredr",        # FRED API
  "tidyverse",    # data + ggplot2
  "lubridate",    # date utilities
  "randomForest", # Random Forest
  "strucchange",  # Chow test, Bai-Perron, recursive CUSUM
  "sandwich",     # HAC covariance
  "lmtest",       # coeftest
  "forecast",     # auto.arima, forecast, dm.test
  "Metrics",      # rmse(), mae()
  "scales"        # axis labels
)

new_pkgs <- setdiff(required_packages, rownames(installed.packages()))
if (length(new_pkgs)) install.packages(new_pkgs, dependencies = TRUE)
suppressPackageStartupMessages(
  invisible(lapply(required_packages, library, character.only = TRUE))
)

# -- 1. FRED download ----------------------------------------------------------
# Obtain a free API key at https://fred.stlouisfed.org/docs/api/api_key.html
# then run:  Sys.setenv(FRED_API_KEY = "your_key_here")
# OR paste your key directly in place of Sys.getenv("FRED_API_KEY") below.

FRED_KEY <- "dd3c25f5fe826af4918778cbbf2ed37a"
fredr_set_key(FRED_KEY)

pull <- function(id, start = "2009-01-01", end = "2023-12-01") {
  fredr(series_id = id,
        observation_start = as.Date(start),
        observation_end   = as.Date(end),
        frequency = "m") |>
    dplyr::select(date, value) |>
    dplyr::rename(!!id := value)
}

cat("Downloading from FRED...\n")
retail <- pull("MRTSSM44000USS")
cpi    <- pull("CPIAUCSL")
unemp  <- pull("UNRATE")
oil    <- pull("DCOILWTICO")

raw <- retail |>
  left_join(cpi,   by = "date") |>
  left_join(unemp, by = "date") |>
  left_join(oil,   by = "date") |>
  arrange(date) |>
  drop_na()

cat(sprintf("Raw data: %d obs, %s to %s\n",
            nrow(raw), min(raw$date), max(raw$date)))

# -- 2. Deflation and log transforms ------------------------------------------
# Deflate nominal sales by CPI; rebase to January 2015
cpi_base <- raw$CPIAUCSL[raw$date == as.Date("2015-01-01")]
stopifnot(length(cpi_base) == 1)

raw <- raw |>
  mutate(
    real_sales = MRTSSM44000USS / CPIAUCSL * cpi_base,
    log_sales  = log(real_sales),
    log_cpi    = log(CPIAUCSL),
    log_oil    = log(DCOILWTICO)
  )

# -- 3. Lags and regime labels (no seasonal index yet — computed in-loop) ------
raw <- raw |>
  mutate(
    lag1  = lag(log_sales, 1),
    lag2  = lag(log_sales, 2),
    lag3  = lag(log_sales, 3),
    lag12 = lag(log_sales, 12),
    yoy   = (real_sales - lag(real_sales, 12)) / lag(real_sales, 12) * 100,
    covid_dummy = as.integer(date >= as.Date("2020-03-01") &
                               date <= as.Date("2021-12-01")),
    regime = case_when(
      date <= as.Date("2012-12-01") ~ "Post-GFC Recovery",
      date <= as.Date("2019-12-01") ~ "Stable Growth",
      date <= as.Date("2021-12-01") ~ "COVID-19 Disruption",
      TRUE                          ~ "Post-COVID Rebound"
    ) |> factor(levels = c("Post-GFC Recovery", "Stable Growth",
                            "COVID-19 Disruption", "Post-COVID Rebound"))
  ) |>
  drop_na(lag1, lag2, lag3, lag12)   # removes first 12 obs

cat(sprintf("After lag construction: %d obs\n", nrow(raw)))

# -- 4. Full-sample OLS + HAC (for reporting; NOT used in walk-forward) --------
#
# NOTE: The seasonal index in the full-sample model is computed on ALL data,
# which is acceptable here because this model is for descriptive reporting only.
# In the walk-forward evaluation the seasonal index is always constructed
# from the training window only, avoiding look-ahead bias.

raw_full <- raw |>
  mutate(
    seasonal_idx = ave(real_sales, month(date),
                       FUN = function(x) mean(x, na.rm = TRUE)) /
                   mean(real_sales, na.rm = TRUE)
  )

ols_formula <- log_sales ~ seasonal_idx + lag1 + lag2 + lag3 + lag12 +
  log_cpi + log_oil + UNRATE + covid_dummy + covid_dummy:log_cpi

ols_full  <- lm(ols_formula, data = raw_full)
hac_table <- coeftest(ols_full, vcov = NeweyWest(ols_full, prewhite = FALSE))

cat("\n=== Full-sample OLS with Newey-West HAC SEs ===\n")
print(hac_table)

# -- 5. Full-sample SARIMA (descriptive) ---------------------------------------
ts_full <- ts(raw_full$log_sales,
              start = c(year(min(raw_full$date)), month(min(raw_full$date))),
              frequency = 12)

sarima_full <- auto.arima(ts_full, seasonal = TRUE, ic = "bic",
                          stepwise = FALSE, approximation = FALSE)
cat("\n=== Full-sample auto-selected SARIMA ===\n")
print(sarima_full)

# -- 6. Full-sample Random Forest (for variable importance) --------------------
set.seed(42)

# Nine predictors entering the Random Forest
rf_features <- c("seasonal_idx", "lag1", "lag2", "lag3", "lag12",
                 "log_cpi", "log_oil", "UNRATE", "covid_dummy")

rf_full <- randomForest(
  log_sales ~ .,
  data       = raw_full[, c("log_sales", rf_features)] |> drop_na(),
  ntree      = 1000,
  mtry       = floor(sqrt(length(rf_features))),
  importance = TRUE
)

imp_df <- importance(rf_full) |>
  as.data.frame() |>
  tibble::rownames_to_column("var") |>
  arrange(desc(`%IncMSE`)) |>
  mutate(
    # Use pmax(0) so negative values (uninformative predictors) don't
    # distort percentage shares — standard practice in RF importance reporting
    imp_pos = pmax(`%IncMSE`, 0),
    pct = imp_pos / sum(imp_pos) * 100,
    cumulative = cumsum(pct),
    label = recode(var,
                   seasonal_idx = "Seasonal Index",
                   lag1 = "AR(1)", lag2 = "AR(2)", lag3 = "AR(3)",
                   lag12 = "AR(12)", log_cpi = "Log CPI",
                   log_oil = "Log Oil Price", UNRATE = "Unemployment",
                   covid_dummy = "COVID Dummy")
  )

cat("\n=== Variable importance (% Inc MSE) ===\n")
print(imp_df[, c("label", "pct", "cumulative")])

seasonal_pct <- imp_df$pct[imp_df$var == "seasonal_idx"]
ar_short_pct <- sum(imp_df$pct[imp_df$var %in% c("lag1","lag2","lag3")])
cat(sprintf("\nSeasonal index: %.1f%%   Short AR lags: %.1f%%   Ratio: %.1f:1\n",
            seasonal_pct, ar_short_pct, seasonal_pct / ar_short_pct))

# -- 7. Structural break analysis ----------------------------------------------
covid_row <- which(raw_full$date == as.Date("2020-03-01"))
chow      <- sctest(ols_formula, data = raw_full, type = "Chow",
                    point = covid_row)
cat(sprintf("\n=== Chow test at 2020-M3: F = %.2f,  p = %.4f ===\n",
            chow$statistic, chow$p.value))

bp <- breakpoints(log_sales ~ log_cpi, data = raw_full, breaks = 5)
cat("\n=== Bai-Perron breakpoints ===\n")
print(summary(bp))
cat("Break dates:\n")
print(raw_full$date[bp$breakpoints])

# Price elasticity split
pre_df  <- raw_full |> filter(date < as.Date("2020-03-01"))
cov_df  <- raw_full |> filter(covid_dummy == 1)
ols_pre <- lm(log_sales ~ log_cpi + lag1 + lag12 + UNRATE, pre_df)
ols_cov <- lm(log_sales ~ log_cpi + lag1 + lag12 + UNRATE, cov_df)

b_pre <- coef(ols_pre)["log_cpi"]
b_cov <- coef(ols_cov)["log_cpi"]
cat(sprintf(
  "\nPrice elasticity:  pre-COVID = %.3f   COVID = %.3f   change = %.1f%%\n",
  b_pre, b_cov, (b_cov - b_pre) / abs(b_pre) * 100
))

# -- 8. Rolling walk-forward evaluation ----------------------------------------
#
# Design: fixed 24-month training window, 1-step-ahead forecast.
# The seasonal index is constructed WITHIN each training window
# to prevent look-ahead bias (important: this differs from the
# full-sample seasonal_idx computed above in raw_full).
#
# Window slides forward one month at a time:
#   Origin 1:  train Jan 2010-Dec 2011, forecast Jan 2012
#   Origin 2:  train Feb 2010-Jan 2012, forecast Feb 2012
#   ...
#   Origin 144: train Jan 2022-Dec 2023 — last forecast is Dec 2023

# Walk-forward OLS uses a stable formula without the covid interaction term
# The full ols_formula with covid_dummy:log_cpi causes near-singular matrices
# in windows where covid_dummy has near-zero variance (Post-COVID period)
ols_wf_formula <- log_sales ~ seasonal_idx + lag1 + lag2 + lag3 + lag12 +
  log_cpi + log_oil + UNRATE + covid_dummy
df_base <- raw    # raw without the full-sample seasonal_idx
n       <- nrow(df_base)
origins <- (WINDOW + 1):n   # indices of forecast targets

out <- tibble(
  date       = df_base$date[origins],
  actual     = df_base$log_sales[origins],
  ols_pred   = NA_real_,
  sar_pred   = NA_real_,
  rf_pred    = NA_real_,
  chow_p     = NA_real_,   # recursive Chow test p-value at this origin
  regime     = df_base$regime[origins]
)

cat(sprintf("\nWalk-forward: %d origins, window = %d months\n",
            length(origins), WINDOW))

set.seed(42)   # for RF reproducibility across all iterations

for (i in seq_along(origins)) {

  tgt   <- origins[i]
  train <- df_base[(tgt - WINDOW):(tgt - 1), ]
  test  <- df_base[tgt, , drop = FALSE]

  # --- seasonal index from training window only (no look-ahead) ---
  month_means <- tapply(train$real_sales, month(train$date), mean, na.rm = TRUE)
  overall_mean <- mean(train$real_sales, na.rm = TRUE)
  si_lookup <- month_means / overall_mean   # named vector length 12

  train <- train |>
    mutate(seasonal_idx = as.numeric(si_lookup[as.character(month(date))]))
  test <- test |>
    mutate(seasonal_idx = as.numeric(si_lookup[as.character(month(test$date))]))

  # OLS (uses simpler walk-forward formula to avoid near-singular matrices
  # in windows where covid_dummy has near-zero variance)
  m_ols <- tryCatch(lm(ols_wf_formula, data = train), error = function(e) NULL)
  if (!is.null(m_ols)) {
    pred_ols <- tryCatch(predict(m_ols, newdata = test), error = function(e) NA_real_)
    # Cap extreme predictions (more than 3 SD from training mean)
    train_mean <- mean(train$log_sales, na.rm = TRUE)
    train_sd   <- sd(train$log_sales, na.rm = TRUE)
    if (!is.na(pred_ols) && abs(pred_ols - train_mean) < 3 * train_sd) {
      out$ols_pred[i] <- pred_ols
    }
    # Recursive Chow test
    ols_simple <- tryCatch(
      lm(log_sales ~ lag1 + lag12 + log_cpi + seasonal_idx, data = train),
      error = function(e) NULL
    )
    if (!is.null(ols_simple)) {
      chow_result <- tryCatch({
        efp_obj <- efp(log_sales ~ lag1 + lag12 + log_cpi + seasonal_idx,
                       data = train, type = "Rec-CUSUM")
        sctest(efp_obj)$p.value
      }, error = function(e) NA_real_)
      out$chow_p[i] <- chow_result
    }
  }

  # SARIMA
  ts_tr <- ts(train$log_sales,
              start = c(year(min(train$date)), month(min(train$date))),
              frequency = 12)
  m_sar <- tryCatch(auto.arima(ts_tr, seasonal = TRUE, ic = "bic",
                               stepwise = TRUE),
                    error = function(e) NULL)
  if (!is.null(m_sar))
    out$sar_pred[i] <- tryCatch(as.numeric(forecast(m_sar, h = 1)$mean),
                                error = function(e) NA_real_)

  # Random Forest (smaller forest for speed)
  rf_train <- train[, c("log_sales", rf_features)] |> drop_na()
  rf_test  <- test[,  c("log_sales", rf_features)] |> drop_na()
  m_rf <- tryCatch(
    randomForest(log_sales ~ ., data = rf_train, ntree = 500),
    error = function(e) NULL
  )
  if (!is.null(m_rf) && nrow(rf_test) > 0)
    out$rf_pred[i] <- tryCatch(predict(m_rf, newdata = rf_test),
                               error = function(e) NA_real_)

  if (i %% 24 == 0) cat(sprintf("  completed %d / %d origins\n",
                                 i, length(origins)))
}

out_clean <- out |> drop_na(ols_pred, sar_pred, rf_pred)
cat(sprintf("Usable evaluation origins: %d\n", nrow(out_clean)))

# -- 8a. Regime-switching strategy ----------------------------------------------
# Rule: if recursive Chow p < 0.10 at that origin, use RF; else use SARIMA
out_clean <- out_clean |>
  mutate(
    switch_pred = if_else(!is.na(chow_p) & chow_p < 0.10,
                          rf_pred, sar_pred)
  )

# -- 8b. Diebold-Mariano tests (RF vs SARIMA, by regime) ----------------------
#
# dm.test() from the forecast package.
# Harvey-Leybourne-Newbold correction applied via h=1, alternative="two.sided"
# We report the test statistic and p-value.

dm_test_regime <- function(df, regime_label) {
  sub <- df |> filter(regime == regime_label) |> drop_na(rf_pred, sar_pred, actual)
  if (nrow(sub) < 5) return(tibble(regime = regime_label, DM_stat = NA, DM_p = NA))
  e_rf  <- sub$actual - sub$rf_pred
  e_sar <- sub$actual - sub$sar_pred
  res <- tryCatch(
    dm.test(e_rf, e_sar, alternative = "two.sided", h = 1, power = 2),
    error = function(e) NULL
  )
  if (is.null(res)) return(tibble(regime = regime_label, DM_stat = NA, DM_p = NA))
  tibble(regime = regime_label,
         DM_stat = round(res$statistic, 2),
         DM_p    = round(res$p.value, 3))
}

regime_levels <- c("Post-GFC Recovery", "Stable Growth",
                   "COVID-19 Disruption", "Post-COVID Rebound")

dm_results <- map_dfr(regime_levels, ~dm_test_regime(out_clean, .x))
cat("\n=== Diebold-Mariano tests: RF vs SARIMA (two-sided, squared loss) ===\n")
print(dm_results)

# Overall performance (now includes Switch)
perf_overall <- tibble(
  Model = c("OLS (HAC)", "SARIMA", "Random Forest", "Switch"),
  RMSE  = c(rmse(out_clean$actual, out_clean$ols_pred),
             rmse(out_clean$actual, out_clean$sar_pred),
             rmse(out_clean$actual, out_clean$rf_pred),
             rmse(out_clean$actual, out_clean$switch_pred)),
  MAE   = c(mae(out_clean$actual, out_clean$ols_pred),
             mae(out_clean$actual, out_clean$sar_pred),
             mae(out_clean$actual, out_clean$rf_pred),
             mae(out_clean$actual, out_clean$switch_pred))
)
cat("\n=== Overall walk-forward performance ===\n")
print(perf_overall)

sar_rmse  <- perf_overall$RMSE[perf_overall$Model == "SARIMA"]
sw_rmse   <- perf_overall$RMSE[perf_overall$Model == "Switch"]
pct_gain  <- round((sar_rmse - sw_rmse) / sar_rmse * 100, 1)
cat(sprintf("\nSwitch vs SARIMA RMSE improvement: %.1f%%\n", pct_gain))

# By regime (includes Switch)
perf_regime <- out_clean |>
  filter(!is.na(regime)) |>
  group_by(regime) |>
  summarise(
    n          = n(),
    OLS_RMSE   = rmse(actual, ols_pred),
    SAR_RMSE   = rmse(actual, sar_pred),
    RF_RMSE    = rmse(actual, rf_pred),
    SW_RMSE    = rmse(actual, switch_pred),
    OLS_MAE    = mae(actual, ols_pred),
    SAR_MAE    = mae(actual, sar_pred),
    RF_MAE     = mae(actual, rf_pred),
    SW_MAE     = mae(actual, switch_pred),
    .groups    = "drop"
  ) |>
  left_join(dm_results, by = "regime")

cat("\n=== Performance by regime (with Switch and DM tests) ===\n")
print(perf_regime)

# -- 9. Figures ----------------------------------------------------------------

# Fig 1: Real retail sales with regime shading
fig1 <- ggplot(raw_full, aes(x = date, y = real_sales / 1e6)) +
  annotate("rect", xmin = as.Date("2009-01-01"), xmax = as.Date("2012-12-01"),
           ymin = -Inf, ymax = Inf, fill = "#6B9EC7", alpha = 0.10) +
  annotate("rect", xmin = as.Date("2020-03-01"), xmax = as.Date("2021-12-01"),
           ymin = -Inf, ymax = Inf, fill = "#E8A838", alpha = 0.18) +
  geom_line(colour = "#1B3A6B", linewidth = 0.75) +
  geom_vline(xintercept = as.Date("2020-03-01"),
             linetype = "dashed", colour = "#C0392B", linewidth = 0.6) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  scale_y_continuous(labels = scales::dollar_format(suffix = "T", scale = 1)) +
  labs(
    title   = "U.S. Real Retail Sales, 2009-2023",
    subtitle = "Blue: Post-GFC Recovery  |  Gold: COVID-19 Disruption",
    x = NULL,
    y = "Real sales (trillions, 2015 USD)",
    caption = "Source: FRED MRTSSM44000USS deflated by CPIAUCSL (2015 base)"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title   = element_text(face = "bold", colour = "#1B3A6B"),
        plot.caption = element_text(size = 8, colour = "grey55"))

ggsave("fig1_real_sales.pdf", fig1, width = 9, height = 4.5)
cat("fig1 saved\n")

# Fig 2: Variable importance
fig2 <- imp_df |>
  mutate(label = fct_reorder(label, pct)) |>
  ggplot(aes(x = pct, y = label)) +
  geom_col(fill = "#1B3A6B", alpha = 0.85, width = 0.65) +
  geom_text(aes(label = sprintf("%.1f%%", pct)), hjust = -0.1, size = 3.3) +
  scale_x_continuous(
    limits = c(0, max(imp_df$pct) * 1.18),
    labels = scales::percent_format(scale = 1)
  ) +
  labs(
    title    = "Random Forest variable importance",
    subtitle = "% increase in out-of-bag MSE when predictor is permuted",
    x = "% Increase MSE", y = NULL,
    caption  = "B = 1,000 trees; mtry = 3; set.seed(42)"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", colour = "#1B3A6B"))

ggsave("fig2_variable_importance.pdf", fig2, width = 7, height = 5)
cat("fig2 saved\n")

# Fig 3: Regime RMSE bar chart
fig3_data <- perf_regime |>
  pivot_longer(c(OLS_RMSE, SAR_RMSE, RF_RMSE),
               names_to = "Model", values_to = "RMSE") |>
  mutate(Model = recode(Model,
                        OLS_RMSE = "OLS (HAC)",
                        SAR_RMSE = "SARIMA",
                        RF_RMSE  = "Random Forest"))

fig3 <- ggplot(fig3_data, aes(x = regime, y = RMSE, fill = Model)) +
  geom_col(position = "dodge", alpha = 0.87, width = 0.72) +
  scale_fill_manual(values = c("OLS (HAC)"     = "#1B3A6B",
                                "SARIMA"        = "#6B9EC7",
                                "Random Forest" = "#E8A838")) +
  labs(
    title    = "Walk-forward RMSE by economic regime",
    subtitle = "Fixed 24-month rolling window; one-step-ahead forecast",
    x = NULL, y = "RMSE (log real sales)", fill = NULL,
    caption  = "Lower = better"
  ) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x     = element_text(angle = 18, hjust = 1),
        plot.title      = element_text(face = "bold", colour = "#1B3A6B"),
        legend.position = "top")

ggsave("fig3_regime_rmse.pdf", fig3, width = 8, height = 5)
cat("fig3 saved\n")

# Fig 4: Rolling price elasticity
roll_window <- 24
roll_idx    <- seq(roll_window, nrow(raw_full) - 1)

roll_beta <- vapply(roll_idx, function(i) {
  sub <- raw_full[(i - roll_window + 1):i, ]
  m   <- tryCatch(lm(log_sales ~ log_cpi + lag1 + UNRATE, sub),
                  error = function(e) NULL)
  if (!is.null(m) && "log_cpi" %in% names(coef(m))) coef(m)["log_cpi"]
  else NA_real_
}, numeric(1))

roll_df <- tibble(date = raw_full$date[roll_idx], beta = roll_beta) |>
  drop_na()

fig4 <- ggplot(roll_df, aes(x = date, y = beta)) +
  annotate("rect", xmin = as.Date("2020-03-01"), xmax = as.Date("2021-12-01"),
           ymin = -Inf, ymax = Inf, fill = "#E8A838", alpha = 0.20) +
  geom_line(colour = "#1B3A6B", linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey60") +
  geom_vline(xintercept = as.Date("2020-03-01"),
             linetype = "dashed", colour = "#C0392B", linewidth = 0.5) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  labs(
    title    = "Rolling price elasticity of retail sales (24-month window)",
    subtitle = "Gold: COVID-19 Disruption  |  Dashed red: March 2020",
    x = NULL, y = "OLS coefficient on log(CPI)",
    caption  = "Rolling 24-month OLS window; dependent variable: log real retail sales"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", colour = "#1B3A6B"))

ggsave("fig4_rolling_elasticity.pdf", fig4, width = 9, height = 4.5)
cat("fig4 saved\n")

# Fig 5: Switching strategy — cumulative absolute error vs SARIMA
# Shows how much the Switch rule saves at each point in time
fig5_data <- out_clean |>
  arrange(date) |>
  mutate(
    abs_err_sar    = abs(actual - sar_pred),
    abs_err_switch = abs(actual - switch_pred),
    cum_saving     = cumsum(abs_err_sar - abs_err_switch),
    switch_active  = !is.na(chow_p) & chow_p < 0.10
  )

fig5 <- ggplot(fig5_data, aes(x = date, y = cum_saving)) +
  annotate("rect", xmin = as.Date("2020-03-01"), xmax = as.Date("2021-12-01"),
           ymin = -Inf, ymax = Inf, fill = "#E8A838", alpha = 0.18) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.5) +
  geom_line(colour = "#1B3A6B", linewidth = 0.9) +
  geom_point(data = filter(fig5_data, switch_active),
             aes(x = date, y = cum_saving),
             colour = "#C0392B", size = 1.8, alpha = 0.7) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  labs(
    title    = "Cumulative forecast error saving: Switch vs. SARIMA",
    subtitle = "Gold: COVID-19 Disruption  |  Red dots: Switch triggers (Chow p < 0.10)",
    x = NULL, y = "Cumulative MAE saving (log real sales)",
    caption  = "Positive values indicate Switch outperforms SARIMA. Walk-forward evaluation."
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", colour = "#1B3A6B"),
        plot.caption = element_text(size = 8, colour = "grey55"))

ggsave("fig5_switch_savings.pdf", fig5, width = 9, height = 4.5)
cat("fig5 saved\n")

# -- 10. Save processed data ---------------------------------------------------
write_csv(
  raw_full |> select(date, real_sales, log_sales, log_cpi, log_oil,
                     UNRATE, seasonal_idx, lag1, lag2, lag3, lag12,
                     covid_dummy, regime, yoy),
  "retail_data_processed.csv"
)

# Also save walk-forward predictions (useful for replication)
write_csv(
  out_clean |> select(date, actual, ols_pred, sar_pred, rf_pred,
                      switch_pred, chow_p, regime),
  "walkforward_predictions.csv"
)
cat("retail_data_processed.csv and walkforward_predictions.csv saved\n")

# -- 11. Session info ----------------------------------------------------------
cat("\n=== Session info ===\n")
print(sessionInfo())
cat("\n=== Done ===\n")
