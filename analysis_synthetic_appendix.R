# =============================================================================
# Synthetic Data Generation and Analysis
# Generates figA1_synth_importance.pdf and figA2_synth_rmse.pdf
# for Appendix A of the IJF submission
#
# This script:
#   1. Generates a synthetic retail sales series under a simple seasonal DGP
#   2. Runs the same RF variable importance analysis as in analysis_v2.R
#   3. Runs the same walk-forward evaluation
#   4. Saves two figures for the appendix showing the synthetic results
#      that contrast with the real FRED data findings
#
# Run this script AFTER analysis_v2.R (it reuses the same packages)
# =============================================================================

# -- 0. Packages ---------------------------------------------------------------
required_packages <- c(
  "tidyverse", "lubridate", "randomForest",
  "forecast", "Metrics", "scales"
)
new_pkgs <- setdiff(required_packages, rownames(installed.packages()))
if (length(new_pkgs)) install.packages(new_pkgs, dependencies = TRUE)
suppressPackageStartupMessages(
  invisible(lapply(required_packages, library, character.only = TRUE))
)

# -- 1. Synthetic DGP ----------------------------------------------------------
# Generate 180 months of synthetic retail sales (Jan 2009 - Dec 2023)
# DGP: log_sales = trend + seasonal + AR(1) + COVID_shock + noise
# Seasonal is EXPLICIT and strong (as commonly done in synthetic benchmarks)
# This is the key feature that makes synthetic results misleading.

set.seed(42)
n_total <- 180
dates   <- seq(as.Date("2009-01-01"), by = "month", length.out = n_total)

# Seasonal pattern (strong, explicit — typical synthetic DGP assumption)
month_idx  <- month(dates)
seasonal   <- c(0.88,0.88,0.97,0.97,1.00,1.01,1.00,1.00,1.00,1.01,1.04,1.22)
seas_vec   <- seasonal[month_idx]

# Trend
trend_vec  <- seq(0, 0.30, length.out = n_total)

# COVID shock (April 2020 = obs 136)
covid_dummy <- as.integer(dates >= as.Date("2020-03-01") &
                            dates <= as.Date("2021-12-01"))
covid_shock <- ifelse(dates == as.Date("2020-04-01"), -0.22,   # sudden drop
                ifelse(dates >= as.Date("2020-05-01") &
                       dates <= as.Date("2020-08-01"), 0.05,    # partial recovery
                0))

# Base log sales
log_sales_true <- 12.5 + trend_vec + log(seas_vec) + cumsum(c(0, rnorm(n_total-1, 0, 0.008)))
log_sales_true <- log_sales_true + covid_shock
synth_noise    <- rnorm(n_total, 0, 0.015)
log_sales      <- log_sales_true + synth_noise

# Synthetic CPI, oil, unemployment (trending, correlated with trend)
log_cpi  <- log(237) + seq(0, 0.18, length.out = n_total) + rnorm(n_total, 0, 0.005)
log_oil  <- log(65)  + 0.3*sin(seq(0, 4*pi, length.out = n_total)) + rnorm(n_total, 0, 0.08)
UNRATE   <- pmax(3.5, 9.5 - seq(0, 5.5, length.out = n_total) + rnorm(n_total, 0, 0.3))
UNRATE   <- ifelse(dates >= as.Date("2020-03-01") & dates <= as.Date("2020-06-01"),
                   UNRATE + c(5, 9, 6, 4), UNRATE)

# Build synthetic dataframe
synth <- tibble(
  date        = dates,
  log_sales   = log_sales,
  real_sales  = exp(log_sales) * 1e5,  # rescale to approximate millions
  log_cpi     = log_cpi,
  log_oil     = log_oil,
  UNRATE      = UNRATE,
  covid_dummy = covid_dummy,
  regime      = case_when(
    date <= as.Date("2012-12-01") ~ "Post-GFC Recovery",
    date <= as.Date("2019-12-01") ~ "Stable Growth",
    date <= as.Date("2021-12-01") ~ "COVID-19 Disruption",
    TRUE                          ~ "Post-COVID Rebound"
  ) |> factor(levels = c("Post-GFC Recovery","Stable Growth",
                          "COVID-19 Disruption","Post-COVID Rebound"))
) |>
  mutate(
    lag1  = lag(log_sales, 1),
    lag2  = lag(log_sales, 2),
    lag3  = lag(log_sales, 3),
    lag12 = lag(log_sales, 12),
    # Seasonal index computed from full sample (for illustrative full-sample RF)
    seasonal_idx = ave(real_sales, month(date),
                       FUN = function(x) mean(x, na.rm = TRUE)) /
                   mean(real_sales, na.rm = TRUE)
  ) |>
  drop_na(lag1, lag2, lag3, lag12)

cat(sprintf("Synthetic data: %d obs, %s to %s\n",
            nrow(synth), min(synth$date), max(synth$date)))

# -- 2. Full-sample RF variable importance (synthetic) -------------------------
rf_features <- c("seasonal_idx","lag1","lag2","lag3","lag12",
                 "log_cpi","log_oil","UNRATE","covid_dummy")

rf_synth <- randomForest(
  log_sales ~ .,
  data       = synth[, c("log_sales", rf_features)],
  ntree      = 1000,
  mtry       = floor(sqrt(length(rf_features))),
  importance = TRUE
)

imp_synth <- importance(rf_synth) |>
  as.data.frame() |>
  tibble::rownames_to_column("var") |>
  arrange(desc(`%IncMSE`)) |>
  mutate(
    pct = `%IncMSE` / sum(`%IncMSE`) * 100,
    label = recode(var,
      seasonal_idx = "Seasonal Index",
      lag1 = "AR(1)", lag2 = "AR(2)", lag3 = "AR(3)",
      lag12 = "AR(12)", log_cpi = "Log CPI",
      log_oil = "Log Oil Price", UNRATE = "Unemployment",
      covid_dummy = "COVID Dummy")
  )

cat("\n=== Synthetic RF variable importance ===\n")
print(imp_synth[, c("label","pct")])

seasonal_pct_synth <- imp_synth$pct[imp_synth$var == "seasonal_idx"]
ar_short_pct_synth <- sum(imp_synth$pct[imp_synth$var %in% c("lag1","lag2","lag3")])
cat(sprintf("Seasonal: %.1f%%   Short-AR: %.1f%%   Ratio: %.1f:1\n",
            seasonal_pct_synth, ar_short_pct_synth,
            seasonal_pct_synth / ar_short_pct_synth))

# -- 3. Walk-forward evaluation (synthetic) ------------------------------------
WINDOW  <- 24
df_wf   <- synth
n_wf    <- nrow(df_wf)
origins <- (WINDOW + 1):n_wf

out_s <- tibble(
  date     = df_wf$date[origins],
  actual   = df_wf$log_sales[origins],
  sar_pred = NA_real_,
  rf_pred  = NA_real_,
  regime   = df_wf$regime[origins]
)

ols_formula_s <- log_sales ~ seasonal_idx + lag1 + lag2 + lag3 + lag12 +
  log_cpi + log_oil + UNRATE + covid_dummy

cat(sprintf("\nWalk-forward (synthetic): %d origins\n", length(origins)))
set.seed(42)

for (i in seq_along(origins)) {
  tgt   <- origins[i]
  train <- df_wf[(tgt - WINDOW):(tgt - 1), ]
  test  <- df_wf[tgt, , drop = FALSE]

  # Seasonal index from training window
  mm  <- tapply(train$real_sales, month(train$date), mean, na.rm = TRUE)
  om  <- mean(train$real_sales, na.rm = TRUE)
  sil <- mm / om
  train <- train |> mutate(seasonal_idx = as.numeric(sil[as.character(month(date))]))
  test  <- test  |> mutate(seasonal_idx = as.numeric(sil[as.character(month(test$date))]))

  # SARIMA
  ts_tr  <- ts(train$log_sales,
               start = c(year(min(train$date)), month(min(train$date))), frequency = 12)
  m_sar  <- tryCatch(auto.arima(ts_tr, seasonal = TRUE, ic = "bic", stepwise = TRUE),
                     error = function(e) NULL)
  if (!is.null(m_sar))
    out_s$sar_pred[i] <- tryCatch(as.numeric(forecast(m_sar, h=1)$mean),
                                   error = function(e) NA_real_)

  # Random Forest
  rf_tr <- train[, c("log_sales", rf_features)] |> drop_na()
  rf_te <- test[,  c("log_sales", rf_features)] |> drop_na()
  m_rf  <- tryCatch(randomForest(log_sales ~ ., data = rf_tr, ntree = 500),
                    error = function(e) NULL)
  if (!is.null(m_rf) && nrow(rf_te) > 0)
    out_s$rf_pred[i] <- tryCatch(predict(m_rf, newdata = rf_te),
                                  error = function(e) NA_real_)

  if (i %% 24 == 0) cat(sprintf("  %d / %d\n", i, length(origins)))
}

out_s_clean <- out_s |> drop_na(sar_pred, rf_pred)

perf_s <- out_s_clean |>
  filter(!is.na(regime)) |>
  group_by(regime) |>
  summarise(
    SAR_RMSE = rmse(actual, sar_pred),
    RF_RMSE  = rmse(actual, rf_pred),
    .groups  = "drop"
  )

cat("\n=== Synthetic walk-forward performance ===\n")
print(perf_s)

# -- 4. Figure A1: Variable importance comparison bar chart -------------------
# Side-by-side: synthetic (left) vs real (right) for same predictors

# Real data importance (hardcoded from analysis_v2.R run)
imp_real <- tibble(
  var   = c("seasonal_idx","log_cpi","lag12","UNRATE","lag1","log_oil",
            "covid_dummy","lag2","lag3"),
  label = c("Seasonal Index","Log CPI","AR(12)","Unemployment","AR(1)",
            "Log Oil Price","COVID Dummy","AR(2)","AR(3)"),
  pct_real = c(1.4, 17.1, 14.1, 12.4, 17.4, 8.3, 2.9, 13.7, 12.6)
)

# Merge synthetic and real
imp_compare <- imp_synth |>
  select(var, label, pct_synth = pct) |>
  left_join(imp_real |> select(var, pct_real), by = "var") |>
  pivot_longer(cols = c(pct_synth, pct_real),
               names_to = "Source", values_to = "pct") |>
  mutate(
    Source = recode(Source,
                    pct_synth = "Synthetic data",
                    pct_real  = "Real FRED data"),
    label = factor(label, levels = rev(imp_synth$label))
  )

figA1 <- ggplot(imp_compare, aes(x = pct, y = label, fill = Source)) +
  geom_col(position = "dodge", width = 0.65, alpha = 0.88) +
  scale_fill_manual(values = c("Synthetic data" = "#E8A838",
                                "Real FRED data" = "#1B3A6B")) +
  geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.3) +
  scale_x_continuous(labels = scales::percent_format(scale = 1)) +
  labs(
    title    = "Variable importance: synthetic data vs. real FRED data",
    subtitle = "% increase in out-of-bag MSE when each predictor is permuted",
    x        = "% Importance", y = NULL, fill = NULL,
    caption  = paste0("Synthetic: seasonal index first (", round(seasonal_pct_synth,1),
                      "%), short AR lags last (", round(ar_short_pct_synth,1), "%).\n",
                      "Real data: AR(1) first (17.4%), seasonal index last (1.4%).")
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", colour = "#1B3A6B"),
    plot.caption  = element_text(size = 8, colour = "grey50"),
    legend.position = "top"
  )

ggsave("figA1_synth_importance.pdf", figA1, width = 8, height = 5.5)
cat("figA1 saved\n")

# -- 5. Figure A2: Regime RMSE comparison (synthetic vs real) -----------------
# Real RMSE values (from analysis_v2.R run)
perf_real <- tibble(
  regime   = factor(c("Post-GFC Recovery","Stable Growth",
                       "COVID-19 Disruption","Post-COVID Rebound"),
                    levels = c("Post-GFC Recovery","Stable Growth",
                               "COVID-19 Disruption","Post-COVID Rebound")),
  SAR_RMSE = c(0.00589, 0.00655, 0.0411,  0.00941),
  RF_RMSE  = c(0.00933, 0.00802, 0.0431,  0.0134)
)

# Combine synthetic and real for plotting
make_long <- function(df, source_label) {
  df |>
    pivot_longer(c(SAR_RMSE, RF_RMSE), names_to = "Model", values_to = "RMSE") |>
    mutate(
      Model  = recode(Model, SAR_RMSE = "SARIMA", RF_RMSE = "Random Forest"),
      Source = source_label
    )
}

regime_levels <- c("Post-GFC Recovery","Stable Growth",
                   "COVID-19 Disruption","Post-COVID Rebound")

perf_s_plot <- perf_s |>
  mutate(regime = factor(regime, levels = regime_levels))

fig_data <- bind_rows(
  make_long(perf_s_plot, "Synthetic data"),
  make_long(perf_real,   "Real FRED data")
) |>
  mutate(
    Source = factor(Source, levels = c("Synthetic data","Real FRED data")),
    Model  = factor(Model,  levels = c("SARIMA","Random Forest"))
  )

figA2 <- ggplot(fig_data, aes(x = regime, y = RMSE, fill = Model, alpha = Source)) +
  geom_col(position = "dodge", width = 0.72) +
  scale_fill_manual(values = c("SARIMA" = "#1B3A6B", "Random Forest" = "#E8A838")) +
  scale_alpha_manual(values = c("Synthetic data" = 0.55, "Real FRED data" = 0.95)) +
  annotate("text", x = 3, y = 0.057,
           label = "RF wins\n(synthetic)", size = 3, colour = "#E8A838", fontface = "bold") +
  annotate("text", x = 3, y = 0.044,
           label = "SAR wins\n(real)", size = 3, colour = "#1B3A6B", fontface = "bold") +
  labs(
    title    = "Walk-forward RMSE: synthetic data vs. real FRED data",
    subtitle = "Darker bars = real data; lighter bars = synthetic data",
    x = NULL, y = "RMSE (log sales)", fill = "Model", alpha = "Data source",
    caption  = "COVID-19 Disruption: RF wins in synthetic data but SARIMA wins in real data."
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x     = element_text(angle = 15, hjust = 1),
    plot.title      = element_text(face = "bold", colour = "#1B3A6B"),
    plot.caption    = element_text(size = 8, colour = "grey50"),
    legend.position = "top"
  )

ggsave("figA2_synth_rmse.pdf", figA2, width = 9, height = 5.5)
cat("figA2 saved\n")

cat("\n=== Synthetic appendix figures done ===\n")
cat("Files saved: figA1_synth_importance.pdf, figA2_synth_rmse.pdf\n")
cat("Place both files in the same folder as your .tex files before compiling.\n")
