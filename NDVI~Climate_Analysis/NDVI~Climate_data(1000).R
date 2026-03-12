# ================================
# 0. packages
# ================================
library(data.table)
library(ggplot2)
library(scales)

# ================================
# 1. read data
# ================================
setwd("C:/Users/HUAWEI/Downloads/Dissertation/juneau_ndvi_1000sites/output/lsat_mean_trends")
mean_ndvi <- fread("juneau_ndvi_1000sites_mean_ndvi_timeseries.csv")
setwd("C:/Users/HUAWEI/Downloads/Dissertation/juneau_ndvi_1000sites/output/lsat_site_trends")
site_trend <- fread("juneau_ndvi_1000sites_site_ndvi_trends.csv")

# 只保留 full_period，避免以后如果有别的时期混进去
mean_ndvi <- mean_ndvi[period == "full_period"]
site_trend <- site_trend[period == "full_period"]

# ================================
# 2. Figure 1: mean NDVI greening trend
# ================================

# 计算 ribbon 上下界
mean_ndvi[, ndvi.lower := ndvi.mean - ndvi.sd]
mean_ndvi[, ndvi.upper := ndvi.mean + ndvi.sd]

# 线性拟合，方便把 slope / p-value 写到图上
fit_mean <- lm(ndvi.mean ~ year, data = mean_ndvi)
fit_sum <- summary(fit_mean)

slope_txt <- sprintf("Slope = %.4f NDVI yr^-1", coef(fit_mean)[2])
p_txt <- if (coef(fit_sum)[2,4] < 0.001) {
  "p < 0.001"
} else {
  sprintf("p = %.3f", coef(fit_sum)[2,4])
}
r2_txt <- sprintf("R² = %.3f", fit_sum$r.squared)

fig1 <- ggplot(mean_ndvi, aes(x = year, y = ndvi.mean)) +
  geom_ribbon(aes(ymin = ndvi.lower, ymax = ndvi.upper),
              fill = "grey80", alpha = 0.5) +
  geom_line(linewidth = 0.8, colour = "black") +
  geom_point(size = 2.2, colour = "black") +
  geom_smooth(method = "lm", se = FALSE,
              linetype = "dashed", linewidth = 0.9) +
  annotate(
    "text",
    x = min(mean_ndvi$year) + 2,
    y = max(mean_ndvi$ndvi.upper, na.rm = TRUE),
    label = paste(slope_txt, p_txt, r2_txt, sep = "\n"),
    hjust = 0, vjust = 1, size = 4
  ) +
  labs(
    title = "Mean annual maximum NDVI through time",
    subtitle = "Ribbon shows ±1 SD across sites",
    x = "Year",
    y = "Mean annual maximum NDVI"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  )

print(fig1)

ggsave(
  filename = "figure1_mean_ndvi_greening_trend.png",
  plot = fig1,
  width = 8,
  height = 5,
  dpi = 400
)

# ================================
# 3. Figure 2: site-level trend frequency
# ================================

# 重新整理类别名称和顺序
site_trend[, trend_group := fifelse(
  trend.cat == "browning.sig.p5",  "Browning (p < 0.05)",
  fifelse(trend.cat == "browning.sig.p10", "Browning (p < 0.10)",
          fifelse(trend.cat == "insig",            "Insignificant",
                  fifelse(trend.cat == "greening.sig.p10", "Greening (p < 0.10)",
                          fifelse(trend.cat == "greening.sig.p5",  "Greening (p < 0.05)", NA_character_))))
)]

trend_levels <- c(
  "Browning (p < 0.05)",
  "Browning (p < 0.10)",
  "Insignificant",
  "Greening (p < 0.10)",
  "Greening (p < 0.05)"
)

site_trend[, trend_group := factor(trend_group, levels = trend_levels)]

trend_freq <- site_trend[, .N, by = trend_group][order(trend_group)]
trend_freq[, pct := N / sum(N)]
trend_freq[, label := paste0(N, "\n(", percent(pct, accuracy = 0.1), ")")]

fig2 <- ggplot(trend_freq, aes(x = trend_group, y = N)) +
  geom_col(width = 0.75, fill = "grey55", colour = "black") +
  geom_text(aes(label = label), vjust = -0.3, size = 4) +
  labs(
    title = "Frequency distribution of site-level NDVI trends",
    x = NULL,
    y = "Number of sites"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 20, hjust = 1)
  ) +
  expand_limits(y = max(trend_freq$N) * 1.15)

print(fig2)

ggsave(
  filename = "figure2_site_ndvi_trend_frequency.png",
  plot = fig2,
  width = 8,
  height = 5,
  dpi = 400
)

# ================================
# 4. optional: summary table for writing
# ================================
trend_freq_out <- copy(trend_freq)
trend_freq_out[, pct := round(pct * 100, 1)]
print(trend_freq_out)

fwrite(trend_freq_out, "site_ndvi_trend_frequency_summary.csv")

# ================================
# 5. trends map
# ================================
print(names(site_trend))

# 重命名类别，方便图例更清楚
site_trend[, trend_group := fifelse(
  trend.cat == "browning.sig.p5",  "Browning (p < 0.05)",
  fifelse(trend.cat == "browning.sig.p10", "Browning (p < 0.10)",
          fifelse(trend.cat == "insig",            "Insignificant",
                  fifelse(trend.cat == "greening.sig.p10", "Greening (p < 0.10)",
                          fifelse(trend.cat == "greening.sig.p5",  "Greening (p < 0.05)", NA_character_))))
)]

trend_levels <- c(
  "Browning (p < 0.05)",
  "Browning (p < 0.10)",
  "Insignificant",
  "Greening (p < 0.10)",
  "Greening (p < 0.05)"
)

site_trend[, trend_group := factor(trend_group, levels = trend_levels)]

# 画图
fig3 <- ggplot(site_trend, aes(x = longitude, y = latitude, colour = trend_group)) +
  geom_point(size = 2.2, alpha = 0.9) +
  coord_equal() +
  labs(
    title = "Spatial pattern of site-level NDVI trends",
    x = "Longitude",
    y = "Latitude",
    colour = "Trend category"
  ) +
  scale_colour_manual(values = c(
    "Browning (p < 0.05)" = "#8B0000",
    "Browning (p < 0.10)" = "#E07A5F",
    "Insignificant"       = "grey65",
    "Greening (p < 0.10)" = "#81B29A",
    "Greening (p < 0.05)" = "#1B7F3B"
  )) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  )

print(fig3)

ggsave(
  "figure3_spatial_ndvi_trend_map.png",
  fig3,
  width = 7.5,
  height = 6,
  dpi = 400
)

#===============================================================================
# Part2. NDVI ~ Climate data
#===============================================================================

# =========================================================
# 0. Packages
# =========================================================
library(tidyr)
library(dplyr)
library(broom)
library(patchwork)
library(scales)

# =========================================================
# 1. Read data
# =========================================================
setwd("C:/Users/HUAWEI/Downloads/Dissertation/juneau_ndvi_1000sites/output/lsat_mean_trends")
ndvi <- fread("juneau_ndvi_1000sites_mean_ndvi_timeseries.csv")
ndvi_trend <- fread("juneau_ndvi_1000sites_mean_ndvi_trends.csv")

setwd("C:/Users/HUAWEI/Downloads/Dissertation/Juneau_Climate")
t2m <- fread("Juneau_nonIce_JJA_1984_2025_t2m_mean.csv")
soilT <- fread("Juneau_cliped_JJA_1984_2025_soilTemp_L1toL4.csv")
soilM <- fread("Juneau_cliped_JJA_1984_2025_soilMoistVWC_L1toL4.csv")
snowD <- fread("Juneau_nonIce_JJA_1984_2025_snowDepth_mean.csv")
snowSM <- fread("Juneau_nonIce_JJA_1984_2025_snowfallSum_snowmeltSum_total.csv")

# =========================================================
# 2. Clean and harmonise
# =========================================================
ndvi <- ndvi[period == "full_period", .(year, ndvi.mean, ndvi.sd, n.sites)]
ndvi[, year := as.integer(year)]

t2m <- t2m[, .(year = as.integer(year), temperature_2m)]
t2m[, air_temp_C := temperature_2m - 273.15]

soilT <- soilT[, .(
  year = as.integer(year),
  soil_temperature_level_1,
  soil_temperature_level_2,
  soil_temperature_level_3,
  soil_temperature_level_4
)]
soilT[, soil_temp_mean_C := rowMeans(.SD) - 273.15,
      .SDcols = patterns("^soil_temperature_level_")]
soilT[, soil_temp_min_C := do.call(pmin, .SD) - 273.15,
      .SDcols = patterns("^soil_temperature_level_")]
soilT[, soil_temp_max_C := do.call(pmax, .SD) - 273.15,
      .SDcols = patterns("^soil_temperature_level_")]

soilM <- soilM[, .(
  year = as.integer(year),
  volumetric_soil_water_layer_1,
  volumetric_soil_water_layer_2,
  volumetric_soil_water_layer_3,
  volumetric_soil_water_layer_4
)]
soilM[, soil_moist_mean := rowMeans(.SD),
      .SDcols = patterns("^volumetric_soil_water_layer_")]
soilM[, soil_moist_min := do.call(pmin, .SD),
      .SDcols = patterns("^volumetric_soil_water_layer_")]
soilM[, soil_moist_max := do.call(pmax, .SD),
      .SDcols = patterns("^volumetric_soil_water_layer_")]

snowD <- snowD[, .(year = as.integer(year), snow_depth)]

snowSM <- snowSM[, .(
  year = as.integer(year),
  snowfall_sum,
  snowmelt_sum
)]

# =========================================================
# 3. Merge
# =========================================================
dt <- Reduce(function(x, y) merge(x, y, by = "year", all.x = TRUE),
             list(
               ndvi,
               t2m[, .(year, air_temp_C)],
               soilT[, .(year,
                         soil_temperature_level_1,
                         soil_temperature_level_2,
                         soil_temperature_level_3,
                         soil_temperature_level_4,
                         soil_temp_mean_C,
                         soil_temp_min_C,
                         soil_temp_max_C)],
               soilM[, .(year,
                         volumetric_soil_water_layer_1,
                         volumetric_soil_water_layer_2,
                         volumetric_soil_water_layer_3,
                         volumetric_soil_water_layer_4,
                         soil_moist_mean,
                         soil_moist_min,
                         soil_moist_max)],
               snowD,
               snowSM
             ))

# =========================================================
# 4. NDVI trend summary
# =========================================================
print(ndvi_trend)

# =========================================================
# 5. Simple regression helper
# =========================================================
run_simple_model <- function(data, xvar, yvar = "ndvi.mean") {
  form <- as.formula(paste(yvar, "~", xvar))
  fit <- lm(form, data = data)
  td <- tidy(fit)
  gd <- glance(fit)
  
  out <- data.table(
    predictor = xvar,
    intercept = td$estimate[td$term == "(Intercept)"],
    slope = td$estimate[td$term == xvar],
    p_value = td$p.value[td$term == xvar],
    r_squared = gd$r.squared,
    adj_r_squared = gd$adj.r.squared,
    n = nobs(fit)
  )
  return(out)
}

results_tbl <- rbindlist(list(
  run_simple_model(dt, "air_temp_C"),
  run_simple_model(dt, "soil_temp_mean_C"),
  run_simple_model(dt, "soil_moist_mean"),
  run_simple_model(dt, "snow_depth"),
  run_simple_model(dt, "snowfall_sum"),
  run_simple_model(dt, "snowmelt_sum")
))

# Pearson correlations
cor_tbl <- rbindlist(lapply(
  c("air_temp_C", "soil_temp_mean_C", "soil_moist_mean",
    "snow_depth", "snowfall_sum", "snowmelt_sum"),
  function(v) {
    ct <- cor.test(dt[[v]], dt[["ndvi.mean"]], method = "pearson")
    data.table(
      predictor = v,
      pearson_r = unname(ct$estimate),
      pearson_p = ct$p.value
    )
  }
))

results_tbl <- merge(results_tbl, cor_tbl, by = "predictor")
print(results_tbl)

fwrite(results_tbl, "ndvi_climate_simple_model_results.csv")

# =========================================================
# 6. Plot 1: NDVI + temperature time series
#    top = temperature
#    bottom = NDVI
# =========================================================
p_temp_ts <- ggplot(dt, aes(x = year)) +
  geom_ribbon(aes(ymin = soil_temp_min_C, ymax = soil_temp_max_C),
              fill = "grey75", alpha = 0.5) +
  geom_line(aes(y = soil_temp_mean_C), linewidth = 1) +
  geom_line(aes(y = air_temp_C), linetype = "dashed", linewidth = 0.9) +
  labs(
    title = "Summer temperature through time",
    subtitle = "Solid line = mean soil temperature; ribbon = range across 4 soil layers; dashed line = 2 m air temperature",
    x = NULL,
    y = expression(Temperature~(degree*C))
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  )

p_ndvi_ts <- ggplot(dt, aes(x = year, y = ndvi.mean)) +
  geom_ribbon(aes(ymin = ndvi.mean - ndvi.sd, ymax = ndvi.mean + ndvi.sd),
              fill = "grey80", alpha = 0.5) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  labs(
    title = "Mean annual maximum NDVI through time",
    subtitle = "Ribbon = ±1 SD across sites",
    x = "Year",
    y = "Mean annual maximum NDVI"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  )

fig_temp_ts <- p_temp_ts / p_ndvi_ts + plot_layout(heights = c(1, 1))
ggsave("figure_temp_ndvi_timeseries.png", fig_temp_ts, width = 9, height = 8, dpi = 400)
print(fig_temp_ts)
#==========================================================
# Multi-panel: NDVI vs climate data
#==========================================================
# Publication-style figures (multi-panel + unified theme)
# =========================================================
# -----------------------------
# Common publication theme
# -----------------------------
theme_pub <- theme_classic(base_size = 12) +
  theme(
    plot.title   = element_text(face = "bold", size = 11, hjust = 0.5),
    plot.subtitle= element_text(size = 9, hjust = 0.5),
    axis.title   = element_text(face = "bold", size = 11),
    axis.text    = element_text(size = 10, colour = "black"),
    strip.text   = element_text(face = "bold", size = 10),
    strip.background = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
    plot.margin  = margin(6, 6, 6, 6)
  )

# -----------------------------
# Consistent colours
# -----------------------------
pt_col   <- "#355C7D"   # muted academic blue
line_col <- "#C06C84"   # muted rose for lm line
rib_col  <- "#BFC7D5"   # soft grey-blue ribbon
ts_col   <- "#2F4858"   # dark blue-grey for time series

# -----------------------------
# Helper: extract label text safely
# -----------------------------
get_label <- function(tbl, pred_name) {
  out <- tbl[predictor == pred_name,
             sprintf("r = %.2f\np = %.3f\nR² = %.2f", pearson_r, pearson_p, r_squared)]
  if (length(out) == 0 || is.na(out)) {
    out <- "r = NA\np = NA\nR² = NA"
  }
  out[1]
}

# -----------------------------
# Helper: make scatter plot
# -----------------------------
mk_scatter <- function(data, xvar, xlab, title_txt, label_txt) {
  
  xvals <- data[[xvar]]
  yvals <- data[["ndvi.mean"]]
  
  x_rng <- range(xvals, na.rm = TRUE)
  y_rng <- range(yvals, na.rm = TRUE)
  
  x_pos <- x_rng[1] + 0.03 * diff(x_rng)
  y_pos <- y_rng[2] - 0.05 * diff(y_rng)
  
  ggplot(data, aes(x = .data[[xvar]], y = .data[["ndvi.mean"]])) +
    geom_point(size = 2.1, alpha = 0.9, colour = pt_col) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.8,
                colour = line_col, fill = rib_col) +
    annotate("text",
             x = x_pos,
             y = y_pos,
             label = label_txt,
             hjust = 0, vjust = 1,
             size = 3.5) +
    labs(
      title = title_txt,
      x = xlab,
      y = "Mean annual maximum NDVI"
    ) +
    theme_pub
}

# =========================================================
# 7. Multi-panel scatter plot: NDVI vs all environmental variables
# =========================================================

label_air   <- get_label(results_tbl, "air_temp_C")
label_soilT <- get_label(results_tbl, "soil_temp_mean_C")
label_sm    <- get_label(results_tbl, "soil_moist_mean")
label_sd    <- get_label(results_tbl, "snow_depth")
label_sf    <- get_label(results_tbl, "snowfall_sum")
label_smelt <- get_label(results_tbl, "snowmelt_sum")

p_air <- mk_scatter(
  dt,
  xvar      = "air_temp_C",
  xlab      = expression(Summer~air~temperature~(degree*C)),
  title_txt = "NDVI vs summer 2 m air temperature",
  label_txt = label_air
)

p_soilT <- mk_scatter(
  dt,
  xvar      = "soil_temp_mean_C",
  xlab      = expression(Mean~summer~soil~temperature~(degree*C)),
  title_txt = "NDVI vs mean summer soil temperature",
  label_txt = label_soilT
)

p_sm_scatter <- mk_scatter(
  dt,
  xvar      = "soil_moist_mean",
  xlab      = "Mean summer volumetric soil water content",
  title_txt = "NDVI vs mean summer soil moisture",
  label_txt = label_sm
)

p_sd <- mk_scatter(
  dt,
  xvar      = "snow_depth",
  xlab      = "Summer mean snow depth",
  title_txt = "NDVI vs snow depth",
  label_txt = label_sd
)

p_sf <- mk_scatter(
  dt,
  xvar      = "snowfall_sum",
  xlab      = "Summer snowfall sum",
  title_txt = "NDVI vs snowfall",
  label_txt = label_sf
)

p_smelt <- mk_scatter(
  dt,
  xvar      = "snowmelt_sum",
  xlab      = "Summer snowmelt sum",
  title_txt = "NDVI vs snowmelt",
  label_txt = label_smelt
)

fig_ndvi_all_scatter <- wrap_plots(
  p_air, p_soilT, p_sm_scatter,
  p_sd, p_sf, p_smelt,
  ncol = 2
) +
  plot_annotation(tag_levels = "a") &
  theme(
    legend.position = "none"
  )

ggsave(
  filename = "figure_ndvi_all_environmental_relationships_multipanel.png",
  plot     = fig_ndvi_all_scatter,
  width    = 180,
  height   = 240,
  units    = "mm",
  dpi      = 600,
  bg       = "white"
)

print(fig_ndvi_all_scatter)
# =========================================================
# 7. Plot 2: NDVI vs temperature
# =========================================================
label_air <- results_tbl[predictor == "air_temp_C",
                         sprintf("r = %.2f\np = %.3f\nR² = %.2f", pearson_r, pearson_p, r_squared)]
label_soilT <- results_tbl[predictor == "soil_temp_mean_C",
                           sprintf("r = %.2f\np = %.3f\nR² = %.2f", pearson_r, pearson_p, r_squared)]

p_air <- ggplot(dt, aes(x = air_temp_C, y = ndvi.mean)) +
  geom_point(size = 2.2, alpha = 0.9) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.9) +
  annotate("text",
           x = min(dt$air_temp_C, na.rm = TRUE),
           y = max(dt$ndvi.mean, na.rm = TRUE),
           label = label_air, hjust = 0, vjust = 1, size = 4) +
  labs(
    title = "NDVI vs summer 2 m air temperature",
    x = expression(Summer~air~temperature~(degree*C)),
    y = "Mean annual maximum NDVI"
  ) +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(face = "bold"),
        axis.title = element_text(face = "bold"))

p_soilT <- ggplot(dt, aes(x = soil_temp_mean_C, y = ndvi.mean)) +
  geom_point(size = 2.2, alpha = 0.9) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.9) +
  annotate("text",
           x = min(dt$soil_temp_mean_C, na.rm = TRUE),
           y = max(dt$ndvi.mean, na.rm = TRUE),
           label = label_soilT, hjust = 0, vjust = 1, size = 4) +
  labs(
    title = "NDVI vs mean summer soil temperature",
    x = expression(Mean~summer~soil~temperature~(degree*C)),
    y = "Mean annual maximum NDVI"
  ) +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(face = "bold"),
        axis.title = element_text(face = "bold"))

fig_temp_scatter <- p_air + p_soilT
ggsave("figure_ndvi_temperature_relationships.png", fig_temp_scatter, width = 11, height = 5, dpi = 400)
print()
# =========================================================
# 8. Plot 3: soil moisture time series
# =========================================================
p_sm_ts <- ggplot(dt, aes(x = year)) +
  geom_ribbon(aes(ymin = soil_moist_min, ymax = soil_moist_max),
              fill = "grey75", alpha = 0.5) +
  geom_line(aes(y = soil_moist_mean), linewidth = 1) +
  labs(
    title = "Summer soil moisture through time",
    subtitle = "Solid line = mean soil moisture; ribbon = range across 4 soil layers",
    x = "Year",
    y = "Volumetric soil water content"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  )

ggsave("figure_soil_moisture_timeseries.png", p_sm_ts, width = 9, height = 4.5, dpi = 400)

# =========================================================
# 9. Plot 4: NDVI vs soil moisture
# =========================================================
label_sm <- results_tbl[predictor == "soil_moist_mean",
                        sprintf("r = %.2f\np = %.3f\nR² = %.2f", pearson_r, pearson_p, r_squared)]

p_sm_scatter <- ggplot(dt, aes(x = soil_moist_mean, y = ndvi.mean)) +
  geom_point(size = 2.2, alpha = 0.9) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.9) +
  annotate("text",
           x = min(dt$soil_moist_mean, na.rm = TRUE),
           y = max(dt$ndvi.mean, na.rm = TRUE),
           label = label_sm, hjust = 0, vjust = 1, size = 4) +
  labs(
    title = "NDVI vs mean summer soil moisture",
    x = "Mean summer volumetric soil water content",
    y = "Mean annual maximum NDVI"
  ) +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(face = "bold"),
        axis.title = element_text(face = "bold"))

ggsave("figure_ndvi_soilmoisture_relationship.png", p_sm_scatter, width = 6, height = 5, dpi = 400)

# =========================================================
# 10. Plot 5: snow time series
# =========================================================
snow_long <- data.table::melt(
  dt[, .(year, snow_depth, snowfall_sum, snowmelt_sum)],
  id.vars = "year",
  variable.name = "metric",
  value.name = "value"
)

setDT(snow_long)

snow_long[, metric := factor(
  metric,
  levels = c("snow_depth", "snowfall_sum", "snowmelt_sum"),
  labels = c("Snow depth (mean)", "Snowfall (sum)", "Snowmelt (sum)")
)]

p_snow_ts <- ggplot(snow_long, aes(x = year, y = value)) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ metric, scales = "free_y", ncol = 1) +
  labs(
    title = "Summer snow variables through time",
    x = "Year",
    y = NULL
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold")
  )

ggsave("figure_snow_timeseries.png", p_snow_ts, width = 8, height = 8, dpi = 400)

# =========================================================
# 11. Plot 6: NDVI vs snow metrics
# =========================================================
mk_scatter <- function(xvar, xlab, title_txt, label_txt) {
  ggplot(dt, aes_string(x = xvar, y = "ndvi.mean")) +
    geom_point(size = 2.2, alpha = 0.9) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.9) +
    annotate("text",
             x = min(dt[[xvar]], na.rm = TRUE),
             y = max(dt$ndvi.mean, na.rm = TRUE),
             label = label_txt, hjust = 0, vjust = 1, size = 4) +
    labs(
      title = title_txt,
      x = xlab,
      y = "Mean annual maximum NDVI"
    ) +
    theme_classic(base_size = 13) +
    theme(plot.title = element_text(face = "bold"),
          axis.title = element_text(face = "bold"))
}

label_sd <- results_tbl[predictor == "snow_depth",
                        sprintf("r = %.2f\np = %.3f\nR² = %.2f", pearson_r, pearson_p, r_squared)]
label_sf <- results_tbl[predictor == "snowfall_sum",
                        sprintf("r = %.2f\np = %.3f\nR² = %.2f", pearson_r, pearson_p, r_squared)]
label_smelt <- results_tbl[predictor == "snowmelt_sum",
                           sprintf("r = %.2f\np = %.3f\nR² = %.2f", pearson_r, pearson_p, r_squared)]

p_sd <- mk_scatter("snow_depth", "Summer mean snow depth", "NDVI vs snow depth", label_sd)
p_sf <- mk_scatter("snowfall_sum", "Summer snowfall sum", "NDVI vs snowfall", label_sf)
p_smelt <- mk_scatter("snowmelt_sum", "Summer snowmelt sum", "NDVI vs snowmelt", label_smelt)

fig_snow_scatter <- (p_sd + p_sf + p_smelt)
ggsave("figure_ndvi_snow_relationships.png", fig_snow_scatter, width = 15, height = 5, dpi = 400)



#===============================================================================
# (Secondly repeat) NDVI ~ climate data
#===============================================================================
# =========================================================
# 1. Paths
# =========================================================
ndvi_dir <- "C:/Users/HUAWEI/Downloads/Dissertation/juneau_ndvi_1000sites/output/lsat_mean_trends"
clim_dir <- "C:/Users/HUAWEI/Downloads/Dissertation/Juneau_Climate"
fig_dir  <- file.path(clim_dir, "figures(1000)")

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# =========================================================
# 2. Read data
# =========================================================
ndvi  <- fread(file.path(ndvi_dir, "juneau_ndvi_1000sites_mean_ndvi_timeseries.csv"))
t2m   <- fread(file.path(clim_dir, "Juneau_nonIce_JJA_1984_2025_t2m_mean.csv"))
soilT <- fread(file.path(clim_dir, "Juneau_cliped_JJA_1984_2025_soilTemp_L1toL4.csv"))
soilM <- fread(file.path(clim_dir, "Juneau_cliped_JJA_1984_2025_soilMoistVWC_L1toL4.csv"))
snowD <- fread(file.path(clim_dir, "Juneau_nonIce_JJA_1984_2025_snowDepth_mean.csv"))
snowSM <- fread(file.path(clim_dir, "Juneau_nonIce_JJA_1984_2025_snowfallSum_snowmeltSum_total.csv"))

# =========================================================
# 3. Clean data
# =========================================================
ndvi <- ndvi[period == "full_period", .(year = as.integer(year), ndvi.mean, ndvi.sd)]

t2m <- t2m[, .(
  year = as.integer(year),
  air_temp_C = temperature_2m - 273.15
)]

soilT <- soilT[, .(
  year = as.integer(year),
  st1 = soil_temperature_level_1 - 273.15,
  st2 = soil_temperature_level_2 - 273.15,
  st3 = soil_temperature_level_3 - 273.15,
  st4 = soil_temperature_level_4 - 273.15
)]
soilT[, soil_temp_mean := rowMeans(.SD), .SDcols = c("st1", "st2", "st3", "st4")]
soilT[, soil_temp_min  := do.call(pmin, .SD), .SDcols = c("st1", "st2", "st3", "st4")]
soilT[, soil_temp_max  := do.call(pmax, .SD), .SDcols = c("st1", "st2", "st3", "st4")]

soilM <- soilM[, .(
  year = as.integer(year),
  sm1 = volumetric_soil_water_layer_1,
  sm2 = volumetric_soil_water_layer_2,
  sm3 = volumetric_soil_water_layer_3,
  sm4 = volumetric_soil_water_layer_4
)]
soilM[, soil_moist_mean := rowMeans(.SD), .SDcols = c("sm1", "sm2", "sm3", "sm4")]
soilM[, soil_moist_min  := do.call(pmin, .SD), .SDcols = c("sm1", "sm2", "sm3", "sm4")]
soilM[, soil_moist_max  := do.call(pmax, .SD), .SDcols = c("sm1", "sm2", "sm3", "sm4")]

snowD <- snowD[, .(
  year = as.integer(year),
  snow_depth
)]

snowSM <- snowSM[, .(
  year = as.integer(year),
  snowfall_sum,
  snowmelt_sum
)]

# =========================================================
# 4. Merge all
# =========================================================
dt <- Reduce(function(x, y) merge(x, y, by = "year", all = FALSE),
             list(ndvi, t2m, soilT, soilM, snowD, snowSM))

# =========================================================
# 5. Standardise to z-score
# =========================================================
zfun <- function(x) as.numeric(scale(x))

dt[, ndvi_z       := zfun(ndvi.mean)]
dt[, air_temp_z   := zfun(air_temp_C)]
dt[, soil_temp_z  := zfun(soil_temp_mean)]
dt[, soil_moist_z := zfun(soil_moist_mean)]
dt[, snow_depth_z := zfun(snow_depth)]
dt[, snowfall_z   := zfun(snowfall_sum)]
dt[, snowmelt_z   := zfun(snowmelt_sum)]

# layer ribbons
dt[, st1_z := zfun(st1)]
dt[, st2_z := zfun(st2)]
dt[, st3_z := zfun(st3)]
dt[, st4_z := zfun(st4)]
dt[, soil_temp_layer_min_z := do.call(pmin, .SD), .SDcols = c("st1_z", "st2_z", "st3_z", "st4_z")]
dt[, soil_temp_layer_max_z := do.call(pmax, .SD), .SDcols = c("st1_z", "st2_z", "st3_z", "st4_z")]

dt[, sm1_z := zfun(sm1)]
dt[, sm2_z := zfun(sm2)]
dt[, sm3_z := zfun(sm3)]
dt[, sm4_z := zfun(sm4)]
dt[, soil_moist_layer_min_z := do.call(pmin, .SD), .SDcols = c("sm1_z", "sm2_z", "sm3_z", "sm4_z")]
dt[, soil_moist_layer_max_z := do.call(pmax, .SD), .SDcols = c("sm1_z", "sm2_z", "sm3_z", "sm4_z")]

# =========================================================
# 6. High / low NDVI shading
# =========================================================
q_hi <- quantile(dt$ndvi.mean, 0.85, na.rm = TRUE)
q_lo <- quantile(dt$ndvi.mean, 0.15, na.rm = TRUE)

dt[, ndvi_phase := fifelse(
  ndvi.mean >= q_hi, "High NDVI",
  fifelse(ndvi.mean <= q_lo, "Low NDVI", "Normal")
)]

shade_dt <- dt[ndvi_phase != "Normal", .(year, ndvi_phase)]
shade_dt[, xmin := year - 0.5]
shade_dt[, xmax := year + 0.5]
shade_dt[, ymin := -Inf]
shade_dt[, ymax := Inf]

# =========================================================
# 7. Common theme
# =========================================================
theme_panel <- theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(face = "bold"),
    legend.title = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 10)
  )

shade_scale <- scale_fill_manual(values = c(
  "High NDVI" = "#A5D6A7",   # green
  "Low NDVI"  = "#FFF59D"    # yellow
))

xlims <- range(dt$year, na.rm = TRUE)

# =========================================================
# 8. Panel A: Temperature + NDVI
# =========================================================
p_temp <- ggplot() +
  geom_rect(
    data = shade_dt,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = ndvi_phase),
    inherit.aes = FALSE,
    alpha = 0.22
  ) +
  geom_ribbon(
    data = dt,
    aes(x = year, ymin = soil_temp_layer_min_z, ymax = soil_temp_layer_max_z),
    fill = "#C49A6C", alpha = 0.18
  ) +
  geom_line(data = dt, aes(x = year, y = ndvi_z, colour = "NDVI"), linewidth = 1.2) +
  geom_line(data = dt, aes(x = year, y = air_temp_z, colour = "Air temperature", linetype = "Air temperature"), linewidth = 0.95) +
  geom_line(data = dt, aes(x = year, y = soil_temp_z, colour = "Mean soil temperature"), linewidth = 1.0) +
  scale_colour_manual(values = c(
    "NDVI" = "black",
    "Air temperature" = "#E66101",
    "Mean soil temperature" = "#8C6D31"
  )) +
  scale_linetype_manual(values = c("Air temperature" = "dashed")) +
  shade_scale +
  coord_cartesian(xlim = xlims) +
  labs(
    title = "A. NDVI and summer temperature",
    x = NULL,
    y = "Standardised anomaly (z-score)"
  ) +
  theme_panel +
  guides(linetype = "none")

# =========================================================
# 9. Panel B: Soil moisture + NDVI
# =========================================================
p_moist <- ggplot() +
  geom_rect(
    data = shade_dt,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = ndvi_phase),
    inherit.aes = FALSE,
    alpha = 0.22
  ) +
  geom_ribbon(
    data = dt,
    aes(x = year, ymin = soil_moist_layer_min_z, ymax = soil_moist_layer_max_z),
    fill = "#6BAED6", alpha = 0.18
  ) +
  geom_line(data = dt, aes(x = year, y = ndvi_z, colour = "NDVI"), linewidth = 1.2) +
  geom_line(data = dt, aes(x = year, y = soil_moist_z, colour = "Mean soil moisture"), linewidth = 1.0) +
  scale_colour_manual(values = c(
    "NDVI" = "black",
    "Mean soil moisture" = "#3182BD"
  )) +
  shade_scale +
  coord_cartesian(xlim = xlims) +
  labs(
    title = "B. NDVI and summer soil moisture",
    x = NULL,
    y = "Standardised anomaly (z-score)"
  ) +
  theme_panel

# =========================================================
# 10. Panel C: Snow depth + NDVI
# =========================================================
p_snowdepth <- ggplot() +
  geom_rect(
    data = shade_dt,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = ndvi_phase),
    inherit.aes = FALSE,
    alpha = 0.22
  ) +
  geom_line(data = dt, aes(x = year, y = ndvi_z, colour = "NDVI"), linewidth = 1.2) +
  geom_line(data = dt, aes(x = year, y = snow_depth_z, colour = "Snow depth"), linewidth = 1.0) +
  scale_colour_manual(values = c(
    "NDVI" = "black",
    "Snow depth" = "#756BB1"
  )) +
  shade_scale +
  coord_cartesian(xlim = xlims) +
  labs(
    title = "C. NDVI and summer snow depth",
    x = "Year",
    y = "Standardised anomaly (z-score)"
  ) +
  theme_panel

# =========================================================
# 11. Panel D: Snowfall + snowmelt + NDVI
# =========================================================
p_snowflux <- ggplot() +
  geom_rect(
    data = shade_dt,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = ndvi_phase),
    inherit.aes = FALSE,
    alpha = 0.22
  ) +
  geom_line(data = dt, aes(x = year, y = ndvi_z, colour = "NDVI"), linewidth = 1.2) +
  geom_line(data = dt, aes(x = year, y = snowfall_z, colour = "Snowfall", linetype = "Snowfall"), linewidth = 0.95) +
  geom_line(data = dt, aes(x = year, y = snowmelt_z, colour = "Snowmelt"), linewidth = 1.0) +
  scale_colour_manual(values = c(
    "NDVI" = "black",
    "Snowfall" = "#E78AC3",
    "Snowmelt" = "#D7301F"
  )) +
  scale_linetype_manual(values = c("Snowfall" = "dotdash")) +
  shade_scale +
  coord_cartesian(xlim = xlims) +
  labs(
    title = "D. NDVI, snowfall and snowmelt",
    x = "Year",
    y = "Standardised anomaly (z-score)"
  ) +
  theme_panel +
  guides(linetype = "none")

# =========================================================
# 12. Combine
# =========================================================
fig_multi <- (p_temp + p_moist) / (p_snowdepth + p_snowflux) +
  plot_annotation(
    title = "Interannual variability of NDVI and summer climate factors (1984–2025)",
    subtitle = "Black line shows NDVI in all panels. Green shading indicates high-NDVI years; yellow shading indicates low-NDVI years.\nBrown ribbon = range across 4 soil temperature layers; blue ribbon = range across 4 soil moisture layers.",
    theme = theme(
      plot.title = element_text(face = "bold", size = 18),
      plot.subtitle = element_text(size = 11)
    )
  )

print(fig_multi)

ggsave(
  filename = file.path(fig_dir, "figure_ndvi_climate_multipanel_edited.png"),
  plot = fig_multi,
  width = 14,
  height = 10,
  dpi = 500
)

#===============================================================================
# Part3. Break year
#===============================================================================
# =========================================================
# PART 1. Set paths
# 目的：避免反复 setwd() 导致读错文件
# =========================================================
out_dir  <- file.path(clim_dir, "break_analysis(1000)")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# =========================================================
# PART 2. Read data
# 目的：读取 annual NDVI 与 JJA climate data
# =========================================================
ndvi  <- fread(file.path(ndvi_dir, "juneau_ndvi_1000sites_mean_ndvi_timeseries.csv"))
t2m   <- fread(file.path(clim_dir, "Juneau_nonIce_JJA_1984_2025_t2m_mean.csv"))
soilT <- fread(file.path(clim_dir, "Juneau_cliped_JJA_1984_2025_soilTemp_L1toL4.csv"))
soilM <- fread(file.path(clim_dir, "Juneau_cliped_JJA_1984_2025_soilMoistVWC_L1toL4.csv"))
snowD <- fread(file.path(clim_dir, "Juneau_nonIce_JJA_1984_2025_snowDepth_mean.csv"))
snowSM <- fread(file.path(clim_dir, "Juneau_nonIce_JJA_1984_2025_snowfallSum_snowmeltSum_total.csv"))

# =========================================================
# PART 3. Clean and harmonise
# 目的：统一 year 列，构造我们真正要分析的变量
# =========================================================

# NDVI: 只用 full_period
ndvi <- ndvi[period == "full_period", .(
  year = as.integer(year),
  ndvi = ndvi.mean
)]

# Air temperature: Kelvin -> Celsius
t2m <- t2m[, .(
  year = as.integer(year),
  air_temp_C = temperature_2m - 273.15
)]

# Soil temperature: 4 layers -> mean
soilT <- soilT[, .(
  year = as.integer(year),
  st1 = soil_temperature_level_1 - 273.15,
  st2 = soil_temperature_level_2 - 273.15,
  st3 = soil_temperature_level_3 - 273.15,
  st4 = soil_temperature_level_4 - 273.15
)]
soilT[, soil_temp_mean_C := rowMeans(.SD), .SDcols = c("st1", "st2", "st3", "st4")]

# Soil moisture: 4 layers -> mean
soilM <- soilM[, .(
  year = as.integer(year),
  sm1 = volumetric_soil_water_layer_1,
  sm2 = volumetric_soil_water_layer_2,
  sm3 = volumetric_soil_water_layer_3,
  sm4 = volumetric_soil_water_layer_4
)]
soilM[, soil_moist_mean := rowMeans(.SD), .SDcols = c("sm1", "sm2", "sm3", "sm4")]

# Snow
snowD <- snowD[, .(
  year = as.integer(year),
  snow_depth
)]

snowSM <- snowSM[, .(
  year = as.integer(year),
  snowfall_sum,
  snowmelt_sum
)]

# =========================================================
# PART 4. Merge all variables by year
# 目的：把所有变量整合到一个表里，后面统一分析
# =========================================================
dt <- Reduce(function(x, y) merge(x, y, by = "year", all = FALSE),
             list(ndvi, t2m, soilT[, .(year, soil_temp_mean_C)],
                  soilM[, .(year, soil_moist_mean)],
                  snowD, snowSM))

# 检查时间范围
print(range(dt$year))
print(names(dt))

# =========================================================
# PART 5. Define candidate break years
# 目的：只扫描 hypothesis 最关注的年份范围
# 说明：这里用 1990–2002，可按需要微调
# =========================================================
candidate_breaks <- 1984:2025 #可以改成1985:2025

# 为了避免某一段太短，设置最少样本数
min_n_each_side <- 8

# 我们要分析的 climate variables
clim_vars <- c(
  "air_temp_C",
  "soil_temp_mean_C",
  "soil_moist_mean",
  "snow_depth",
  "snowfall_sum",
  "snowmelt_sum"
)

# =========================================================
# PART 6. Function: fit one linear model and extract results
# 目的：给 early / late period 统一输出 slope, r, p, R2
# =========================================================
fit_simple_stats <- function(data, xvar, yvar = "ndvi") {
  tmp <- data[complete.cases(data[, ..xvar], data[, ..yvar])]
  
  if (nrow(tmp) < 5) {
    return(data.table(
      n = nrow(tmp),
      slope = NA_real_,
      intercept = NA_real_,
      r = NA_real_,
      p = NA_real_,
      r2 = NA_real_
    ))
  }
  
  form <- as.formula(paste(yvar, "~", xvar))
  fit <- lm(form, data = tmp)
  sm  <- summary(fit)
  ct  <- cor.test(tmp[[xvar]], tmp[[yvar]], method = "pearson")
  
  data.table(
    n = nrow(tmp),
    slope = coef(fit)[2],
    intercept = coef(fit)[1],
    r = unname(ct$estimate),
    p = ct$p.value,
    r2 = sm$r.squared
  )
}

# =========================================================
# PART 7. Function: scan break years for one climate variable
# 目的：对每个候选 break year 计算两段回归总 RSS
# =========================================================
scan_break_year <- function(data, xvar, yvar = "ndvi",
                            candidate_years, min_n_each_side = 8) {
  
  out_list <- list()
  
  for (b in candidate_years) {
    early <- data[year <= b & complete.cases(data[, ..xvar], data[, ..yvar])]
    late  <- data[year >  b & complete.cases(data[, ..xvar], data[, ..yvar])]
    
    # 保证前后两段都有足够样本
    if (nrow(early) < min_n_each_side || nrow(late) < min_n_each_side) next
    
    fit_early <- lm(as.formula(paste(yvar, "~", xvar)), data = early)
    fit_late  <- lm(as.formula(paste(yvar, "~", xvar)), data = late)
    
    rss_early <- sum(residuals(fit_early)^2)
    rss_late  <- sum(residuals(fit_late)^2)
    rss_total <- rss_early + rss_late
    
    out_list[[length(out_list) + 1]] <- data.table(
      variable = xvar,
      break_year = b,
      n_early = nrow(early),
      n_late = nrow(late),
      rss_early = rss_early,
      rss_late = rss_late,
      rss_total = rss_total,
      slope_early = coef(fit_early)[2],
      slope_late = coef(fit_late)[2]
    )
  }
  
  rbindlist(out_list, fill = TRUE)
}

# =========================================================
# PART 8. Run break scan for all variables
# 目的：自动找到每个变量最优 break year
# =========================================================
scan_results <- rbindlist(lapply(clim_vars, function(v) {
  scan_break_year(
    data = dt,
    xvar = v,
    yvar = "ndvi",
    candidate_years = candidate_breaks,
    min_n_each_side = min_n_each_side
  )
}), fill = TRUE)

# 保存完整扫描结果
fwrite(scan_results, file.path(out_dir, "break_year_scan_all_variables.csv"))

# 取每个变量 RSS 最小的年份 = 最优 break year
best_breaks <- scan_results[, .SD[which.min(rss_total)], by = variable]

# 增加方向解释
best_breaks[, direction_early := fifelse(slope_early > 0, "positive", "negative")]
best_breaks[, direction_late  := fifelse(slope_late  > 0, "positive", "negative")]

print(best_breaks)
fwrite(best_breaks, file.path(out_dir, "best_break_years.csv"))

# =========================================================
# PART 9. Build early vs late correlation table
# 目的：用最优 break year 正式输出前后两段统计量
# =========================================================
early_late_list <- list()

for (i in 1:nrow(best_breaks)) {
  v <- best_breaks$variable[i]
  b <- best_breaks$break_year[i]
  
  early <- dt[year <= b]
  late  <- dt[year >  b]
  
  stats_early <- fit_simple_stats(early, xvar = v, yvar = "ndvi")
  stats_late  <- fit_simple_stats(late,  xvar = v, yvar = "ndvi")
  
  stats_early[, `:=`(
    variable = v,
    break_year = b,
    period = paste0("Early (<= ", b, ")")
  )]
  
  stats_late[, `:=`(
    variable = v,
    break_year = b,
    period = paste0("Late (> ", b, ")")
  )]
  
  early_late_list[[length(early_late_list) + 1]] <- stats_early
  early_late_list[[length(early_late_list) + 1]] <- stats_late
}

early_late_tbl <- rbindlist(early_late_list, fill = TRUE)

# 调整列顺序
setcolorder(early_late_tbl, c(
  "variable", "break_year", "period", "n",
  "slope", "intercept", "r", "p", "r2"
))

print(early_late_tbl)
fwrite(early_late_tbl, file.path(out_dir, "early_late_correlation_table.csv"))

# =========================================================
# PART 10. Make break-year scan figure
# 目的：直观看每个变量的 RSS 在哪一年最低
# =========================================================
nice_names <- c(
  air_temp_C = "Air temperature",
  soil_temp_mean_C = "Mean soil temperature",
  soil_moist_mean = "Mean soil moisture",
  snow_depth = "Snow depth",
  snowfall_sum = "Snowfall",
  snowmelt_sum = "Snowmelt"
)

scan_results[, variable_label := factor(nice_names[variable], levels = nice_names)]
best_breaks[, variable_label := factor(nice_names[variable], levels = nice_names)]

fig_break_scan <- ggplot(scan_results, aes(x = break_year, y = rss_total)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  geom_vline(data = best_breaks, aes(xintercept = break_year),
             linetype = "dashed", colour = "red") +
  facet_wrap(~ variable_label, scales = "free_y", ncol = 2) +
  labs(
    title = "Automatic detection of break year in NDVI–climate relationships",
    subtitle = "For each candidate break year, NDVI ~ climate was fitted separately for early and late periods.\nThe selected break year is the one that minimises total RSS.",
    x = "Candidate break year",
    y = "Total residual sum of squares (RSS)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 10.5),
    axis.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold")
  )

ggsave(
  file.path(out_dir, "figure_break_year_scan.png"),
  fig_break_scan,
  width = 11,
  height = 9,
  dpi = 500
)
print(fig_break_scan)
# =========================================================
# PART 11. Make early vs late scatter plots
# 目的：直观看关系是否翻转
# =========================================================
plot_one_variable <- function(data, xvar, break_year, label_txt) {
  tmp <- copy(data)
  tmp[, period_group := ifelse(year <= break_year,
                               paste0("Early (<= ", break_year, ")"),
                               paste0("Late (> ", break_year, ")"))]
  
  ggplot(tmp, aes_string(x = xvar, y = "ndvi", colour = "period_group")) +
    geom_point(size = 2.2, alpha = 0.9) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
    labs(
      title = label_txt,
      x = label_txt,
      y = "Mean annual maximum NDVI",
      colour = NULL
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 13),
      axis.title = element_text(face = "bold"),
      legend.position = "bottom"
    )
}

scatter_list <- list()

for (i in 1:nrow(best_breaks)) {
  v <- best_breaks$variable[i]
  b <- best_breaks$break_year[i]
  lab <- nice_names[[v]]
  
  scatter_list[[i]] <- plot_one_variable(
    data = dt,
    xvar = v,
    break_year = b,
    label_txt = lab
  )
}

fig_scatter <- wrap_plots(scatter_list, ncol = 2) +
  plot_annotation(
    title = "Early vs late NDVI–climate relationships using the automatically detected break year",
    subtitle = "Each panel is split using the break year that minimises total RSS."
  )

ggsave(
  file.path(out_dir, "figure_early_late_scatter.png"),
  fig_scatter,
  width = 12,
  height = 13,
  dpi = 500
)

print(fig_scatter)

# =========================================================
# (with stats in plot) PART 11. Make early vs late scatter plots with in-panel stats
# 目的：图内直接显示 slope, r, p, R2
# =========================================================
# =========================================================
# PART 11. Make early vs late scatter plots with in-panel stats
# 目的：图内直接显示 slope, r, p, R2
# =========================================================
fmt_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}

plot_one_variable <- function(data, stats_tbl, xvar, b_year, label_txt) {
  
  tmp <- copy(data)
  tmp[, period_group := ifelse(
    year <= b_year,
    paste0("Early (<= ", b_year, ")"),
    paste0("Late (> ", b_year, ")")
  )]
  
  early_stats <- stats_tbl[
    variable == xvar & period == paste0("Early (<= ", b_year, ")")
  ]
  late_stats <- stats_tbl[
    variable == xvar & period == paste0("Late (> ", b_year, ")")
  ]
  
  early_lab <- sprintf(
    "Early (≤ %d)\nr = %.2f\np = %s\nR² = %.2f",
    b_year,
   # early_stats$slope[1],
    early_stats$r[1],
    fmt_p(early_stats$p[1]),
    early_stats$r2[1]
  )
  
  late_lab <- sprintf(
    "Late (> %d)\nr = %.2f\np = %s\nR² = %.2f",
    b_year,
    #late_stats$slope[1],
    late_stats$r[1],
    fmt_p(late_stats$p[1]),
    late_stats$r2[1]
  )
  
  x_rng <- range(tmp[[xvar]], na.rm = TRUE)
  y_rng <- range(tmp$ndvi, na.rm = TRUE)
  
  x_left  <- x_rng[1] + 0.03 * diff(x_rng)
  x_right <- x_rng[2] - 0.03 * diff(x_rng)
  y_top   <- y_rng[2] - 0.02 * diff(y_rng)
  
  ggplot(tmp, aes(x = .data[[xvar]], y = ndvi, colour = period_group)) +
    geom_point(size = 2.2, alpha = 0.9) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
    annotate(
      "text",
      x = x_left, y = y_top,
      label = early_lab,
      hjust = 0, vjust = 1,
      colour = "#F8766D", size = 3.3
    ) +
    annotate(
      "text",
      x = x_right, y = y_top,
      label = late_lab,
      hjust = 1, vjust = 1,
      colour = "#00BFC4", size = 3.3
    ) +
    scale_colour_manual(
      values = setNames(
        c("#F8766D", "#00BFC4"),
        c(
          paste0("Early (<= ", b_year, ")"),
          paste0("Late (> ", b_year, ")")
        )
      )
    ) +
    labs(
      title = label_txt,
      x = label_txt,
      y = "Mean annual maximum NDVI",
      colour = NULL
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 13),
      axis.title = element_text(face = "bold"),
      legend.position = "bottom",
      legend.text = element_text(size = 9)
    )
}

scatter_list <- list()

for (i in 1:nrow(best_breaks)) {
  v <- best_breaks$variable[i]
  b <- best_breaks$break_year[i]
  lab <- nice_names[[v]]
  
  scatter_list[[i]] <- plot_one_variable(
    data = dt,
    stats_tbl = early_late_tbl,
    xvar = v,
    b_year = b,
    label_txt = lab
  )
}

fig_scatter <- wrap_plots(scatter_list, ncol = 2) +
  plot_annotation(
    title = "Early vs late NDVI–climate relationships using the automatically detected break year",
    subtitle = "Each panel is split using the break year that minimises total RSS."
  )

ggsave(
  file.path(out_dir, "figure_early_late_scatter_simple.png"),
  fig_scatter,
  width = 13,
  height = 14,
  dpi = 500
)

print(fig_scatter)
