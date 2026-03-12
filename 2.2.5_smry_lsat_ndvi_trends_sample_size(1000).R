# ABOUT THIS SCRIPT =========================================================================
# Summarize sample-size sensitivity analysis for NDVI trends

rm(list = ls())
library(data.table)
library(reshape2)
library(ggplot2)
library(ggpubr)

base_dir <- "C:/Users/HUAWEI/Downloads/Dissertation/Juneau_NDVI"
setwd(base_dir)

source("0.1_fun_lsat_tools.R")

# =========================
# USER SETTINGS
# =========================
dataset_tag <- "juneau_ndvi_1000sites"
run_dir <- "C:/Users/HUAWEI/Downloads/Dissertation/Juneau_NDVI/runs/juneau_ndvi_1000sites"
overwrite_files <- FALSE

output_dir <- file.path(run_dir, "output")
fig_dir <- file.path(run_dir, "figures")
safe_mkdir(file.path(output_dir, "lsat_sample_size_test"), overwrite = TRUE)
safe_mkdir(fig_dir, overwrite = TRUE)

biome_file <- file.path(output_dir, "lsat_sample_size_test", "trend_mag_mcreps", paste0(dataset_tag, "_biome_trend.csv"))
site_file  <- file.path(output_dir, "lsat_sample_size_test", "trend_freq_mcreps", paste0(dataset_tag, "_site_trend_freq.csv"))
ndvi_file  <- file.path(output_dir, "ndvi_max_timeseries", paste0(dataset_tag, "_site_lsat_ndvi_max_timeseries.csv"))

outfile_smry <- file.path(output_dir, "lsat_sample_size_test", paste0(dataset_tag, "_sample_size_analysis_smry.csv"))
outfile_fig <- file.path(fig_dir, paste0(dataset_tag, "_Landsat_NDVI_trend_sample_size.jpg"))

if (!file.exists(biome_file)) stop(paste0("Biome trend MC file not found:\n", biome_file))
if (!file.exists(site_file))  stop(paste0("Site trend MC file not found:\n", site_file))
if (!file.exists(ndvi_file))  stop(paste0("NDVI max time series file not found:\n", ndvi_file))

biome.trnd.dt <- fread(biome_file)
site.trnd.dt  <- fread(site_file)

year_count <- fread(ndvi_file)[, uniqueN(year)]
biome.trnd.dt[, total.change.pcnt := (slope * year_count) / int * 100]

biome.trnd.smry.dt <- biome.trnd.dt[, .(
  slope = median(slope, na.rm = TRUE),
  slope.q025 = quantile(slope, 0.025, na.rm = TRUE),
  slope.q975 = quantile(slope, 0.975, na.rm = TRUE),
  total.change.pcnt = median(total.change.pcnt, na.rm = TRUE),
  total.change.pcnt.q025 = quantile(total.change.pcnt, 0.025, na.rm = TRUE),
  total.change.pcnt.q975 = quantile(total.change.pcnt, 0.975, na.rm = TRUE)
), by = sample.size]

biome.trnd.smry.dt[, CI95 := total.change.pcnt.q975 - total.change.pcnt.q025]

site.trnd.smry.dt <- site.trnd.dt[, .(
  pcnt.sites = median(pcnt.sites, na.rm = TRUE),
  pcnt.sites.q025 = quantile(pcnt.sites, 0.025, na.rm = TRUE),
  pcnt.sites.q975 = quantile(pcnt.sites, 0.975, na.rm = TRUE)
), by = c("sample.size","trend.cat")]

site.trnd.smry.dt[, CI95 := pcnt.sites.q975 - pcnt.sites.q025]

# summary table
biome.fancy <- copy(biome.trnd.smry.dt)
biome.fancy[, total.change.pcnt := paste0(
  sprintf("%.2f", total.change.pcnt), " [",
  sprintf("%.2f", total.change.pcnt.q025), ", ",
  sprintf("%.2f", total.change.pcnt.q975), "]"
)]
biome.fancy <- biome.fancy[, .(sample.size, total.change.pcnt)]

site.fancy <- copy(site.trnd.smry.dt)
site.fancy[, pcnt.sites := paste0(
  sprintf("%.2f", pcnt.sites), " [",
  sprintf("%.2f", pcnt.sites.q025), ", ",
  sprintf("%.2f", pcnt.sites.q975), "]"
)]
site.fancy <- dcast(site.fancy[, .(sample.size, trend.cat, pcnt.sites)],
                    sample.size ~ trend.cat, value.var = "pcnt.sites")

smry.table <- biome.fancy[site.fancy, on = "sample.size"]
safe_fwrite(smry.table, outfile_smry, overwrite = overwrite_files)

# figures
biome.trnd.fig <- ggplot(biome.trnd.smry.dt, aes(sample.size, total.change.pcnt)) +
  geom_ribbon(aes(ymin = total.change.pcnt.q025, ymax = total.change.pcnt.q975), alpha = 0.4) +
  geom_line() +
  theme_bw() +
  labs(y = "Change in mean NDVI (%)", x = "Sample size")

site.trnd.smry.plot <- site.trnd.smry.dt[trend.cat != "insig"]

site.trnd.fig <- ggplot(site.trnd.smry.plot, aes(sample.size, pcnt.sites, colour = trend.cat, fill = trend.cat)) +
  geom_ribbon(aes(ymin = pcnt.sites.q025, ymax = pcnt.sites.q975), alpha = 0.3, colour = NA) +
  geom_line() +
  theme_bw() +
  labs(y = "Percent of sites", x = "Sample size")

combo.fig <- ggarrange(biome.trnd.fig, site.trnd.fig, labels = c("(a)", "(b)"), ncol = 2, nrow = 1)

if (file.exists(outfile_fig) && !overwrite_files) {
  stop(paste0("Figure already exists and overwrite_files = FALSE:\n", outfile_fig))
}
jpeg(outfile_fig, width = 8, height = 4, units = "in", res = 400)
print(combo.fig)
dev.off()

cat("Summary finished.\n")
