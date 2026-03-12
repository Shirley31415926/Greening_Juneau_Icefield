# ABOUT THIS SCRIPT =========================================================================
# Compute NDVI trends from annual site-level NDVImax time series
# Simplified local version: site-level + overall mean only

rm(list = ls())
library(data.table)
library(dplyr)
library(tidyr)
library(zyp)
library(R.utils)

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
safe_mkdir(file.path(output_dir, "lsat_site_trends"), overwrite = TRUE)
safe_mkdir(file.path(output_dir, "lsat_mean_trends"), overwrite = TRUE)

ndvi_file <- file.path(output_dir, "ndvi_max_timeseries", paste0(dataset_tag, "_site_lsat_ndvi_max_timeseries.csv"))
outfile_mean_ts <- file.path(output_dir, "lsat_mean_trends", paste0(dataset_tag, "_mean_ndvi_timeseries.csv"))
outfile_mean_trend <- file.path(output_dir, "lsat_mean_trends", paste0(dataset_tag, "_mean_ndvi_trends.csv"))
outfile_site_trend <- file.path(output_dir, "lsat_site_trends", paste0(dataset_tag, "_site_ndvi_trends.csv"))
outfile_trend_freq <- file.path(output_dir, "lsat_site_trends", paste0(dataset_tag, "_site_ndvi_trend_frequency.csv"))

calc.trends <- function(x, y){
  keep <- complete.cases(x, y)
  x <- x[keep]
  y <- y[keep]
  
  if (length(unique(x)) < 5 || length(y) < 5) {
    return(data.table(slope = NA_real_, int = NA_real_, tau = NA_real_, pval = NA_real_))
  }
  
  xx <- zyp.yuepilon(y, x)
  return(data.table(
    slope = round(unname(xx["trend"]), 5),
    int   = round(unname(xx["intercept"]), 5),
    tau   = round(unname(xx["tau"]), 3),
    pval  = round(unname(xx["sig"]), 4)
  ))
}

# read annual NDVImax
if (!file.exists(ndvi_file)) {
  stop(paste0("NDVI max file not found:\n", ndvi_file))
}

ndvi.ts <- fread(ndvi_file, fill = TRUE)

if ("ndvi.max" %in% colnames(ndvi.ts)) {
  setnames(ndvi.ts, "ndvi.max", "ndvi")
}

ndvi.ts <- ndvi.ts[!is.na(ndvi)]

yr_min <- min(ndvi.ts$year, na.rm = TRUE)
yr_max <- max(ndvi.ts$year, na.rm = TRUE)

cat("Time range:", yr_min, "-", yr_max, "\n")
cat("Sites:", uniqueN(ndvi.ts$site), "\n")

# define periods dynamically
period_list <- list()
period_list[["full_period"]] <- copy(ndvi.ts)

if (yr_max >= 2005) {
  period_list[["post_2000"]] <- copy(ndvi.ts[year >= 2000])
}
if (yr_max >= 2015) {
  period_list[["post_2013"]] <- copy(ndvi.ts[year >= 2013])
}

ndvi.periods <- rbindlist(lapply(names(period_list), function(p) {
  tmp <- copy(period_list[[p]])
  tmp[, period := p]
  tmp[, year.rsc := year - min(year, na.rm = TRUE)]
  tmp
}), fill = TRUE)

# overall mean trend
mean.ts <- ndvi.periods[, .(
  ndvi.mean = mean(ndvi, na.rm = TRUE),
  ndvi.sd = sd(ndvi, na.rm = TRUE),
  n.sites = uniqueN(site)
), by = c("period","year","year.rsc")]

mean.trend <- mean.ts %>%
  group_by(period) %>%
  do(out = calc.trends(x = .$year.rsc, y = .$ndvi.mean)) %>%
  tidyr::unnest(cols = c(out)) %>%
  data.table()

mean.trend <- mean.trend[!is.na(slope)]
n.yrs.dt <- mean.ts[, .(n.yrs = uniqueN(year)), by = period]
mean.trend <- n.yrs.dt[mean.trend, on = "period"]
mean.trend[, total.change := slope * n.yrs]
mean.trend[, total.change.pcnt := total.change / int * 100]

safe_fwrite(mean.ts, outfile_mean_ts, overwrite = overwrite_files)
safe_fwrite(mean.trend, outfile_mean_trend, overwrite = overwrite_files)

# site-level trend
site.trend <- ndvi.periods %>%
  group_by(period, site) %>%
  do(out = calc.trends(x = .$year.rsc, y = .$ndvi)) %>%
  tidyr::unnest(cols = c(out)) %>%
  data.table()

site.trend <- site.trend[!is.na(slope)]

coord.cols <- intersect(c("site","latitude","longitude"), colnames(ndvi.ts))
if (all(c("site","latitude","longitude") %in% coord.cols)) {
  coords <- unique(ndvi.ts[, .(site, latitude, longitude)])
  site.trend <- coords[site.trend, on = "site"]
}

site.nyrs <- ndvi.periods[, .(n.yrs = uniqueN(year)), by = .(period, site)]
site.trend <- site.nyrs[site.trend, on = c("period","site")]
site.trend[, total.change := slope * n.yrs]
site.trend[, total.change.pcnt := total.change / int * 100]
site.trend[, sig := cut(pval, c(-Inf, 0.05, 0.10, Inf), c("sig.p5", "sig.p10", "insig"))]
site.trend[, slope.cat := cut(slope, c(-Inf, 0, Inf), c("browning", "greening"))]
site.trend[, trend.cat := paste(slope.cat, sig, sep = ".")]
site.trend[trend.cat %in% c("greening.insig","browning.insig"), trend.cat := "insig"]

safe_fwrite(site.trend, outfile_site_trend, overwrite = overwrite_files)

trend.freq <- site.trend[, .(n.sites = .N), by = c("period","trend.cat")]
trend.freq[, n.total := sum(n.sites), by = period]
trend.freq[, pcnt.sites := n.sites / n.total * 100]

safe_fwrite(trend.freq, outfile_trend_freq, overwrite = overwrite_files)

cat("Trend analysis finished.\n")

