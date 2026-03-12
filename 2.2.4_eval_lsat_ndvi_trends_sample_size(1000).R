# ABOUT THIS SCRIPT =========================================================================
# Assess sample-size effects on NDVI trend estimation
# Simplified for local single-file workflow

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
set.seed(123)

output_dir <- file.path(run_dir, "output")
safe_mkdir(file.path(output_dir, "lsat_sample_size_test", "trend_freq_mcreps"), overwrite = TRUE)
safe_mkdir(file.path(output_dir, "lsat_sample_size_test", "trend_mag_mcreps"), overwrite = TRUE)

ndvi_file <- file.path(output_dir, "ndvi_max_timeseries", paste0(dataset_tag, "_site_lsat_ndvi_max_timeseries.csv"))
outfile_site_freq <- file.path(output_dir, "lsat_sample_size_test", "trend_freq_mcreps", paste0(dataset_tag, "_site_trend_freq.csv"))
outfile_biome_trend <- file.path(output_dir, "lsat_sample_size_test", "trend_mag_mcreps", paste0(dataset_tag, "_biome_trend.csv"))

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

if (!file.exists(ndvi_file)) {
  stop(paste0("NDVI max time series file not found:\n", ndvi_file))
}

site.ts.dt <- fread(ndvi_file, fill = TRUE)

if ("ndvi.max" %in% colnames(site.ts.dt)) {
  setnames(site.ts.dt, "ndvi.max", "ndvi")
}

site.ts.dt <- site.ts.dt[!is.na(ndvi)]
site.ts.dt[, year.rsc := year - min(year, na.rm = TRUE)]

sites <- unique(site.ts.dt$site)
n.site.total <- length(sites)

cat("Total sites available:", n.site.total, "\n")

if (n.site.total < 10) stop("Too few sites for sample-size evaluation.")

sample.sizes <- unique(round(seq(5, n.site.total, length.out = min(10, n.site.total))))
sample.sizes <- sample.sizes[sample.sizes < n.site.total]
sample.sizes <- unique(c(sample.sizes, n.site.total))

reps <- 20

biome.trnd.list <- list()
site.trnd.list <- list()
cnt <- 1

for (j in sample.sizes) {
  for (r in 1:reps) {
    samp.sites <- sample(sites, j, replace = FALSE)
    site.ts.ss.dt <- site.ts.dt[site %in% samp.sites]
    
    biome.ts.ss.dt <- site.ts.ss.dt[, .(
      ndvi.avg = mean(ndvi, na.rm = TRUE),
      ndvi.sd = sd(ndvi, na.rm = TRUE),
      n.sites = uniqueN(site)
    ), by = c("year","year.rsc")]
    
    biome.trnd.ss.dt <- biome.ts.ss.dt %>%
      group_by() %>%
      do(out = calc.trends(x = .$year.rsc, y = .$ndvi.avg)) %>%
      tidyr::unnest(cols = c(out)) %>%
      data.table()
    
    biome.trnd.ss.dt[, `:=`(rep = r, sample.size = j)]
    biome.trnd.list[[cnt]] <- biome.trnd.ss.dt
    
    site.trnd.ss.dt <- site.ts.ss.dt %>%
      group_by(site) %>%
      do(out = calc.trends(x = .$year.rsc, y = .$ndvi)) %>%
      tidyr::unnest(cols = c(out)) %>%
      data.table()
    
    site.trnd.ss.dt <- site.trnd.ss.dt[!is.na(slope)]
    site.trnd.ss.dt[, trend.cat := cut(slope, c(-Inf, 0, Inf), c("browning","greening"))]
    site.trnd.ss.dt[pval >= 0.10001, trend.cat := "insig"]
    
    site.trnd.freq.ss.dt <- site.trnd.ss.dt[, .(n.sites = .N), by = trend.cat]
    site.trnd.freq.ss.dt[, n.sites.tot := sum(n.sites)]
    site.trnd.freq.ss.dt[, pcnt.sites := (n.sites / n.sites.tot) * 100]
    site.trnd.freq.ss.dt[, `:=`(rep = r, sample.size = j)]
    
    site.trnd.list[[cnt]] <- site.trnd.freq.ss.dt
    cnt <- cnt + 1
  }
  
  cat("Finished sample size", j, "\n")
}

biome.trnd.smry.dt <- rbindlist(biome.trnd.list, fill = TRUE)
site.trnd.smry.dt  <- rbindlist(site.trnd.list, fill = TRUE)

safe_fwrite(site.trnd.smry.dt, outfile_site_freq, overwrite = overwrite_files)
safe_fwrite(biome.trnd.smry.dt, outfile_biome_trend, overwrite = overwrite_files)

cat("Sample-size evaluation done.\n")

