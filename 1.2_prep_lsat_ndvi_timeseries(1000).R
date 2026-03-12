# ABOUT THIS SCRIPT =========================================================================
# Prepare Landsat NDVI time series from one cleaned local CSV
# Includes optional cross-calibration and phenology correction
# Safe version using run-specific folder

rm(list = ls())
library(data.table)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(ranger)
library(R.utils)

base_dir <- "C:/Users/HUAWEI/Downloads/Dissertation/Juneau_NDVI"
setwd(base_dir)

source("0.1_fun_lsat_tools.R")

# =========================
# USER SETTINGS
# =========================
dataset_tag <- "juneau_ndvi_1000sites"

# 把这里改成你 1.1 跑完后打印出来的 run_dir
run_dir <- "C:/Users/HUAWEI/Downloads/Dissertation/Juneau_NDVI/runs/juneau_ndvi_1000sites"

overwrite_files <- FALSE
use_sensor_uncertainty <- FALSE

clean_file <- file.path(run_dir, "cleaned", paste0(dataset_tag, "_cleaned_lsat.csv"))
output_dir <- file.path(run_dir, "output")

safe_mkdir(file.path(output_dir, "xcal_ndvi"), overwrite = TRUE)
safe_mkdir(file.path(output_dir, "pheno_curves"), overwrite = TRUE)
safe_mkdir(file.path(output_dir, "pheno_timeseries"), overwrite = TRUE)
safe_mkdir(file.path(output_dir, "pheno_max_eval"), overwrite = TRUE)
safe_mkdir(file.path(output_dir, "ndvi_max_timeseries"), overwrite = TRUE)

if (!file.exists(clean_file)) {
  stop(paste0("Cleaned input file not found:\n", clean_file))
}

# =========================
# READ CLEANED DATA
# =========================
lsat <- fread(clean_file)

lsat <- lsat[!is.na(ndvi)]
lsat <- lsat[ndvi > -1 & ndvi < 1]

cat("Number of sites:", uniqueN(lsat$site), "\n")
cat("Satellites present:", paste(sort(unique(lsat$satellite)), collapse = ", "), "\n")

# =========================
# OPTIONAL SENSOR UNCERTAINTY
# =========================
if (use_sensor_uncertainty) {
  set.seed(123)
  
  if ("LT05" %in% lsat$satellite) {
    LT05.red.scalers <- runif(nrow(lsat[satellite == "LT05"]), -0.07, 0.07)
    LT05.nir.scalers <- runif(nrow(lsat[satellite == "LT05"]), -0.07, 0.07)
    lsat[satellite == "LT05", nir := nir + nir * LT05.nir.scalers]
    lsat[satellite == "LT05", red := red + red * LT05.red.scalers]
  }
  if ("LE07" %in% lsat$satellite) {
    LE07.red.scalers <- runif(nrow(lsat[satellite == "LE07"]), -0.05, 0.05)
    LE07.nir.scalers <- runif(nrow(lsat[satellite == "LE07"]), -0.05, 0.05)
    lsat[satellite == "LE07", nir := nir + nir * LE07.nir.scalers]
    lsat[satellite == "LE07", red := red + red * LE07.red.scalers]
  }
  if ("LC08" %in% lsat$satellite) {
    LC08.red.scalers <- runif(nrow(lsat[satellite == "LC08"]), -0.03, 0.03)
    LC08.nir.scalers <- runif(nrow(lsat[satellite == "LC08"]), -0.03, 0.03)
    lsat[satellite == "LC08", nir := nir + nir * LC08.nir.scalers]
    lsat[satellite == "LC08", red := red + red * LC08.red.scalers]
  }
  
  lsat <- lsat_spec_index(lsat, "ndvi")
}

# =========================
# CROSS-CALIBRATION
# =========================
cat("Starting cross-calibration...\n")
lsat <- lsat_xcal_rf(
  dt = lsat,
  band = "ndvi",
  doy.rng = 150:250,
  min.obs = 5,
  frac.eval = 0.33,
  outfile.prefix = dataset_tag,
  outdir = file.path(output_dir, "xcal_ndvi")
)

# =========================
# LIGHT SITE FILTERING
# =========================
lsat[, ndvi.xcal.avg := mean(ndvi.xcal, na.rm = TRUE), by = site]
lsat <- lsat[ndvi.xcal.avg >= 0.05]
cat("After low-NDVI filter, sites:", uniqueN(lsat$site), "\n")

lsat[, n.yrs := uniqueN(year), by = site]
lsat <- lsat[n.yrs >= 5]
cat("After record-length filter, sites:", uniqueN(lsat$site), "\n")

lsat[, n.obs := .N, by = site]
lsat <- lsat[n.obs >= 15]
cat("After min observation filter, sites:", uniqueN(lsat$site), "\n")

if (nrow(lsat) == 0) {
  stop("No observations remain after filtering.")
}

# =========================
# PHENOLOGY CORRECTION
# =========================
cat("Starting phenology correction...\n")

year_span <- uniqueN(lsat$year)
window_yrs_use <- min(15, max(5, ifelse(year_span %% 2 == 0, year_span - 1, year_span)))

spl.outfile <- file.path(output_dir, "pheno_curves", paste0(dataset_tag, "_pheno_curves.csv"))

lsat.pheno <- lsat_pheno(
  dt = lsat,
  vi = "ndvi.xcal",
  window.yrs = window_yrs_use,
  window.min.obs = 20,
  spar = 0.7,
  spl.fit.outfile = spl.outfile
)

pheno_ts_file <- file.path(output_dir, "pheno_timeseries", paste0(dataset_tag, "_pheno_corrected_timeseries.csv"))
safe_fwrite(lsat.pheno, pheno_ts_file, overwrite = overwrite_files)

# =========================
# ANNUAL MAX
# =========================
lsat.max <- lsat_pheno_max(lsat.pheno, vi = "ndvi.xcal", min.frac.of.max = 0.75)

pheno.eval <- tryCatch({
  lsat_pheno_max_eval(
    dt = lsat.pheno,
    vi = "ndvi.xcal",
    min.frac.of.max = 0.75,
    min.obs = 11,
    reps = 10,
    outdir = file.path(output_dir, "pheno_max_eval"),
    outfile.suffix = dataset_tag
  )
}, error = function(e) {
  message("Phenology max evaluation skipped: ", e$message)
  NULL
})

# rename columns carefully
nm <- colnames(lsat.max)
nm <- gsub("ndvi.xcal.gs.med", "ndvi.gs.med", nm, fixed = TRUE)
nm <- gsub("ndvi.xcal.gs.q90", "ndvi.gs.q90", nm, fixed = TRUE)
nm <- gsub("ndvi.xcal.max.pred.min", "ndvi.max.lower", nm, fixed = TRUE)
nm <- gsub("ndvi.xcal.max.pred.max", "ndvi.max.upper", nm, fixed = TRUE)
nm <- gsub("ndvi.xcal.max.pred", "ndvi.max", nm, fixed = TRUE)
colnames(lsat.max) <- nm

# site record metadata
lsat.max[, first.yr := min(year, na.rm = TRUE), by = site]
lsat.max[, n.yr.obs := uniqueN(year), by = site]

# expand to full site-year combinations
yr_min <- min(lsat.max$year, na.rm = TRUE)
yr_max <- max(lsat.max$year, na.rm = TRUE)
full.fac <- CJ(site = unique(lsat.max$site), year = yr_min:yr_max)
lsat.max <- lsat.max[full.fac, on = c("site", "year")]

ndvi_max_file <- file.path(output_dir, "ndvi_max_timeseries", paste0(dataset_tag, "_site_lsat_ndvi_max_timeseries.csv"))
safe_fwrite(lsat.max, ndvi_max_file, overwrite = overwrite_files)

cat("All done.\n")
cat("Run folder:\n", run_dir, "\n")
cat("NDVI max file:\n", ndvi_max_file, "\n")
