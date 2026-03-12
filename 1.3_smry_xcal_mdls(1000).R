# ABOUT THIS SCRIPT =========================================================================
# Summarize Landsat cross-sensor calibration model evaluation

rm(list = ls())
library(data.table)

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
xcal_dir <- file.path(output_dir, "xcal_ndvi")
outfile_smry <- file.path(output_dir, paste0(dataset_tag, "_lsat_ndvi_xcal_smry.csv"))

if (!dir.exists(xcal_dir)) {
  stop(paste0("Cross-calibration directory not found:\n", xcal_dir))
}

files <- list.files(xcal_dir, pattern = "xcal_rf_eval\\.csv$", full.names = TRUE)

if (length(files) == 0) {
  smry <- data.table(
    note = "No cross-calibration evaluation files found. Calibration may have been skipped."
  )
  print(smry)
  safe_fwrite(smry, outfile_smry, overwrite = overwrite_files)
} else {
  xcal.dt <- rbindlist(lapply(files, fread), fill = TRUE)
  
  if ("note" %in% colnames(xcal.dt)) {
    print(xcal.dt)
    safe_fwrite(xcal.dt, outfile_smry, overwrite = overwrite_files)
  } else {
    xcal.smry.dt <- xcal.dt[, .(
      rf.r2 = median(rf.r2, na.rm = TRUE),
      rf.rmse = median(rf.rmse, na.rm = TRUE),
      rf.n = median(rf.n, na.rm = TRUE),
      xval.r2 = median(xval.r2, na.rm = TRUE),
      xval.rmse = median(xval.rmse, na.rm = TRUE),
      xval.bias = median(xval.bias, na.rm = TRUE),
      xval.n = median(xval.n, na.rm = TRUE)
    ), by = "sat"]
    
    print(xcal.smry.dt)
    safe_fwrite(xcal.smry.dt, outfile_smry, overwrite = overwrite_files)
  }
}

