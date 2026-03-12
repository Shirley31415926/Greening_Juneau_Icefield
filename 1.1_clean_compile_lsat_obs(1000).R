# ABOUT THIS SCRIPT =========================================================================
# Clean one local Landsat CSV exported from GEE and compute NDVI-ready observations
# Safe version with unique run folder and overwrite protection

rm(list = ls())
library(data.table)
library(R.utils)

base_dir <- "C:/Users/HUAWEI/Downloads/Dissertation/Juneau_NDVI"
setwd(base_dir)

source("0.1_fun_lsat_tools.R")

# =========================
# USER SETTINGS
# =========================
infile <- "C:/Users/HUAWEI/Downloads/Dissertation/Juneau_NDVI/lsat_tundra_samples_1000.csv"

dataset_tag <- "juneau_ndvi_1000sites"
run_id <- make_run_id()   # 也可以自己写，例如 "test01"
overwrite_existing_run <- FALSE
overwrite_files <- FALSE

# =========================
# CREATE SAFE RUN FOLDER
# =========================
run_dir <- make_run_dir(
  base_dir = base_dir,
  dataset_tag = dataset_tag,
  run_id = run_id,
  overwrite = overwrite_existing_run
)

clean_dir  <- file.path(run_dir, "cleaned")
output_dir <- file.path(run_dir, "output")
log_dir    <- file.path(run_dir, "logs")

outfile_clean <- file.path(clean_dir,  paste0(dataset_tag, "_cleaned_lsat.csv"))
outfile_smry  <- file.path(output_dir, paste0(dataset_tag, "_cleaning_summary.csv"))
# =========================
# READ + VALIDATE RAW CSV
# =========================
if (!file.exists(infile)) {
  stop(paste0("Input file not found:\n", infile))
}

validated_dir <- file.path(run_dir, "validated")
safe_mkdir(validated_dir, overwrite = TRUE)

validated_file <- file.path(validated_dir, paste0(dataset_tag, "_validated_raw.csv"))
badline_file   <- file.path(log_dir, paste0(dataset_tag, "_bad_lines.csv"))

lines <- readLines(infile, warn = FALSE)
if (length(lines) < 2) {
  stop("Input file appears empty or has no data rows.")
}

header <- lines[1]
header_n <- length(strsplit(header, ",", fixed = TRUE)[[1]])

field_n <- sapply(lines, function(x) length(strsplit(x, ",", fixed = TRUE)[[1]]))
good_idx <- which(field_n == header_n)
bad_idx  <- which(field_n != header_n)

cat("Total lines in raw file:", length(lines), "\n")
cat("Expected fields per line:", header_n, "\n")
cat("Good lines:", length(good_idx), "\n")
cat("Bad lines:", length(bad_idx), "\n")

if (1L %notin% good_idx) {
  stop("Header line is malformed. Please check the raw CSV.")
}

writeLines(lines[good_idx], validated_file)

if (length(bad_idx) > 0) {
  bad_dt <- data.table(
    line_number = bad_idx,
    field_count = field_n[bad_idx],
    raw_text = lines[bad_idx]
  )
  safe_fwrite(bad_dt, badline_file, overwrite = overwrite_files)
} else {
  bad_dt <- data.table(note = "No malformed lines found.")
  safe_fwrite(bad_dt, badline_file, overwrite = overwrite_files)
}

cat("Validated raw CSV written to:\n", validated_file, "\n")
cat("Bad line log written to:\n", badline_file, "\n")

lsat.samples <- fread(validated_file)

# remove GEE null records written as 0
if ("LANDSAT_PRODUCT_ID" %in% colnames(lsat.samples)) {
  lsat.samples <- lsat.samples[LANDSAT_PRODUCT_ID != 0 & !is.na(LANDSAT_PRODUCT_ID)]
}

# =========================
# PREP
# =========================
lsat.samples <- lsat_general_prep(lsat.samples)

# year and seasonal filters
lsat.samples <- lsat.samples[year >= 1984]
lsat.samples <- lsat.samples[doy >= 150 & doy <= 250]

n.sites.all <- uniqueN(lsat.samples$site)
n.obs.all   <- nrow(lsat.samples)

# remove water first
if ("pixel.qa" %in% colnames(lsat.samples)) {
  lsat.samples[, water := vapply(pixel.qa, water_flag, integer(1))]
  lsat.samples <- lsat.samples[water == 0]
}
if ("jrc.water" %in% colnames(lsat.samples)) {
  lsat.samples[, jrc.water := as.numeric(jrc.water)]
  lsat.samples <- lsat.samples[jrc.water == 0]
}

n.sites.land <- uniqueN(lsat.samples$site)
n.obs.land   <- nrow(lsat.samples)

# QAQC
lsat.samples <- lsat_qaqc_flags(lsat.samples, filter.water = FALSE)

n.sites.land.clear <- uniqueN(lsat.samples$site)
n.obs.land.clear   <- nrow(lsat.samples)

# neighborhood mean
lsat.samples <- lsat_ngb_mean(lsat.samples)
n.obs.land.clear.ngb <- nrow(lsat.samples)

# spectral index
lsat.samples <- lsat_spec_index(lsat.samples, "ndvi")

# keep only sensible NDVI
lsat.samples <- lsat.samples[!is.na(ndvi)]
lsat.samples <- lsat.samples[ndvi > -1 & ndvi < 1]

# =========================
# SAVE OUTPUTS SAFELY
# =========================
safe_fwrite(lsat.samples, outfile_clean, overwrite = overwrite_files)

clean.smry <- data.table(
  dataset_tag = dataset_tag,
  run_id = run_id,
  input_file = infile,
  output_clean_file = outfile_clean,
  n.sites.all = n.sites.all,
  n.obs.all = n.obs.all,
  n.sites.land = n.sites.land,
  n.obs.land = n.obs.land,
  n.sites.land.clear = n.sites.land.clear,
  n.obs.land.clear = n.obs.land.clear,
  n.obs.land.clear.ngb = n.obs.land.clear.ngb,
  n.obs.final = nrow(lsat.samples),
  n.sites.final = uniqueN(lsat.samples$site)
)

print(clean.smry)
safe_fwrite(clean.smry, outfile_smry, overwrite = overwrite_files)

write_run_manifest(
  run_dir,
  list(
    dataset_tag = dataset_tag,
    run_id = run_id,
    infile = infile,
    outfile_clean = outfile_clean,
    outfile_summary = outfile_smry
  ),
  overwrite = overwrite_files
)

cat("Finished cleaning.\n")
cat("Run folder:\n", run_dir, "\n")
cat("Cleaned file written to:\n", outfile_clean, "\n")

