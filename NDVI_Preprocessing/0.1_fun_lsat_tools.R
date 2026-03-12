# ABOUT THIS SCRIPT ================================================================================
# Custom functions for processing site-level Landsat data exported from Google Earth Engine
# Rewritten for single local CSV workflow
# Date: 2026-03-11

# =========================
# SAFE PATH HELPERS
# =========================
safe_mkdir <- function(path, overwrite = FALSE) {
  if (dir.exists(path)) {
    if (!overwrite) {
      stop(paste0("Directory already exists and overwrite = FALSE:\n", path))
    }
  } else {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  return(normalizePath(path, winslash = "/", mustWork = FALSE))
}

safe_parent_mkdir <- function(path) {
  parent <- dirname(path)
  if (!dir.exists(parent)) dir.create(parent, recursive = TRUE, showWarnings = FALSE)
}

safe_fwrite <- function(x, file, overwrite = FALSE, ...) {
  safe_parent_mkdir(file)
  if (file.exists(file) && !overwrite) {
    stop(paste0("File already exists and overwrite = FALSE:\n", file))
  }
  data.table::fwrite(x, file = file, ...)
}

make_run_id <- function() {
  format(Sys.time(), "%Y%m%d_%H%M%S")
}

make_run_dir <- function(base_dir, dataset_tag, run_id = NULL, overwrite = FALSE) {
  if (is.null(run_id)) run_id <- make_run_id()
  run_dir <- file.path(base_dir, "runs", paste0(dataset_tag, "_", run_id))
  safe_mkdir(run_dir, overwrite = overwrite)
  safe_mkdir(file.path(run_dir, "cleaned"), overwrite = TRUE)
  safe_mkdir(file.path(run_dir, "output"), overwrite = TRUE)
  safe_mkdir(file.path(run_dir, "figures"), overwrite = TRUE)
  safe_mkdir(file.path(run_dir, "logs"), overwrite = TRUE)
  return(normalizePath(run_dir, winslash = "/", mustWork = FALSE))
}

write_run_manifest <- function(run_dir, manifest_list, overwrite = FALSE) {
  manifest_path <- file.path(run_dir, "logs", "run_manifest.txt")
  if (file.exists(manifest_path) && !overwrite) {
    stop(paste0("Manifest already exists:\n", manifest_path))
  }
  lines <- paste(names(manifest_list), unlist(manifest_list), sep = ": ")
  writeLines(lines, manifest_path)
  invisible(manifest_path)
}

# =========================
# LANDSAT GENERAL PREP
# =========================
lsat_general_prep <- function(dt){
  require(data.table)
  
  dt <- data.table(dt)
  
  # standardize names
  colnames(dt) <- tolower(colnames(dt))
  colnames(dt) <- gsub("_", ".", colnames(dt))
  
  # rename key fields
  if ("landsat.product.id" %in% colnames(dt)) setnames(dt, "landsat.product.id", "landsat.id")
  if ("qa.pixel" %in% colnames(dt)) setnames(dt, "qa.pixel", "pixel.qa")
  if ("qa.radsat" %in% colnames(dt)) setnames(dt, "qa.radsat", "radsat.qa")
  if ("max.extent" %in% colnames(dt)) setnames(dt, "max.extent", "jrc.water")
  
  # create solar zenith angle
  if (!"solar.zenith.angle" %in% colnames(dt) && "sun.elevation" %in% colnames(dt)) {
    dt[, solar.zenith.angle := 90 - as.numeric(sun.elevation)]
  }
  
  # must have landsat id
  if (!"landsat.id" %in% colnames(dt)) {
    stop("landsat.id not found. Check LANDSAT_PRODUCT_ID in the CSV.")
  }
  
  # parse satellite
  dt[, satellite := substr(as.character(landsat.id), 1, 4)]
  dt[, collection := "C02"]
  
  # parse date from LANDSAT_PRODUCT_ID
  date_str <- regmatches(as.character(dt$landsat.id), regexpr("[0-9]{8}", as.character(dt$landsat.id)))
  bad_date <- is.na(date_str) | date_str == ""
  if (any(bad_date)) {
    dt <- dt[!bad_date]
    date_str <- date_str[!bad_date]
  }
  
  dates <- as.POSIXlt(date_str, format = "%Y%m%d", tz = "UTC")
  dt[, year := dates$year + 1900]
  dt[, doy := dates$yday + 1]
  
  # split by satellite
  lsat57.dt <- copy(dt[satellite %in% c("LT05", "LE07")])
  lsat8.dt  <- copy(dt[satellite %in% c("LC08", "LC09")])
  
  # rename LT05 / LE07 bands
  if (nrow(lsat57.dt) > 0) {
    if ("sr.b1" %in% colnames(lsat57.dt)) setnames(lsat57.dt, "sr.b1", "blue")
    if ("sr.b2" %in% colnames(lsat57.dt)) setnames(lsat57.dt, "sr.b2", "green")
    if ("sr.b3" %in% colnames(lsat57.dt)) setnames(lsat57.dt, "sr.b3", "red")
    if ("sr.b4" %in% colnames(lsat57.dt)) setnames(lsat57.dt, "sr.b4", "nir")
    if ("sr.b5" %in% colnames(lsat57.dt)) setnames(lsat57.dt, "sr.b5", "swir1")
    if ("sr.b7" %in% colnames(lsat57.dt)) setnames(lsat57.dt, "sr.b7", "swir2")
  }
  
  # rename LC08 / LC09 bands
  if (nrow(lsat8.dt) > 0) {
    if ("sr.b1" %in% colnames(lsat8.dt)) setnames(lsat8.dt, "sr.b1", "ublue")
    if ("sr.b2" %in% colnames(lsat8.dt)) setnames(lsat8.dt, "sr.b2", "blue")
    if ("sr.b3" %in% colnames(lsat8.dt)) setnames(lsat8.dt, "sr.b3", "green")
    if ("sr.b4" %in% colnames(lsat8.dt)) setnames(lsat8.dt, "sr.b4", "red")
    if ("sr.b5" %in% colnames(lsat8.dt)) setnames(lsat8.dt, "sr.b5", "nir")
    if ("sr.b6" %in% colnames(lsat8.dt)) setnames(lsat8.dt, "sr.b6", "swir1")
    if ("sr.b7" %in% colnames(lsat8.dt)) setnames(lsat8.dt, "sr.b7", "swir2")
  }
  
  # C2 L2 SR scaling
  scale_sr <- function(x) as.numeric(x) * 0.0000275 - 0.2
  
  sr_cols_57 <- intersect(c("blue","green","red","nir","swir1","swir2"), colnames(lsat57.dt))
  sr_cols_8  <- intersect(c("ublue","blue","green","red","nir","swir1","swir2"), colnames(lsat8.dt))
  
  if (length(sr_cols_57) > 0) {
    lsat57.dt[, (sr_cols_57) := lapply(.SD, scale_sr), .SDcols = sr_cols_57]
  }
  if (length(sr_cols_8) > 0) {
    lsat8.dt[, (sr_cols_8) := lapply(.SD, scale_sr), .SDcols = sr_cols_8]
  }
  
  dt <- rbind(lsat57.dt, lsat8.dt, fill = TRUE)
  
  # drop unneeded columns
  drop_cols <- intersect(c("landsat.id", "id", "time"), colnames(dt))
  if (length(drop_cols) > 0) dt[, (drop_cols) := NULL]
  
  wanted_cols <- c(
    "site","cid","orig.fid","latitude","longitude","jrc.water",
    "satellite","year","doy","collection",
    "solar.zenith.angle","sun.azimuth",
    "pixel.qa","radsat.qa","cloud.cover",
    "ublue","blue","green","red","nir","swir1","swir2"
  )
  wanted_cols <- wanted_cols[wanted_cols %in% colnames(dt)]
  setcolorder(dt, wanted_cols)
  
  ord_cols <- intersect(c("site","year","doy","satellite"), colnames(dt))
  if (length(ord_cols) > 0) setorderv(dt, ord_cols)
  
  return(dt)
}

# =========================
# QA_PIXEL bit helpers
# Landsat C2 L2:
# bit 0 fill
# bit 1 dilated cloud
# bit 2 cirrus
# bit 3 cloud
# bit 4 cloud shadow
# bit 5 snow
# bit 6 clear
# bit 7 water
# =========================
is_bit_set <- function(x, bit) {
  if (is.na(x)) return(FALSE)
  as.logical(bitwAnd(as.integer(x), bitwShiftL(1L, bit)))
}

clear_value <- function(x) {
  if (is.na(x)) return(0L)
  
  fill          <- is_bit_set(x, 0)
  dilated_cloud <- is_bit_set(x, 1)
  cirrus        <- is_bit_set(x, 2)
  cloud         <- is_bit_set(x, 3)
  cloud_shadow  <- is_bit_set(x, 4)
  
  if (fill || dilated_cloud || cirrus || cloud || cloud_shadow) 0L else 1L
}

snow_flag <- function(x) {
  if (is.na(x)) return(0L)
  if (is_bit_set(x, 5)) 1L else 0L
}

water_flag <- function(x) {
  if (is.na(x)) return(0L)
  if (is_bit_set(x, 7)) 1L else 0L
}

# =========================
# QAQC FLAGS
# =========================
lsat_qaqc_flags <- function(dt, cloud.max = 80, geom.max = 30, sza.max = 60,
                            filter.snow = TRUE, filter.water = TRUE){
  require(data.table)
  dt <- data.table(dt)
  n.orig <- nrow(dt)
  
  if ("pixel.qa" %in% colnames(dt)) {
    dt[, clear := vapply(pixel.qa, clear_value, integer(1))]
    dt <- dt[clear == 1]
  }
  
  if (filter.snow && "pixel.qa" %in% colnames(dt)) {
    dt[, snow := vapply(pixel.qa, snow_flag, integer(1))]
    dt <- dt[snow == 0]
  }
  
  if (filter.water) {
    if ("pixel.qa" %in% colnames(dt)) {
      dt[, water := vapply(pixel.qa, water_flag, integer(1))]
      dt <- dt[water == 0]
    }
    if ("jrc.water" %in% colnames(dt)) {
      dt[, jrc.water := as.numeric(jrc.water)]
      dt <- dt[jrc.water == 0]
    }
  }
  
  if ("cloud.cover" %in% colnames(dt)) {
    dt <- dt[as.numeric(cloud.cover) <= cloud.max]
  }
  if ("solar.zenith.angle" %in% colnames(dt)) {
    dt <- dt[as.numeric(solar.zenith.angle) <= sza.max]
  }
  if ("geometric.rmse.model" %in% colnames(dt)) {
    dt <- dt[as.numeric(geometric.rmse.model) <= geom.max]
  }
  if ("radsat.qa" %in% colnames(dt)) {
    dt <- dt[as.numeric(radsat.qa) == 0]
  }
  
  if ("blue" %in% colnames(dt))  dt <- dt[blue  > 0.005 & blue  < 1]
  if ("green" %in% colnames(dt)) dt <- dt[green > 0.005 & green < 1]
  if ("red" %in% colnames(dt))   dt <- dt[red   > 0.005 & red   < 1]
  if ("nir" %in% colnames(dt))   dt <- dt[nir   > 0.005 & nir   < 1]
  
  n.final <- nrow(dt)
  n.removed <- n.orig - n.final
  print(paste("removed", n.removed, "of", n.orig, "observations (", round(n.removed / n.orig * 100, 2), "%)"))
  
  return(dt)
}

# =========================
# NEIGHBOR OUTLIER QAQC
# =========================
zscore <- function(x){
  s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

lsat_ngb_qaqc <- function(dt, zscore.lim = 2){
  require(data.table)
  dt <- data.table(dt)
  n.orig <- nrow(dt)
  
  dt[, blue_z  := abs(zscore(blue)),  by = c("site","year","doy","satellite")]
  dt[, green_z := abs(zscore(green)), by = c("site","year","doy","satellite")]
  dt[, red_z   := abs(zscore(red)),   by = c("site","year","doy","satellite")]
  dt[, nir_z   := abs(zscore(nir)),   by = c("site","year","doy","satellite")]
  
  dt <- dt[blue_z <= zscore.lim][green_z <= zscore.lim][red_z <= zscore.lim][nir_z <= zscore.lim]
  dt[, c("blue_z","green_z","red_z","nir_z") := NULL]
  
  n.final <- nrow(dt)
  n.removed <- n.orig - n.final
  print(paste("removed", n.removed, "of", n.orig, "observations (", round(n.removed / n.orig * 100, 2), "%)"))
  
  return(dt)
}

# =========================
# NEIGHBORHOOD MEAN
# =========================
lsat_ngb_mean <- function(dt){
  require(data.table)
  dt <- data.table(dt)
  
  agg_cols <- c("latitude","longitude","ublue","blue","green","red","nir","swir1","swir2","tir")
  agg_cols <- agg_cols[agg_cols %in% colnames(dt)]
  
  dt <- dt[, lapply(.SD, mean, na.rm = TRUE),
           by = c("site","year","doy","satellite"),
           .SDcols = agg_cols]
  
  return(dt)
}

# =========================
# SPECTRAL INDICES
# =========================
lsat_spec_index <- function(dt, si){
  require(data.table)
  dt <- data.table(dt)
  
  if (si == "ndvi")  dt[, ndvi  := (nir - red) / (nir + red)]
  if (si == "nirv")  dt[, nirv  := (nir * (nir - red)) / (nir + red)]
  if (si == "evi")   dt[, evi   := 2.5 * (nir - red) / (nir + 6 * red - 7.5 * blue + 1)]
  if (si == "evi2")  dt[, evi2  := 2.5 * (nir - red) / (nir + 2.5 * red + 1)]
  if (si == "nbr" && all(c("swir1","swir2") %in% colnames(dt))) {
    dt[, nbr := (swir1 - swir2) / (swir1 + swir2)]
  }
  if (si == "ndwi")  dt[, ndwi  := (green - nir) / (green + nir)]
  if (si == "psri")  dt[, psri  := (red - blue) / nir]
  if (si == "msi" && "swir1" %in% colnames(dt)) dt[, msi := swir1 / nir]
  if (si == "ndii" && "swir1" %in% colnames(dt)) dt[, ndii := (nir - swir1) / (nir + swir1)]
  if (si == "ndvsi" && "swir1" %in% colnames(dt)) dt[, ndvsi := (swir1 - red) / (swir1 + red)]
  if (si == "satvi" && all(c("swir1","swir2") %in% colnames(dt))) {
    dt[, satvi := 1.5 * (swir1 - red) / (swir1 + red + 0.5) - swir2 / 2]
  }
  
  return(dt)
}

# =========================
# CROSS CALIBRATION
# skip if only one sensor or no LE07 or insufficient overlap
# =========================
lsat_xcal_rf <- function(dt, band, doy.rng, min.obs = 5, frac.eval = 0.33,
                         outfile.prefix = "ndvi", outdir = "output/xcal_ndvi"){
  require(data.table)
  require(R.utils)
  mkdirs(outdir)
  
  dt <- data.table(dt)
  sats <- sort(unique(dt$satellite))
  
  if (length(sats) <= 1 || !"LE07" %in% sats) {
    message("Cross-calibration skipped: only one sensor present or LE07 absent.")
    dt[, xcal := get(band)]
    setnames(dt, "xcal", paste0(band, ".xcal"))
    
    smry <- data.table(
      band = band,
      sat = paste(sats, collapse = ";"),
      note = "Cross-calibration skipped: only one sensor present or LE07 absent."
    )
    fwrite(smry, file.path(outdir, paste0(outfile.prefix, "_xcal_rf_eval.csv")))
    return(dt)
  }
  
  require(ggplot2)
  require(ggpubr)
  require(ranger)
  require(zoo)
  
  sats_to_cal <- sats[sats != "LE07"]
  dt[, xcal := NA_real_]
  
  model.eval.df <- data.frame(matrix(NA, nrow = length(sats_to_cal), ncol = 9))
  colnames(model.eval.df) <- c("band","sat","rf.r2","rf.rmse","rf.n","xval.r2","xval.rmse","xval.n","xval.bias")
  model.eval.df$sat <- sats_to_cal
  
  for (i in sats_to_cal) {
    xcal.dt <- dt[doy %in% doy.rng]
    xcal.dt <- xcal.dt[satellite %in% c(i, "LE07")]
    
    site.yr.dt <- xcal.dt[, .(year = unique(year)), by = .(site, satellite)]
    site.yr.dt <- dcast.data.table(site.yr.dt, site + year ~ satellite, value.var = "year")
    site.yr.dt <- na.omit(site.yr.dt)
    
    if (nrow(site.yr.dt) == 0) {
      message(paste("Skipping", i, "- no overlapping site-years with LE07."))
      next
    }
    
    site.yr.dt <- melt.data.table(site.yr.dt, id = "site", measure = c("LE07", i),
                                  variable.name = "satellite", value.name = "year")
    xcal.dt <- xcal.dt[site.yr.dt, on = c("site","year","satellite")]
    
    site.doy.dt <- xcal.dt[, .(n.obs = .N), by = c("site","satellite","doy")]
    full.fac <- data.table(expand.grid(site = unique(xcal.dt$site),
                                       satellite = unique(xcal.dt$satellite),
                                       doy = doy.rng))
    site.doy.dt <- site.doy.dt[full.fac, on = c("site","satellite","doy")]
    setorderv(site.doy.dt, c("site","satellite","doy"))
    site.doy.dt[is.na(n.obs), n.obs := 0]
    
    site.doy.dt <- site.doy.dt[, n.obs.15days :=
                                 zoo::rollapplyr(n.obs, FUN = sum, align = "center",
                                                 width = 15, partial = TRUE),
                               by = c("site","satellite")]
    
    site.doy.dt <- dcast.data.table(site.doy.dt, site + doy ~ satellite, value.var = "n.obs.15days")
    site.doy.dt <- site.doy.dt[get("LE07") >= min.obs & get(i) >= min.obs]
    
    if (nrow(site.doy.dt) == 0) {
      message(paste("Skipping", i, "- insufficient overlapping 15-day windows."))
      next
    }
    
    site.doy.dt <- site.doy.dt[, .SD[sample(.N, 1)], by = "site"]
    site.doy.dt[, c("LE07", i) := NULL]
    setnames(site.doy.dt, "doy", "focal.doy")
    
    site.doy.win.dt <- data.table(
      site = site.doy.dt$site,
      matrix(unlist(lapply(site.doy.dt$focal.doy, function(x) x - seq(-7, 7, 1))),
             ncol = 15, byrow = TRUE)
    )
    site.doy.win.dt <- melt.data.table(site.doy.win.dt, id.vars = "site", value.name = "doy")
    site.doy.win.dt[, variable := NULL]
    
    xcal.dt <- site.doy.win.dt[xcal.dt, on = c("site","doy"), nomatch = 0]
    xcal.dt <- xcal.dt[site.doy.dt, on = "site"]
    
    coord.dt <- xcal.dt[, .(latitude = mean(latitude, na.rm = TRUE),
                            longitude = mean(longitude, na.rm = TRUE)),
                        by = "site"]
    
    rf.dat <- xcal.dt[, .(mov.med = median(get(band), na.rm = TRUE)),
                      by = c("site","satellite","focal.doy")]
    setnames(rf.dat, "focal.doy", "doy")
    rf.dat <- dcast.data.table(rf.dat, site + doy ~ satellite, value.var = "mov.med")
    rf.dat <- rf.dat[coord.dt, on = "site"]
    setnames(rf.dat, c("LE07", i), c(paste0("LE07.", band), band))
    
    rf.dat <- na.omit(rf.dat)
    
    if (nrow(rf.dat) < 10) {
      message(paste("Skipping", i, "- too few matched samples for RF."))
      next
    }
    
    sites <- unique(rf.dat$site)
    if (length(sites) < 4) {
      message(paste("Skipping", i, "- too few sites for train/eval split."))
      next
    }
    
    sites.train <- sample(sites, max(2, round(length(sites) * (1 - frac.eval))), replace = FALSE)
    sites.eval <- sites[!sites %in% sites.train]
    
    rf.dat.train <- rf.dat[site %in% sites.train]
    rf.dat.eval  <- rf.dat[site %in% sites.eval]
    
    if (nrow(rf.dat.train) < 5) {
      message(paste("Skipping", i, "- training data too small."))
      next
    }
    
    rf.form <- formula(paste0("LE07.", band, " ~ ", band, " + doy + latitude + longitude"))
    rf.xcal <- ranger(rf.form, rf.dat.train, importance = "impurity")
    
    sat.dt <- dt[satellite == i]
    rf.pred <- predict(rf.xcal, sat.dt)
    dt[satellite == i, xcal := rf.pred$predictions]
    
    saveRDS(rf.xcal, file.path(outdir, paste0(outfile.prefix, "_", i, "_xcal_rf.Rds")))
    
    if (nrow(rf.dat.eval) > 0) {
      rf.dat.eval[, (paste0("LE07.", band, ".pred")) := predict(rf.xcal, rf.dat.eval)$predictions]
      fwrite(rf.dat.eval, file.path(outdir, paste0(outfile.prefix, "_", i, "_xcal_rf_eval_data.csv")))
      
      lm.form <- formula(paste0("LE07.", band, " ~ LE07.", band, ".pred"))
      xval.lm.smry <- summary(lm(lm.form, rf.dat.eval))
      xval.rmse <- sqrt(mean((rf.dat.eval[[paste0("LE07.", band)]] -
                                rf.dat.eval[[paste0("LE07.", band, ".pred")]])^2, na.rm = TRUE))
      xval.bias <- mean(rf.dat.eval[[paste0("LE07.", band)]] -
                          rf.dat.eval[[paste0("LE07.", band, ".pred")]], na.rm = TRUE)
      
      model.eval.df$band[model.eval.df$sat == i] <- band
      model.eval.df$rf.r2[model.eval.df$sat == i] <- round(rf.xcal$r.squared, 3)
      model.eval.df$rf.rmse[model.eval.df$sat == i] <- round(sqrt(rf.xcal$prediction.error), 5)
      model.eval.df$rf.n[model.eval.df$sat == i] <- rf.xcal$num.samples
      model.eval.df$xval.r2[model.eval.df$sat == i] <- round(xval.lm.smry$r.squared, 3)
      model.eval.df$xval.rmse[model.eval.df$sat == i] <- round(xval.rmse, 5)
      model.eval.df$xval.bias[model.eval.df$sat == i] <- round(xval.bias, 5)
      model.eval.df$xval.n[model.eval.df$sat == i] <- nrow(rf.dat.eval)
    }
  }
  
  dt[satellite == "LE07", xcal := get(band)]
  dt[is.na(xcal), xcal := get(band)]
  
  setnames(dt, "xcal", paste0(band, ".xcal"))
  fwrite(as.data.table(model.eval.df), file.path(outdir, paste0(outfile.prefix, "_xcal_rf_eval.csv")))
  
  return(dt)
}

# =========================
# PHENOLOGY
# =========================
lsat_pheno <- function(dt, vi, window.yrs = 15, window.min.obs = 30, spar = 0.7, spl.fit.outfile = NA){
  require(data.table)
  
  dt <- copy(dt)
  dt <- dt[, c("site","latitude","longitude","year","doy",vi), with = FALSE]
  setnames(dt, vi, "vi")
  setorderv(dt, c("site","doy"))
  
  all.yrs <- sort(unique(dt$year))
  
  if (length(all.yrs) < window.yrs) {
    warning("Not enough years for moving phenology window. Using all years as one window.")
    
    dt[, n.obs.focal.win := .N, by = "site"]
    dt <- dt[n.obs.focal.win >= window.min.obs]
    if (nrow(dt) == 0) stop("No sites meet window.min.obs requirement.")
    
    doy.rng <- min(dt$doy):max(dt$doy)
    
    splines.dt <- dt[, .(spl.fit = list(smooth.spline(doy, vi, spar = spar))), by = "site"]
    spline.fits.dt <- splines.dt[, .(
      spl.fit = unlist(Map(function(mod, doyvals) predict(mod, x = doyvals)$y,
                           spl.fit, MoreArgs = list(doyvals = doy.rng)))
    ), by = "site"]
    spline.fits.dt[, doy := doy.rng, by = "site"]
    spline.fits.dt[, spl.fit.max := max(spl.fit), by = "site"]
    spline.fits.dt[, vi.adjustment := abs(spl.fit - spl.fit.max)]
    spline.fits.dt[, spl.frac.max := spl.fit / spl.fit.max]
    spline.fits.dt[, pos.doy.avg := doy[which.max(spl.fit)], by = "site"]
    spline.fits.dt[, focal.yr := mean(all.yrs)]
    
    dt <- spline.fits.dt[dt, on = c("site","doy")]
    dt[, vi.max.pred := vi + vi.adjustment]
    
    if (!is.na(spl.fit.outfile)) fwrite(spline.fits.dt, spl.fit.outfile)
    setnames(dt, "vi", vi)
    return(dt)
  }
  
  focal.yrs <- all.yrs[-c(1:(round(window.yrs/2)), (length(all.yrs)-round(window.yrs/2)+1):length(all.yrs))]
  n.focal.yrs <- length(focal.yrs)
  
  data.list <- list()
  splines.list <- list()
  
  for (ii in seq_along(focal.yrs)) {
    focal.yr <- focal.yrs[ii]
    focal.win <- seq(focal.yr - round(window.yrs/2), focal.yr + round(window.yrs/2))
    focal.dt <- dt[year %in% focal.win]
    setorderv(focal.dt, c("site","doy"))
    
    focal.dt[, n.obs.focal.win := .N, by = "site"]
    focal.dt <- focal.dt[n.obs.focal.win >= window.min.obs]
    if (nrow(focal.dt) == 0) next
    
    doy.rng <- min(focal.dt$doy):max(focal.dt$doy)
    
    splines.dt <- focal.dt[, .(spl.fit = list(smooth.spline(doy, vi, spar = spar))), by = "site"]
    spline.fits.dt <- splines.dt[, .(
      spl.fit = unlist(Map(function(mod, doyvals) predict(mod, x = doyvals)$y,
                           spl.fit, MoreArgs = list(doyvals = doy.rng)))
    ), by = "site"]
    spline.fits.dt[, doy := doy.rng, by = "site"]
    
    refitting <- TRUE
    while (refitting) {
      focal.dt <- spline.fits.dt[focal.dt, on = c("site","doy")]
      focal.dt[, abs.pcnt.dif := abs((vi - spl.fit) / ((vi + spl.fit) / 2) * 100)]
      refit.sites <- unique(focal.dt[abs.pcnt.dif > 100]$site)
      focal.dt <- focal.dt[abs.pcnt.dif <= 100]
      focal.dt[, c("spl.fit","abs.pcnt.dif") := NULL]
      
      if (length(refit.sites) > 0) {
        refit.dt <- focal.dt[site %in% refit.sites]
        spline.fits.dt <- spline.fits.dt[!site %in% refit.sites]
        
        spline.refits.dt <- refit.dt[, .(spl.fit = list(smooth.spline(doy, vi, spar = spar))), by = "site"]
        spline.refits.dt <- spline.refits.dt[, .(
          spl.fit = unlist(Map(function(mod, doyvals) predict(mod, x = doyvals)$y,
                               spl.fit, MoreArgs = list(doyvals = doy.rng)))
        ), by = "site"]
        spline.refits.dt[, doy := doy.rng, by = "site"]
        
        spline.fits.dt <- rbind(spline.fits.dt, spline.refits.dt)
      } else {
        refitting <- FALSE
      }
    }
    
    site.doy.smry <- focal.dt[, .(min.doy = min(doy), max.doy = max(doy)), by = "site"]
    spline.fits.dt <- spline.fits.dt[site.doy.smry, on = "site"]
    spline.fits.dt <- spline.fits.dt[doy >= min.doy & doy <= max.doy]
    spline.fits.dt[, spl.fit.max := max(spl.fit), by = "site"]
    spline.fits.dt[, vi.adjustment := abs(spl.fit - spl.fit.max)]
    spline.fits.dt[, spl.frac.max := spl.fit / spl.fit.max]
    spline.fits.dt[, pos.doy.avg := doy[which.max(spl.fit)], by = "site"]
    spline.fits.dt[, focal.yr := focal.yr]
    
    splines.list[[length(splines.list) + 1]] <- spline.fits.dt
    
    focal.dt <- spline.fits.dt[focal.dt, on = c("site","doy")]
    focal.dt[, vi.max.pred := vi + vi.adjustment]
    
    if (ii == 1) {
      yr.win <- c((focal.yr - round(window.yrs/2)):focal.yr)
      data.list[[length(data.list) + 1]] <- focal.dt[year %in% yr.win]
    } else if (ii == n.focal.yrs) {
      yr.win <- c(focal.yr:(focal.yr + round(window.yrs/2)))
      data.list[[length(data.list) + 1]] <- focal.dt[year %in% yr.win]
    } else {
      data.list[[length(data.list) + 1]] <- focal.dt[year == focal.yr]
    }
    
    print(ii / n.focal.yrs)
  }
  
  spline.dt <- rbindlist(splines.list, fill = TRUE)
  if (!is.na(spl.fit.outfile)) fwrite(spline.dt, spl.fit.outfile)
  
  dt.out <- rbindlist(data.list, fill = TRUE)
  setorderv(dt.out, c("site","year","doy"))
  colnames(dt.out) <- gsub("vi", vi, colnames(dt.out))
  
  return(dt.out)
}

# =========================
# PHENO MAX
# =========================
lsat_pheno_max <- function(dt, vi, min.frac.of.max = 0.75){
  require(data.table)
  
  dt <- copy(dt)
  colnames(dt) <- gsub(vi, "vi", colnames(dt), fixed = TRUE)
  dt <- dt[spl.frac.max >= min.frac.of.max]
  
  dt[, `:=`(
    avg = mean(vi.max.pred, na.rm = TRUE),
    sd  = sd(vi.max.pred, na.rm = TRUE),
    n   = .N
  ), by = c("site","year")]
  
  dt[, abs.zscore := abs((vi.max.pred - avg) / sd)]
  dt <- dt[is.na(abs.zscore) | abs.zscore < 2]
  
  dt.smry <- dt[, .(
    latitude        = first(latitude),
    longitude       = first(longitude),
    n.obs           = .N,
    vi.gs.med       = median(vi, na.rm = TRUE),
    vi.gs.q90       = quantile(vi, 0.9, na.rm = TRUE),
    vi.max.pred     = mean(vi.max.pred, na.rm = TRUE),
    vi.max.pred.min = min(vi.max.pred, na.rm = TRUE),
    vi.max.pred.max = max(vi.max.pred, na.rm = TRUE),
    pos.doy.avg     = mean(pos.doy.avg, na.rm = TRUE)
  ), by = c("site","year")]
  
  dt.smry[vi.max.pred < vi.gs.q90, vi.max.pred := vi.gs.q90]
  colnames(dt.smry) <- gsub("vi", vi, colnames(dt.smry), fixed = TRUE)
  
  return(dt.smry)
}

# =========================
# PHENO MAX EVALUATION
# =========================
lsat_pheno_max_eval <- function(dt, vi, min.frac.of.max = 0.75, min.obs = 10, reps = 10,
                                outdir = "output/pheno_max_eval", outfile.suffix = "run"){
  require(data.table)
  require(ggplot2)
  require(R.utils)
  
  colnames(dt) <- gsub(vi, "vi", colnames(dt), fixed = TRUE)
  dt <- dt[spl.frac.max >= min.frac.of.max]
  dt[, n.obs.gs := .N, by = c("site","year")]
  dt <- dt[n.obs.gs > min.obs]
  
  if (nrow(dt) == 0) {
    warning("No site-years meet min.obs for pheno_max_eval.")
    return(NULL)
  }
  
  dt[, `:=`(
    avg = mean(vi.max.pred, na.rm = TRUE),
    sd  = sd(vi.max.pred, na.rm = TRUE),
    n   = .N
  ), by = c("site","year")]
  
  dt[, abs.zscore := abs((vi.max.pred - avg) / sd)]
  dt <- dt[is.na(abs.zscore) | abs.zscore < 2]
  dt[, vi.max.obs := quantile(vi, 0.90, na.rm = TRUE), by = c("site","year")]
  
  out.list <- list()
  cnt <- 1
  
  for (i in 1:(min.obs - 1)) {
    for (j in 1:reps) {
      rep.dt <- dt[, .SD[sample(.N, i)], by = c("site","year")]
      rep.dt[, `:=`(n.obs = i, rep = j)]
      
      rep.dt <- rep.dt[, .(
        vi.max.obs   = first(vi.max.obs),
        vi.max.uncor = quantile(vi, 0.9, na.rm = TRUE),
        vi.max.cor   = median(vi.max.pred, na.rm = TRUE)
      ), by = c("n.obs","rep","site","year")]
      
      rep.dt[vi.max.cor < vi.max.uncor, vi.max.cor := vi.max.uncor]
      rep.dt[, vi.uncor.pcntdif := (vi.max.uncor - vi.max.obs) / vi.max.obs * 100]
      rep.dt[, vi.cor.pcntdif   := (vi.max.cor   - vi.max.obs) / vi.max.obs * 100]
      
      out.list[[cnt]] <- rep.dt
      cnt <- cnt + 1
    }
  }
  
  eval.dt <- rbindlist(out.list)
  eval.dt <- na.omit(eval.dt)
  
  eval.smry.dt <- eval.dt[, .(
    vi.uncor.pcntdif.med = median(vi.uncor.pcntdif, na.rm = TRUE),
    vi.cor.pcntdif.med   = median(vi.cor.pcntdif,   na.rm = TRUE)
  ), by = c("site","year","n.obs")]
  
  eval.smry.dt <- melt.data.table(
    eval.smry.dt,
    id.vars = c("site","year","n.obs"),
    value.name = "pcnt.dif",
    variable.name = "correction"
  )
  
  eval.smry.dt[, n.obs.fac := as.factor(n.obs)]
  eval.smry.dt[, correction := factor(correction, labels = c("Raw","Corrected"))]
  
  ylab.pcntdif <- paste0("Difference from observed ", toupper(gsub(".xcal", "", vi, fixed = TRUE)), "max (%)")
  
  fig <- ggplot(eval.smry.dt, aes(n.obs.fac, pcnt.dif, fill = correction)) +
    geom_boxplot(outlier.size = 0.7, outlier.color = "gray") +
    theme_bw() +
    labs(y = ylab.pcntdif, x = "Number of observations") +
    theme(
      legend.position = "right",
      axis.text  = element_text(size = 12),
      axis.title = element_text(size = 14, face = "bold")
    )
  
  colnames(eval.smry.dt) <- gsub("vi", vi, colnames(eval.smry.dt), fixed = TRUE)
  mkdirs(outdir)
  
  fig.outname <- file.path(outdir, paste0("pheno_max_eval_", outfile.suffix, ".jpg"))
  jpeg(fig.outname, width = 6, height = 4, res = 400, units = "in")
  print(fig)
  dev.off()
  
  data.outname <- file.path(outdir, paste0("pheno_max_eval_", outfile.suffix, ".csv"))
  fwrite(eval.smry.dt, data.outname)
  
  return(eval.smry.dt)
}
