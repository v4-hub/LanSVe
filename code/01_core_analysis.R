#!/usr/bin/env Rscript
# =============================================================================
# 01_core_analysis.R — Core breakpoint detection and bootstrap confidence
#                       intervals for the vegetation-landslide relationship
#
# Paper: "Landslide susceptibility peaks at intermediate vegetation density:
#         a global satellite analysis of critical NDVI thresholds"
#
# Reproduces the two primary results:
#   1. Two-breakpoint segmented regression on pre-failure NDVI (BP1, BP2)
#      identifying the Critical Transition Zone (CTZ) — this is the paper's
#      central claim.
#   2. BCa bootstrap 95% CIs (10,000 resamples) for both breakpoints.
#   3. Zone distribution counts (IDZ / CTZ / SDZ).
#
# Outputs (written to models_rev/ for downstream scripts):
#   models_rev/ndvi_seg.rds            segmented-regression object
#   models_rev/bootstrap_ci_ndvi.csv   bootstrap BCa CIs for BP1 and BP2
#   models_rev/bootstrap_ndvi.rds      full boot object (for audit)
#   models_rev/zone_distribution.csv   event counts per zone
#
# Also writes:
#   supplementary/Table_S1.csv         breakpoint estimates + CIs (final table)
#
# Run:
#   Rscript code/01_core_analysis.R
#
# Requires: segmented, boot, dplyr
# =============================================================================
suppressMessages({library(segmented); library(boot); library(dplyr)})
set.seed(123)

DATA <- "data/landslide_data_with_terrain.csv"
d <- read.csv(DATA, stringsAsFactors = FALSE)
cat(sprintf("Loaded %d events\n", nrow(d)))

dir.create("models_rev", showWarnings = FALSE)
dir.create("supplementary", showWarnings = FALSE)

# ── Segmented-regression helper (histogram + Tukey 3R/S/3R smoothing) ─────
# This is the exact method described in the Methods section.
advanced_segmented_analysis <- function(data, var_name, breaks = 2) {
  classify_distribution <- function(y_values) {
    n <- length(y_values); q1 <- floor(n/4); q3 <- floor(3*n/4)
    sm <- mean(y_values[1:q1], na.rm=TRUE)
    mm <- mean(y_values[(q1+1):(q3-1)], na.rm=TRUE)
    em <- mean(y_values[q3:n], na.rm=TRUE)
    if (mm > sm*1.1 && mm > em*1.1) return("peaked")
    else if (mm < sm*0.9 && mm < em*0.9) return("valley")
    else return("other")
  }
  valid_data <- data[!is.na(data[[var_name]]), ]
  x_vals <- valid_data[[var_name]]
  n_bins <- min(150, floor(length(x_vals)/20))
  x_range <- range(x_vals, na.rm=TRUE)
  breaks_seq <- seq(x_range[1], x_range[2], length.out=n_bins+1)
  h <- hist(x_vals, breaks=breaks_seq, plot=FALSE)
  counts <- h$counts; mid_points <- h$mids
  s1 <- stats::smooth(counts, kind="3R")
  s2 <- stats::smooth(s1, kind="S")
  counts_final <- stats::smooth(s2, kind="3R")
  df <- data.frame(x=mid_points, y=as.numeric(counts_final))
  lin_mod <- lm(y ~ x, data=df)
  shape <- classify_distribution(df$y)
  if (shape=="peaked") {
    pk <- df$x[which.max(df$y)]; q10<-quantile(df$x,.10); q90<-quantile(df$x,.90)
    initial_psi <- c((q10+pk)/2,(pk+q90)/2)
  } else if (shape=="valley") {
    vl <- df$x[which.min(df$y)]; q25<-quantile(df$x,.25); q75<-quantile(df$x,.75)
    initial_psi <- c((q25+vl)/2,(vl+q75)/2)
  } else initial_psi <- quantile(df$x, probs=seq(.25,.75,length.out=breaks))
  initial_psi <- pmax(x_range[1]+0.05*diff(x_range),
                      pmin(x_range[2]-0.05*diff(x_range), initial_psi))
  seg_mod <- segmented::segmented(lin_mod, seg.Z=~x, psi=list(x=initial_psi))
  bp <- seg_mod$psi[,"Est."]
  ci <- confint(seg_mod); ci_rows <- grep("psi", rownames(ci))
  thr <- ci[ci_rows,,drop=FALSE]
  list(breakpoints=bp, confidence_intervals=thr, data=df,
       r_squared=summary(seg_mod)$r.squared, shape=shape, n=length(x_vals))
}

# ── Run segmented regression for NDVI, EVI, LAI ───────────────────────────
for (v in c("NDVI_1","EVI_1","LAI_1")) {
  r <- tryCatch(advanced_segmented_analysis(d, v),
                error=function(e){cat("ERR",v,e$message,"\n"); NULL})
  if (!is.null(r)) {
    cat(sprintf("\n=== %s (n=%d, shape=%s) ===\n", v, r$n, r$shape))
    cat(sprintf("  R^2 = %.4f\n", r$r_squared))
    cat(sprintf("  BP1 = %.4f  CI[%.4f, %.4f]\n",
        r$breakpoints[1], r$confidence_intervals[1,2], r$confidence_intervals[1,3]))
    cat(sprintf("  BP2 = %.4f  CI[%.4f, %.4f]\n",
        r$breakpoints[2], r$confidence_intervals[2,2], r$confidence_intervals[2,3]))
    if (v=="NDVI_1") saveRDS(r, "models_rev/ndvi_seg.rds")
  }
}

# ── Bootstrap BCa CIs for NDVI breakpoints (10,000 resamples) ─────────────
ndvi <- readRDS("models_rev/ndvi_seg.rds")
bdata <- ndvi$data

boot_fn <- function(data, indices) {
  dd <- data[indices,]
  tryCatch({
    lm0 <- lm(y~x, data=dd)
    sm  <- segmented::segmented(lm0, seg.Z=~x, npsi=2)
    sort(sm$psi[,"Est."])
  }, error=function(e) rep(NA,2))
}

set.seed(42)
ncpu <- min(parallel::detectCores()-1, 4)
cat("\nRunning 10,000 bootstrap resamples (this takes ~2 min)...\n")
br <- tryCatch(
  boot::boot(bdata, boot_fn, R=10000, parallel="multicore", ncpus=ncpu),
  error=function(e) boot::boot(bdata, boot_fn, R=10000)
)

# BCa intervals
ci1 <- tryCatch(boot::boot.ci(br, type="bca", index=1)$bca[4:5], error=function(e) c(NA,NA))
ci2 <- tryCatch(boot::boot.ci(br, type="bca", index=2)$bca[4:5], error=function(e) c(NA,NA))

boot_ci <- data.frame(
  Breakpoint = 1:2,
  Estimate   = ndvi$breakpoints,
  CI_Lower   = c(ci1[1], ci2[1]),
  CI_Upper   = c(ci1[2], ci2[2]),
  Boot_SE    = apply(br$t, 2, sd, na.rm=TRUE),
  Success_Rate = apply(br$t, 2, function(x) mean(!is.na(x)))
)

cat("\nBootstrap BCa CIs:\n")
print(boot_ci, row.names=FALSE)

saveRDS(br, "models_rev/bootstrap_ndvi.rds")
write.csv(boot_ci, "models_rev/bootstrap_ci_ndvi.csv", row.names=FALSE)
write.csv(boot_ci, "supplementary/Table_S1.csv",        row.names=FALSE)
cat("-> supplementary/Table_S1.csv  (breakpoints + bootstrap CIs)\n")

# ── Zone distribution ──────────────────────────────────────────────────────
BP1 <- boot_ci$Estimate[1]; BP2 <- boot_ci$Estimate[2]
nv  <- d$NDVI_1[!is.na(d$NDVI_1)]
N_analyzed <- length(nv); N_full <- nrow(d)

zone_dist <- data.frame(
  Zone        = c("IDZ","CTZ","SDZ","Missing_NDVI"),
  N_events    = c(sum(nv<BP1), sum(nv>=BP1&nv<=BP2), sum(nv>BP2), N_full-N_analyzed),
  Pct_analyzed= c(round(100*sum(nv<BP1)/N_analyzed,1),
                  round(100*sum(nv>=BP1&nv<=BP2)/N_analyzed,1),
                  round(100*sum(nv>BP2)/N_analyzed,1), NA),
  Pct_full    = c(round(100*sum(nv<BP1)/N_full,1),
                  round(100*sum(nv>=BP1&nv<=BP2)/N_full,1),
                  round(100*sum(nv>BP2)/N_full,1),
                  round(100*(N_full-N_analyzed)/N_full,1))
)
cat(sprintf("\nZone distribution (N_analyzed=%d, N_full=%d):\n", N_analyzed, N_full))
print(zone_dist, row.names=FALSE)
write.csv(zone_dist, "models_rev/zone_distribution.csv", row.names=FALSE)

cat(sprintf("\n=== CORE ANALYSIS COMPLETE ===\n"))
cat(sprintf("  BP1 = %.3f  [%.3f, %.3f] BCa 95%%\n", BP1, ci1[1], ci1[2]))
cat(sprintf("  BP2 = %.3f  [%.3f, %.3f] BCa 95%%\n", BP2, ci2[1], ci2[2]))
cat(sprintf("  CTZ: %.1f%% of analysed events (%.1f%% of full inventory)\n",
    zone_dist$Pct_analyzed[2], zone_dist$Pct_full[2]))
cat("DONE. Authoritative numbers saved to models_rev/\n")
