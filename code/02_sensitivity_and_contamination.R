#!/usr/bin/env Rscript
# =============================================================================
# 02_sensitivity_and_contamination.R — Robustness checks and contamination test
#
# Paper: "Landslide susceptibility peaks at intermediate vegetation density:
#         a global satellite analysis of critical NDVI thresholds"
#
# This script provides two classes of validation:
#
# PART A: Breakpoint sensitivity analyses (Reviewer-requested robustness checks)
#   A2  Cloud/QA sensitivity + Sentinel-2 recovery of missing MODIS NDVI
#   A3  Binning/bandwidth sensitivity (histogram bins vs. KDE methods)
#   A4  MODIS NDVI sensor-noise propagation (sigma = 0.02–0.05)
#   A5  Phenology control (within-event trajectory baseline comparison)
#   A6  Pixel RF: train vs. CV R^2 (overfitting transparency)
#
# PART B: Post-failure compositing contamination test
#   Shows that the apparent pre-failure NDVI "early-warning" decline is an
#   artefact of the terminal 16-day MOD13Q1 composite straddling the failure
#   date (median 9 days post-event). With clean pre-failure windows the
#   decline is not significant (IDZ p=0.19, CTZ p=0.34).
#   The central non-monotonic breakpoint result is unaffected.
#
# Must run AFTER 01_core_analysis.R (reads models_rev/ndvi_seg.rds).
#
# Outputs:
#   supplementary/Table_S12.csv   pixel RF train vs. CV R^2
#   supplementary/Table_S4.csv   binning/bandwidth sensitivity
#   supplementary/Table_S5.csv   sensor-noise propagation
#   supplementary/Table_S6.csv   compositing contamination test (Part B)
#   models_rev/                   intermediate tables for 04_supplementary_figures.R
#
# Note: Table_S3 (subset breakpoints) is written by 04_supplementary_figures.R.
#
# Run:
#   Rscript code/02_sensitivity_and_contamination.R
#
# Requires: segmented, randomForest, dplyr, scales
# =============================================================================
suppressMessages({library(segmented); library(randomForest); library(dplyr)})
set.seed(123)

d <- read.csv("data/landslide_data_with_terrain.csv", stringsAsFactors=FALSE)
dir.create("models_rev", showWarnings=FALSE)
dir.create("supplementary", showWarnings=FALSE)
out <- function(...) cat(sprintf(...))

# Load authoritative breakpoints from 01_core_analysis.R
bp_auth <- read.csv("models_rev/bootstrap_ci_ndvi.csv")
BP1 <- bp_auth$Estimate[1]; BP2 <- bp_auth$Estimate[2]
cat(sprintf("Authoritative breakpoints: BP1=%.4f BP2=%.4f\n", BP1, BP2))

zone_of <- function(v) ifelse(v<BP1,"IDZ", ifelse(v<=BP2,"CTZ","SDZ"))

# ── Shared segmented helper (histogram + Tukey 3R/S/3R) ───────────────────
seg_bp <- function(x_vals, n_bins=NULL) {
  x_vals <- x_vals[!is.na(x_vals)]
  if (is.null(n_bins)) n_bins <- min(150, floor(length(x_vals)/20))
  rng <- range(x_vals); bs <- seq(rng[1], rng[2], length.out=n_bins+1)
  h <- hist(x_vals, breaks=bs, plot=FALSE)
  cf <- stats::smooth(stats::smooth(stats::smooth(h$counts,"3R"),"S"),"3R")
  df <- data.frame(x=h$mids, y=as.numeric(cf)); lm0 <- lm(y~x, data=df)
  psi0 <- quantile(df$x, c(.4,.8))
  psi0 <- pmax(rng[1]+0.05*diff(rng), pmin(rng[2]-0.05*diff(rng), psi0))
  sm <- tryCatch(segmented(lm0, seg.Z=~x, psi=list(x=psi0)), error=function(e) NULL)
  if (is.null(sm)) return(c(NA,NA,NA))
  bp <- sort(sm$psi[,"Est."]); c(bp[1], bp[2], summary(sm)$r.squared)
}

B0 <- seg_bp(d$NDVI_1)
out("Baseline NDVI check: BP1=%.4f BP2=%.4f R2=%.4f\n", B0[1], B0[2], B0[3])


# =============================================================================
# PART A: SENSITIVITY ANALYSES
# =============================================================================

# ─────────────────────────────────────────────────────────────────────────────
# A2. Cloud/QA sensitivity + Sentinel-2 recovery
# ─────────────────────────────────────────────────────────────────────────────
out("\n===== A2: CLOUD-QA SENSITIVITY =====\n")
qa_tab <- table(d$QA_1[!is.na(d$NDVI_1)], useNA="ifany"); print(qa_tab)

s14a <- data.frame()
for (lab in c("QA0_strict","QA01_reliable_current","all_nonNA")) {
  sub <- switch(lab,
    QA0_strict             = d$NDVI_1[!is.na(d$NDVI_1) & d$QA_1==0],
    QA01_reliable_current  = d$NDVI_1[!is.na(d$NDVI_1) & d$QA_1 %in% c(0,1)],
    all_nonNA              = d$NDVI_1[!is.na(d$NDVI_1)])
  b <- seg_bp(sub)
  s14a <- rbind(s14a, data.frame(Subset=lab, N=length(sub), BP1=b[1], BP2=b[2], R2=b[3]))
  out("  %-22s n=%-6d BP1=%.4f BP2=%.4f R2=%.4f\n", lab, length(sub), b[1], b[2], b[3])
}

both <- d[!is.na(d$NDVI_1) & !is.na(d$S2_NDVI_1) & d$S2_NDVI_1>=-1 & d$S2_NDVI_1<=1, ]
cal <- lm(NDVI_1 ~ S2_NDVI_1, data=both)
out("  S2->MODIS calibration: MODIS = %.3f + %.3f*S2  (n=%d, R2=%.3f)\n",
    coef(cal)[1], coef(cal)[2], nrow(both), summary(cal)$r.squared)
miss <- d[is.na(d$NDVI_1) & !is.na(d$S2_NDVI_1) & d$S2_NDVI_1>=-1 & d$S2_NDVI_1<=1, ]
miss$NDVI_recovered <- pmin(1, pmax(0, predict(cal, miss)))
rec_zone  <- zone_of(miss$NDVI_recovered)
main_zone <- zone_of(d$NDVI_1[!is.na(d$NDVI_1)])
s14b <- data.frame(
  Zone=c("IDZ","CTZ","SDZ"),
  Main_pct      = as.numeric(prop.table(table(factor(main_zone,c("IDZ","CTZ","SDZ"))))*100),
  Recovered_pct = as.numeric(prop.table(table(factor(rec_zone,c("IDZ","CTZ","SDZ"))))*100))
out("  Recovered MODIS-missing events via S2 (n=%d):\n", nrow(miss))
print(s14b, row.names=FALSE)
ma <- miss$Country_Name %in% c("Bangladesh","Vietnam","Myanmar","India","Nepal",
                                "Philippines","Indonesia")
out("  ...of which %d in monsoon Asia; their CTZ+high share (NDVI>=BP1) = %.1f%%\n",
    sum(ma), 100*mean(miss$NDVI_recovered[ma] >= BP1))
write.csv(s14a, "models_rev/TableS14a_QA_sensitivity.csv", row.names=FALSE)
write.csv(s14b, "models_rev/TableS14b_S2_recovery.csv",    row.names=FALSE)

# ─────────────────────────────────────────────────────────────────────────────
# A3. Binning / bandwidth sensitivity
# ─────────────────────────────────────────────────────────────────────────────
out("\n===== A3: BINNING / BANDWIDTH SENSITIVITY =====\n")
s15 <- data.frame()
for (nb in c(75,100,125,150,175,200)) {
  b <- seg_bp(d$NDVI_1, n_bins=nb)
  s15 <- rbind(s15, data.frame(Method=sprintf("hist_%dbins",nb), BP1=b[1],BP2=b[2],R2=b[3]))
}
nv <- d$NDVI_1[!is.na(d$NDVI_1)]
for (bw in c("SJ","nrd0")) for (mult in c(0.5,1,2)) {
  k <- density(nv, from=0, to=1, n=512, bw=bw, adjust=mult)
  df2 <- data.frame(x=k$x, y=k$y); lm0 <- lm(y~x, data=df2)
  sm <- tryCatch(segmented(lm0, seg.Z=~x, psi=list(x=quantile(df2$x,c(.4,.8)))),
                 error=function(e) NULL)
  if (!is.null(sm)) {
    bp <- sort(sm$psi[,"Est."])
    s15 <- rbind(s15, data.frame(Method=sprintf("KDE_%s_x%.1f",bw,mult),
                                  BP1=bp[1], BP2=bp[2], R2=summary(sm)$r.squared))
  }
}
print(s15, row.names=FALSE)
out("  BP1 range [%.3f, %.3f]  BP2 range [%.3f, %.3f]\n",
    min(s15$BP1,na.rm=TRUE), max(s15$BP1,na.rm=TRUE),
    min(s15$BP2,na.rm=TRUE), max(s15$BP2,na.rm=TRUE))
write.csv(s15, "models_rev/TableS15_bandwidth_sensitivity.csv", row.names=FALSE)
write.csv(s15, "supplementary/Table_S4.csv",                   row.names=FALSE)
cat("-> supplementary/Table_S4.csv\n")

# ─────────────────────────────────────────────────────────────────────────────
# A4. MODIS NDVI sensor-noise propagation (sigma = 0.00 / 0.02 / 0.03 / 0.05)
# ─────────────────────────────────────────────────────────────────────────────
out("\n===== A4: SENSOR-NOISE PROPAGATION =====\n")
s16 <- data.frame()
for (sg in c(0.00,0.02,0.03,0.05)) {
  bp1s <- c(); bp2s <- c()
  ntrial <- if (sg==0) 1 else 200
  for (t in 1:ntrial) {
    xn <- nv + rnorm(length(nv), 0, sg)
    b  <- seg_bp(pmin(1,pmax(0,xn)))
    bp1s <- c(bp1s,b[1]); bp2s <- c(bp2s,b[2])
  }
  s16 <- rbind(s16, data.frame(Sigma=sg,
    BP1_mean=mean(bp1s,na.rm=TRUE), BP1_sd=sd(bp1s,na.rm=TRUE),
    BP2_mean=mean(bp2s,na.rm=TRUE), BP2_sd=sd(bp2s,na.rm=TRUE)))
  out("  sigma=%.2f  BP1=%.4f+/-%.4f  BP2=%.4f+/-%.4f\n",
      sg, mean(bp1s,na.rm=TRUE), sd(bp1s,na.rm=TRUE),
      mean(bp2s,na.rm=TRUE), sd(bp2s,na.rm=TRUE))
}
write.csv(s16, "models_rev/TableS16_noise_propagation.csv", row.names=FALSE)
write.csv(s16, "supplementary/Table_S5.csv",               row.names=FALSE)
cat("-> supplementary/Table_S5.csv\n")

# ─────────────────────────────────────────────────────────────────────────────
# A5. Phenology control (within-event baseline anomaly)
# ─────────────────────────────────────────────────────────────────────────────
out("\n===== A5: PHENOLOGY CONTROL =====\n")
dd <- d[!is.na(d$NDVI_1), ]; dd$Zone <- zone_of(dd$NDVI_1)
dd$term_slope  <- dd$NDVI_1 - dd$NDVI_2
dd$early_slope <- (dd$NDVI_3 - dd$NDVI_5)/2
s17_baseline <- dd %>% group_by(Zone) %>% summarise(
  n=sum(!is.na(term_slope)&!is.na(early_slope)),
  terminal_slope=mean(term_slope,na.rm=TRUE),
  early_baseline_slope=mean(early_slope,na.rm=TRUE),
  anomaly=mean(term_slope-early_slope,na.rm=TRUE),
  p_wilcox=tryCatch(wilcox.test(term_slope, early_slope, paired=TRUE)$p.value,
                    error=function(e) NA),
  .groups="drop")
print(as.data.frame(s17_baseline), row.names=FALSE)
write.csv(as.data.frame(s17_baseline), "models_rev/TableS17_baseline_anomaly.csv", row.names=FALSE)

# ─────────────────────────────────────────────────────────────────────────────
# A6. Pixel RF: train vs. spatial-CV R^2 (overfitting check)
# ─────────────────────────────────────────────────────────────────────────────
out("\n===== A6: PIXEL RF TRAIN vs SPATIAL-CV R^2 =====\n")
k <- density(nv, from=0, to=1, n=101, bw="SJ")
dens <- data.frame(g=round(k$x,2), dt=scales::rescale(k$y))
feat <- c("EVI_1","LAI_1","MODIS_IGBP_Simplified","Climate_Zone","Season",
          "NDVI_change_1_to_2","EVI_change_1_to_2","LAI_change_1_to_2",
          "Hansen_Tree_Cover_2000_Percent")
md <- d %>% filter(!is.na(NDVI_1)) %>% mutate(g=round(NDVI_1,2)) %>%
  left_join(dens, by="g") %>% filter(!is.na(dt)) %>%
  dplyr::select(all_of(feat), dt, Latitude) %>% na.omit() %>%
  mutate(across(where(is.character), as.factor))
set.seed(123)
ql <- unique(quantile(md$Latitude, 0:5/5))
md$Fold <- as.numeric(cut(md$Latitude, ql, include.lowest=TRUE))
tr <- c(); te <- c()
for (kf in 1:5) {
  trn <- md[md$Fold!=kf,]; tst <- md[md$Fold==kf,]
  rf  <- randomForest(as.formula(paste("dt ~", paste(feat,collapse="+"))),
                      data=trn, ntree=200)
  tr <- c(tr, cor(trn$dt, predict(rf,trn))^2)
  te <- c(te, cor(tst$dt, predict(rf,tst))^2)
}
out("  Pixel RF: TRAIN R2=%.3f+/-%.3f   CV R2=%.3f+/-%.3f   gap=%.3f\n",
    mean(tr),sd(tr),mean(te),sd(te),mean(tr)-mean(te))
s13 <- data.frame(metric=c("train_R2","spatial_cv_R2"),
                  mean=c(mean(tr),mean(te)), sd=c(sd(tr),sd(te)))
write.csv(s13, "models_rev/TableS13_trainvtest_pixel.csv", row.names=FALSE)
write.csv(s13, "supplementary/Table_S12.csv",              row.names=FALSE)
cat("-> supplementary/Table_S12.csv\n")


# =============================================================================
# PART B: POST-FAILURE COMPOSITING CONTAMINATION TEST
# =============================================================================
# The terminal MODIS window ("0–16 d before failure") is a 16-day MOD13Q1
# composite. For 96% of events it extends past the failure date (median 9 days
# post-event) and captures the fresh landslide scar, inflating the apparent
# zone-ordered NDVI decline (-5.4 / -2.1 / +2.1 %). Restricted to windows
# that are unambiguously pre-failure (NDVI_2 – NDVI_3), the decline
# disappears (IDZ p=0.19, CTZ p=0.34). The breakpoints are robust to window.
# =============================================================================
out("\n\n===== PART B: COMPOSITING CONTAMINATION TEST =====\n")

# B1. How far does the terminal composite straddle the failure date?
db1  <- d$Days_Before_1[!is.na(d$Days_Before_1)]
db2  <- d$Days_Before_2[!is.na(d$Days_Before_2)]
post <- pmax(0, 16 - db1)
out("  Days_Before_1: median=%.1f  mean=%.1f\n", median(db1), mean(db1))
out("  Post-event days in window-1 composite: median=%.1f  frac>0=%.1f%%\n",
    median(post), 100*mean(post>0))
out("  Window 2 (16-32 d): %.1f%% fully pre-event\n", 100*mean(db2>=16))

# B2. Apparent vs. clean-window terminal NDVI step
m <- d[!is.na(d$NDVI_1),]; m$Zone <- zone_of(m$NDVI_1)
tab <- data.frame()
out("\n  Zone | apparent step (NDVI1-NDVI2) | clean step (NDVI2-NDVI3) | p\n")
for (z in c("IDZ","CTZ","SDZ")) {
  s <- m[m$Zone==z,]
  app_pct  <- 100*mean((s$NDVI_1-s$NDVI_2)/pmax(0.01,abs(s$NDVI_1)), na.rm=TRUE)
  app_slp  <- mean(s$NDVI_1-s$NDVI_2, na.rm=TRUE)
  term23   <- s$NDVI_2 - s$NDVI_3
  early45  <- s$NDVI_4 - s$NDVI_5
  ok <- !is.na(term23) & !is.na(early45)
  clean_slp  <- mean(term23[ok])
  clean_anom <- mean(term23[ok] - early45[ok])
  pwil <- tryCatch(wilcox.test(term23[ok], early45[ok], paired=TRUE)$p.value,
                   error=function(e) NA)
  tab <- rbind(tab, data.frame(
    Zone=z, n=nrow(s),
    Apparent_terminal_step_dNDVI_1_2 = round(app_slp,4),
    Post_event_days_in_window1_median = round(median(post),1),
    Clean_terminal_step_dNDVI_2_3    = round(clean_slp,4),
    Clean_anomaly_vs_NDVI_4_5        = round(clean_anom,4),
    Clean_p_wilcoxon                 = signif(pwil,3)))
  out("  %s: apparent %+0.2f%% (step %+.4f) -> clean step %+.4f  p=%.2g\n",
      z, app_pct, app_slp, clean_slp, pwil)
}
print(tab, row.names=FALSE)
write.csv(tab, "supplementary/Table_S6.csv",             row.names=FALSE)
write.csv(tab, "models_rev/TableS17_contamination.csv",   row.names=FALSE)
cat("-> supplementary/Table_S6.csv\n")

# B3. Breakpoints are robust to window choice (central result unaffected)
out("\n  Breakpoints by window (robustness):\n")
for (col in c("NDVI_1","NDVI_2","NDVI_3","NDVI_5")) {
  b <- seg_bp(d[[col]]); out("  %-7s BP1=%.3f  BP2=%.3f\n", col, b[1], b[2])
}
both2 <- d[!is.na(d$NDVI_1)&!is.na(d$NDVI_3),]
out("  Zone-label churn NDVI_1 vs NDVI_3 = %.1f%%\n",
    100*mean(zone_of(both2$NDVI_1)!=zone_of(both2$NDVI_3)))

cat("\n=== PART A+B COMPLETE ===\n")
cat("Tables S13, S15, S16, S17 written to supplementary/\n")
