#!/usr/bin/env Rscript
# =============================================================================
# 04_supplementary_figures.R — Reproduce supplementary Figures S1–S10
#
# Paper: "Landslide susceptibility peaks at intermediate vegetation density:
#         a global satellite analysis of critical NDVI thresholds"
#
# (The former global inventory map is now panel A of main Figure 1; see
#  03_main_figures.R. Supplementary figures are numbered to match the paper.)
#
# Figure S1 (1-panel): Bootstrap distribution of BP1 and BP2 (n = 10,000).
# Figure S2 (2-panel): EVI and LAI breakpoints (index comparison).
# Figure S3 (2-panel): GAM vs. segmented fit and first derivative — NDVI.
# Figure S4 (2-panel): GAM vs. segmented fit and first derivative — EVI.
# Figure S5 (2-panel): GAM vs. segmented fit and first derivative — LAI.
# Figure S6 (1-panel): Cross-platform land cover Sankey/alluvial diagram
#   (MODIS → Copernicus → Hansen → Unified). Requires ggalluvial.
# Figure S7 (2-panel): NDVI–landslide breakpoints within climate zones.
# Figure S8 (2-panel): Breakpoints across QA/geographic subsets + missingness.
#   Writes supplementary/Table_S3.csv.
# Figure S9 (1-panel): Breakpoint stability under MODIS NDVI sensor noise.
# Figure S10 (2-panel): Compositing-contamination diagnostic.
#
# (A Normalized Frequency Ratio figure was removed in revision; see the note in
#  the body. The exposure question is addressed by the grid-based RF, 05_*.py.)
#
# Must run AFTER:
#   01_core_analysis.R          (→ models_rev/bootstrap_ci_ndvi.csv,
#                                   models_rev/bootstrap_ndvi.rds)
#   02_sensitivity_and_contamination.R
#                               (→ models_rev/TableS16_noise_propagation.csv,
#                                   supplementary/Table_S6.csv)
#
# Outputs: supplementary/Figure_S1.png … Figure_S10.png  + Table_S3.csv
#
# Run:
#   Rscript code/04_supplementary_figures.R
#
# Requires: segmented, mgcv, ggplot2, dplyr, patchwork, tidyr, scales, ggalluvial
# =============================================================================
suppressMessages({
  library(segmented); library(mgcv)
  library(ggplot2); library(dplyr)
  library(patchwork); library(tidyr); library(scales)
})
set.seed(123)

d    <- read.csv("data/landslide_data_with_terrain.csv", stringsAsFactors=FALSE)
boot <- read.csv("models_rev/bootstrap_ci_ndvi.csv")
s16  <- read.csv("models_rev/TableS16_noise_propagation.csv")
dir.create("supplementary", showWarnings=FALSE)

th <- theme_bw(base_size=11) + theme(
  panel.grid.minor=element_blank(),
  plot.title=element_text(face="bold", size=12),
  axis.title=element_text(size=11), legend.position="none")
col_idz<-"#2C7FB8"; col_ctz<-"#C51B7D"; col_sdz<-"#1B9E77"

BP1   <- boot$Estimate[1]; BP2 <- boot$Estimate[2]
bp1ci <- c(boot$CI_Lower[1], boot$CI_Upper[1])
bp2ci <- c(boot$CI_Lower[2], boot$CI_Upper[2])

# ── Shared segmented helper ────────────────────────────────────────────────
seg_full <- function(x_vals) {
  x_vals <- x_vals[!is.na(x_vals)]
  n_bins <- min(150, floor(length(x_vals)/20)); rng<-range(x_vals)
  bs <- seq(rng[1],rng[2],length.out=n_bins+1); h<-hist(x_vals,breaks=bs,plot=FALSE)
  cf <- stats::smooth(stats::smooth(stats::smooth(h$counts,"3R"),"S"),"3R")
  df <- data.frame(x=h$mids,y=as.numeric(cf)); lm0<-lm(y~x,data=df)
  psi0<-quantile(df$x,c(.4,.8))
  psi0<-pmax(rng[1]+.05*diff(rng),pmin(rng[2]-.05*diff(rng),psi0))
  sm<-segmented(lm0,seg.Z=~x,psi=list(x=psi0))
  pr<-data.frame(x=seq(min(df$x),max(df$x),length.out=300)); pr$y<-predict(sm,pr)
  list(df=df, pred=pr, bp=sort(sm$psi[,"Est."]), r2=summary(sm)$r.squared)
}


# =============================================================================
# FIGURE S2 — EVI and LAI breakpoints (index comparison)
# =============================================================================
mk_panel <- function(var, lab, color, tag) {
  r <- seg_full(d[[var]]); ymx<-max(r$df$y)
  ggplot(r$df, aes(x,y)) +
    geom_point(alpha=0.4, color="grey55", size=0.8) +
    geom_line(data=r$pred, aes(x,y), color=color, linewidth=1.2) +
    geom_vline(xintercept=r$bp, linetype="dashed", color="#E31A1C", linewidth=0.7) +
    annotate("text", x=ifelse(var=="LAI_1",max(r$df$x)*0.95,0.02), y=ymx*0.96,
             hjust=ifelse(var=="LAI_1",1,0), fontface="bold", size=3.2,
             label=sprintf("R² = %.3f", r$r2)) +
    annotate("label", x=r$bp[1], y=ymx*0.5, size=2.5, color="#E31A1C", fontface="bold",
             label=sprintf("BP1=%.3f",r$bp[1])) +
    annotate("label", x=r$bp[2], y=ymx*0.75, size=2.5, color="#E31A1C", fontface="bold",
             label=sprintf("BP2=%.3f",r$bp[2])) +
    labs(title=tag, x=lab, y="Landslide frequency (smoothed)") + th
}
pS9 <- mk_panel("EVI_1","Pre-failure EVI","#D95F02","(A) EVI") +
       mk_panel("LAI_1","Pre-failure LAI","#7570B3","(B) LAI")
ggsave("supplementary/Figure_S2.png", pS9, width=11, height=4.4, dpi=300)
cat("-> supplementary/Figure_S2.png\n")


# =============================================================================
# FIGURE S8 — Breakpoints across QA and geographic subsets + missingness map
# Writes Table_S3 (subset breakpoint stability table)
# =============================================================================
monsoon <- c("Bangladesh","Vietnam","Myanmar","India","Nepal","Philippines",
             "Indonesia","Thailand","Sri Lanka")
subsets <- list(
  "QA=0 (strict)"          = d$NDVI_1[!is.na(d$NDVI_1) & d$QA_1==0],
  "QA<=1 (full reliable)"  = d$NDVI_1[!is.na(d$NDVI_1)],
  "Excl. monsoon Asia"     = d$NDVI_1[!is.na(d$NDVI_1) & !(d$Country_Name %in% monsoon)],
  "Asia/Oceania"           = d$NDVI_1[!is.na(d$NDVI_1) & d$Continent=="Asia/Oceania"],
  "Americas"               = d$NDVI_1[!is.na(d$NDVI_1) & d$Continent %in% c("South America","North America")],
  "Europe"                 = d$NDVI_1[!is.na(d$NDVI_1) & d$Continent=="Europe"])
s10 <- data.frame()
for (nm in names(subsets)) {
  x <- subsets[[nm]]
  if (length(x)>=300) {
    b <- seg_full(x)
    s10 <- rbind(s10, data.frame(Subset=nm, n=length(x), BP1=b$bp[1], BP2=b$bp[2]))
  }
}
s10$Subset <- factor(s10$Subset, levels=rev(s10$Subset))
write.csv(s10, "supplementary/Table_S3.csv", row.names=FALSE)
cat("-> supplementary/Table_S3.csv\n")

s10l <- tidyr::pivot_longer(s10, c(BP1,BP2), names_to="BP", values_to="val")
pS10a <- ggplot(s10l, aes(val, Subset, color=BP)) +
  annotate("rect", xmin=bp1ci[1], xmax=bp1ci[2], ymin=-Inf, ymax=Inf, fill=col_ctz, alpha=0.10) +
  annotate("rect", xmin=bp2ci[1], xmax=bp2ci[2], ymin=-Inf, ymax=Inf, fill=col_sdz, alpha=0.10) +
  geom_vline(xintercept=c(BP1,BP2), linetype="dotted", color="grey55") +
  geom_point(size=4.0) +
  geom_text(aes(label=sprintf("%.3f",val)), vjust=-0.9, size=4.0, show.legend=FALSE) +
  geom_text(data=s10, aes(x=0.695, y=Subset, label=sprintf("n=%s",format(n,big.mark=","))),
            inherit.aes=FALSE, hjust=0, size=3.6, color="grey40") +
  scale_color_manual(values=c(BP1=col_ctz, BP2=col_sdz)) +
  coord_cartesian(xlim=c(0.69,0.90)) +
  theme_bw(base_size=16) + theme(panel.grid.minor=element_blank(),
    plot.title=element_text(face="bold", size=rel(0.95)), legend.position="top",
    legend.title=element_blank(), legend.text=element_text(size=rel(0.95))) +
  labs(title="(A) Stability across QA & geographic subsets",
       x="NDVI breakpoint (shaded = full-sample 95% bootstrap CI)", y=NULL)

miss_by <- d %>% group_by(Country_Name) %>%
  summarise(total=n(), missing=sum(is.na(NDVI_1)), .groups="drop") %>%
  filter(total>=80) %>% mutate(pct=100*missing/total) %>% arrange(desc(pct)) %>% head(12)
pS10b <- ggplot(miss_by, aes(reorder(Country_Name,pct), pct)) +
  geom_col(fill=col_ctz, alpha=0.85) + coord_flip() +
  labs(title="(B) NDVI-missing fraction by country",
       x=NULL, y="% events missing reliable NDVI") +
  theme_bw(base_size=16) + theme(panel.grid.minor=element_blank(),
    plot.title=element_text(face="bold", size=rel(0.95)), legend.position="none")

ggsave("supplementary/Figure_S8.png",
       pS10a+pS10b+patchwork::plot_layout(widths=c(1.25,1)), width=12, height=5.0, dpi=300)
cat(sprintf("-> supplementary/Figure_S8.png  (BP1 range [%.3f,%.3f], BP2 range [%.3f,%.3f])\n",
    min(s10$BP1),max(s10$BP1),min(s10$BP2),max(s10$BP2)))


# =============================================================================
# FIGURE S9 — NDVI sensor-noise propagation
# =============================================================================
s16l <- s16 %>% select(Sigma,BP1_mean,BP1_sd,BP2_mean,BP2_sd) %>%
  pivot_longer(-Sigma, names_to=c("BP",".value"), names_sep="_")
pS11 <- ggplot(s16l, aes(Sigma, mean, color=BP)) +
  geom_ribbon(aes(ymin=mean-sd, ymax=mean+sd, fill=BP), alpha=0.2, color=NA) +
  geom_line(linewidth=1) + geom_point(size=2.5) +
  scale_color_manual(values=c(BP1=col_ctz, BP2=col_sdz)) +
  scale_fill_manual(values=c(BP1=col_ctz,  BP2=col_sdz)) +
  annotate("rect", xmin=0.02, xmax=0.05, ymin=-Inf, ymax=Inf, fill="grey70", alpha=0.15) +
  annotate("text", x=0.035, y=0.80, label="MODIS NDVI\nuncertainty\n(±0.02–0.05)", size=2.8) +
  theme_bw(base_size=11) + theme(panel.grid.minor=element_blank(),
    plot.title=element_text(face="bold"), legend.position=c(0.88,0.5)) +
  labs(title="Breakpoint stability under added MODIS NDVI noise",
       x="Added Gaussian noise σ (NDVI units)", y="Estimated breakpoint", color="", fill="")
ggsave("supplementary/Figure_S9.png", pS11, width=7, height=4.4, dpi=300)
cat("-> supplementary/Figure_S9.png\n")


# =============================================================================
# FIGURE S10 — Compositing contamination diagnostic
# =============================================================================
col_app   <- "#E31A1C"; col_clean <- "#1B9E77"
postd     <- pmax(0, 16 - d$Days_Before_1[!is.na(d$Days_Before_1)])

pS12a <- ggplot(data.frame(x=postd), aes(x)) +
  geom_histogram(binwidth=1, fill=col_app, alpha=0.85, color="white", linewidth=0.2) +
  geom_vline(xintercept=median(postd), linetype="dashed", linewidth=0.9) +
  annotate("text", x=median(postd)+0.4, y=Inf, vjust=1.8, hjust=0, size=4.3,
           label=sprintf("median = %.0f d after failure\n96%% of events > 0", median(postd))) +
  theme_bw(base_size=16) + theme(panel.grid.minor=element_blank(),
    plot.title=element_text(face="bold", size=rel(0.95))) +
  labs(title="(A) The terminal '0-16 d' composite straddles failure",
       x="Days of the terminal 16-day MODIS composite falling AFTER failure",
       y="Number of events")

s17 <- read.csv("supplementary/Table_S6.csv")
bdf <- rbind(
  data.frame(Zone=s17$Zone, Phase="Apparent (straddles failure)",
             val=s17$Apparent_terminal_step_dNDVI_1_2),
  data.frame(Zone=s17$Zone, Phase="Clean (pre-failure)",
             val=s17$Clean_terminal_step_dNDVI_2_3))
bdf$Zone  <- factor(bdf$Zone, levels=c("IDZ","CTZ","SDZ"))
bdf$Phase <- factor(bdf$Phase,
  levels=c("Apparent (straddles failure)","Clean (pre-failure)"))
lab     <- s17
lab$txt <- ifelse(lab$Clean_p_wilcoxon>=0.05, "n.s.", sprintf("p=%.2g", lab$Clean_p_wilcoxon))
pS12b <- ggplot(bdf, aes(Zone, val, fill=Phase)) +
  geom_col(position=position_dodge(0.7), width=0.62, alpha=0.9) +
  geom_hline(yintercept=0, color="grey40") +
  scale_fill_manual(values=setNames(c(col_app,col_clean), levels(bdf$Phase))) +
  geom_text(data=lab, aes(x=Zone, y=Clean_terminal_step_dNDVI_2_3, label=txt),
            inherit.aes=FALSE, vjust=-0.6, hjust=-0.15, size=4.3, color=col_clean) +
  theme_bw(base_size=16) + theme(panel.grid.minor=element_blank(),
    plot.title=element_text(face="bold", size=rel(0.95)), legend.position="top",
    legend.title=element_blank(), legend.text=element_text(size=rel(0.9))) +
  labs(title="(B) Decline vanishes with clean windows",
       x="Biophysical zone", y="Terminal NDVI step (per 16 d)")

ggsave("supplementary/Figure_S10.png",
       pS12a+pS12b+patchwork::plot_layout(widths=c(1.1,1)), width=12, height=5.0, dpi=300)
cat("-> supplementary/Figure_S10.png\n")


# NOTE: The global inventory map (formerly a standalone supplementary panel)
# is now panel (A) of main-text Figure 1, generated by 03_main_figures.R.


# =============================================================================
# FIGURE S6 — Sankey / alluvial cross-platform land cover consistency
# =============================================================================
if (requireNamespace("ggalluvial", quietly=TRUE)) {
  suppressMessages(library(ggalluvial))

  unified_class <- function(modis, cop, han) {
    dplyr::case_when(
      grepl("Forest",    modis)|grepl("Forest",    cop)|grepl("Forest",    han) ~ "Forest",
      grepl("Woody",     modis)|grepl("Woody",     cop)|grepl("Woody",     han) ~ "Woody Vegetation",
      grepl("Grassland", modis)|grepl("Grassland", cop)                         ~ "Grassland",
      grepl("Cropland",  modis)|grepl("Cropland",  cop)                         ~ "Cropland/Vegetation",
      grepl("Urban",     modis)|grepl("Urban",     cop)                         ~ "Urban/Built-up",
      grepl("Barren",    modis)|grepl("Barren",    cop)                         ~ "Barren/Sparse",
      grepl("Water",     modis)|grepl("Water",     cop)                         ~ "Water/Wetland/Ice",
      han == "Non-forest"                                                        ~ "Non-forest",
      TRUE                                                                       ~ "Unknown")
  }

  cls <- c("Forest","Woody Vegetation","Grassland","Cropland/Vegetation",
           "Urban/Built-up","Barren/Sparse","Water/Wetland/Ice","Non-forest")
  pal <- c(Forest="#1a9850","Woody Vegetation"="#91cf60",Grassland="#d9ef8b",
           "Cropland/Vegetation"="#fee08b","Urban/Built-up"="#999999",
           "Barren/Sparse"="#d0d0d0","Water/Wetland/Ice"="#4393c3",
           "Non-forest"="#bababa",Unknown="#525252")

  sk <- d %>%
    filter(!is.na(MODIS_IGBP_Simplified) & !is.na(Copernicus_LC_Simplified) &
           !is.na(Hansen_Simplified) &
           MODIS_IGBP_Simplified!="Unknown" & Copernicus_LC_Simplified!="Unknown") %>%
    group_by(MODIS_IGBP_Simplified, Copernicus_LC_Simplified, Hansen_Simplified) %>%
    summarise(n=n(), .groups="drop") %>%
    mutate(alluvium=row_number(),
           Unified=unified_class(MODIS_IGBP_Simplified, Copernicus_LC_Simplified, Hansen_Simplified))

  sk_long <- sk %>%
    pivot_longer(c(MODIS_IGBP_Simplified, Copernicus_LC_Simplified, Hansen_Simplified),
                 names_to="Source_", values_to="stratum") %>%
    mutate(Source_=gsub("_Simplified","",Source_))

  sk_uni <- sk %>%
    transmute(alluvium, n, Source_="Unified", stratum=Unified)

  pd <- bind_rows(sk_long %>% select(alluvium,n,Source_,stratum), sk_uni) %>%
    mutate(stratum=factor(stratum, levels=cls),
           Source_=factor(Source_, levels=c("MODIS_IGBP","Copernicus_LC","Hansen","Unified")))

  p_s2 <- ggplot(pd, aes(x=Source_, stratum=stratum, alluvium=alluvium, y=n, fill=stratum)) +
    geom_flow(stat="alluvium", lode.guidance="forward", color="darkgray", alpha=0.6) +
    geom_stratum(alpha=0.9) +
    scale_fill_manual(values=pal, name="Land Cover Classification",
                      breaks=cls, labels=cls) +
    scale_x_discrete(limits=c("MODIS_IGBP","Copernicus_LC","Hansen","Unified"),
                     labels=c("MODIS","Copernicus","Hansen","Unified"),
                     name="Data Source / Classification",
                     expand=expansion(mult=c(0.02,0.02))) +
    scale_y_continuous(expand=expansion(mult=c(0,0))) +
    labs(title="Cross-Platform Land Cover Classification Consistency",
         subtitle="Flow of landslide event locations between MODIS, Copernicus, and Hansen classifications") +
    theme_minimal(base_family="sans") +
    theme(legend.position="top", legend.box="horizontal",
          legend.title=element_text(size=9), legend.text=element_text(size=8),
          legend.key.size=unit(0.5,"cm"),
          axis.text.y=element_blank(), axis.title.y=element_blank(),
          panel.grid=element_blank(),
          plot.margin=ggplot2::margin(5,5,5,5))

  ggsave("supplementary/Figure_S6.png", p_s2, width=8.3, height=6.8, dpi=300)
  cat("-> supplementary/Figure_S6.png\n")
} else {
  cat("SKIP Figure_S6 — install ggalluvial\n")
}


# NOTE: A Normalized Frequency Ratio (NFR) figure was removed in revision. The
# NFR requires a global per-month NDVI *area* distribution as background; a
# rigorous monthly-matched MOD13C2 background was not available, and proxy
# backgrounds make the NFR's post-peak shape background-dependent. The
# "elevated susceptibility is not a land-area exposure artefact" point is made
# instead by the grid-based Random Forest (05_grid_rf_validation.py), which
# predicts landslide density per fixed-area 0.5-degree cell.


# =============================================================================
# FIGURE S7 — Breakpoints by climate zone (Temperate and Tropical subsets)
# =============================================================================
cz_panels <- list(
  Temperate = d$NDVI_1[!is.na(d$NDVI_1) & d$Climate_Zone=="Temperate"],
  Tropical  = d$NDVI_1[!is.na(d$NDVI_1) & d$Climate_Zone=="Tropical"])

cz_colors <- c(Temperate="#2C7FB8", Tropical="#1B9E77")
ps3_list  <- list()
for (nm in names(cz_panels)) {
  r <- seg_full(cz_panels[[nm]])
  ymax_cz <- max(r$df$y)
  ps3_list[[nm]] <- ggplot(r$df, aes(x,y)) +
    geom_point(alpha=0.35, color="grey60", size=1.0) +
    geom_line(data=r$pred, aes(x,y), color=cz_colors[nm], linewidth=1.4) +
    geom_vline(xintercept=r$bp, linetype="dashed", color="#E31A1C", linewidth=0.9) +
    annotate("label", x=r$bp[1], y=ymax_cz*0.55, size=4.3, color="#E31A1C", fontface="bold",
             label=sprintf("BP1=%.3f", r$bp[1])) +
    annotate("label", x=r$bp[2], y=ymax_cz*0.80, size=4.3, color="#E31A1C", fontface="bold",
             label=sprintf("BP2=%.3f", r$bp[2])) +
    annotate("text", x=0.02, y=ymax_cz*0.97, hjust=0, size=4.6,
             label=sprintf("R²=%.3f  n=%s", r$r2, format(length(cz_panels[[nm]]), big.mark=","))) +
    labs(title=nm, x="NDVI", y="Smoothed frequency") +
    theme_bw(base_size=16) +
    theme(panel.grid.minor=element_blank(), plot.title=element_text(face="bold"))
}
pS3 <- ps3_list$Temperate | ps3_list$Tropical    # no stray plot title
ggsave("supplementary/Figure_S7.png", pS3, width=11, height=4.6, dpi=300)
cat("-> supplementary/Figure_S7.png\n")


# =============================================================================
# FIGURE S1 — Bootstrap distribution of BP1 and BP2 (n = 10,000)
# =============================================================================
br_path <- "models_rev/bootstrap_ndvi.rds"
if (file.exists(br_path)) {
  br <- readRDS(br_path)
  bp1_boot <- br$t[,1]; bp2_boot <- br$t[,2]
  bp_df <- rbind(
    data.frame(Breakpoint="BP1", val=bp1_boot),
    data.frame(Breakpoint="BP2", val=bp2_boot))
  bp_df$Breakpoint <- factor(bp_df$Breakpoint)

  pS5 <- ggplot(bp_df, aes(val, fill=Breakpoint)) +
    geom_histogram(binwidth=0.004, alpha=0.8, color="white", linewidth=0.2,
                   position="identity") +
    geom_vline(xintercept=c(BP1, BP2), color="red", linewidth=0.9) +
    geom_vline(xintercept=c(boot$CI_Lower[1], boot$CI_Upper[1],
                            boot$CI_Lower[2], boot$CI_Upper[2]),
               color="red", linetype="dashed", linewidth=0.7) +
    scale_fill_manual(values=c(BP1="#2C7FB8", BP2="#E6AB02")) +
    scale_y_continuous(labels=scales::comma) +
    labs(title=sprintf("Bootstrap Distribution of Breakpoints (n = %s)", format(nrow(br$t), big.mark=",")),
         x="Breakpoint Value", y="Count", fill="Breakpoint") +
    theme_bw(base_size=11) +
    theme(panel.grid.minor=element_blank(), plot.title=element_text(face="bold"),
          legend.position=c(0.88,0.75), legend.background=element_rect(fill="white",color="grey85"))
  ggsave("supplementary/Figure_S1.png", pS5, width=10, height=5, dpi=300)
  cat("-> supplementary/Figure_S1.png\n")
} else {
  cat("SKIP Figure_S1 — models_rev/bootstrap_ndvi.rds not found (run 01_core_analysis.R first)\n")
}


# =============================================================================
# FIGURES S6–S8 — GAM vs. segmented regression + first derivative
# =============================================================================
gam_panel <- function(var, lab, line_col) {
  x_vals <- d[[var]][!is.na(d[[var]])]
  n_bins <- min(150, floor(length(x_vals)/20)); rng <- range(x_vals)
  bs <- seq(rng[1], rng[2], length.out=n_bins+1)
  h  <- hist(x_vals, breaks=bs, plot=FALSE)
  cf <- stats::smooth(stats::smooth(stats::smooth(h$counts,"3R"),"S"),"3R")
  bin_df <- data.frame(x=h$mids, y=as.numeric(cf))

  # Segmented fit
  lm0  <- lm(y~x, data=bin_df)
  psi0 <- quantile(bin_df$x, c(.4,.8))
  psi0 <- pmax(rng[1]+.05*diff(rng), pmin(rng[2]-.05*diff(rng), psi0))
  sm   <- segmented(lm0, seg.Z=~x, psi=list(x=psi0))
  seg_bp <- sort(sm$psi[,"Est."])

  # GAM fit
  gm <- mgcv::gam(y ~ s(x, k=20, bs="tp"), data=bin_df, method="REML")
  r2_gam <- summary(gm)$r.sq
  pred_x <- seq(rng[1], rng[2], length.out=400)
  pred_y <- as.numeric(predict(gm, newdata=data.frame(x=pred_x)))
  dy <- diff(pred_y); dx <- diff(pred_x)
  deriv_x <- pred_x[-1] - dx/2
  deriv_y <- dy/dx
  sign_chg <- which(diff(sign(deriv_y)) != 0)
  gam_tp   <- deriv_x[sign_chg]
  gam_df   <- data.frame(x=pred_x, y=pred_y)
  deriv_df <- data.frame(x=deriv_x, d=deriv_y)

  pa <- ggplot(bin_df, aes(x,y)) +
    geom_point(alpha=0.35, color="grey60", size=0.8) +
    geom_line(data=gam_df, aes(x,y), color=line_col, linewidth=1.2) +
    geom_vline(xintercept=seg_bp, linetype="dashed", color="#E31A1C", linewidth=0.7) +
    {if(length(gam_tp)>0) geom_vline(xintercept=gam_tp, linetype="dotted",
                                      color=line_col, linewidth=0.7)} +
    annotate("text", x=rng[1]+0.02*diff(rng), y=max(bin_df$y)*0.95, hjust=0, size=2.8,
             label=sprintf("GAM R² = %.3f", r2_gam), color=line_col, fontface="bold") +
    labs(title=sprintf("(a) GAM vs Segmented: %s", var),
         subtitle="Red dashed = segmented BPs, Green dotted = GAM turning points",
         x=lab, y="Landslide Frequency (Smoothed)") +
    theme_bw(base_size=10) + theme(panel.grid.minor=element_blank())

  pb <- ggplot(deriv_df, aes(x,d)) +
    geom_line(color=line_col, linewidth=1.0) +
    geom_hline(yintercept=0, linetype="dashed", color="grey40") +
    geom_vline(xintercept=seg_bp, linetype="dashed", color="#E31A1C", linewidth=0.7) +
    labs(title=sprintf("(b) GAM First Derivative: %s", var),
         x=lab, y="dFrequency/dNDVI") +
    theme_bw(base_size=10) + theme(panel.grid.minor=element_blank())

  list(pa=pa, pb=pb)
}

panels_ndvi <- gam_panel("NDVI_1", "NDVI_1 Value", "#1B7837")
ggsave("supplementary/Figure_S3.png",
       (panels_ndvi$pa / panels_ndvi$pb) + plot_layout(heights=c(1.1,0.9)),
       width=7, height=8, dpi=300)
cat("-> supplementary/Figure_S3.png\n")

panels_evi  <- gam_panel("EVI_1",  "EVI_1 Value",  "#D95F02")
ggsave("supplementary/Figure_S4.png",
       (panels_evi$pa / panels_evi$pb) + plot_layout(heights=c(1.1,0.9)),
       width=7, height=8, dpi=300)
cat("-> supplementary/Figure_S4.png\n")

panels_lai  <- gam_panel("LAI_1",  "LAI_1 Value",  "#7570B3")
ggsave("supplementary/Figure_S5.png",
       (panels_lai$pa / panels_lai$pb) + plot_layout(heights=c(1.1,0.9)),
       width=7, height=8, dpi=300)
cat("-> supplementary/Figure_S5.png\n")


cat("\n=== SUPPLEMENTARY FIGURES COMPLETE: S1–S12 ===\n")
cat("(S1 skipped if sf/rnaturalearth not installed; S2 skipped if ggalluvial not installed)\n")
