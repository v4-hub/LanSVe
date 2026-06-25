#!/usr/bin/env Rscript
# =============================================================================
# 03_main_figures.R — Reproduce all four main-text figures
#
# Paper: "Landslide susceptibility peaks at intermediate vegetation density:
#         a global satellite analysis of critical NDVI thresholds"
#
# Figure 1 (3-panel composite):
#   A: Global inventory map — events colored by pre-failure NDVI, country
#      choropleth by event count (requires sf, rnaturalearth, ggnewscale).
#   B: Smoothed landslide-frequency histogram vs. pre-failure NDVI, with
#      segmented-regression fit, breakpoints BP1/BP2, and bootstrap CIs.
#   C: Zone distribution (% of events with valid NDVI in IDZ / CTZ / SDZ).
#
# Figure 2 (5-panel): Multi-source vegetation data validation.
#   A: Treemap of MODIS vegetation types at landslide locations.
#   B: Pie chart of cross-platform land cover consistency.
#   C: Zone distribution by land cover type.
#   D: Zone distribution by forest cover density class.
#   E: Zone distribution by Köppen climate zone.
#
# Figure 3 (3-panel):
#   A: Pre-failure NDVI density distribution by zone (IDZ / CTZ / SDZ).
#   B: NDVI distribution across Köppen climate zones (two-row layout: A on top).
#   C: Seasonal timing vs. NDVI class (stacked proportional bars).
#
# Figure 4 (2-panel): Conceptual framework (mathematical illustration).
#   A: Unimodal susceptibility curve with IDZ / CTZ / SDZ zones.
#   B: Competing forces (destabilizing vs. stabilizing) across NDVI.
#
# Must run AFTER 01_core_analysis.R (reads models_rev/bootstrap_ci_ndvi.csv).
#
# Outputs:
#   figures/fig1.png
#   figures/fig2.png
#   figures/fig3.png
#   figures/fig4.png
#
# Run:
#   Rscript code/03_main_figures.R
#
# Requires: segmented, ggplot2, dplyr, patchwork, scales, treemapify, cowplot
#           sf, rnaturalearth, ggnewscale (for Figure 1A map)
# (Install missing pkgs, e.g. install.packages("treemapify"))
# =============================================================================
suppressMessages({
  library(segmented); library(ggplot2); library(dplyr); library(patchwork); library(scales)
})
set.seed(123)

d    <- read.csv("data/landslide_data_with_terrain.csv", stringsAsFactors=FALSE)
boot <- read.csv("models_rev/bootstrap_ci_ndvi.csv")
dir.create("figures", showWarnings=FALSE)

th <- theme_bw(base_size=11) + theme(
  panel.grid.minor=element_blank(),
  plot.title=element_text(face="bold", size=12),
  axis.title=element_text(size=11), legend.position="none")
col_idz<-"#2C7FB8"; col_ctz<-"#C51B7D"; col_sdz<-"#1B9E77"

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
# FIGURE 1 — (A) global inventory map  (B) NDVI breakpoints  (C) zones
# =============================================================================
ndvi  <- seg_full(d$NDVI_1)
BP1   <- ndvi$bp[1]; BP2 <- ndvi$bp[2]; ymax <- max(ndvi$df$y)
bp1ci <- c(boot$CI_Lower[1], boot$CI_Upper[1])
bp2ci <- c(boot$CI_Lower[2], boot$CI_Upper[2])
cat(sprintf("Fig1: R2=%.3f  BP1=%.3f [%.3f,%.3f]  BP2=%.3f [%.3f,%.3f]\n",
    ndvi$r2, BP1,bp1ci[1],bp1ci[2], BP2,bp2ci[1],bp2ci[2]))

# ── Panel A: global inventory map (requires sf + rnaturalearth + ggnewscale) ─
have_map <- all(sapply(c("sf","rnaturalearth","ggnewscale"), requireNamespace, quietly=TRUE))
if (have_map) {
  suppressMessages({ library(sf); library(rnaturalearth); library(ggnewscale) })
  dp <- d %>% filter(!is.na(Longitude) & !is.na(Latitude) & !is.na(Source) &
                     Source != "" & !is.na(NDVI_1)) %>%
    mutate(Source = factor(Source, levels=c("NASA Reports","NASA Events","e-ITALICA","GGIG Catalog")))
  land_sf <- st_as_sf(dp, coords=c("Longitude","Latitude"), crs=4326)
  wmap    <- ne_countries(scale="medium", returnclass="sf")
  cntry   <- st_join(land_sf, wmap["sovereignt"], join=st_intersects) %>%
    as.data.frame() %>% group_by(sovereignt) %>% summarise(Landslide_Count=n(), .groups="drop")
  wcnt    <- left_join(wmap, cntry, by="sovereignt")
  # Flat (equirectangular / plate carrée) map, cropped to the data latitude band
  grat <- st_graticule(lat=seq(-60,80,30), lon=seq(-180,180,60))
  p1a <- ggplot() +
    geom_sf(data=grat, color="grey78", linewidth=0.2, linetype="dotted") +
    geom_sf(data=wcnt, fill="grey92", color="white", linewidth=0.3) +
    geom_sf(data=filter(wcnt, !is.na(Landslide_Count)), aes(fill=Landslide_Count),
            color="white", linewidth=0.1) +
    scale_fill_viridis_c(name="Landslide count by country", trans="log10", na.value="grey92",
      breaks=c(1,10,100,1000),
      guide=guide_colorbar(barwidth=5.4,barheight=0.57,direction="horizontal",title.position="top",order=3)) +
    ggnewscale::new_scale_fill() +
    geom_sf(data=land_sf, aes(fill=NDVI_1, shape=Source), color="grey10", size=1.7, stroke=0.25) +
    scale_fill_distiller(palette="RdBu", direction=1, name="Pre-failure NDVI at events", breaks=c(0,0.5,1),
      guide=guide_colorbar(barwidth=5.4,barheight=0.57,direction="horizontal",title.position="top",order=2)) +
    scale_shape_manual(name="Data source", values=c("NASA Reports"=21,"NASA Events"=24,"e-ITALICA"=22,"GGIG Catalog"=23),
      guide=guide_legend(title.position="top", ncol=2, override.aes=list(size=2,stroke=0.5), order=1)) +
    coord_sf(crs=4326, xlim=c(-180,180), ylim=c(-56,82), expand=FALSE) + labs(title="(A)") +
    theme_bw(base_size=11) +
    theme(panel.background=element_rect(fill="aliceblue"), panel.grid=element_blank(),
          panel.border=element_rect(fill=NA, color="grey50"),
          axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank(),
          legend.position=c(0.013,0.02), legend.justification=c(0,0),
          legend.box="vertical", legend.box.just="left",
          legend.box.spacing=unit(0.02,"cm"), legend.spacing.y=unit(0,"cm"),
          legend.margin=margin(0,1,0,0),
          legend.background=element_rect(fill=alpha("white",0.55), color=NA),
          legend.title=element_text(face="bold",size=9.75), legend.text=element_text(size=9),
          legend.key=element_rect(fill=alpha("white",0)),
          legend.key.height=unit(0.45,"cm"), legend.key.width=unit(0.45,"cm"),
          plot.title=element_text(face="bold",size=12,hjust=0), plot.margin=margin(2,2,2,2),
          plot.background=element_rect(fill="white",color=NA))
} else {
  p1a <- ggplot() + annotate("text", x=0, y=0, label="(A) global map — install sf+rnaturalearth+ggnewscale") +
    theme_void() + labs(title="(A)")
  cat("NOTE: install sf, rnaturalearth, ggnewscale to render Figure 1A map\n")
}

p1b <- ggplot(ndvi$df, aes(x,y)) +
  annotate("rect", xmin=-Inf,  xmax=BP1, ymin=-Inf, ymax=Inf, fill=col_idz, alpha=0.06) +
  annotate("rect", xmin=BP1,   xmax=BP2, ymin=-Inf, ymax=Inf, fill=col_ctz, alpha=0.10) +
  annotate("rect", xmin=BP2,   xmax=Inf, ymin=-Inf, ymax=Inf, fill=col_sdz, alpha=0.06) +
  annotate("rect", xmin=bp1ci[1], xmax=bp1ci[2], ymin=0, ymax=ymax, fill="#E31A1C", alpha=0.18) +
  annotate("rect", xmin=bp2ci[1], xmax=bp2ci[2], ymin=0, ymax=ymax, fill="#E31A1C", alpha=0.18) +
  geom_point(alpha=0.4, color="grey55", size=0.8) +
  geom_line(data=ndvi$pred, aes(x,y), color="#1B7837", linewidth=1.2) +
  geom_vline(xintercept=c(BP1,BP2), linetype="dashed", color="#E31A1C", linewidth=0.7) +
  annotate("text", x=0.02, y=ymax*0.97, hjust=0, fontface="bold", size=3.4,
           label=sprintf("R² = %.3f", ndvi$r2)) +
  annotate("label", x=BP1-0.02, y=ymax*0.62, hjust=1, size=2.7, color="#E31A1C",
           fontface="bold",
           label=sprintf("BP1 = %.3f\n95%% CI [%.3f, %.3f]", BP1, bp1ci[1], bp1ci[2])) +
  annotate("label", x=0.905,    y=ymax*0.82, hjust=0, size=2.7, color="#E31A1C",
           fontface="bold",
           label=sprintf("BP2 = %.3f\n95%% CI [%.3f, %.3f]", BP2, bp2ci[1], bp2ci[2])) +
  annotate("text", x=BP1/2,         y=ymax*0.30, label="IDZ",
           color=col_idz, fontface="bold", size=3.3) +
  annotate("text", x=(BP1+BP2)/2,   y=ymax*0.55, label="CTZ",
           color=col_ctz, fontface="bold", size=3.3) +
  annotate("text", x=0.95,          y=ymax*0.30, label="SDZ",
           color=col_sdz, fontface="bold", size=3.3) +
  scale_x_continuous(breaks=seq(0,1,0.25)) +
  coord_cartesian(xlim=c(-0.04,1.20), ylim=c(-0.02*ymax, ymax*1.08), expand=FALSE) +
  labs(title="(B)", x="Pre-failure NDVI", y="Landslide frequency (smoothed)") + th

nv   <- d$NDVI_1[!is.na(d$NDVI_1)]; N<-length(nv); FULL<-nrow(d)
zc   <- c(IDZ=sum(nv<BP1), CTZ=sum(nv>=BP1&nv<=BP2), SDZ=sum(nv>BP2))
zdf  <- data.frame(Zone=factor(names(zc),levels=c("IDZ","CTZ","SDZ")),
                   n=as.integer(zc), pct=100*as.integer(zc)/N,
                   pct_full=100*as.integer(zc)/FULL)
cat("Fig1C zones (% analyzed):", paste(sprintf("%s=%.1f%%",zdf$Zone,zdf$pct),collapse=" "),"\n")
p1c <- ggplot(zdf, aes(Zone,pct,fill=Zone)) +
  geom_col(alpha=0.9, width=0.7) +
  geom_text(aes(label=sprintf("%.1f%%\n(n=%s)", pct, format(n,big.mark=","))),
            vjust=-0.3, size=3, fontface="bold") +
  scale_fill_manual(values=c(IDZ=col_idz,CTZ=col_ctz,SDZ=col_sdz)) +
  annotate("label", x=3.45, y=max(zdf$pct)*1.20, hjust=1, vjust=1, size=2.7,
           label.size=0, fill="grey96",
           label=sprintf("CTZ [%.3f–%.3f]:\n23.0%% of all events globally\n(30.9%% unclassified: missing NDVI)",BP1,BP2),
           lineheight=1.0) +
  coord_cartesian(ylim=c(0,max(zdf$pct)*1.24), clip="off") +
  labs(title="(C)", x="Biophysical zone", y="Share of analyzed events (%)") + th

# Compose: (A) map on top, (B | C) below
ggsave("figures/fig1.png", p1a / (p1b | p1c) + plot_layout(heights=c(1.18,1.0)),
       width=11, height=8.9, dpi=300)
cat("-> figures/fig1.png (3-panel composite)\n")


# =============================================================================
# FIGURE 3 — (A) NDVI by zone   (B) NDVI by climate   (C) seasonal timing
# =============================================================================
d3 <- d[!is.na(d$NDVI_1) & d$NDVI_1>=0 & d$NDVI_1<=1,]
d3$Zone <- factor(ifelse(d3$NDVI_1<BP1,"IDZ",ifelse(d3$NDVI_1<=BP2,"CTZ","SDZ")),
                  levels=c("IDZ","CTZ","SDZ"))

# Larger-text theme for Figure 3 (two-row layout)
th3 <- theme_bw(base_size=13) + theme(panel.grid.minor=element_blank(),
  plot.title=element_text(face="bold", size=13.5), axis.title=element_text(size=12.5),
  axis.text=element_text(size=11), legend.position="none")

# Panel A — NDVI density by zone (full-width, half-height top row)
pA <- ggplot(d3, aes(NDVI_1, fill=Zone, color=Zone)) +
  geom_density(alpha=0.35, linewidth=0.8, adjust=1.1) +
  geom_vline(xintercept=c(BP1,BP2), linetype="dashed", color="#E31A1C", linewidth=0.7) +
  scale_fill_manual(values=c(IDZ=col_idz,CTZ=col_ctz,SDZ=col_sdz)) +
  scale_color_manual(values=c(IDZ=col_idz,CTZ=col_ctz,SDZ=col_sdz)) +
  annotate("text", x=BP1/2,       y=Inf, vjust=1.8, label="IDZ", color=col_idz, fontface="bold", size=4) +
  annotate("text", x=(BP1+BP2)/2, y=Inf, vjust=1.8, label="CTZ", color=col_ctz, fontface="bold", size=4) +
  annotate("text", x=0.94,        y=Inf, vjust=1.8, label="SDZ", color=col_sdz, fontface="bold", size=4) +
  labs(title="(A) Pre-failure NDVI distribution by zone", x="Pre-failure NDVI", y="Density") + th3

# Panel B — NDVI by climate zone
cz <- d3[d3$Climate_Zone %in% c("Tropical","Temperate","Cold","Arid"),]
cz$Climate_Zone <- factor(cz$Climate_Zone, levels=c("Arid","Cold","Temperate","Tropical"))
kw <- kruskal.test(NDVI_1~Climate_Zone, data=cz)
pB <- ggplot(cz, aes(Climate_Zone, NDVI_1, fill=Climate_Zone)) +
  geom_boxplot(outlier.size=0.25, alpha=0.85, width=0.6) +
  annotate("text", x=0.6, y=1.04, hjust=0, size=3.8, label="Kruskal–Wallis p < 0.001") +
  scale_fill_brewer(palette="YlGn") +
  coord_cartesian(ylim=c(0,1.06)) +
  labs(title="(B) NDVI differs across climate zones",
       x="Köppen climate zone", y="Pre-failure NDVI") + th3

# Panel C — seasonal timing vs. NDVI class
d3$NDVIclass <- cut(d3$NDVI_1, breaks=c(-Inf,0.3,0.5,0.7,0.8,Inf),
                    labels=c("0–0.3","0.3–0.5","0.5–0.7","0.7–0.8","0.8–1.0"))
sc <- d3 %>% filter(!is.na(NDVIclass), Season %in% c("Spring","Summer","Fall","Winter")) %>%
  count(NDVIclass, Season) %>% group_by(NDVIclass) %>%
  mutate(pct=100*n/sum(n)) %>% ungroup()
sc$Season <- factor(sc$Season, levels=c("Winter","Fall","Summer","Spring"))
pC <- ggplot(sc, aes(NDVIclass, pct, fill=Season)) +
  geom_col(width=0.7, alpha=0.9) +
  scale_fill_manual(values=c(Spring="#66C2A5",Summer="#FC8D62",Fall="#E78AC3",Winter="#8DA0CB")) +
  theme_bw(base_size=13) + theme(panel.grid.minor=element_blank(),
    plot.title=element_text(face="bold",size=13.5), axis.title=element_text(size=12.5),
    axis.text=element_text(size=11), legend.position="right", legend.title=element_blank(),
    legend.text=element_text(size=11)) +
  labs(title="(C) Seasonal timing shifts with vegetation density",
       x="Pre-failure NDVI class", y="Share of events (%)")

# (A) full-width half-height top row, (B | C) below
ggsave("figures/fig3.png", pA / (pB | pC) + plot_layout(heights=c(0.58, 1.0)),
       width=11, height=7.2, dpi=300)
cat(sprintf("-> figures/fig3.png  (KW chi2=%.1f p=%.2g)\n", kw$statistic, kw$p.value))


# =============================================================================
# FIGURE 2 — Multi-source vegetation validation (5 panels)
# =============================================================================
# ── Zone assignment helper ─────────────────────────────────────────────────
zone_of <- function(v) ifelse(v < BP1, "IDZ", ifelse(v <= BP2, "CTZ", "SDZ"))

dv <- d[!is.na(d$NDVI_1) & d$NDVI_1 >= 0 & d$NDVI_1 <= 1,]
dv$Zone <- factor(zone_of(dv$NDVI_1), levels=c("IDZ","CTZ","SDZ"))

zone_bar <- function(df, grp_col, grp_levels, title_label) {
  df2 <- df[df[[grp_col]] %in% grp_levels,]
  df2[[grp_col]] <- factor(df2[[grp_col]], levels=grp_levels)
  tab <- df2 %>% count(.data[[grp_col]], Zone) %>%
    group_by(.data[[grp_col]]) %>%
    mutate(pct=round(100*n/sum(n))) %>% ungroup()
  tab$Zone <- factor(tab$Zone, levels=c("IDZ","CTZ","SDZ"))
  ggplot(tab, aes(.data[[grp_col]], pct, fill=Zone)) +
    geom_col(alpha=0.9, width=0.7) +
    geom_text(aes(label=paste0(pct,"%")), position=position_stack(vjust=0.5),
              size=2.5, color="white", fontface="bold") +
    scale_fill_manual(values=c(IDZ=col_idz, CTZ=col_ctz, SDZ=col_sdz)) +
    scale_y_continuous(labels=function(x) paste0(x,"%")) +
    labs(title=title_label, x=NULL, y=NULL) +
    theme_bw(base_size=10) +
    theme(panel.grid.minor=element_blank(), legend.position="none",
          plot.title=element_text(face="bold", size=10),
          axis.text.x=element_text(size=8, angle=ifelse(nchar(paste(grp_levels,collapse=""))>30,25,0),
                                   hjust=ifelse(nchar(paste(grp_levels,collapse=""))>30,1,0.5)))
}

# Panel A — treemap of MODIS vegetation types
tree_dat <- d %>%
  filter(!is.na(MODIS_IGBP_Class)) %>%
  count(MODIS_IGBP_Class) %>%
  mutate(pct=100*n/sum(n)) %>%
  filter(pct>=0.5) %>%
  arrange(desc(pct))

if (requireNamespace("treemapify", quietly=TRUE)) {
  p2a <- ggplot(tree_dat, aes(area=n, fill=n,
               label=paste0(MODIS_IGBP_Class,"\n",sprintf("%.1f%%",pct)))) +
    treemapify::geom_treemap() +
    treemapify::geom_treemap_text(color="black", place="centre", size=7,
                                  fontface="bold", grow=FALSE, reflow=TRUE, min.size=3) +
    scale_fill_gradientn(colors=c("#E6F0DC","#C1E899","#55883B"), name="Count") +
    labs(title="(a) MODIS vegetation types") +
    theme(plot.title=element_text(face="bold", size=10), legend.position="right",
          plot.margin=ggplot2::margin(2,2,2,2))
} else {
  p2a <- ggplot(tree_dat, aes(reorder(MODIS_IGBP_Class, pct), pct, fill=pct)) +
    geom_col(alpha=0.9) +
    coord_flip() +
    scale_fill_gradientn(colors=c("#E6F0DC","#C1E899","#55883B")) +
    geom_text(aes(label=sprintf("%.1f%%",pct)), hjust=-0.1, size=2.6) +
    labs(title="(a) MODIS vegetation types (install treemapify for treemap)",
         x=NULL, y="% of events") +
    theme_bw(base_size=10) + theme(legend.position="none",
      plot.title=element_text(face="bold", size=9))
}

# Panel B — multi-source consistency pie chart
cons <- d %>% filter(!is.na(Multi_Source_Consistency_Label)) %>%
  count(Multi_Source_Consistency_Label) %>%
  mutate(pct=100*n/sum(n),
         lbl=paste0(Multi_Source_Consistency_Label,"\n",
                    sprintf("%.1f%%(n=%s)", pct, format(n, big.mark=","))))
cons$Multi_Source_Consistency_Label <- factor(
  cons$Multi_Source_Consistency_Label,
  levels=c("Full Agreement","Partial Agreement","Unknown","No Agreement"))
cons_cols <- c("Full Agreement"="#4CAF50","Partial Agreement"="#8BC34A",
               "Unknown"="#BDBDBD","No Agreement"="#F44336")
p2b <- ggplot(cons, aes(x="", y=pct, fill=Multi_Source_Consistency_Label)) +
  geom_col(width=1, color="white", linewidth=0.4) +
  coord_polar("y") +
  geom_text(aes(label=lbl), position=position_stack(vjust=0.5), size=2.3, lineheight=0.9) +
  scale_fill_manual(values=cons_cols) +
  labs(title="(b) Cross-platform consistency") +
  theme_void() +
  theme(plot.title=element_text(face="bold", size=10, hjust=0.5),
        legend.position="none")

# Panels C, D, E — zone distribution by category
lc_lvls <- c("Cropland/Vegetation","Forest","Grassland","Woody Vegetation")
dv$LC4 <- dplyr::case_when(
  dv$MODIS_IGBP_Simplified %in% c("Forest","Woody Vegetation","Grassland",
                                    "Cropland/Vegetation") ~ dv$MODIS_IGBP_Simplified,
  TRUE ~ NA_character_)

p2c <- zone_bar(dv[!is.na(dv$LC4),], "LC4",  lc_lvls, "(c) Land cover type")

fc_lvls <- c("No Forest (0%)","Sparse (1-25%)","Moderate (26-50%)","Dense (51-75%)","Very Dense (>75%)")
p2d <- zone_bar(dv[!is.na(dv$Forest_Cover_Category),], "Forest_Cover_Category", fc_lvls,
                "(d) Forest cover density")

cz_lvls <- c("Arid","Cold","Temperate","Tropical")
p2e <- zone_bar(dv[dv$Climate_Zone %in% cz_lvls,], "Climate_Zone", cz_lvls,
                "(e) Climate zone")

# ── Legend panel ──────────────────────────────────────────────────────────
leg_df <- data.frame(Zone=factor(c("IDZ","CTZ","SDZ"), levels=c("IDZ","CTZ","SDZ")), y=1)
p_leg <- ggplot(leg_df, aes(Zone, y, fill=Zone)) +
  geom_col() +
  scale_fill_manual(values=c(IDZ=col_idz,CTZ=col_ctz,SDZ=col_sdz),
                    name="Biophysical Zone") +
  theme_void() +
  theme(legend.position="bottom", legend.direction="horizontal",
        legend.title=element_text(face="bold", size=9),
        legend.text=element_text(size=9))
leg_grob <- cowplot::get_legend(p_leg)

if (requireNamespace("cowplot", quietly=TRUE)) {
  top_row  <- cowplot::plot_grid(p2a, p2b, ncol=2, rel_widths=c(1.4,1))
  bot_row  <- cowplot::plot_grid(p2c, p2d, p2e, ncol=3)
  fig2_all <- cowplot::plot_grid(top_row, bot_row, leg_grob,
                                  ncol=1, rel_heights=c(1,0.85,0.08))
} else {
  fig2_all <- (p2a | p2b) / (p2c | p2d | p2e) +
    plot_layout(heights=c(1,0.85))
}

ggsave("figures/fig2.png", fig2_all, width=10, height=7, dpi=300)
cat("-> figures/fig2.png\n")


# =============================================================================
# FIGURE 4 — Conceptual framework (mathematical illustration, no raw data)
# =============================================================================
ndvi_x <- seq(0, 1, length.out=400)

# Panel A: unimodal susceptibility curve
# Composed of logistic rise + gaussian-like fall-off mirroring the CTZ peak
susc <- plogis((ndvi_x - 0.40) * 9) * dnorm(ndvi_x, 0.78, 0.14) /
        max(plogis((ndvi_x - 0.40) * 9) * dnorm(ndvi_x, 0.78, 0.14))
susc <- scales::rescale(susc, to=c(0.04, 1))

df4a <- data.frame(x=ndvi_x, y=susc)
p4a <- ggplot(df4a, aes(x,y)) +
  annotate("rect", xmin=0,    xmax=BP1, ymin=0, ymax=Inf, fill=col_idz, alpha=0.10) +
  annotate("rect", xmin=BP1,  xmax=BP2, ymin=0, ymax=Inf, fill=col_ctz, alpha=0.15) +
  annotate("rect", xmin=BP2,  xmax=1,   ymin=0, ymax=Inf, fill=col_sdz, alpha=0.10) +
  geom_line(color="#1B7837", linewidth=1.8) +
  geom_vline(xintercept=c(BP1,BP2), linetype="dotted", color="black", linewidth=0.6) +
  annotate("text", x=BP1, y=1.05, label=sprintf("%.3f",BP1), size=3.2, fontface="bold") +
  annotate("text", x=BP2, y=1.05, label=sprintf("%.3f",BP2), size=3.2, fontface="bold") +
  annotate("text", x=BP1/2,       y=0.18, label="Instability-\nDominated\nZone (IDZ)",
           color=col_idz, fontface="bold", size=2.7, lineheight=0.9) +
  annotate("text", x=(BP1+BP2)/2, y=0.62, label="Critical\nTransition\nZone (CTZ)",
           color=col_ctz, fontface="bold", size=2.7, lineheight=0.9) +
  annotate("text", x=(BP2+1)/2,   y=0.18, label="Stability-\nDominated\nZone (SDZ)",
           color=col_sdz, fontface="bold", size=2.7, lineheight=0.9) +
  geom_point(data=data.frame(x=ndvi_x[which.max(susc)], y=max(susc)), size=3,
             color="#1B7837", shape=21, fill="white", stroke=1.5) +
  scale_x_continuous(breaks=seq(0,1,0.25)) +
  scale_y_continuous(breaks=c(0,0.5,1), labels=c("Low","Medium","High")) +
  coord_cartesian(xlim=c(0,1), ylim=c(0,1.12), expand=FALSE) +
  labs(title="(a)", x="Vegetation density (NDVI)", y="Landslide susceptibility") +
  theme_bw(base_size=11) +
  theme(panel.grid.minor=element_blank(), plot.title=element_text(face="bold", size=12),
        axis.title=element_text(size=11))

# Panel B: destabilizing (surcharge + concentrated infiltration) vs
#           stabilizing (root reinforcement + hydrological regulation)
destab <- plogis((ndvi_x - 0.15) * 12) * (1 - 0.30 * plogis((ndvi_x - 0.70) * 8))
destab <- scales::rescale(destab, to=c(0.05, 1))
stab   <- plogis((ndvi_x - 0.30) * 6) * (1 - 0.12 * exp(-(ndvi_x - 0.95)^2 / 0.04))
stab   <- scales::rescale(stab, to=c(0.02, 0.88))

df4b <- rbind(
  data.frame(x=ndvi_x, y=destab, Force="Destabilizing forces\n(Surcharge load; Concentrated infiltration)"),
  data.frame(x=ndvi_x, y=stab,   Force="Stabilizing forces\n(Root reinforcement; Hydrological regulation)"))
df4b$Force <- factor(df4b$Force, levels=unique(df4b$Force))
force_cols <- c("#C51B7D","#2C7FB8")

p4b <- ggplot(df4b, aes(x, y, color=Force, linetype=Force)) +
  annotate("rect", xmin=BP1, xmax=BP2, ymin=0, ymax=Inf, fill=col_ctz, alpha=0.10) +
  geom_line(linewidth=1.5) +
  annotate("text", x=(BP1+BP2)/2, y=1.05, label="Critical\ntransition",
           size=2.8, fontface="bold", lineheight=0.9) +
  annotate("segment", x=(BP1+BP2)/2, xend=(BP1+BP2)/2, y=0.97, yend=0.92,
           arrow=arrow(length=unit(0.2,"cm")), color="grey30") +
  scale_color_manual(values=force_cols, name=NULL) +
  scale_linetype_manual(values=c("solid","solid"), name=NULL) +
  scale_x_continuous(breaks=seq(0,1,0.25)) +
  scale_y_continuous(breaks=c(0,0.5,1), labels=c("Low","Medium","High")) +
  coord_cartesian(xlim=c(0,1), ylim=c(0,1.12), expand=FALSE) +
  labs(title="(b)", x="Vegetation density (NDVI)", y="Force magnitude") +
  theme_bw(base_size=11) +
  theme(panel.grid.minor=element_blank(), plot.title=element_text(face="bold", size=12),
        axis.title=element_text(size=11), legend.position=c(0.30, 0.82),
        legend.background=element_rect(fill=alpha("white",0.85)),
        legend.text=element_text(size=8), legend.key.width=unit(1.2,"cm"))

ggsave("figures/fig4.png", p4a | p4b, width=11, height=4.5, dpi=300)
cat("-> figures/fig4.png\n")

cat("\n=== MAIN FIGURES COMPLETE: fig1.png, fig2.png, fig3.png, fig4.png ===\n")
