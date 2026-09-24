###############################################################################
## GSS3318 SensiScope — Skin Microbiome Alpha & Beta Diversity
## -----------------------------------------------------------------------------
## Study design
##   * 66 samples = 22 subjects x 3 body sites (Back / Hand / Leg)
##   * Sensitivity (Sensitive / NotSensitive) is a BETWEEN-subject factor
##     (11 Sensitive, 11 NotSensitive subjects)
##   * BodySite is a WITHIN-subject (paired) factor: every subject (SID) is
##     sampled at all three sites
##   * SiteSensitivity = BodySite x Sensitivity (6 groups, 11 subjects each)
##
## What this script does
##   1. Loads the cleaned count table, filters low-abundance taxa and computes
##      per-sample alpha-diversity indices.
##   2. Alpha diversity: association with the SS-10 sensitivity score,
##      Kruskal-Wallis omnibus tests, Sensitive-vs-NotSensitive comparisons
##      within each site, paired cross-site comparisons, fold-changes, plus
##      requested t-test and Pearson-correlation figures.
##   3. Beta diversity: Bray-Curtis on relative abundances, NMDS/PCoA ordination,
##      reference-style SiteSensitivity/SID hull figures, PERMANOVA (adonis2,
##      preserving the original analysis), multivariate dispersion (betadisper)
##      and pairwise PERMANOVA.
##   4. Writes publication-ready figures (PNG only) and an Excel workbook into
##      <input_filename>_result/ (e.g. count7_cleandata_result/). No
##      PowerPoint is created.
##
## Lines 1-221 of the original script (data loading + the intended alpha
## statistics for BodySite / Sensitivity) are preserved in spirit; the join is
## fixed so the grouping columns stay clean, and everything downstream is
## rewritten for this design.
###############################################################################

rm(list = ls())

suppressPackageStartupMessages({
  library(vegan)      # diversity, vegdist, metaMDS, adonis2, betadisper
  library(permute)    # restricted permutation designs (how / blocks)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(ggpubr)     # publication-ready ggplot helpers + stat_compare_means
  library(rstatix)
  library(scales)
  library(openxlsx)
})

set.seed(1234)  # reproducible ordination / permutation tests

## ---------------------------------------------------------------------------
## 0. Configuration, helpers and a shared publication theme
## ---------------------------------------------------------------------------
args     <- commandArgs(trailingOnly = TRUE)
filename <- if (length(args) >= 1) args[1] else "count7_cleandata"
meta_file <- "GSS3318_readsCount_meta.txt"
## Output directory = <input data file name>_result (was a fixed "results"
## folder before), so outputs from different input files land in their own,
## clearly-named directory instead of overwriting one shared "results/".
outdir   <- paste0(basename(filename), "_result")
dir.create(outdir, showWarnings = FALSE)
out <- function(...) file.path(outdir, paste0(...))

## Signed fold-change: values in (0,1) are expressed as negative reciprocals so
## that "2x higher" and "2x lower" are symmetric (kept from the original script).
truefc <- function(ratio) {
  if (!is.finite(ratio) || ratio == 0) return(NA_real_)
  if (ratio < 1) return(-1 / ratio)
  ratio
}

## p-value helpers that never abort the pipeline on degenerate input.
safe_p <- function(expr) tryCatch(expr, error = function(e) NA_real_)
tt_p    <- function(x, y, paired = FALSE)
  safe_p(t.test(x, y, paired = paired)$p.value)
wilcox_p <- function(x, y, paired = FALSE)
  safe_p(suppressWarnings(wilcox.test(x, y, paired = paired, exact = FALSE)$p.value))
kw_p <- function(data, metric, group)
  safe_p(kruskal.test(reformulate(group, response = metric), data = data)$p.value)

## Pair two body sites by subject (SID) so paired tests are correctly aligned.
paired_by_sid <- function(data, metric, siteA, siteB, sens = NULL) {
  d <- data
  if (!is.null(sens)) d <- d[d$Sensitivity == sens, ]
  a <- d[d$BodySite == siteA, c("SID", metric)]
  b <- d[d$BodySite == siteB, c("SID", metric)]
  m <- merge(a, b, by = "SID", suffixes = c(".A", ".B"))
  list(x = m[[paste0(metric, ".A")]], y = m[[paste0(metric, ".B")]])
}

## Colour palettes (colour-blind friendly, consistent across every figure).
pal_sens <- c(NotSensitive = "#4C72B0", Sensitive = "#C44E52")
pal_site <- c(Back = "#55A868", Hand = "#8172B3", Leg = "#DD8452")
pal_ss   <- c(BackNotSensitive = "#9ecae1", BackSensitive = "#08519c",
              HandNotSensitive = "#bcbddc", HandSensitive = "#54278f",
              LegNotSensitive  = "#fdae6b", LegSensitive  = "#a63603")

theme_pub <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = "grey92"),
          strip.background = element_rect(fill = "grey95", colour = NA),
          strip.text       = element_text(face = "bold"),
          plot.title       = element_text(face = "bold", hjust = 0),
          legend.position  = "top")
}

save_fig <- function(plot, name, width, height) {
  ggsave(out(filename, ".", name, ".png"), plot, width = width, height = height,
         dpi = 300, bg = "white")
}

# ## Nice axis / panel labels for the diversity metrics.
# metric_labels <- c(
#   observed  = "Observed richness",
#   shannon   = "Shannon index (H')")
## Nice axis / panel labels for the diversity metrics.
metric_labels <- c(
  observed  = "Observed richness",
  shannon   = "Shannon index",
  invsimp   = "Inverse Simpson",
  simpson   = "Simpson index",
  pielou    = "Pielou's evenness",
  menhinick = "Menhinick index",
  margalef  = "Margalef index")

## ---------------------------------------------------------------------------
## 1. Load counts, filter taxa and compute alpha diversity  (per original 1-85)
## ---------------------------------------------------------------------------
test <- read.table(filename, header = TRUE, row.names = 1, sep = "\t")
test[is.na(test)] <- 0
d <- dim(test)

## Keep taxa whose MEAN count per sample is >= 10 (original filter: X / d[2] >= 10)
keep_taxa   <- (rowSums(test) / d[2]) >= 10
test_filter <- test[keep_taxa, ]
test_d      <- as.data.frame(t(test_filter))   # samples x taxa

## Parse sample names  "BodySite.Sensitivity.SID"
parts       <- strsplit(rownames(test_d), "[.]")
BodySite    <- vapply(parts, `[`, character(1), 1)
Sensitivity <- vapply(parts, `[`, character(1), 2)
SID         <- vapply(parts, `[`, character(1), 3)

## Alpha-diversity indices (richness, diversity, evenness).
N        <- rowSums(test_d)
observed <- rowSums(test_d > 0)
shannon  <- diversity(test_d, index = "shannon")
simpson  <- diversity(test_d, index = "simpson")
invsimp  <- diversity(test_d, index = "invsimpson")
pielou   <- shannon / log(observed)                 # evenness
menhinick <- observed / sqrt(N)
margalef  <- (observed - 1) / log(N)


mydata0 <- data.frame(
  shortID = rownames(test_d),
  BodySite, Sensitivity, SID,
  SiteSensitivity = paste0(BodySite, Sensitivity),
  observed, shannon,simpson, invsimp, pielou, menhinick, margalef,
  stringsAsFactors = FALSE)

## Merge metadata.  The metadata file also carries BodySite/Sensitivity/SID, so
## drop those duplicates before joining to keep single, clean grouping columns
## (this is the fix that lets the intended alpha statistics run).
meta <- read.table(meta_file, sep = "\t", header = TRUE, check.names = TRUE)
meta2 <- meta %>% select(-any_of(c("BodySite", "Sensitivity", "SID")))
mydata <- inner_join(mydata0, meta2, by = c("shortID" = "NewID"))

## Ordered factors for consistent plotting / modelling.
mydata$BodySite    <- factor(mydata$BodySite, levels = c("Back", "Hand", "Leg"))
mydata$Sensitivity <- factor(mydata$Sensitivity, levels = c("NotSensitive", "Sensitive"))
mydata$SID         <- factor(mydata$SID)
mydata$SiteSensitivity <- factor(mydata$SiteSensitivity,
  levels = c("BackNotSensitive", "BackSensitive",
             "HandNotSensitive", "HandSensitive",
             "LegNotSensitive",  "LegSensitive"))

SS10 <- "All.SS10.Total_score_of_all_SS.10_questions"   # sensitivity questionnaire score
alpha_metrics <- c("observed", "shannon", "invsimp", "simpson",
                   "pielou", "menhinick", "margalef")
#alpha_metrics <- c("observed", "shannon")
core_metrics  <- c("observed", "shannon")

## Per-sample alpha table (tidy export).
alpha_tbl <- mydata %>%
  select(shortID, SID, BodySite, Sensitivity, SiteSensitivity,
         all_of(alpha_metrics), all_of(SS10)) %>%
  rename(SS10_total = all_of(SS10))
write.table(alpha_tbl, file = out(filename, ".alphadiversity.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

## ---------------------------------------------------------------------------
## 2. Alpha-diversity statistics (one comprehensive row per metric)
## ---------------------------------------------------------------------------
sites <- c("Back", "Hand", "Leg")

alpha_stat_row <- function(metric) {
  y_all <- mydata[[metric]]

  ## (a) Association with the SS-10 sensitivity score, within each body site
  cors <- lapply(sites, function(s) {
    dd <- mydata[mydata$BodySite == s, ]
    sp <- tryCatch(cor.test(dd[[SS10]], dd[[metric]], method = "spearman", exact = FALSE),
                   error = function(e) NULL)
    pe <- tryCatch(cor.test(dd[[SS10]], dd[[metric]], method = "pearson",  exact = FALSE),
                   error = function(e) NULL)
    c(sp_rho = if (is.null(sp)) NA else unname(sp$estimate),
      sp_p   = if (is.null(sp)) NA else sp$p.value,
      pe_r   = if (is.null(pe)) NA else unname(pe$estimate),
      pe_p   = if (is.null(pe)) NA else pe$p.value)
  })
  names(cors) <- sites

  ## (b) Kruskal-Wallis omnibus tests
  KP_BodySite        <- kw_p(mydata, metric, "BodySite")
  KP_Sensitivity     <- kw_p(mydata, metric, "Sensitivity")
  KP_SiteSensitivity <- kw_p(mydata, metric, "SiteSensitivity")

  ## (c) Sensitive vs NotSensitive within each site + overall (between-subject)
  gS <- y_all[mydata$Sensitivity == "Sensitive"]
  gN <- y_all[mydata$Sensitivity == "NotSensitive"]
  sv <- lapply(sites, function(s) {
    xs <- y_all[mydata$BodySite == s & mydata$Sensitivity == "Sensitive"]
    xn <- y_all[mydata$BodySite == s & mydata$Sensitivity == "NotSensitive"]
    c(wilcox = wilcox_p(xs, xn), t = tt_p(xs, xn),
      meanS = mean(xs, na.rm = TRUE), meanN = mean(xn, na.rm = TRUE),
      fc = truefc(mean(xs, na.rm = TRUE) / mean(xn, na.rm = TRUE)))
  })
  names(sv) <- sites

  ## (d) Paired cross-site comparisons (aligned by SID) for all/Sensitive/NotSensitive
  site_pairs <- list(c("Back", "Hand"), c("Back", "Leg"), c("Hand", "Leg"))
  paired_block <- function(sens) {
    vapply(site_pairs, function(pr) {
      pv <- paired_by_sid(mydata, metric, pr[1], pr[2], sens = sens)
      c(w = wilcox_p(pv$x, pv$y, paired = TRUE),
        t = tt_p(pv$x, pv$y, paired = TRUE),
        fc = truefc(mean(pv$x, na.rm = TRUE) / mean(pv$y, na.rm = TRUE)))
    }, numeric(3))
  }
  pb_all <- paired_block(NULL)
  pb_S   <- paired_block("Sensitive")
  pb_N   <- paired_block("NotSensitive")
  colnames(pb_all) <- colnames(pb_S) <- colnames(pb_N) <-
    c("Back_Hand", "Back_Leg", "Hand_Leg")

  tibble(
    metric = metric,
    KP_BodySite, KP_Sensitivity, KP_SiteSensitivity,
    ## SS-10 correlations
    spearman_rho_Back = cors$Back["sp_rho"], spearman_p_Back = cors$Back["sp_p"],
    spearman_rho_Hand = cors$Hand["sp_rho"], spearman_p_Hand = cors$Hand["sp_p"],
    spearman_rho_Leg  = cors$Leg["sp_rho"],  spearman_p_Leg  = cors$Leg["sp_p"],
    pearson_r_Back = cors$Back["pe_r"], pearson_p_Back = cors$Back["pe_p"],
    pearson_r_Hand = cors$Hand["pe_r"], pearson_p_Hand = cors$Hand["pe_p"],
    pearson_r_Leg  = cors$Leg["pe_r"],  pearson_p_Leg  = cors$Leg["pe_p"],
    ## Sensitive vs NotSensitive
    SvN_Back_wilcox = sv$Back["wilcox"], SvN_Back_t = sv$Back["t"], SvN_Back_fc = sv$Back["fc"],
    SvN_Hand_wilcox = sv$Hand["wilcox"], SvN_Hand_t = sv$Hand["t"], SvN_Hand_fc = sv$Hand["fc"],
    SvN_Leg_wilcox  = sv$Leg["wilcox"],  SvN_Leg_t  = sv$Leg["t"],  SvN_Leg_fc  = sv$Leg["fc"],
    SvN_All_wilcox  = wilcox_p(gS, gN),  SvN_All_t  = tt_p(gS, gN),
    SvN_All_fc = truefc(mean(gS, na.rm = TRUE) / mean(gN, na.rm = TRUE)),
    ## paired cross-site (all subjects)
    All_BackHand_w = pb_all["w", 1], All_BackLeg_w = pb_all["w", 2], All_HandLeg_w = pb_all["w", 3],
    All_BackHand_fc = pb_all["fc", 1], All_BackLeg_fc = pb_all["fc", 2], All_HandLeg_fc = pb_all["fc", 3],
    ## paired cross-site (Sensitive only / NotSensitive only)
    S_BackHand_w = pb_S["w", 1], S_BackLeg_w = pb_S["w", 2], S_HandLeg_w = pb_S["w", 3],
    N_BackHand_w = pb_N["w", 1], N_BackLeg_w = pb_N["w", 2], N_HandLeg_w = pb_N["w", 3]
  )
}

alpha_stats <- bind_rows(lapply(alpha_metrics, alpha_stat_row))

## Group-level summary (n, mean, sd, se, median) per SiteSensitivity per metric.
alpha_summary <- mydata %>%
  select(SiteSensitivity, BodySite, Sensitivity, all_of(alpha_metrics)) %>%
  pivot_longer(all_of(alpha_metrics), names_to = "metric", values_to = "value") %>%
  group_by(metric, BodySite, Sensitivity, SiteSensitivity) %>%
  summarise(n = sum(!is.na(value)), mean = mean(value, na.rm = TRUE),
            sd = sd(value, na.rm = TRUE), se = sd / sqrt(n),
            median = median(value, na.rm = TRUE), .groups = "drop")

## Manuscript-friendly long-form SS-10 correlation table. Spearman is the safer
## primary correlation for ordinal/questionnaire data; Pearson is included and
## plotted because it was specifically requested.
alpha_correlations <- bind_rows(lapply(alpha_metrics, function(metric) {
  bind_rows(lapply(sites, function(s) {
    dd <- mydata %>%
      filter(BodySite == s) %>%
      transmute(score = .data[[SS10]], value = .data[[metric]]) %>%
      filter(is.finite(score), is.finite(value))

    bind_rows(lapply(c("spearman", "pearson"), function(method) {
      ct <- tryCatch({
        if (method == "spearman") {
          cor.test(dd$score, dd$value, method = method, exact = FALSE)
        } else {
          cor.test(dd$score, dd$value, method = method)
        }
      }, error = function(e) NULL)
      tibble(metric = metric,
             BodySite = s,
             method = method,
             n = nrow(dd),
             estimate = if (is.null(ct)) NA_real_ else unname(ct$estimate),
             conf_low = if (!is.null(ct) && !is.null(ct$conf.int)) ct$conf.int[1] else NA_real_,
             conf_high = if (!is.null(ct) && !is.null(ct$conf.int)) ct$conf.int[2] else NA_real_,
             p_value = if (is.null(ct)) NA_real_ else ct$p.value)
    }))
  }))
}))

## Stats paired with the two new figures below (3c2 / 3c3):
##  (i)  Sensitive vs NotSensitive pooled across all three body sites — this is
##       the same pooled Sensitive/NotSensitive contrast already summarised in
##       alpha_stats (SvN_All_*); pulled out here as its own tidy table so it is
##       reported alongside the dedicated "all sites pooled" figure.
##  (ii) BodySite (Kruskal-Wallis) computed SEPARATELY within the Sensitive and
##       within the NotSensitive subjects — this is new: alpha_stats only had
##       the BodySite Kruskal test pooled across sensitivity.
alpha_overall_sens_stats <- alpha_stats %>%
  transmute(metric,
            wilcox_p = SvN_All_wilcox, t_p = SvN_All_t, fold_change = SvN_All_fc)

alpha_bodysite_by_sens_stats <- bind_rows(lapply(alpha_metrics, function(metric) {
  tibble(
    metric = metric,
    KW_BodySite_Sensitive    = kw_p(mydata[mydata$Sensitivity == "Sensitive", ],
                                     metric, "BodySite"),
    KW_BodySite_NotSensitive = kw_p(mydata[mydata$Sensitivity == "NotSensitive", ],
                                     metric, "BodySite"))
}))

## ---------------------------------------------------------------------------
## 3. Alpha-diversity figures (publication ready)
## ---------------------------------------------------------------------------
## 3a. Per-metric boxplot: BodySite on x, Sensitivity by colour, with the
##     within-site Sensitive-vs-NotSensitive p-value.
alpha_box <- function(metric, test_method = "wilcox.test",
                      test_subtitle = "Sensitive:NotSensitive (Wilcoxon)") {
  ggplot(mydata, aes(BodySite, .data[[metric]], fill = Sensitivity)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.75, width = 0.65,
                 position = position_dodge(0.75)) +
    geom_point(aes(colour = Sensitivity),
               position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.75),
               size = 1.6, alpha = 0.85, show.legend = FALSE) +
    stat_compare_means(aes(group = Sensitivity), method = test_method,
                       label = "p.format", size = 3.4,
                       label.y = max(mydata[[metric]], na.rm = TRUE) * 1.04) +
    scale_fill_manual(values = pal_sens) +
    scale_colour_manual(values = pal_sens) +
    labs(x = NULL, y = metric_labels[[metric]],
         title = metric_labels[[metric]],
         subtitle = test_subtitle) +
    theme_pub()
}
for (m in alpha_metrics) save_fig(alpha_box(m), paste0("alpha.", m), 6.5, 5)

## Requested parametric t-test version. Kept separate from the Wilcoxon figure
## because alpha-diversity distributions are often non-normal in small cohorts.
alpha_ttest_box <- function(metric) {
  alpha_box(metric, test_method = "t.test",
            test_subtitle = "Sensitive:NotSensitive (two-sample t-test)")
}
for (m in alpha_metrics) save_fig(alpha_ttest_box(m), paste0("alpha.", m, ".ttest"), 6.5, 5)

## 3b. Combined 4-panel figure of the core metrics.
alpha_panel <- ggarrange(plotlist = lapply(core_metrics, alpha_box),
                         ncol = 2, nrow = 1, common.legend = TRUE, legend = "top")
save_fig(alpha_panel, "alpha.core_panel", 7, 5)

alpha_ttest_panel <- ggarrange(plotlist = lapply(core_metrics, alpha_ttest_box),
                               ncol = 2, nrow = 1, common.legend = TRUE, legend = "top")
save_fig(alpha_ttest_panel, "alpha.core_ttest_panel", 7, 5)

## 3c. Site-level view (all subjects) with Kruskal-Wallis across the 3 sites.
alpha_site_box <- function(metric) {
  ggplot(mydata, aes(BodySite, .data[[metric]], fill = BodySite)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.6) +
    geom_jitter(width = 0.12, size = 1.6, alpha = 0.7) +
    stat_compare_means(method = "kruskal.test", label = "p.format", size = 3.6,
                       label.y = max(mydata[[metric]], na.rm = TRUE) * 1.08) +
    scale_fill_manual(values = pal_site) +
    labs(x = NULL, y = metric_labels[[metric]], title = metric_labels[[metric]], subtitle = "Bodysite (Kruskal)") +
    theme_pub() + theme(legend.position = "none")
}
alpha_site_panel <- ggarrange(plotlist = lapply(core_metrics, alpha_site_box),
                              ncol = 2, nrow = 1)
save_fig(alpha_site_panel, "alpha.bySite_panel", 6, 4)

## 3c2. All sites pooled: Sensitive vs NotSensitive (ignores BodySite).
## Complements 3a (which is stratified per body site) with the simple
## between-subject contrast; stats are in alpha_overall_sens_stats.
alpha_overall_sens_box <- function(metric, test_method = "wilcox.test",
                                   test_subtitle = "All sites pooled (Wilcoxon)") {
  ggplot(mydata, aes(Sensitivity, .data[[metric]], fill = Sensitivity)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.55) +
    geom_jitter(width = 0.12, size = 1.6, alpha = 0.75, show.legend = FALSE) +
    stat_compare_means(method = test_method, label = "p.format", size = 3.6,
                       label.y = max(mydata[[metric]], na.rm = TRUE) * 1.08) +
    scale_fill_manual(values = pal_sens) +
    labs(x = NULL, y = metric_labels[[metric]], title = metric_labels[[metric]],
         subtitle = test_subtitle) +
    theme_pub() + theme(legend.position = "none")
}
for (m in alpha_metrics)
  save_fig(alpha_overall_sens_box(m), paste0("alpha.", m, ".overall_SensVsNotSens"), 5, 5)

alpha_overall_sens_panel <- ggarrange(plotlist = lapply(core_metrics, alpha_overall_sens_box),
                                      ncol = 2, nrow = 1)
save_fig(alpha_overall_sens_panel, "alpha.core_overall_SensVsNotSens_panel", 7, 5)

## 3c3. BodySite differences shown separately for Sensitive and NotSensitive
## subjects (two panels, one per Sensitivity group), with the Kruskal-Wallis
## p-value for the BodySite effect displayed within each panel; stats are in
## alpha_bodysite_by_sens_stats.
alpha_site_by_sens_box <- function(metric) {
  ggplot(mydata, aes(BodySite, .data[[metric]], fill = BodySite)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.6) +
    geom_jitter(width = 0.12, size = 1.6, alpha = 0.7) +
    stat_compare_means(method = "kruskal.test", label = "p.format", size = 3.6,
                       label.y = max(mydata[[metric]], na.rm = TRUE) * 1.08) +
    facet_wrap(~ Sensitivity) +
    scale_fill_manual(values = pal_site) +
    labs(x = NULL, y = metric_labels[[metric]], title = metric_labels[[metric]],
         subtitle = "Bodysite within each Sensitivity group (Kruskal)") +
    theme_pub() + theme(legend.position = "none")
}
for (m in alpha_metrics)
  save_fig(alpha_site_by_sens_box(m), paste0("alpha.", m, ".bySite_withinSensitivity"), 8, 5)

## Compact multi-metric version (core metrics) of the same 2-panel comparison,
## laid out as metric (rows) x Sensitivity (columns) for a quick overview.
alpha_site_by_sens_df <- mydata %>%
  select(Sensitivity, BodySite, all_of(core_metrics)) %>%
  pivot_longer(all_of(core_metrics), names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric, levels = core_metrics, labels = metric_labels[core_metrics]))
alpha_site_by_sens_compact <- ggplot(alpha_site_by_sens_df, aes(BodySite, value, fill = BodySite)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.6) +
  geom_jitter(width = 0.12, size = 1.4, alpha = 0.7) +
  stat_compare_means(method = "kruskal.test", label = "p.format", size = 3.2) +
  facet_grid(metric ~ Sensitivity, scales = "free_y", switch = "y") +
  scale_fill_manual(values = pal_site) +
  labs(x = NULL, y = NULL,
       title = "Body site differences within each Sensitivity group",
       subtitle = "Kruskal-Wallis p-value per panel") +
  theme_pub() + theme(legend.position = "none", strip.placement = "outside")
save_fig(alpha_site_by_sens_compact, "alpha.core_bySite_withinSensitivity_panel", 7, 5.5)

## 3d. Association between alpha diversity and the SS-10 sensitivity score.
ss_df <- mydata %>%
  select(BodySite, all_of(core_metrics), all_of(SS10)) %>%
  rename(SS10 = all_of(SS10)) %>%
  pivot_longer(all_of(core_metrics), names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric, levels = core_metrics,
                         labels = metric_labels[core_metrics]))
ss10_plot <- ggplot(ss_df, aes(SS10, value, colour = BodySite, fill = BodySite)) +
  geom_point(size = 1.8, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.15, linewidth = 0.7) +
  stat_cor(method = "spearman", size = 3, label.x.npc = "left", show.legend = FALSE,
           p.accuracy = 0.001, r.accuracy = 0.01) +
  facet_grid(metric ~ BodySite, scales = "free_y", switch = "y") +
  scale_colour_manual(values = pal_site) + scale_fill_manual(values = pal_site) +
  labs(x = "SS-10 total sensitivity score", y = NULL,
       title = "Alpha diversity vs SS-10 sensitivity score",
       subtitle = "Spearman correlation per body site") +
  theme_pub() + theme(legend.position = "none",
                      strip.placement = "outside")
save_fig(ss10_plot, "alpha.SS10_correlation", 9, 5)

ss10_plot_pearson <- ggplot(ss_df, aes(SS10, value, colour = BodySite, fill = BodySite)) +
  geom_point(size = 1.8, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.15, linewidth = 0.7) +
  stat_cor(method = "pearson", size = 3, label.x.npc = "left", show.legend = FALSE,
           p.accuracy = 0.001, r.accuracy = 0.01) +
  facet_grid(metric ~ BodySite, scales = "free_y", switch = "y") +
  scale_colour_manual(values = pal_site) + scale_fill_manual(values = pal_site) +
  labs(x = "SS-10 total sensitivity score", y = NULL,
       title = "Alpha diversity vs SS-10 sensitivity score",
       subtitle = "Pearson correlation per body site (requested parametric analysis)") +
  theme_pub() + theme(legend.position = "none",
                      strip.placement = "outside")
save_fig(ss10_plot_pearson, "alpha.SS10_pearson_correlation", 9, 5)

## ---------------------------------------------------------------------------
## 4. Beta diversity
## ---------------------------------------------------------------------------
## Bray-Curtis on relative abundances (essential: sequencing depth spans
## 813 -> 1.5e7 reads across samples).
comm     <- test_d[rowSums(test_d) > 0, , drop = FALSE]
meta_b   <- mydata[match(rownames(comm), mydata$shortID), ]
comm_ra  <- decostand(comm, method = "total")
bray     <- vegdist(comm_ra, method = "bray")

## 4a. NMDS ordination
nmds   <- metaMDS(bray, k = 2, trymax = 200, trace = 0)
stress <- nmds$stress
ord <- as.data.frame(scores(nmds, display = "sites")) %>%
  rownames_to_column("shortID") %>%
  left_join(mydata, by = "shortID")
stress_lab <- sprintf("NMDS (Bray-Curtis) — stress = %.3f", stress)

## Convex hulls (stat_chull) are used instead of parametric stat_ellipse
## throughout: with only ~11 subjects per group per site, a normal-theory
## ellipse can be a poor / misleading approximation, whereas a hull is a
## purely descriptive, distribution-free envelope of the observed points.
nmds_site <- ggplot(ord, aes(NMDS1, NMDS2, colour = BodySite)) +
  stat_chull(aes(fill = BodySite), geom = "polygon", alpha = 0.12, colour = NA) +
  geom_point(aes(shape = Sensitivity), size = 2.8, alpha = 0.9) +
  scale_colour_manual(values = pal_site) + scale_fill_manual(values = pal_site) +
  labs(title = "Community structure by body site", subtitle = stress_lab) +
  theme_pub()
save_fig(nmds_site, "beta.NMDS_bySite", 7, 6)

nmds_sens <- ggplot(ord, aes(NMDS1, NMDS2, colour = Sensitivity)) +
  stat_chull(aes(fill = Sensitivity), geom = "polygon", alpha = 0.12, colour = NA) +
  geom_point(size = 2.8, alpha = 0.9) +
  scale_colour_manual(values = pal_sens) + scale_fill_manual(values = pal_sens) +
  labs(title = "Community structure by sensitivity", subtitle = stress_lab) +
  theme_pub()
save_fig(nmds_sens, "beta.NMDS_bySensitivity", 7, 6)

nmds_facet <- nmds_site + facet_wrap(~ Sensitivity) +
  labs(title = "Community structure by body site, split by sensitivity")
save_fig(nmds_facet, "beta.NMDS_site_by_sensitivity", 10, 5.5)

## Reference-inspired beta-diversity views from GSS3318Diversity.R, upgraded for
## manuscript use: ordination is still Bray-Curtis on relative abundance rather
## than pseudocount raw counts, and convex hulls are descriptive overlays only.
nmds_sitesens <- ggplot(ord, aes(NMDS1, NMDS2,
                                 colour = SiteSensitivity,
                                 fill = SiteSensitivity)) +
  stat_chull(alpha = 0.18, geom = "polygon", linewidth = 0.35) +
  geom_point(aes(shape = Sensitivity), size = 2.4, alpha = 0.9) +
  scale_colour_manual(values = pal_ss) + scale_fill_manual(values = pal_ss) +
  labs(title = "Community structure by body site and sensitivity",
       subtitle = stress_lab) +
  theme_pub()
save_fig(nmds_sitesens, "beta.NMDS_SiteSensitivity", 7, 5.5)

nmds_sitesens_by_site <- nmds_sitesens +
  facet_wrap(~ BodySite) +
  labs(title = "Community structure by sensitivity within each body site")
save_fig(nmds_sitesens_by_site, "beta.NMDS_SiteSensitivity_bySite", 9, 5.5)

nmds_sitesens_by_sensitivity <- nmds_sitesens +
  facet_wrap(~ Sensitivity) +
  labs(title = "Community structure by body site within each sensitivity group")
save_fig(nmds_sitesens_by_sensitivity, "beta.NMDS_SiteSensitivity_bySensitivity", 9, 5.5)

nmds_sid <- ggplot(ord, aes(NMDS1, NMDS2, colour = SID, fill = SID)) +
  stat_chull(alpha = 0.10, geom = "polygon", linewidth = 0.25, show.legend = FALSE) +
  geom_point(size = 2.1, alpha = 0.85, show.legend = FALSE) +
  facet_wrap(~ Sensitivity) +
  labs(title = "Within-subject sample triplets across body sites",
       subtitle = paste(stress_lab, "| each hull connects one subject's three sites")) +
  theme_pub() + theme(legend.position = "none")
save_fig(nmds_sid, "beta.NMDS_bySID", 10, 5.5)

## 4b. PCoA (classical MDS) with variance explained on the axes.
pcoa    <- cmdscale(bray, k = 2, eig = TRUE)
eig     <- pcoa$eig
var_expl <- 100 * eig[1:2] / sum(eig[eig > 0])
pcoa_df <- as.data.frame(pcoa$points) %>%
  setNames(c("PCoA1", "PCoA2")) %>%
  rownames_to_column("shortID") %>%
  left_join(mydata, by = "shortID")
pcoa_plot <- ggplot(pcoa_df, aes(PCoA1, PCoA2, colour = BodySite)) +
  stat_chull(aes(fill = BodySite), geom = "polygon", alpha = 0.12, colour = NA) +
  geom_point(aes(shape = Sensitivity), size = 2.8, alpha = 0.9) +
  scale_colour_manual(values = pal_site) + scale_fill_manual(values = pal_site) +
  labs(title = "PCoA (Bray-Curtis)",
       x = sprintf("PCoA1 (%.1f%%)", var_expl[1]),
       y = sprintf("PCoA2 (%.1f%%)", var_expl[2])) +
  theme_pub()
save_fig(pcoa_plot, "beta.PCoA_bySite", 7, 6)

## 4c. PERMANOVA (adonis2)
## Original PERMANOVA analyses retained, with explicit design notes:
##  * BodySite is within-subject  -> restrict permutations within SID (blocks)
##  * Sensitivity is between-subject -> standard permutation
##  * BodySite * Sensitivity and SiteSensitivity are useful omnibus/descriptive
##    models, but for manuscript inference they should be interpreted together
##    with the subject-restricted BodySite test, within-site sensitivity tests
##    and dispersion diagnostics because this is a repeated-measures design.
ctrl_site <- how(blocks = meta_b$SID, nperm = 999)
perm_site <- adonis2(bray ~ BodySite, data = meta_b, permutations = ctrl_site)
perm_sens <- adonis2(bray ~ Sensitivity, data = meta_b, permutations = 999)
perm_int  <- adonis2(bray ~ BodySite * Sensitivity, data = meta_b,
                     by = "margin", permutations = 999)
perm_ss   <- adonis2(bray ~ SiteSensitivity, data = meta_b, permutations = 999)
perm_sid  <- adonis2(bray ~ SID, data = meta_b, permutations = 999)

tidy_adonis <- function(fit, model) {
  as.data.frame(fit) %>%
    rownames_to_column("term") %>%
    filter(!term %in% c("Residual", "Total")) %>%
    transmute(model, term,
              Df, R2 = round(R2, 4), F = round(`F`, 3),
              p_value = `Pr(>F)`)
}
permanova_tbl <- bind_rows(
  tidy_adonis(perm_site, "BodySite (subject-restricted)"),
  tidy_adonis(perm_sens, "Sensitivity"),
  tidy_adonis(perm_int,  "BodySite * Sensitivity (marginal)"),
  tidy_adonis(perm_ss,   "SiteSensitivity"),
  tidy_adonis(perm_sid,  "SID (subject)"))

## 4d. Multivariate homogeneity of dispersion (interpretation aid for PERMANOVA)
disp_tbl <- bind_rows(lapply(
  list(BodySite = meta_b$BodySite,
       Sensitivity = meta_b$Sensitivity,
       SiteSensitivity = meta_b$SiteSensitivity),
  function(g) {
    bd <- betadisper(bray, g)
    pt <- permutest(bd, permutations = 999)
    data.frame(F = round(pt$tab$F[1], 3), p_value = pt$tab$`Pr(>F)`[1])
  }), .id = "grouping")

## 4e. Pairwise PERMANOVA
pairwise_adonis <- function(groups, block = NULL, label) {
  lvls <- levels(droplevels(as.factor(groups)))
  cmb  <- combn(lvls, 2)
  res  <- lapply(seq_len(ncol(cmb)), function(i) {
    keep <- groups %in% cmb[, i]
    sub_d <- as.dist(as.matrix(bray)[keep, keep])
    meta_sub <- meta_b[keep, , drop = FALSE]
    ctrl <- if (is.null(block)) 999 else how(blocks = droplevels(meta_sub[[block]]), nperm = 999)
    fit <- adonis2(sub_d ~ groups[keep], permutations = ctrl)
    data.frame(comparison = paste(cmb[1, i], "vs", cmb[2, i]),
               R2 = round(fit$R2[1], 4), F = round(fit$`F`[1], 3),
               p_value = fit$`Pr(>F)`[1])
  })
  bind_rows(res) %>%
    mutate(p_adj_BH = round(p.adjust(p_value, "BH"), 4), grouping = label)
}
## body sites: paired -> restrict permutations within subject
pw_site <- pairwise_adonis(meta_b$BodySite, block = "SID", label = "BodySite (paired)")
## Sensitive vs NotSensitive within each site (between-subject)
pw_sens_site <- bind_rows(lapply(sites, function(s) {
  keep <- meta_b$BodySite == s
  sub_d <- as.dist(as.matrix(bray)[keep, keep])
  fit <- adonis2(sub_d ~ meta_b$Sensitivity[keep], permutations = 999)
  data.frame(comparison = paste0(s, ": Sensitive vs NotSensitive"),
             R2 = round(fit$R2[1], 4), F = round(fit$`F`[1], 3),
             p_value = fit$`Pr(>F)`[1])
})) %>% mutate(p_adj_BH = round(p.adjust(p_value, "BH"), 4),
               grouping = "Sensitivity within site")
## All sites pooled: Sensitive vs NotSensitive (between-subject; all samples,
## ignoring BodySite). This is the same overall Sensitivity contrast as
## perm_sens in permanova_tbl, added here as its own row so the "all sites
## pooled" Sensitive-vs-NotSensitive result is easy to find alongside the
## other pairwise comparisons.
pw_sens_all <- pairwise_adonis(meta_b$Sensitivity, block = NULL,
                               label = "Sensitivity (all sites pooled)") %>%
  mutate(comparison = "All sites: Sensitive vs NotSensitive")
## all 6 SiteSensitivity groups
pw_ss <- pairwise_adonis(meta_b$SiteSensitivity, block = NULL, label = "SiteSensitivity")

beta_pairwise <- bind_rows(pw_site, pw_sens_site, pw_sens_all, pw_ss)

## ---------------------------------------------------------------------------
## 5. Export: Excel workbook; figures were already saved above
## ---------------------------------------------------------------------------
wb <- createWorkbook()
addWorksheet(wb, "alpha_perSample"); writeData(wb, "alpha_perSample", alpha_tbl)
addWorksheet(wb, "alpha_summary");   writeData(wb, "alpha_summary", alpha_summary)
addWorksheet(wb, "alpha_stats");     writeData(wb, "alpha_stats", alpha_stats)
addWorksheet(wb, "alpha_correlations"); writeData(wb, "alpha_correlations", alpha_correlations)
addWorksheet(wb, "alpha_overall_SensVsNotSens"); writeData(wb, "alpha_overall_SensVsNotSens", alpha_overall_sens_stats)
addWorksheet(wb, "alpha_bySite_withinSens");      writeData(wb, "alpha_bySite_withinSens", alpha_bodysite_by_sens_stats)
addWorksheet(wb, "beta_PERMANOVA");  writeData(wb, "beta_PERMANOVA", permanova_tbl)
addWorksheet(wb, "beta_dispersion"); writeData(wb, "beta_dispersion", disp_tbl)
addWorksheet(wb, "beta_pairwise");   writeData(wb, "beta_pairwise", beta_pairwise)
saveWorkbook(wb, out(filename, ".Diversity.xlsx"), overwrite = TRUE)

## ---------------------------------------------------------------------------
## Console summary
## ---------------------------------------------------------------------------
cat("\n==================== DONE ====================\n")
cat(sprintf("Samples: %d | Taxa kept: %d (of %d)\n", nrow(test_d), ncol(test_d), d[1]))
cat(sprintf("NMDS stress: %.3f\n", stress))
cat("\nPERMANOVA:\n"); print(permanova_tbl, row.names = FALSE)
cat("\nBeta dispersion:\n"); print(disp_tbl, row.names = FALSE)
cat("\nOutputs written to: ", normalizePath(outdir), "\n")
