###############################################################################
## GSS3318 SensiScope — genus-level abundance analysis
## -----------------------------------------------------------------------------
## New script that leaves get_paired_wilcox_site.R unchanged.
##
## Purpose
##   Analyze the genus abundance table in the same spirit as the alpha-diversity
##   workflow: grouped abundance summaries, Sensitivity contrasts, paired
##   BodySite contrasts, SS-10 correlations, figures, and Excel statistics.
##
## Study design encoded by sample names
##   Sample columns must look like: BodySite.Sensitivity.SID
##   Example: Back.NotSensitive.1
##
## Outputs
##   A timestamped folder:
##     <genus_filename>_genus_result_<YYYYmmdd_HHMMSS>/
##   containing:
##     * PNG figures
##     * TSV statistics tables
##     * Excel workbook with all statistics
##
## Usage
##   Rscript get_paired_wilcox_site_genus_analysis.R [genus_file] [metadata_file] [top_n] [outdir] [plot_threshold] [threshold_type]
##
##   threshold_type can be "p" or "fdr".  The default is raw p-value filtering
##   with threshold 0.10.  If the 4th argument is numeric and no outdir is
##   supplied, it is interpreted as plot_threshold and the timestamped output
##   folder is used.
##
## Examples
##   Rscript get_paired_wilcox_site_genus_analysis.R
##   Rscript get_paired_wilcox_site_genus_analysis.R genus_cleandata GSS3318_readsCount_meta.txt 25
##   Rscript get_paired_wilcox_site_genus_analysis.R genus_cleandata GSS3318_readsCount_meta.txt 20 0.1
##   Rscript get_paired_wilcox_site_genus_analysis.R genus_cleandata GSS3318_readsCount_meta.txt 20 my_results 0.05
##   Rscript get_paired_wilcox_site_genus_analysis.R genus_cleandata GSS3318_readsCount_meta.txt 20 my_results_fdr 0.10 fdr
###############################################################################

rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(openxlsx)
  library(scales)
})

set.seed(1234)

## ---------------------------------------------------------------------------
## 0. User options and shared helper functions
## ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

find_file <- function(path) {
  ## Try the current directory first and then the parent directory.  This makes
  ## the script usable from either Stat/ or Stat/code1/.
  candidates <- unique(c(path, file.path(".", path), file.path("..", path)))
  hit <- candidates[file.exists(candidates)]
  if (length(hit) == 0) {
    stop("Cannot find file: ", path, "\nChecked: ",
         paste(candidates, collapse = ", "), call. = FALSE)
  }
  hit[1]
}

genus_file <- find_file(if (length(args) >= 1) args[1] else "genus_cleandata")
meta_file  <- find_file(if (length(args) >= 2) args[2] else "GSS3318_readsCount_meta.txt")
top_n      <- if (length(args) >= 3) suppressWarnings(as.integer(args[3])) else 20L
if (is.na(top_n) || top_n < 1) top_n <- 20L

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
default_outdir <- paste0(basename(genus_file), "_genus_result_", timestamp)
is_numeric_arg <- function(x) !is.na(suppressWarnings(as.numeric(x)))
if (length(args) >= 4 && is_numeric_arg(args[4]) && length(args) < 5) {
  ## Backward-compatible shortcut: make the threshold easy to change without
  ## forcing the user to provide a custom output directory.
  outdir <- default_outdir
  plot_threshold <- as.numeric(args[4])
  threshold_type <- "p"
} else {
  outdir <- if (length(args) >= 4) args[4] else default_outdir
  plot_threshold <- if (length(args) >= 5) suppressWarnings(as.numeric(args[5])) else 0.10
  threshold_type <- if (length(args) >= 6) tolower(args[6]) else "p"
}
if (!is.finite(plot_threshold) || plot_threshold <= 0 || plot_threshold > 1) {
  stop("plot_threshold must be a number between 0 and 1.", call. = FALSE)
}
if (!threshold_type %in% c("p", "fdr")) {
  stop("threshold_type must be either 'p' or 'fdr'.", call. = FALSE)
}
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
out <- function(...) file.path(outdir, paste0(...))

## Metadata score used for the requested correlation analysis.  Because
## read.table(check.names = TRUE) sanitizes column names, "SS-10" becomes
## "SS.10" in R.
ss10_col <- "All.SS10.Total_score_of_all_SS.10_questions"

sites <- c("Back", "Hand", "Leg")
sens_levels <- c("NotSensitive", "Sensitive")
site_sens_levels <- c("BackNotSensitive", "BackSensitive",
                      "HandNotSensitive", "HandSensitive",
                      "LegNotSensitive", "LegSensitive")

pal_sens <- c(NotSensitive = "#4C72B0", Sensitive = "#C44E52")
pal_site <- c(Back = "#55A868", Hand = "#8172B3", Leg = "#DD8452")

theme_pub <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = "grey92"),
          strip.background = element_rect(fill = "grey95", colour = NA),
          strip.text = element_text(face = "bold"),
          plot.title = element_text(face = "bold", hjust = 0),
          legend.position = "top")
}

save_fig <- function(plot, name, width, height) {
  ggsave(out(basename(genus_file), ".", name, ".png"), plot,
         width = width, height = height, dpi = 300, bg = "white")
}

safe_p <- function(expr) tryCatch(expr, error = function(e) NA_real_)

wilcox_p <- function(x, y, paired = FALSE) {
  ## For paired tests, x and y must already be aligned by subject.  For unpaired
  ## tests, finite values are kept independently.
  if (paired) {
    keep <- is.finite(x) & is.finite(y)
    x <- x[keep]
    y <- y[keep]
  } else {
    x <- x[is.finite(x)]
    y <- y[is.finite(y)]
  }
  if (length(x) < 2 || length(y) < 2) return(NA_real_)
  safe_p(suppressWarnings(wilcox.test(x, y, paired = paired, exact = FALSE)$p.value))
}

tt_p <- function(x, y, paired = FALSE) {
  if (paired) {
    keep <- is.finite(x) & is.finite(y)
    x <- x[keep]
    y <- y[keep]
  } else {
    x <- x[is.finite(x)]
    y <- y[is.finite(y)]
  }
  if (length(x) < 2 || length(y) < 2) return(NA_real_)
  safe_p(t.test(x, y, paired = paired)$p.value)
}

kw_p <- function(data, value_col, group_col) {
  dd <- data %>%
    filter(is.finite(.data[[value_col]]), !is.na(.data[[group_col]])) %>%
    mutate(group_tmp = droplevels(as.factor(.data[[group_col]])))
  if (nrow(dd) < 3 || n_distinct(dd$group_tmp) < 2) return(NA_real_)
  safe_p(kruskal.test(reformulate("group_tmp", response = value_col), data = dd)$p.value)
}

cor_stats <- function(x, y, method) {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  if (length(x) < 3 || length(unique(x)) < 2 || length(unique(y)) < 2) {
    return(tibble(n = length(x), estimate = NA_real_, conf_low = NA_real_,
                  conf_high = NA_real_, p_value = NA_real_))
  }
  ct <- tryCatch({
    if (method == "spearman") {
      suppressWarnings(cor.test(x, y, method = method, exact = FALSE))
    } else {
      suppressWarnings(cor.test(x, y, method = method))
    }
  }, error = function(e) NULL)
  tibble(n = length(x),
         estimate = if (is.null(ct)) NA_real_ else unname(ct$estimate),
         conf_low = if (!is.null(ct) && !is.null(ct$conf.int)) ct$conf.int[1] else NA_real_,
         conf_high = if (!is.null(ct) && !is.null(ct$conf.int)) ct$conf.int[2] else NA_real_,
         p_value = if (is.null(ct)) NA_real_ else ct$p.value)
}

truefc <- function(ratio) {
  ## Signed fold-change: 2x higher = +2; 2x lower = -2.
  if (!is.finite(ratio) || ratio == 0) return(NA_real_)
  if (ratio < 1) return(-1 / ratio)
  ratio
}

sig_mark <- function(p_value, fdr_value) {
  ## Four-level annotation requested for figure heatmaps:
  ##   p < 0.10   = *
  ##   p < 0.05   = **
  ##   FDR < 0.10 = ***
  ##   FDR < 0.05 = ****
  ## FDR thresholds take priority over raw p-value thresholds.
  case_when(
    !is.na(fdr_value) & fdr_value < 0.05 ~ "****",
    !is.na(fdr_value) & fdr_value < 0.10 ~ "***",
    !is.na(p_value) & p_value < 0.05 ~ "**",
    !is.na(p_value) & p_value < 0.10 ~ "*",
    TRUE ~ ""
  )
}

clean_genus_name <- function(x) {
  x <- gsub("^[a-z]__", "", x)
  x <- gsub("_", " ", x)
  ifelse(is.na(x) | x == "", "Unclassified", x)
}

safe_filename <- function(x) {
  ## File names are built from genus names, groups and p-values.  Keep them
  ## portable and readable across operating systems.
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  substr(x, 1, 80)
}

safe_feature_filename <- function(x) {
  ## EC input rows often look like "EC_1.1.1.1: Alcohol dehydrogenase".
  ## Use only the EC number in plot filenames.  For non-EC features, fall back
  ## to a sanitized, length-limited feature ID.
  x <- as.character(x)
  ec_match <- regexpr("^EC_[0-9]+(\\.[0-9]+)*", x)
  ifelse(ec_match > 0, regmatches(x, ec_match), safe_filename(x))
}

p_label <- function(p) {
  ifelse(is.na(p), "NA", formatC(p, format = "e", digits = 2))
}

passes_threshold <- function(p_value, fdr_value) {
  ## Choose whether individual plots are selected by raw p-value or by FDR.
  ## This affects only which individual plots are made; all p-values and FDRs
  ## are still exported in the statistics tables.
  value <- if (threshold_type == "fdr") fdr_value else p_value
  is.finite(value) & value < plot_threshold
}

threshold_label <- function() {
  paste0(ifelse(threshold_type == "fdr", "FDR", "p"), " < ", plot_threshold)
}

loo_cor_diagnostics <- function(x, y, method, full_estimate) {
  ## Leave-one-out (LOO) correlation diagnostics help flag associations that are
  ## driven by one or two samples.  A stable signal should usually keep the same
  ## direction after each single sample is removed; if LOO p-values or estimates
  ## change dramatically, inspect the scatter plot before interpreting the hit.
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  n <- length(x)
  if (n < 4 || !is.finite(full_estimate) || full_estimate == 0) {
    return(tibble(loo_n = n, loo_same_direction_all = NA,
                  loo_min_estimate = NA_real_, loo_max_estimate = NA_real_,
                  loo_min_abs_estimate = NA_real_, loo_max_p = NA_real_))
  }
  loo <- lapply(seq_len(n), function(i) {
    cor_stats(x[-i], y[-i], method) %>%
      transmute(estimate, p_value)
  }) %>% bind_rows()
  tibble(loo_n = n,
         loo_same_direction_all = all(sign(loo$estimate) == sign(full_estimate), na.rm = TRUE),
         loo_min_estimate = min(loo$estimate, na.rm = TRUE),
         loo_max_estimate = max(loo$estimate, na.rm = TRUE),
         loo_min_abs_estimate = min(abs(loo$estimate), na.rm = TRUE),
         loo_max_p = max(loo$p_value, na.rm = TRUE))
}

## ---------------------------------------------------------------------------
## 1. Load genus data and metadata
## ---------------------------------------------------------------------------
genus_wide <- read.table(genus_file, header = TRUE, row.names = 1, sep = "\t",
                         quote = "", comment.char = "",
                         check.names = FALSE, stringsAsFactors = FALSE)

## Force abundance columns to numeric and treat missing/non-numeric abundance as
## zero.  The input genus table is expected to be a numeric abundance matrix.
genus_wide[] <- lapply(genus_wide, function(x) {
  z <- suppressWarnings(as.numeric(x))
  z[is.na(z)] <- 0
  z
})

sample_parts <- strsplit(colnames(genus_wide), "[.]")
bad_samples <- vapply(sample_parts, length, integer(1)) < 3
if (any(bad_samples)) {
  stop("Sample names must look like BodySite.Sensitivity.SID. Bad names: ",
       paste(colnames(genus_wide)[bad_samples], collapse = ", "), call. = FALSE)
}

sample_map <- tibble(
  sample_id = colnames(genus_wide),
  BodySite = vapply(sample_parts, `[`, character(1), 1),
  Sensitivity = vapply(sample_parts, `[`, character(1), 2),
  SID = vapply(sample_parts, `[`, character(1), 3),
  SiteSensitivity = paste0(BodySite, Sensitivity)
)

meta <- read.table(meta_file, header = TRUE, sep = "\t", check.names = TRUE,
                   stringsAsFactors = FALSE)
if (!ss10_col %in% colnames(meta)) {
  stop("Required metadata column is missing after R name cleanup: ", ss10_col,
       call. = FALSE)
}

meta2 <- meta %>%
  select(-any_of(c("BodySite", "Sensitivity", "SID", "SiteSensitivity"))) %>%
  rename(sample_id = NewID, SS10_total = all_of(ss10_col))

positive_min <- min(as.matrix(genus_wide)[as.matrix(genus_wide) > 0], na.rm = TRUE)
if (!is.finite(positive_min)) positive_min <- 1
pseudocount <- positive_min / 2

genus_long <- genus_wide %>%
  rownames_to_column("genus_id") %>%
  mutate(genus = clean_genus_name(genus_id)) %>%
  pivot_longer(-c(genus_id, genus), names_to = "sample_id",
               values_to = "abundance") %>%
  left_join(sample_map, by = "sample_id") %>%
  left_join(meta2, by = "sample_id") %>%
  mutate(BodySite = factor(BodySite, levels = sites),
         Sensitivity = factor(Sensitivity, levels = sens_levels),
         SiteSensitivity = factor(SiteSensitivity, levels = site_sens_levels),
         SID = factor(SID),
         abundance_log10 = log10(abundance + pseudocount))

missing_meta <- genus_long %>% filter(is.na(SS10_total)) %>% distinct(sample_id)
if (nrow(missing_meta) > 0) {
  warning("Metadata was not found for sample(s): ",
          paste(missing_meta$sample_id, collapse = ", "))
}

taxa_rank <- genus_long %>%
  group_by(genus_id, genus) %>%
  summarise(mean_abundance = mean(abundance, na.rm = TRUE),
            median_abundance = median(abundance, na.rm = TRUE),
            max_abundance = max(abundance, na.rm = TRUE),
            prevalence = mean(abundance > 0, na.rm = TRUE),
            percent_zero = mean(abundance == 0 | is.na(abundance)),
            .groups = "drop") %>%
  arrange(desc(mean_abundance), desc(prevalence), genus)

top_taxa <- head(taxa_rank$genus_id, top_n)

## ---------------------------------------------------------------------------
## 2. Statistical analysis
## ---------------------------------------------------------------------------
genus_summary <- genus_long %>%
  group_by(genus_id, genus, BodySite, Sensitivity, SiteSensitivity) %>%
  summarise(n = sum(is.finite(abundance)),
            mean = mean(abundance, na.rm = TRUE),
            sd = sd(abundance, na.rm = TRUE),
            se = sd / sqrt(n),
            median = median(abundance, na.rm = TRUE),
            prevalence = mean(abundance > 0, na.rm = TRUE),
            .groups = "drop")

paired_site_values <- function(data, site_a, site_b, sens = NULL) {
  ## Align BodySite values by SID before paired tests.
  dd <- data
  if (!is.null(sens)) dd <- dd %>% filter(Sensitivity == sens)
  a <- dd %>% filter(BodySite == site_a) %>% select(SID, abundance)
  b <- dd %>% filter(BodySite == site_b) %>% select(SID, abundance)
  inner_join(a, b, by = "SID", suffix = c("_a", "_b"))
}

genus_stat_one <- function(gid) {
  dd <- genus_long %>% filter(genus_id == gid)
  gname <- unique(dd$genus)[1]

  omnibus <- tibble(
    genus_id = gid,
    genus = gname,
    KW_BodySite = kw_p(dd, "abundance", "BodySite"),
    KW_Sensitivity = kw_p(dd, "abundance", "Sensitivity"),
    KW_SiteSensitivity = kw_p(dd, "abundance", "SiteSensitivity"),
    percent_zero = mean(dd$abundance == 0 | is.na(dd$abundance)),
    mean_abundance = mean(dd$abundance, na.rm = TRUE),
    max_abundance = max(dd$abundance, na.rm = TRUE)
  )

  sensitivity_tests <- bind_rows(lapply(c(sites, "All"), function(site) {
    xx <- if (site == "All") dd else dd %>% filter(BodySite == site)
    s <- xx$abundance[xx$Sensitivity == "Sensitive"]
    n <- xx$abundance[xx$Sensitivity == "NotSensitive"]
    tibble(genus_id = gid, genus = gname, scope = site,
           test_family = "Sensitivity_within_site",
           group_a = "Sensitive", group_b = "NotSensitive",
           paired = FALSE,
           n_a = sum(is.finite(s)), n_b = sum(is.finite(n)),
           mean_a = mean(s, na.rm = TRUE), mean_b = mean(n, na.rm = TRUE),
           median_a = median(s, na.rm = TRUE), median_b = median(n, na.rm = TRUE),
           fold_change = truefc(mean(s, na.rm = TRUE) / mean(n, na.rm = TRUE)),
           log2_fc = log2((mean(s, na.rm = TRUE) + 1e-9) /
                          (mean(n, na.rm = TRUE) + 1e-9)),
           wilcox_p = wilcox_p(s, n, paired = FALSE),
           t_p = tt_p(s, n, paired = FALSE))
  }))

  site_pairs <- list(c("Back", "Hand"), c("Back", "Leg"), c("Hand", "Leg"))
  paired_tests <- bind_rows(lapply(c("All", sens_levels), function(scope) {
    bind_rows(lapply(site_pairs, function(pr) {
      m <- paired_site_values(dd, pr[1], pr[2],
                              sens = if (scope == "All") NULL else scope)
      x <- m$abundance_a
      y <- m$abundance_b
      tibble(genus_id = gid, genus = gname, scope = scope,
             test_family = "Paired_BodySite",
             group_a = pr[1], group_b = pr[2],
             paired = TRUE,
             n_a = sum(is.finite(x)), n_b = sum(is.finite(y)),
             mean_a = mean(x, na.rm = TRUE), mean_b = mean(y, na.rm = TRUE),
             median_a = median(x, na.rm = TRUE), median_b = median(y, na.rm = TRUE),
             fold_change = truefc(mean(x, na.rm = TRUE) / mean(y, na.rm = TRUE)),
             log2_fc = log2((mean(x, na.rm = TRUE) + 1e-9) /
                            (mean(y, na.rm = TRUE) + 1e-9)),
             wilcox_p = wilcox_p(x, y, paired = TRUE),
             t_p = tt_p(x, y, paired = TRUE))
    }))
  }))

  list(omnibus = omnibus, contrasts = bind_rows(sensitivity_tests, paired_tests))
}

stat_list <- lapply(taxa_rank$genus_id, genus_stat_one)

genus_omnibus <- bind_rows(lapply(stat_list, `[[`, "omnibus")) %>%
  mutate(KW_BodySite_FDR = p.adjust(KW_BodySite, method = "BH"),
         KW_Sensitivity_FDR = p.adjust(KW_Sensitivity, method = "BH"),
         KW_SiteSensitivity_FDR = p.adjust(KW_SiteSensitivity, method = "BH"))

genus_contrasts <- bind_rows(lapply(stat_list, `[[`, "contrasts")) %>%
  group_by(test_family, scope, group_a, group_b) %>%
  mutate(wilcox_FDR = p.adjust(wilcox_p, method = "BH"),
         t_FDR = p.adjust(t_p, method = "BH")) %>%
  ungroup()

genus_correlations <- genus_long %>%
  group_by(genus_id, genus, BodySite) %>%
  group_modify(~ bind_rows(
    cor_stats(.x$SS10_total, .x$abundance, "spearman") %>% mutate(method = "spearman"),
    cor_stats(.x$SS10_total, .x$abundance, "pearson") %>% mutate(method = "pearson")
  )) %>%
  ungroup() %>%
  select(genus_id, genus, BodySite, method, n, estimate, conf_low, conf_high, p_value) %>%
  group_by(BodySite, method) %>%
  mutate(p_FDR = p.adjust(p_value, method = "BH")) %>%
  ungroup()

## ---------------------------------------------------------------------------
## 3. Figures
## ---------------------------------------------------------------------------
top_label_levels <- taxa_rank$genus[match(top_taxa, taxa_rank$genus_id)]

plot_df <- genus_long %>%
  filter(genus_id %in% top_taxa) %>%
  mutate(genus = factor(genus, levels = top_label_levels))

top_box <- ggplot(plot_df, aes(BodySite, abundance, fill = Sensitivity)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.78, width = 0.65,
               position = position_dodge(0.75)) +
  geom_point(aes(colour = Sensitivity),
             position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.75),
             size = 1.3, alpha = 0.8, show.legend = FALSE) +
  facet_wrap(~ genus, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = pal_sens) +
  scale_colour_manual(values = pal_sens) +
  labs(x = NULL, y = "Relative abundance",
       title = paste0("Top ", length(top_taxa), " genera by body site and sensitivity"),
       subtitle = "Boxplots show genus abundance; points are individual samples") +
  theme_pub()
save_fig(top_box, "top_genera.boxplot_bodySite_sensitivity", 12, 8)

heat_df <- genus_summary %>%
  filter(genus_id %in% top_taxa) %>%
  mutate(genus = factor(genus, levels = top_label_levels))

mean_heat <- ggplot(heat_df, aes(SiteSensitivity, genus, fill = log10(mean + 1e-6))) +
  geom_tile(colour = "white", linewidth = 0.25) +
  scale_fill_gradientn(name = "log10(mean + 1e-6)",
                       colours = c("#F7FBFF", "#6BAED6", "#08306B")) +
  labs(x = NULL, y = NULL, title = "Mean abundance of top genera",
       subtitle = "Rows are ordered by overall mean abundance") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_fig(mean_heat, "top_genera.mean_abundance_heatmap", 8, 7)

sens_effect <- genus_contrasts %>%
  filter(test_family == "Sensitivity_within_site", genus_id %in% top_taxa) %>%
  mutate(genus = factor(genus, levels = top_label_levels),
         scope = factor(scope, levels = c(sites, "All")),
         significance = sig_mark(wilcox_p, wilcox_FDR))

effect_heat <- ggplot(sens_effect, aes(scope, genus, fill = log2_fc)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  geom_text(aes(label = significance), size = 4) +
  scale_fill_gradient2(name = "log2 FC\nSensitive/Not",
                       low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, oob = squish) +
  labs(x = NULL, y = NULL, title = "Sensitivity-associated genus shifts",
       subtitle = "* p<0.10; ** p<0.05; *** FDR<0.10; **** FDR<0.05") +
  theme_pub()
save_fig(effect_heat, "top_genera.sensitivity_log2FC_heatmap", 7, 7)

cor_plot_df <- genus_correlations %>%
  filter(method == "spearman", genus_id %in% top_taxa) %>%
  mutate(genus = factor(genus, levels = top_label_levels),
         significance = sig_mark(p_value, p_FDR))

cor_heat <- ggplot(cor_plot_df, aes(BodySite, genus, fill = estimate)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  geom_text(aes(label = significance), size = 4) +
  scale_fill_gradient2(name = "Spearman rho",
                       low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-1, 1), oob = squish) +
  labs(x = NULL, y = NULL, title = "Genus abundance vs SS-10 score",
       subtitle = "Spearman correlation by body site; * p<0.10; ** p<0.05; *** FDR<0.10; **** FDR<0.05") +
  theme_pub()
save_fig(cor_heat, "top_genera.SS10_spearman_heatmap", 6.5, 7)

pearson_plot_df <- genus_correlations %>%
  filter(method == "pearson", genus_id %in% top_taxa) %>%
  mutate(genus = factor(genus, levels = top_label_levels),
         significance = sig_mark(p_value, p_FDR))

pearson_heat <- ggplot(pearson_plot_df, aes(BodySite, genus, fill = estimate)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  geom_text(aes(label = significance), size = 4) +
  scale_fill_gradient2(name = "Pearson r",
                       low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-1, 1), oob = squish) +
  labs(x = NULL, y = NULL, title = "Genus abundance vs SS-10 score",
       subtitle = "Pearson correlation by body site; * p<0.10; ** p<0.05; *** FDR<0.10; **** FDR<0.05") +
  theme_pub()
save_fig(pearson_heat, "top_genera.SS10_pearson_heatmap", 6.5, 7)

scatter_taxa <- head(top_taxa, min(6, length(top_taxa)))
scatter_levels <- taxa_rank$genus[match(scatter_taxa, taxa_rank$genus_id)]

scatter_df <- genus_long %>%
  filter(genus_id %in% scatter_taxa) %>%
  mutate(genus = factor(genus, levels = scatter_levels))

ss10_scatter <- ggplot(scatter_df, aes(SS10_total, abundance, colour = BodySite, fill = BodySite)) +
  geom_point(size = 1.7, alpha = 0.85) +
  geom_smooth(aes(group = 1), method = "lm", se = TRUE, alpha = 0.15,
              linewidth = 0.65, colour = "black", fill = "grey70") +
  facet_grid(genus ~ BodySite, scales = "free_y") +
  scale_colour_manual(values = pal_site) +
  scale_fill_manual(values = pal_site) +
  labs(x = "SS-10 total sensitivity score", y = "Relative abundance",
       title = "Top genus abundance vs SS-10 score",
       subtitle = "Linear trend shown for visualization; exact Pearson/Spearman tests are in Excel") +
  theme_pub() +
  theme(legend.position = "none")
save_fig(ss10_scatter, "top6_genera.SS10_scatter", 10, 8)

pearson_scatter_df <- genus_long %>%
  filter(genus_id %in% top_taxa) %>%
  mutate(genus = factor(genus, levels = top_label_levels))

pearson_summary_scatter <- ggplot(
  pearson_scatter_df,
  aes(SS10_total, abundance, colour = BodySite, fill = BodySite)
) +
  geom_point(size = 1.4, alpha = 0.78) +
  geom_smooth(aes(group = BodySite), method = "lm", se = TRUE,
              alpha = 0.12, linewidth = 0.55) +
  facet_wrap(~ genus, scales = "free_y", ncol = 4) +
  scale_colour_manual(values = pal_site) +
  scale_fill_manual(values = pal_site) +
  labs(x = "SS-10 total sensitivity score", y = "Relative abundance",
       title = "Summary Pearson view: genus abundance vs SS-10 score",
       subtitle = paste0("Top ", length(top_taxa),
                         " genera; site-specific linear trends shown for visualization")) +
  theme_pub()
save_fig(pearson_summary_scatter, "top_genera.SS10_pearson_summary_scatter", 12, 9)

## ---------------------------------------------------------------------------
## 4. Individual plots for any raw p-value below the adjustable threshold
## ---------------------------------------------------------------------------
## These plots are meant for follow-up review of every result that "shows up"
## at the chosen exploratory threshold.  The default is raw p-value p < 0.10;
## pass threshold_type="fdr" to select by FDR instead.
##   Rscript get_paired_wilcox_site_genus_analysis.R genus_cleandata meta.txt 20 0.05
##   Rscript get_paired_wilcox_site_genus_analysis.R genus_cleandata meta.txt 20 out_fdr 0.10 fdr
individual_root <- file.path(outdir, "individual_plots")
comparison_root <- file.path(individual_root, "comparisons")
sensitivity_comparison_dir <- file.path(comparison_root, "sensitive_vs_not_sensitive")
bodysite_comparison_dir <- file.path(comparison_root, "body_site_differences")
correlation_dir <- file.path(individual_root, "correlations")
dir.create(sensitivity_comparison_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(bodysite_comparison_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(correlation_dir, recursive = TRUE, showWarnings = FALSE)

individual_plot_records <- list()
record_plot <- function(type, genus_id, genus, scope, test, p_value, fdr_value, file) {
  individual_plot_records[[length(individual_plot_records) + 1]] <<- tibble(
    plot_type = type,
    genus_id = genus_id,
    genus = genus,
    scope = scope,
    test = test,
    p_value = p_value,
    fdr_value = fdr_value,
    file = file
  )
}

## 4a. Individual comparison plots.
comparison_hits <- genus_contrasts %>%
  filter(passes_threshold(wilcox_p, wilcox_FDR)) %>%
  arrange(wilcox_p, genus, test_family, scope, group_a, group_b)

plotted_comparison_keys <- character()
for (ii in seq_len(nrow(comparison_hits))) {
  hit <- comparison_hits[ii, ]
  comparison_key <- paste(hit$test_family, hit$genus_id, sep = "::")
  if (comparison_key %in% plotted_comparison_keys) next
  plotted_comparison_keys <- c(plotted_comparison_keys, comparison_key)

  dd <- genus_long %>% filter(genus_id == hit$genus_id)

  if (hit$test_family == "Sensitivity_within_site") {
    ## Unpaired Sensitive-vs-NotSensitive comparison.  Every individual plot
    ## shows the three body sites plus an "All" panel with all samples pooled,
    ## so a site-specific hit can be compared with the overall pattern.
    plot_data <- bind_rows(
      dd %>% mutate(Panel = as.character(BodySite)),
      dd %>% mutate(Panel = "All")
    ) %>%
      mutate(Panel = factor(Panel, levels = c(sites, "All")))

    panel_stats <- genus_contrasts %>%
      filter(genus_id == hit$genus_id,
             test_family == "Sensitivity_within_site",
             scope %in% c(sites, "All")) %>%
      transmute(Panel = factor(scope, levels = c(sites, "All")),
                label = paste0("Wilcoxon p=", p_label(wilcox_p)))

    panel_y <- plot_data %>%
      group_by(Panel) %>%
      summarise(y = max(abundance, na.rm = TRUE), .groups = "drop") %>%
      mutate(y = ifelse(is.finite(y) & y > 0, y * 1.08, 0.05))
    panel_stats <- panel_stats %>%
      left_join(panel_y, by = "Panel") %>%
      mutate(x = 1.5)

    subtitle_text <- paste0(
      "One plot per genus; generated because at least one Sensitive-vs-NotSensitive contrast met ",
      threshold_label()
    )
    p <- ggplot(plot_data, aes(Sensitivity, abundance, fill = Sensitivity)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.78, width = 0.55) +
      geom_jitter(aes(colour = Sensitivity), width = 0.12,
                  size = 1.8, alpha = 0.85, show.legend = FALSE) +
      scale_fill_manual(values = pal_sens) +
      scale_colour_manual(values = pal_sens) +
      labs(x = NULL, y = "Relative abundance",
           title = paste0(hit$genus, ": Sensitive vs NotSensitive"),
           subtitle = subtitle_text) +
      geom_text(data = panel_stats, aes(x = x, y = y, label = label),
                inherit.aes = FALSE, size = 3.2) +
      facet_wrap(~ Panel, scales = "free_y", nrow = 1) +
      theme_pub() +
      theme(legend.position = "none")
  } else {
    ## Paired body-site comparison.  Once any paired site comparison passes the
    ## threshold for this genus, show all three body sites so the full subject
    ## trajectory is visible.  Separate panels clarify whether the threshold-
    ## passing result came from all subjects, Sensitive-only, or NotSensitive-only.
    plot_data <- bind_rows(
      dd %>% select(SID, Sensitivity, BodySite, abundance) %>% mutate(Panel = "All"),
      dd %>% filter(Sensitivity == "Sensitive") %>%
        select(SID, Sensitivity, BodySite, abundance) %>% mutate(Panel = "Sensitive"),
      dd %>% filter(Sensitivity == "NotSensitive") %>%
        select(SID, Sensitivity, BodySite, abundance) %>% mutate(Panel = "NotSensitive")
    ) %>%
      mutate(Panel = factor(Panel, levels = c("All", "Sensitive", "NotSensitive")),
             BodySite = factor(BodySite, levels = sites))

    threshold_stats <- genus_contrasts %>%
      filter(genus_id == hit$genus_id,
             test_family == "Paired_BodySite",
             passes_threshold(wilcox_p, wilcox_FDR)) %>%
      mutate(Panel = factor(scope, levels = c("All", "Sensitive", "NotSensitive")),
             comparison_label = paste0(group_a, " vs ", group_b,
                                       ": p=", p_label(wilcox_p),
                                       ", FDR=", p_label(wilcox_FDR))) %>%
      group_by(Panel) %>%
      arrange(wilcox_p, .by_group = TRUE) %>%
      mutate(label_rank = row_number()) %>%
      ungroup()

    panel_y <- plot_data %>%
      group_by(Panel) %>%
      summarise(y_base = max(abundance, na.rm = TRUE), .groups = "drop") %>%
      mutate(y_base = ifelse(is.finite(y_base) & y_base > 0, y_base, 0.05))
    threshold_stats <- threshold_stats %>%
      left_join(panel_y, by = "Panel") %>%
      mutate(x = 2,
             y = y_base * (1.08 + 0.12 * (label_rank - 1)))

    kw_hit <- genus_omnibus %>% filter(genus_id == hit$genus_id) %>% slice(1)
    subtitle_text <- paste0(
      "Kruskal-Wallis BodySite p=", p_label(kw_hit$KW_BodySite),
      ", FDR=", p_label(kw_hit$KW_BodySite_FDR),
      " | SiteSensitivity p=", p_label(kw_hit$KW_SiteSensitivity),
      ", FDR=", p_label(kw_hit$KW_SiteSensitivity_FDR),
      "\nPanel labels show paired Wilcoxon site comparisons with ",
      threshold_label()
    )
    p <- ggplot(plot_data, aes(BodySite, abundance, group = SID)) +
      geom_line(aes(colour = Sensitivity), alpha = 0.42, linewidth = 0.45) +
      geom_point(aes(colour = Sensitivity), size = 1.9, alpha = 0.9) +
      geom_boxplot(data = plot_data, mapping = aes(BodySite, abundance, group = BodySite),
                   inherit.aes = FALSE, width = 0.42, alpha = 0.16,
                   outlier.shape = NA) +
      geom_text(data = threshold_stats,
                aes(x = x, y = y, label = comparison_label),
                inherit.aes = FALSE, size = 3.0, lineheight = 0.9) +
      facet_wrap(~ Panel, scales = "free_y", nrow = 1) +
      scale_colour_manual(values = pal_sens) +
      labs(x = NULL, y = "Relative abundance",
           title = paste0(hit$genus, ": paired body-site differences"),
           subtitle = subtitle_text) +
      theme_pub()
  }

  file_name <- paste0(
    ifelse(hit$test_family == "Sensitivity_within_site",
           "sensitive_vs_not_sensitive.",
           "body_site_differences."),
    safe_feature_filename(hit$genus_id),
    ".png"
  )
  file_path <- file.path(
    if (hit$test_family == "Sensitivity_within_site") {
      sensitivity_comparison_dir
    } else {
      bodysite_comparison_dir
    },
    file_name
  )
  ggsave(file_path, p, width = ifelse(hit$test_family == "Sensitivity_within_site", 9, 11),
         height = ifelse(hit$test_family == "Sensitivity_within_site", 4.8, 5.5),
         dpi = 300, bg = "white")
  record_plot("comparison", hit$genus_id, hit$genus, hit$scope,
              paste(hit$test_family, "summary plot"),
              hit$wilcox_p, hit$wilcox_FDR, file_path)
}

## 4b. Individual SS-10 correlation plots.
## Plot only correlations supported by BOTH Pearson and Spearman at the chosen
## p-value/FDR threshold and with the same sign.  This is stricter than plotting
## either method alone and reduces visually misleading hits driven by scale or
## rank artifacts.
correlation_hits <- genus_correlations %>%
  select(genus_id, genus, BodySite, method, estimate, p_value, p_FDR) %>%
  pivot_wider(names_from = method,
              values_from = c(estimate, p_value, p_FDR),
              names_sep = "_") %>%
  filter(passes_threshold(p_value_pearson, p_FDR_pearson),
         passes_threshold(p_value_spearman, p_FDR_spearman),
         sign(estimate_pearson) == sign(estimate_spearman)) %>%
  arrange(pmax(p_value_pearson, p_value_spearman), genus, BodySite)

correlation_plot_diagnostics <- list()
plotted_correlation_keys <- character()
for (ii in seq_len(nrow(correlation_hits))) {
  hit <- correlation_hits[ii, ]
  if (hit$genus_id %in% plotted_correlation_keys) next
  plotted_correlation_keys <- c(plotted_correlation_keys, hit$genus_id)

  genus_correlation_hits <- correlation_hits %>%
    filter(genus_id == hit$genus_id) %>%
    mutate(BodySite = factor(BodySite, levels = sites))

  plot_data <- genus_long %>%
    filter(genus_id == hit$genus_id) %>%
    filter(is.finite(SS10_total), is.finite(abundance))

  panel_y <- plot_data %>%
    group_by(BodySite) %>%
    summarise(y_base = max(abundance, na.rm = TRUE), .groups = "drop") %>%
    mutate(y_base = ifelse(is.finite(y_base) & y_base > 0, y_base, 0.05))
  panel_stats <- genus_correlation_hits %>%
    left_join(panel_y, by = "BodySite") %>%
    mutate(
      x = min(plot_data$SS10_total, na.rm = TRUE),
      y = y_base * 1.08,
      label = paste0(
        "Pearson r=", round(estimate_pearson, 3),
        ", p=", p_label(p_value_pearson),
        ", FDR=", p_label(p_FDR_pearson),
        "\nSpearman rho=", round(estimate_spearman, 3),
        ", p=", p_label(p_value_spearman),
        ", FDR=", p_label(p_FDR_spearman)
      )
    )

  subtitle_text <- paste0(
    "One plot per genus; panels are labeled only where Pearson and Spearman both meet ",
    threshold_label(), " and have the same direction"
  )
  p <- ggplot(plot_data, aes(SS10_total, abundance, colour = Sensitivity, fill = Sensitivity)) +
    geom_point(size = 2, alpha = 0.9) +
    geom_smooth(aes(group = 1), method = "lm", se = TRUE, alpha = 0.15,
                linewidth = 0.7, colour = "black", fill = "grey70") +
    geom_text(data = panel_stats, aes(x = x, y = y, label = label),
              inherit.aes = FALSE, hjust = 0, size = 3.0, lineheight = 0.9) +
    facet_wrap(~ BodySite, scales = "free_y", nrow = 1) +
    scale_colour_manual(values = pal_sens) +
    scale_fill_manual(values = pal_sens) +
    labs(x = "SS-10 total sensitivity score", y = "Relative abundance",
         title = paste0(hit$genus, " vs SS-10"),
         subtitle = subtitle_text,
         caption = "Black regression line is fitted to Sensitive and NotSensitive samples together.") +
    theme_pub()

  file_name <- paste0(
    "correlation.",
    "pearson_spearman_concordant.",
    safe_feature_filename(hit$genus_id),
    ".png"
  )
  file_path <- file.path(correlation_dir, file_name)
  ggsave(file_path, p, width = 10, height = 4.8, dpi = 300, bg = "white")
  record_plot("correlation", hit$genus_id, hit$genus, "all_body_sites",
              "Pearson and Spearman concordant SS10 correlation summary plot",
              min(pmax(genus_correlation_hits$p_value_pearson,
                       genus_correlation_hits$p_value_spearman), na.rm = TRUE),
              min(pmax(genus_correlation_hits$p_FDR_pearson,
                       genus_correlation_hits$p_FDR_spearman), na.rm = TRUE),
              file_path)

  for (jj in seq_len(nrow(genus_correlation_hits))) {
    dhit <- genus_correlation_hits[jj, ]
    diag_data <- plot_data %>% filter(BodySite == dhit$BodySite)
    loo_pearson <- loo_cor_diagnostics(diag_data$SS10_total, diag_data$abundance,
                                       "pearson", dhit$estimate_pearson)
    loo_spearman <- loo_cor_diagnostics(diag_data$SS10_total, diag_data$abundance,
                                        "spearman", dhit$estimate_spearman)
    correlation_plot_diagnostics[[length(correlation_plot_diagnostics) + 1]] <- tibble(
      genus_id = dhit$genus_id,
      genus = dhit$genus,
      BodySite = as.character(dhit$BodySite),
      n = nrow(diag_data),
      pearson_r = dhit$estimate_pearson,
      pearson_p = dhit$p_value_pearson,
      pearson_FDR = dhit$p_FDR_pearson,
      pearson_LOO_same_direction_all = loo_pearson$loo_same_direction_all,
      pearson_LOO_min_r = loo_pearson$loo_min_estimate,
      pearson_LOO_max_r = loo_pearson$loo_max_estimate,
      pearson_LOO_min_abs_r = loo_pearson$loo_min_abs_estimate,
      pearson_LOO_max_p = loo_pearson$loo_max_p,
      spearman_rho = dhit$estimate_spearman,
      spearman_p = dhit$p_value_spearman,
      spearman_FDR = dhit$p_FDR_spearman,
      spearman_LOO_same_direction_all = loo_spearman$loo_same_direction_all,
      spearman_LOO_min_rho = loo_spearman$loo_min_estimate,
      spearman_LOO_max_rho = loo_spearman$loo_max_estimate,
      spearman_LOO_min_abs_rho = loo_spearman$loo_min_abs_estimate,
      spearman_LOO_max_p = loo_spearman$loo_max_p,
      file = file_path
    )
  }
}

individual_plot_index <- if (length(individual_plot_records) == 0) {
  tibble(plot_type = character(), genus_id = character(), genus = character(),
         scope = character(), test = character(), p_value = numeric(),
         fdr_value = numeric(), file = character())
} else {
  bind_rows(individual_plot_records) %>%
    arrange(plot_type, p_value, genus, scope, test)
}

correlation_plot_diagnostics <- if (length(correlation_plot_diagnostics) == 0) {
  tibble(genus_id = character(), genus = character(), BodySite = character(),
         n = integer(), pearson_r = numeric(), pearson_p = numeric(),
         pearson_FDR = numeric(), pearson_LOO_same_direction_all = logical(),
         pearson_LOO_min_r = numeric(), pearson_LOO_max_r = numeric(),
         pearson_LOO_min_abs_r = numeric(), pearson_LOO_max_p = numeric(),
         spearman_rho = numeric(), spearman_p = numeric(), spearman_FDR = numeric(),
         spearman_LOO_same_direction_all = logical(),
         spearman_LOO_min_rho = numeric(), spearman_LOO_max_rho = numeric(),
         spearman_LOO_min_abs_rho = numeric(), spearman_LOO_max_p = numeric(),
         file = character())
} else {
  bind_rows(correlation_plot_diagnostics) %>%
    arrange(pmax(pearson_p, spearman_p), genus, BodySite)
}

## ---------------------------------------------------------------------------
## 5. Export TSV and Excel results
## ---------------------------------------------------------------------------
write.table(taxa_rank, out(basename(genus_file), ".taxa_rank.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(genus_summary, out(basename(genus_file), ".group_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(genus_omnibus, out(basename(genus_file), ".omnibus_stats.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(genus_contrasts, out(basename(genus_file), ".contrast_stats.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(genus_correlations, out(basename(genus_file), ".SS10_correlations.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(individual_plot_index, out(basename(genus_file), ".individual_plot_index.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(correlation_plot_diagnostics,
            out(basename(genus_file), ".correlation_plot_diagnostics.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

wb <- createWorkbook()
addWorksheet(wb, "taxa_rank");           writeData(wb, "taxa_rank", taxa_rank)
addWorksheet(wb, "group_summary");       writeData(wb, "group_summary", genus_summary)
addWorksheet(wb, "omnibus_stats");       writeData(wb, "omnibus_stats", genus_omnibus)
addWorksheet(wb, "contrast_stats");      writeData(wb, "contrast_stats", genus_contrasts)
addWorksheet(wb, "SS10_correlations");   writeData(wb, "SS10_correlations", genus_correlations)
addWorksheet(wb, "topN_long_abundance"); writeData(wb, "topN_long_abundance", plot_df)
addWorksheet(wb, "individual_plot_index"); writeData(wb, "individual_plot_index", individual_plot_index)
addWorksheet(wb, "correlation_plot_diagnostics"); writeData(wb, "correlation_plot_diagnostics", correlation_plot_diagnostics)
saveWorkbook(wb, out(basename(genus_file), ".genus_stats.xlsx"), overwrite = TRUE)

## ---------------------------------------------------------------------------
## 6. Console summary
## ---------------------------------------------------------------------------
cat("\n==================== DONE ====================\n")
cat("Genus abundance file: ", normalizePath(genus_file), "\n", sep = "")
cat("Metadata file:        ", normalizePath(meta_file), "\n", sep = "")
cat("Samples: ", n_distinct(genus_long$sample_id),
    " | Genera: ", n_distinct(genus_long$genus_id),
    " | Top genera plotted: ", length(top_taxa), "\n", sep = "")
cat("Individual plot threshold: ", threshold_label(),
    " | Individual plots: ", nrow(individual_plot_index), "\n", sep = "")
cat("Output folder:        ", normalizePath(outdir), "\n", sep = "")
cat("\nMost abundant genera:\n")
print(head(taxa_rank %>% select(genus_id, genus, mean_abundance, prevalence), 10),
      row.names = FALSE)
