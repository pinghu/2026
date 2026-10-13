###############################################################################
## GSS3318 SensiScope — feature-level abundance analysis
## -----------------------------------------------------------------------------
## New script that leaves get_paired_wilcox_site.R unchanged.
##
## Purpose
##   Analyze the feature abundance table in the same spirit as the alpha-diversity
##   workflow: grouped abundance summaries, Sensitivity contrasts, paired
##   BodySite contrasts, SS-10 correlations, figures, and Excel statistics.
##
## Study design encoded by sample names
##   Sample columns must look like: BodySite.Sensitivity.SID
##   Example: Back.NotSensitive.1
##
## Outputs
##   A timestamped folder:
##     <data_filename>_result_<YYYYmmdd_HHMMSS>/
##   containing:
##     * PNG figures
##     * TSV statistics tables (incl. <data_file>.skin_measure_correlations.tsv:
##       pairwise Pearson/Spearman r/rho + p for Dryness, TotalDNAYield, Redness,
##       Skin smoothness, Tewameter and Corneometer, per BodySite and All)
##     * <data_file>.master_summary.tsv: one row per feature holding n/mean/SE for
##       12 sample groups, Kruskal-Wallis p (Sensitivity, BodySite, SID,
##       BodySite x Sensitivity), Sensitive-vs-NotSensitive fold change and p for
##       All/Leg/Back/Hand, prevalence, and Pearson/Spearman correlations with
##       SS-10, Dryness, Smoothness, Redness, Tewameter, Corneometer and
##       TotalDNAYield, plus a "solid" flag when both methods agree at p <= 0.05.
##       The same table is written to the colour-highlighted Excel sheet
##       "master_summary".
##     * <data_file>.consistent_SvsNS_SS10cor_summary.tsv (and Excel sheet
##       "consistent_summary"): one row per feature with Sensitive vs
##       NotSensitive p<=0.1 in any of All/Leg/Hand/Back.  For each group:
##       S vs NS log2 fold change (* p<0.1, ** p<0.05, *** FDR<0.1,
##       **** FDR<0.05; red up / blue down) and the 7 master-style correlation
##       direction columns.  The feature ID is coloured when the fold change
##       agrees with the SS-10 correlation direction.
##       Numeric details: <data_file>.consistent_SvsNS_SS10cor_details.tsv
##     * Individual plots named <feature>.S_vs_NS.png,
##       <feature>.BodySiteDiff.png and <feature>.SS10PearsonSpearmanCor.png
##     * Excel workbook with master_summary and consistent_summary sheets
##
## Usage
##   Rscript get_paired_wilcox_site_GSS3318_analysis4.R [data_file] [metadata_file] [top_n] [outdir] [plot_threshold] [threshold_type]
##
##   threshold_type can be "p" or "fdr".  The default is raw p-value filtering
##   with threshold 0.10.  If the 4th argument is numeric and no outdir is
##   supplied, it is interpreted as plot_threshold and the timestamped output
##   folder is used.
##
##   Individual paired body-site difference figures are NOT produced by
##   default.  Add the flag --bodysite (anywhere on the command line) to
##   turn them on.
##
## Examples
##   Rscript get_paired_wilcox_site_GSS3318_analysis4.R
##   Rscript get_paired_wilcox_site_GSS3318_analysis4.R genus_cleandata GSS3318_readsCount_meta.txt 25
##   Rscript get_paired_wilcox_site_GSS3318_analysis4.R species60_cleandata GSS3318_readsCount_meta.txt 25
##   Rscript get_paired_wilcox_site_GSS3318_analysis4.R genus_cleandata GSS3318_readsCount_meta.txt 20 0.1
##   Rscript get_paired_wilcox_site_GSS3318_analysis4.R genus_cleandata GSS3318_readsCount_meta.txt 20 my_results 0.05
##   Rscript get_paired_wilcox_site_GSS3318_analysis4.R genus_cleandata GSS3318_readsCount_meta.txt 20 my_results_fdr 0.10 fdr
##   Rscript get_paired_wilcox_site_GSS3318_analysis4.R genus_cleandata GSS3318_readsCount_meta.txt 20 0.1 --bodysite
###############################################################################

rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(openxlsx)
  library(scales)
  library(cowplot)
})

set.seed(1234)

## ---------------------------------------------------------------------------
## 0. User options and shared helper functions
## ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
## Optional flag (any position): --bodysite turns on the individual paired
## body-site difference figures, which are skipped by default.
bodysite_flags <- c("--bodysite", "--plot-bodysite", "--bodysite=TRUE", "bodysite=TRUE")
plot_bodysite <- any(tolower(args) %in% tolower(bodysite_flags))
args <- args[!tolower(args) %in% tolower(bodysite_flags)]

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

data_file <- find_file(if (length(args) >= 1) args[1] else "genus_cleandata")
meta_file  <- find_file(if (length(args) >= 2) args[2] else "GSS3318_readsCount_meta.txt")
top_n      <- if (length(args) >= 3) suppressWarnings(as.integer(args[3])) else 20L
if (is.na(top_n) || top_n < 1) top_n <- 20L

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
default_outdir <- paste0(basename(data_file), "_result_", timestamp)
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
  ggsave(out(basename(data_file), ".", name, ".png"), plot,
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

clean_feature_name <- function(x) {
  x <- gsub("^[a-z]__", "", x)
  x <- gsub("_", " ", x)
  ifelse(is.na(x) | x == "", "Unclassified", x)
}

safe_filename <- function(x) {
  ## File names are built from feature names, groups and p-values.  Keep them
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
## 1. Load feature data and metadata
## ---------------------------------------------------------------------------
feature_wide <- read.table(data_file, header = TRUE, row.names = 1, sep = "\t",
                         quote = "", comment.char = "",
                         check.names = FALSE, stringsAsFactors = FALSE)

## Force abundance columns to numeric and treat missing/non-numeric abundance as
## zero.  The input feature table is expected to be a numeric abundance matrix.
feature_wide[] <- lapply(feature_wide, function(x) {
  z <- suppressWarnings(as.numeric(x))
  z[is.na(z)] <- 0
  z
})

sample_parts <- strsplit(colnames(feature_wide), "[.]")
bad_samples <- vapply(sample_parts, length, integer(1)) < 3
if (any(bad_samples)) {
  stop("Sample names must look like BodySite.Sensitivity.SID. Bad names: ",
       paste(colnames(feature_wide)[bad_samples], collapse = ", "), call. = FALSE)
}

sample_map <- tibble(
  sample_id = colnames(feature_wide),
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

positive_min <- min(as.matrix(feature_wide)[as.matrix(feature_wide) > 0], na.rm = TRUE)
if (!is.finite(positive_min)) positive_min <- 1
pseudocount <- positive_min / 2

feature_long <- feature_wide %>%
  rownames_to_column("feature_id") %>%
  mutate(feature = clean_feature_name(feature_id)) %>%
  pivot_longer(-c(feature_id, feature), names_to = "sample_id",
               values_to = "abundance") %>%
  left_join(sample_map, by = "sample_id") %>%
  left_join(meta2, by = "sample_id") %>%
  mutate(BodySite = factor(BodySite, levels = sites),
         Sensitivity = factor(Sensitivity, levels = sens_levels),
         SiteSensitivity = factor(SiteSensitivity, levels = site_sens_levels),
         SID = factor(SID),
         abundance_log10 = log10(abundance + pseudocount))

missing_meta <- feature_long %>% filter(is.na(SS10_total)) %>% distinct(sample_id)
if (nrow(missing_meta) > 0) {
  warning("Metadata was not found for sample(s): ",
          paste(missing_meta$sample_id, collapse = ", "))
}

feature_ranking <- feature_long %>%
  group_by(feature_id, feature) %>%
  summarise(mean_abundance = mean(abundance, na.rm = TRUE),
            median_abundance = median(abundance, na.rm = TRUE),
            max_abundance = max(abundance, na.rm = TRUE),
            prevalence = mean(abundance > 0, na.rm = TRUE),
            percent_zero = mean(abundance == 0 | is.na(abundance)),
            .groups = "drop") %>%
  arrange(desc(mean_abundance), desc(prevalence), feature)

top_features <- head(feature_ranking$feature_id, top_n)

## ---------------------------------------------------------------------------
## 1b. Pairwise Pearson / Spearman correlations among skin measurements
## ---------------------------------------------------------------------------
## Computed per BodySite and for all samples together.
skin_vars <- c(Dryness = "Dryness",
               TotalDNAYield = "TotalDNAYield_ng",
               Redness = "Redness",
               Skin_smoothness = "Skin_smoothness_visual",
               Tewameter = "Tewameter_Measurements",
               Corneometer = "Corneometer_Measurements")
missing_skin <- setdiff(skin_vars, colnames(meta))
if (length(missing_skin) > 0) {
  stop("Missing skin measurement column(s) in metadata: ",
       paste(missing_skin, collapse = ", "), call. = FALSE)
}

skin_df <- meta %>%
  filter(NewID %in% sample_map$sample_id) %>%
  mutate(across(all_of(unname(skin_vars)), ~ suppressWarnings(as.numeric(.x))))

skin_pairs <- t(combn(names(skin_vars), 2))
skin_scopes <- c(sites, "All")

skin_correlations <- lapply(skin_scopes, function(scope) {
  dd <- if (scope == "All") skin_df else skin_df %>% filter(BodySite == scope)
  lapply(seq_len(nrow(skin_pairs)), function(i) {
    v1 <- skin_pairs[i, 1]
    v2 <- skin_pairs[i, 2]
    x <- dd[[skin_vars[[v1]]]]
    y <- dd[[skin_vars[[v2]]]]
    pr <- cor_stats(x, y, "pearson")
    sp <- cor_stats(x, y, "spearman")
    tibble(BodySite = scope, variable_1 = v1, variable_2 = v2, n = pr$n,
           pearson_r = pr$estimate, pearson_p = pr$p_value,
           spearman_rho = sp$estimate, spearman_p = sp$p_value)
  }) %>% bind_rows()
}) %>% bind_rows() %>%
  group_by(BodySite) %>%
  mutate(pearson_FDR = p.adjust(pearson_p, method = "BH"),
         spearman_FDR = p.adjust(spearman_p, method = "BH")) %>%
  ungroup() %>%
  select(BodySite, variable_1, variable_2, n,
         pearson_r, pearson_p, pearson_FDR,
         spearman_rho, spearman_p, spearman_FDR) %>%
  mutate(BodySite = factor(BodySite, levels = skin_scopes)) %>%
  arrange(BodySite) %>%
  mutate(BodySite = as.character(BodySite))

## ---------------------------------------------------------------------------
## 2. Statistical analysis
## ---------------------------------------------------------------------------
feature_group_summary <- feature_long %>%
  group_by(feature_id, feature, BodySite, Sensitivity, SiteSensitivity) %>%
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

feature_stat_one <- function(id) {
  dd <- feature_long %>% filter(feature_id == id)
  feature_name <- unique(dd$feature)[1]

  omnibus <- tibble(
    feature_id = id,
    feature = feature_name,
    KW_BodySite = kw_p(dd, "abundance", "BodySite"),
    KW_BodySite_Sensitive = kw_p(dd %>% filter(Sensitivity == "Sensitive"),
                                 "abundance", "BodySite"),
    KW_BodySite_NotSensitive = kw_p(dd %>% filter(Sensitivity == "NotSensitive"),
                                     "abundance", "BodySite"),
    KW_Sensitivity = kw_p(dd, "abundance", "Sensitivity"),
    KW_SiteSensitivity = kw_p(dd, "abundance", "SiteSensitivity"),
    KW_SID = kw_p(dd, "abundance", "SID"),
    percent_zero = mean(dd$abundance == 0 | is.na(dd$abundance)),
    mean_abundance = mean(dd$abundance, na.rm = TRUE),
    max_abundance = max(dd$abundance, na.rm = TRUE)
  )

  sensitivity_tests <- bind_rows(lapply(c(sites, "All"), function(site) {
    xx <- if (site == "All") dd else dd %>% filter(BodySite == site)
    s <- xx$abundance[xx$Sensitivity == "Sensitive"]
    n <- xx$abundance[xx$Sensitivity == "NotSensitive"]
    tibble(feature_id = id, feature = feature_name, scope = site,
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
      tibble(feature_id = id, feature = feature_name, scope = scope,
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

stat_list <- lapply(feature_ranking$feature_id, feature_stat_one)

feature_omnibus <- bind_rows(lapply(stat_list, `[[`, "omnibus")) %>%
  mutate(KW_BodySite_FDR = p.adjust(KW_BodySite, method = "BH"),
         KW_BodySite_Sensitive_FDR = p.adjust(KW_BodySite_Sensitive, method = "BH"),
         KW_BodySite_NotSensitive_FDR = p.adjust(KW_BodySite_NotSensitive, method = "BH"),
         KW_Sensitivity_FDR = p.adjust(KW_Sensitivity, method = "BH"),
         KW_SiteSensitivity_FDR = p.adjust(KW_SiteSensitivity, method = "BH"),
         KW_SID_FDR = p.adjust(KW_SID, method = "BH"))

feature_contrasts <- bind_rows(lapply(stat_list, `[[`, "contrasts")) %>%
  group_by(test_family, scope, group_a, group_b) %>%
  mutate(wilcox_FDR = p.adjust(wilcox_p, method = "BH"),
         t_FDR = p.adjust(t_p, method = "BH")) %>%
  ungroup()

feature_correlations <- feature_long %>%
  group_by(feature_id, feature, BodySite) %>%
  group_modify(~ bind_rows(
    cor_stats(.x$SS10_total, .x$abundance, "spearman") %>% mutate(method = "spearman"),
    cor_stats(.x$SS10_total, .x$abundance, "pearson") %>% mutate(method = "pearson")
  )) %>%
  ungroup() %>%
  select(feature_id, feature, BodySite, method, n, estimate, conf_low, conf_high, p_value) %>%
  group_by(BodySite, method) %>%
  mutate(p_FDR = p.adjust(p_value, method = "BH")) %>%
  ungroup()

## ---------------------------------------------------------------------------
## 2b. Master table: one row per feature with every statistic side by side
## ---------------------------------------------------------------------------
## Columns are grouped as:
##   n / mean / SE for 12 sample groups
##   Kruskal-Wallis p (Sensitivity, BodySite, SID, BodySite x Sensitivity)
##   Sensitive vs NotSensitive fold change + p for All / Leg / Back / Hand
##   prevalence over all samples
##   Pearson and Spearman correlations with the clinical measurements
master_groups <- tibble(
  label = c("All", "All_Sensitive", "All_NotSensitive",
            "Leg", "Back", "Hand",
            "Leg_Sensitive", "Leg_NotSensitive",
            "Back_Sensitive", "Back_NotSensitive",
            "Hand_Sensitive", "Hand_NotSensitive"),
  site = c(NA, NA, NA, "Leg", "Back", "Hand",
           "Leg", "Leg", "Back", "Back", "Hand", "Hand"),
  sens = c(NA, "Sensitive", "NotSensitive", NA, NA, NA,
           "Sensitive", "NotSensitive", "Sensitive", "NotSensitive",
           "Sensitive", "NotSensitive")
)

group_block <- lapply(seq_len(nrow(master_groups)), function(i) {
  g <- master_groups[i, ]
  dd <- feature_long
  if (!is.na(g$site)) dd <- dd %>% filter(BodySite == g$site)
  if (!is.na(g$sens)) dd <- dd %>% filter(Sensitivity == g$sens)
  dd %>%
    group_by(feature_id, feature) %>%
    summarise(n = sum(is.finite(abundance)),
              mean = mean(abundance, na.rm = TRUE),
              SE = sd(abundance, na.rm = TRUE) / sqrt(sum(is.finite(abundance))),
              .groups = "drop") %>%
    rename_with(~ paste0(g$label, "_", .x), all_of(c("n", "mean", "SE")))
})
group_block <- Reduce(function(a, b) left_join(a, b, by = c("feature_id", "feature")),
                      group_block)

kw_block <- feature_omnibus %>%
  select(feature_id, feature,
         KW_Sensitivity, KW_Sensitivity_FDR,
         KW_BodySite, KW_BodySite_FDR,
         KW_SID, KW_SID_FDR,
         KW_SiteSensitivity, KW_SiteSensitivity_FDR)

## Sensitive vs NotSensitive fold change and p-value, one column pair per scope.
sens_block <- feature_contrasts %>%
  filter(test_family == "Sensitivity_within_site") %>%
  mutate(scope = factor(scope, levels = c("All", sites))) %>%
  arrange(scope) %>%
  select(feature_id, feature, scope, fold_change, wilcox_p, wilcox_FDR) %>%
  pivot_wider(names_from = scope,
              values_from = c(fold_change, wilcox_p, wilcox_FDR),
              names_glue = "SvsNS_{scope}_{.value}")
sens_block <- sens_block %>%
  select(feature_id, feature,
         unlist(lapply(c("All", sites), function(s)
           paste0("SvsNS_", s, "_", c("fold_change", "wilcox_p", "wilcox_FDR")))))

prev_block <- feature_ranking %>%
  transmute(feature_id, feature, prevalence_all_pct = prevalence * 100)

## Correlations between feature abundance and the clinical measurements, using all
## samples.  A correlation is called "solid" when Pearson and Spearman are both
## significant at p <= 0.05 and point in the same direction.
cor_vars <- c(SS10 = "SS10_total",
              Dryness = "Dryness",
              Smoothness = "Skin_smoothness_visual",
              Redness = "Redness",
              Tewameter = "Tewameter_Measurements",
              Corneometer = "Corneometer_Measurements",
              TotalDNAYield = "TotalDNAYield_ng")
missing_cor_vars <- setdiff(cor_vars, colnames(feature_long))
if (length(missing_cor_vars) > 0) {
  stop("Missing metadata column(s) for correlations: ",
       paste(missing_cor_vars, collapse = ", "), call. = FALSE)
}

## Short display labels used on the summary heatmap panel.
cor_var_labels <- c(SS10 = "SS-10", Dryness = "Dryness", Smoothness = "Smoothness",
                    Redness = "Redness", Tewameter = "TEWL", Corneometer = "Hydration",
                    TotalDNAYield = "DNA yield")

clinical_cor_long <- lapply(names(cor_vars), function(vn) {
  col <- cor_vars[[vn]]
  feature_long %>%
    group_by(feature_id, feature) %>%
    group_modify(~ {
      pr <- cor_stats(.x[[col]], .x$abundance, "pearson")
      sp <- cor_stats(.x[[col]], .x$abundance, "spearman")
      tibble(n = pr$n,
             pearson_r = pr$estimate, pearson_p = pr$p_value,
             spearman_rho = sp$estimate, spearman_p = sp$p_value)
    }) %>%
    ungroup() %>%
    mutate(variable = vn,
           pearson_FDR = p.adjust(pearson_p, method = "BH"),
           spearman_FDR = p.adjust(spearman_p, method = "BH"),
           solid = !is.na(pearson_p) & !is.na(spearman_p) &
                   pearson_p <= 0.05 & spearman_p <= 0.05 &
                   sign(pearson_r) == sign(spearman_rho),
           direction = ifelse(!solid, "", ifelse(pearson_r > 0, "positive", "negative")))
}) %>%
  bind_rows() %>%
  select(feature_id, feature, variable, n, pearson_r, pearson_p, pearson_FDR,
         spearman_rho, spearman_p, spearman_FDR, solid, direction)

cor_block <- clinical_cor_long %>%
  pivot_wider(id_cols = c(feature_id, feature), names_from = variable,
              values_from = c(n, pearson_r, pearson_p, pearson_FDR,
                              spearman_rho, spearman_p, spearman_FDR,
                              solid, direction),
              names_glue = "cor_{variable}_{.value}") %>%
  select(feature_id, feature,
         unlist(lapply(names(cor_vars), function(v)
           paste0("cor_", v, "_", c("n", "pearson_r", "pearson_p", "pearson_FDR",
                                    "spearman_rho", "spearman_p", "spearman_FDR",
                                    "solid", "direction")))))

feature_master <- feature_ranking %>%
  select(feature_id, feature, mean_abundance, max_abundance) %>%
  left_join(prev_block, by = c("feature_id", "feature")) %>%
  left_join(group_block, by = c("feature_id", "feature")) %>%
  left_join(kw_block, by = c("feature_id", "feature")) %>%
  left_join(sens_block, by = c("feature_id", "feature")) %>%
  left_join(cor_block, by = c("feature_id", "feature")) %>%
  arrange(desc(mean_abundance))

## ---------------------------------------------------------------------------
## 2c. Consistent summary: Sensitive vs NotSensitive + clinical correlations
## ---------------------------------------------------------------------------
## For each group (All samples, Leg, Hand, Back):
##   * <group>_SvsNS_log2FC: Sensitive vs NotSensitive log2 fold change with
##     * p<0.1, ** p<0.05, *** FDR<0.1, **** FDR<0.05; blank when p > 0.1.
##   * <group>_cor_<variable>_direction: the same "direction" rule used in the
##     master sheet (Pearson and Spearman both p <= solid_p with the same sign),
##     computed within the group for the 7 clinical variables.
## A feature is kept when its S vs NS result is shown in at least one group.
## <group>s where the shown fold change agrees with the SS-10 direction
## (S>NS & positive, or S<NS & negative) are listed in SS10_agree_groups and the
## feature ID cell is coloured like the fold change in the Excel sheet.
consistent_p <- 0.10
solid_p <- 0.05
consistent_scopes <- c("All", "Leg", "Hand", "Back")

scope_data <- bind_rows(
  feature_long %>% mutate(scope = "All"),
  feature_long %>% mutate(scope = as.character(BodySite))
)

scope_cor_direction <- lapply(names(cor_vars), function(vn) {
  col <- cor_vars[[vn]]
  scope_data %>%
    group_by(feature_id, feature, scope) %>%
    group_modify(~ {
      pr <- cor_stats(.x[[col]], .x$abundance, "pearson")
      sp <- cor_stats(.x[[col]], .x$abundance, "spearman")
      tibble(pearson_r = pr$estimate, pearson_p = pr$p_value,
             spearman_rho = sp$estimate, spearman_p = sp$p_value)
    }) %>%
    ungroup() %>%
    mutate(variable = vn,
           solid = !is.na(pearson_p) & !is.na(spearman_p) &
                   pearson_p <= solid_p & spearman_p <= solid_p &
                   sign(pearson_r) == sign(spearman_rho),
           direction = ifelse(!solid, "", ifelse(pearson_r > 0, "positive", "negative")))
}) %>%
  bind_rows()

svsns_scope <- feature_contrasts %>%
  filter(test_family == "Sensitivity_within_site", scope %in% consistent_scopes) %>%
  transmute(feature_id, feature, scope, log2_fc,
            SvsNS_p = wilcox_p, SvsNS_FDR = wilcox_FDR,
            shown = !is.na(SvsNS_p) & SvsNS_p <= consistent_p &
                    is.finite(log2_fc) & log2_fc != 0,
            SvsNS_log2FC = ifelse(shown, paste0(formatC(log2_fc, format = "f", digits = 2),
                                                sig_mark(SvsNS_p, SvsNS_FDR)), "")) %>%
  left_join(scope_cor_direction %>% filter(variable == "SS10") %>%
              select(feature_id, scope, SS10_direction = direction),
            by = c("feature_id", "scope")) %>%
  mutate(SS10_agree = shown & !is.na(SS10_direction) &
           ((log2_fc > 0 & SS10_direction == "positive") |
            (log2_fc < 0 & SS10_direction == "negative")))

consistent_details <- svsns_scope %>%
  filter(shown) %>%
  left_join(scope_cor_direction %>%
              select(feature_id, scope, variable, direction) %>%
              pivot_wider(names_from = variable, values_from = direction,
                          names_glue = "cor_{variable}_direction"),
            by = c("feature_id", "scope")) %>%
  mutate(scope = factor(scope, levels = consistent_scopes)) %>%
  arrange(feature_id, scope) %>%
  select(-shown, -SS10_direction)

consistent_block_cols <- c("SvsNS_log2FC", paste0("cor_", names(cor_vars), "_direction"))
consistent_cols <- unlist(lapply(consistent_scopes, function(s)
  paste0(s, "_", consistent_block_cols)))

kept_ids <- unique(svsns_scope$feature_id[svsns_scope$shown])

consistent_info <- svsns_scope %>%
  filter(feature_id %in% kept_ids) %>%
  mutate(scope = factor(scope, levels = consistent_scopes)) %>%
  arrange(scope) %>%
  group_by(feature_id, feature) %>%
  summarise(SS10_agree_groups = paste(as.character(scope[SS10_agree]), collapse = ";"),
            SS10_agree_direction = case_when(
              !any(SS10_agree) ~ "",
              all(log2_fc[SS10_agree] > 0) ~ "up",
              all(log2_fc[SS10_agree] < 0) ~ "down",
              TRUE ~ "mixed"),
            n_SvsNS_groups = sum(shown),
            best_p = min(SvsNS_p[shown]),
            .groups = "drop")

consistent_wide <- bind_rows(
  svsns_scope %>%
    filter(feature_id %in% kept_ids) %>%
    transmute(feature_id, scope, column = "SvsNS_log2FC", value = SvsNS_log2FC),
  scope_cor_direction %>%
    filter(feature_id %in% kept_ids) %>%
    transmute(feature_id, scope, column = paste0("cor_", variable, "_direction"),
              value = direction)
) %>%
  mutate(column = paste0(scope, "_", column)) %>%
  select(-scope) %>%
  pivot_wider(names_from = column, values_from = value, values_fill = "")
for (cl in setdiff(consistent_cols, colnames(consistent_wide))) {
  consistent_wide[[cl]] <- rep("", nrow(consistent_wide))
}

consistent_summary <- consistent_info %>%
  left_join(consistent_wide, by = "feature_id") %>%
  arrange(SS10_agree_direction == "", desc(nchar(SS10_agree_groups)),
          desc(n_SvsNS_groups), best_p) %>%
  select(feature_id, feature, SS10_agree_groups, SS10_agree_direction,
         n_SvsNS_groups, all_of(consistent_cols))

## ---------------------------------------------------------------------------
## 3. Figures
## ---------------------------------------------------------------------------
top_label_levels <- feature_ranking$feature[match(top_features, feature_ranking$feature_id)]

plot_df <- feature_long %>%
  filter(feature_id %in% top_features) %>%
  mutate(feature = factor(feature, levels = top_label_levels))

top_box <- ggplot(plot_df, aes(BodySite, abundance, fill = Sensitivity)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.78, width = 0.65,
               position = position_dodge(0.75)) +
  geom_point(aes(colour = Sensitivity),
             position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.75),
             size = 1.3, alpha = 0.8, show.legend = FALSE) +
  facet_wrap(~ feature, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = pal_sens) +
  scale_colour_manual(values = pal_sens) +
  labs(x = NULL, y = "Relative abundance",
       title = paste0("Top ", length(top_features), " features by body site and sensitivity"),
       subtitle = "Boxplots show feature abundance; points are individual samples") +
  theme_pub()
save_fig(top_box, "top_features.boxplot_bodySite_sensitivity", 12, 8)

heat_df <- feature_group_summary %>%
  filter(feature_id %in% top_features) %>%
  mutate(feature = factor(feature, levels = top_label_levels))

mean_heat_main <- ggplot(heat_df, aes(SiteSensitivity, feature, fill = log10(mean + 1e-6))) +
  geom_tile(colour = "white", linewidth = 0.25) +
  scale_fill_gradientn(name = "log10(mean + 1e-6)",
                       colours = c("#F7FBFF", "#6BAED6", "#08306B")) +
  labs(x = NULL, y = NULL, title = "Mean abundance of top features",
       subtitle = "Rows are ordered by overall mean abundance") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## Companion columns showing the Kruskal-Wallis BodySite p-value computed on
## (1) all subjects, (2) sensitive subjects only, and (3) not-sensitive
## subjects only.  These are drawn as a separate aligned panel so the
## continuous fill scale of the abundance heatmap is not disturbed.
kw_bodysite_cols <- c("All" = "KW_BodySite",
                      "Sensitive" = "KW_BodySite_Sensitive",
                      "Not-sensitive" = "KW_BodySite_NotSensitive")

kw_pval_long <- function(col_map, suffix) {
  cols <- paste0(col_map, suffix)
  feature_omnibus %>%
    filter(feature_id %in% top_features) %>%
    mutate(feature = factor(feature, levels = top_label_levels)) %>%
    select(feature, all_of(cols)) %>%
    setNames(c("feature", names(col_map))) %>%
    pivot_longer(-feature, names_to = "test", values_to = "value")
}

kw_pval_df <- kw_pval_long(kw_bodysite_cols, "") %>%
  rename(p_value = value) %>%
  left_join(kw_pval_long(kw_bodysite_cols, "_FDR") %>% rename(fdr_value = value),
            by = c("feature", "test")) %>%
  mutate(test = factor(test, levels = names(kw_bodysite_cols)),
         p_text = paste0(p_label(p_value), sig_mark(p_value, fdr_value)))

mean_heat_pvals <- ggplot(kw_pval_df, aes(test, feature, fill = -log10(p_value))) +
  geom_tile(colour = "white", linewidth = 0.25) +
  geom_text(aes(label = p_text), size = 2.4) +
  scale_fill_gradientn(name = "-log10(p)\nKW BodySite",
                       colours = c("#FFF5F0", "#FC9272", "#99000D"),
                       na.value = "grey90") +
  labs(x = NULL, y = NULL, title = "Kruskal-Wallis BodySite p-value",
       subtitle = "* p<0.10  ** p<0.05  *** FDR<0.10  **** FDR<0.05") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

mean_heat <- plot_grid(mean_heat_main, mean_heat_pvals, ncol = 2,
                       align = "h", axis = "tb", rel_widths = c(1, 0.85))
save_fig(mean_heat, "top_features.mean_abundance_heatmap", 13, 7)

sens_effect <- feature_contrasts %>%
  filter(test_family == "Sensitivity_within_site", feature_id %in% top_features) %>%
  mutate(feature = factor(feature, levels = top_label_levels),
         scope = factor(scope, levels = c(sites, "All")),
         significance = sig_mark(wilcox_p, wilcox_FDR))

effect_heat <- ggplot(sens_effect, aes(scope, feature, fill = log2_fc)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  geom_text(aes(label = significance), size = 4) +
  scale_fill_gradient2(name = "log2 FC",
                       low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, oob = squish) +
  labs(x = NULL, y = NULL, title = "A. Sensitive vs\nNot-sensitive") +
  theme_pub()

cor_plot_df <- feature_correlations %>%
  filter(method == "spearman", feature_id %in% top_features) %>%
  mutate(feature = factor(feature, levels = top_label_levels),
         significance = sig_mark(p_value, p_FDR))

cor_heat <- ggplot(cor_plot_df, aes(BodySite, feature, fill = estimate)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  geom_text(aes(label = significance), size = 4) +
  scale_fill_gradient2(name = "Spearman rho",
                       low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-1, 1), oob = squish) +
  labs(x = NULL, y = NULL, title = "B. SS-10 by\nbody site") +
  theme_pub()

## Panel C summarizes the Spearman correlation between feature abundance and each
## clinical measurement across all samples.  A dot marks correlations where
## Pearson and Spearman agree (both p <= 0.05, same direction).
clinical_heat_df <- clinical_cor_long %>%
  filter(feature_id %in% top_features) %>%
  mutate(feature = factor(feature, levels = top_label_levels),
         variable = factor(cor_var_labels[variable],
                           levels = unname(cor_var_labels)),
         significance = sig_mark(spearman_p, spearman_FDR))

clinical_heat <- ggplot(clinical_heat_df, aes(variable, feature, fill = spearman_rho)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  geom_tile(data = clinical_heat_df %>% filter(solid),
            fill = NA, colour = "black", linewidth = 0.6) +
  geom_text(aes(label = significance), size = 4) +
  scale_fill_gradient2(name = "Spearman rho",
                       low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-1, 1), oob = squish) +
  labs(x = NULL, y = NULL, title = "C. Clinical measures,\nall samples") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## All three panels share the same feature rows, so only panel A keeps the y-axis
## labels and the panels can be read across as one figure.
drop_y <- theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

## Every tile column gets the same absolute width (wide enough for "****"), so
## panel A (log2 FC) widens and panels B/C (Spearman rho) shrink to match.
heat_col_in <- 0.55
fix_heat_cols <- function(p, n_cols) {
  p + theme(panel.widths = unit(n_cols * heat_col_in, "in")) +
    guides(fill = guide_colourbar(title.position = "top",
                                  barwidth = unit(min(1.4, n_cols * heat_col_in), "in")))
}
heat_panels <- list(
  fix_heat_cols(effect_heat, nlevels(droplevels(sens_effect$scope))),
  fix_heat_cols(cor_heat + drop_y, length(unique(cor_plot_df$BodySite))),
  fix_heat_cols(clinical_heat + drop_y, nlevels(droplevels(clinical_heat_df$variable)))
)
heat_widths_in <- vapply(heat_panels, function(p) {
  sum(grid::convertWidth(ggplotGrob(p)$widths, "in", valueOnly = TRUE))
}, numeric(1))
sens_cor_heat_row <- plot_grid(plotlist = heat_panels,
                               ncol = 3, align = "h", axis = "tb",
                               rel_widths = heat_widths_in)
sens_cor_heat <- plot_grid(
  sens_cor_heat_row,
  ggdraw() + draw_label(
    paste("Colour: red = higher / positive, blue = lower / negative.",
          "* p<0.10   ** p<0.05   *** FDR<0.10   **** FDR<0.05.\n",
          "Black outline in C: Pearson and Spearman agree at p<=0.05."),
    size = 9, hjust = 0.5, x = 0.5),
  ncol = 1, rel_heights = c(1, 0.07))
save_fig(sens_cor_heat, "top_features.sensitivity_and_correlation_heatmap",
         max(9, sum(heat_widths_in) + 0.3), 7.5)


scatter_features <- head(top_features, min(6, length(top_features)))
scatter_levels <- feature_ranking$feature[match(scatter_features, feature_ranking$feature_id)]

scatter_df <- feature_long %>%
  filter(feature_id %in% scatter_features) %>%
  mutate(feature = factor(feature, levels = scatter_levels))

ss10_scatter <- ggplot(scatter_df, aes(SS10_total, abundance, colour = BodySite, fill = BodySite)) +
  geom_point(size = 1.7, alpha = 0.85) +
  geom_smooth(aes(group = 1), method = "lm", se = TRUE, alpha = 0.15,
              linewidth = 0.65, colour = "black", fill = "grey70") +
  facet_grid(feature ~ BodySite, scales = "free_y") +
  scale_colour_manual(values = pal_site) +
  scale_fill_manual(values = pal_site) +
  labs(x = "SS-10 total sensitivity score", y = "Relative abundance",
       title = "Top feature abundance vs SS-10 score",
       subtitle = "Linear trend shown for visualization; exact Pearson/Spearman tests are in Excel") +
  theme_pub() +
  theme(legend.position = "none")
save_fig(ss10_scatter, "top6_features.SS10_scatter", 10, 8)

pearson_scatter_df <- feature_long %>%
  filter(feature_id %in% top_features) %>%
  mutate(feature = factor(feature, levels = top_label_levels))

pearson_summary_scatter <- ggplot(
  pearson_scatter_df,
  aes(SS10_total, abundance, colour = BodySite, fill = BodySite)
) +
  geom_point(size = 1.4, alpha = 0.78) +
  geom_smooth(aes(group = BodySite), method = "lm", se = TRUE,
              alpha = 0.12, linewidth = 0.55) +
  facet_wrap(~ feature, scales = "free_y", ncol = 4) +
  scale_colour_manual(values = pal_site) +
  scale_fill_manual(values = pal_site) +
  labs(x = "SS-10 total sensitivity score", y = "Relative abundance",
       title = "Summary Pearson view: feature abundance vs SS-10 score",
       subtitle = paste0("Top ", length(top_features),
                         " features; site-specific linear trends shown for visualization")) +
  theme_pub()
save_fig(pearson_summary_scatter, "top_features.SS10_pearson_summary_scatter", 12, 9)

## ---------------------------------------------------------------------------
## 4. Individual plots for any raw p-value below the adjustable threshold
## ---------------------------------------------------------------------------
## These plots are meant for follow-up review of every result that "shows up"
## at the chosen exploratory threshold.  The default is raw p-value p < 0.10;
## pass threshold_type="fdr" to select by FDR instead.
##   Rscript get_paired_wilcox_site_GSS3318_analysis4.R genus_cleandata meta.txt 20 0.05
##   Rscript get_paired_wilcox_site_GSS3318_analysis4.R genus_cleandata meta.txt 20 out_fdr 0.10 fdr
individual_root <- file.path(outdir, "individual_plots")
comparison_root <- file.path(individual_root, "comparisons")
sensitivity_comparison_dir <- file.path(comparison_root, "sensitive_vs_not_sensitive")
bodysite_comparison_dir <- file.path(comparison_root, "body_site_differences")
correlation_dir <- file.path(individual_root, "correlations")
dir.create(sensitivity_comparison_dir, recursive = TRUE, showWarnings = FALSE)
if (plot_bodysite) dir.create(bodysite_comparison_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(correlation_dir, recursive = TRUE, showWarnings = FALSE)

individual_plot_records <- list()
record_plot <- function(type, feature_id, feature, scope, test, p_value, fdr_value, file) {
  individual_plot_records[[length(individual_plot_records) + 1]] <<- tibble(
    plot_type = type,
    feature_id = feature_id,
    feature = feature,
    scope = scope,
    test = test,
    p_value = p_value,
    fdr_value = fdr_value,
    file = file
  )
}

## 4a. Individual comparison plots.
comparison_hits <- feature_contrasts %>%
  filter(passes_threshold(wilcox_p, wilcox_FDR)) %>%
  filter(plot_bodysite | test_family == "Sensitivity_within_site") %>%
  arrange(wilcox_p, feature, test_family, scope, group_a, group_b)

plotted_comparison_keys <- character()
for (ii in seq_len(nrow(comparison_hits))) {
  hit <- comparison_hits[ii, ]
  comparison_key <- paste(hit$test_family, hit$feature_id, sep = "::")
  if (comparison_key %in% plotted_comparison_keys) next
  plotted_comparison_keys <- c(plotted_comparison_keys, comparison_key)

  dd <- feature_long %>% filter(feature_id == hit$feature_id)

  if (hit$test_family == "Sensitivity_within_site") {
    ## Unpaired Sensitive-vs-NotSensitive comparison shown as a single-panel
    ## bar plot (mean +/- SE) with Back, Hand, Leg and All (pooled) on the x
    ## axis; colour separates Sensitive from NotSensitive.  Significance marks
    ## from the Wilcoxon test are drawn above each pair of bars.
    site_levels <- c(sites, "All")
    bar_data <- bind_rows(
      dd %>% mutate(Site = as.character(BodySite)),
      dd %>% mutate(Site = "All")
    ) %>%
      mutate(Site = factor(Site, levels = site_levels),
             Sensitivity = factor(Sensitivity, levels = sens_levels)) %>%
      group_by(Site, Sensitivity) %>%
      summarise(n = sum(is.finite(abundance)),
                mean = mean(abundance, na.rm = TRUE),
                se = ifelse(n > 1, sd(abundance, na.rm = TRUE) / sqrt(n), 0),
                .groups = "drop")

    y_top <- max(bar_data$mean + bar_data$se, na.rm = TRUE)
    if (!is.finite(y_top) || y_top <= 0) y_top <- 0.05

    sig_data <- feature_contrasts %>%
      filter(feature_id == hit$feature_id,
             test_family == "Sensitivity_within_site",
             scope %in% site_levels) %>%
      transmute(Site = factor(scope, levels = site_levels),
                label = sig_mark(wilcox_p, wilcox_FDR)) %>%
      left_join(bar_data %>% group_by(Site) %>%
                  summarise(y_pair = max(mean + se, na.rm = TRUE), .groups = "drop"),
                by = "Site") %>%
      mutate(x = as.numeric(Site),
             y_bracket = y_pair + 0.05 * y_top,
             y_label = y_bracket + 0.015 * y_top) %>%
      filter(label != "")

    dodge_w <- 0.8
    p <- ggplot(bar_data, aes(Site, mean, fill = Sensitivity)) +
      geom_col(position = position_dodge(width = dodge_w), width = 0.72,
               colour = "black", linewidth = 0.3) +
      geom_errorbar(aes(ymin = pmax(mean - se, 0), ymax = mean + se),
                    position = position_dodge(width = dodge_w),
                    width = 0.25, linewidth = 0.4) +
      geom_segment(data = sig_data,
                   aes(x = x - dodge_w / 4, xend = x + dodge_w / 4,
                       y = y_bracket, yend = y_bracket),
                   inherit.aes = FALSE, linewidth = 0.4) +
      geom_text(data = sig_data, aes(x = x, y = y_label, label = label),
                inherit.aes = FALSE, size = 4, vjust = 0) +
      scale_fill_manual(values = pal_sens, name = NULL) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
      labs(x = NULL, y = hit$feature,
           title = hit$feature,
           caption = "* p<0.10  ** p<0.05  *** FDR<0.10  **** FDR<0.05") +
      theme_pub() +
      theme(panel.grid.major.x = element_blank(),
            plot.title = element_text(face = "bold.italic", hjust = 0.5),
            axis.title.y = element_text(face = "bold.italic"),
            plot.caption = element_text(size = 7, hjust = 0.5))
  } else {
    ## Paired body-site comparison.  Once any paired site comparison passes the
    ## threshold for this feature, show all three body sites so the full subject
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

    threshold_stats <- feature_contrasts %>%
      filter(feature_id == hit$feature_id,
             test_family == "Paired_BodySite",
             passes_threshold(wilcox_p, wilcox_FDR)) %>%
      mutate(Panel = factor(scope, levels = c("All", "Sensitive", "NotSensitive")),
             comparison_label = paste0(group_a, " vs ", group_b,
                                       ": p=", p_label(wilcox_p),
                                       sig_mark(wilcox_p, wilcox_FDR),
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

    kw_hit <- feature_omnibus %>% filter(feature_id == hit$feature_id) %>% slice(1)
    subtitle_text <- paste0(
      "Kruskal-Wallis BodySite p=", p_label(kw_hit$KW_BodySite),
      sig_mark(kw_hit$KW_BodySite, kw_hit$KW_BodySite_FDR),
      ", FDR=", p_label(kw_hit$KW_BodySite_FDR),
      " | SiteSensitivity p=", p_label(kw_hit$KW_SiteSensitivity),
      sig_mark(kw_hit$KW_SiteSensitivity, kw_hit$KW_SiteSensitivity_FDR),
      ", FDR=", p_label(kw_hit$KW_SiteSensitivity_FDR),
      "\nPanel labels show paired Wilcoxon site comparisons with ",
      threshold_label(), " (* p<0.10 ** p<0.05 *** FDR<0.10 **** FDR<0.05)"
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
           title = paste0(hit$feature, ": paired body-site differences"),
           subtitle = subtitle_text) +
      theme_pub()
  }

  file_name <- paste0(
    safe_feature_filename(hit$feature_id),
    ifelse(hit$test_family == "Sensitivity_within_site", ".S_vs_NS", ".BodySiteDiff"),
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
  ggsave(file_path, p, width = ifelse(hit$test_family == "Sensitivity_within_site", 3.6, 11),
         height = ifelse(hit$test_family == "Sensitivity_within_site", 4, 5.5),
         dpi = 300, bg = "white")
  record_plot("comparison", hit$feature_id, hit$feature, hit$scope,
              paste(hit$test_family, "summary plot"),
              hit$wilcox_p, hit$wilcox_FDR, file_path)
}

## 4b. Individual SS-10 correlation plots.
## Plot only correlations supported by BOTH Pearson and Spearman at the chosen
## p-value/FDR threshold and with the same sign.  This is stricter than plotting
## either method alone and reduces visually misleading hits driven by scale or
## rank artifacts.
correlation_hits <- feature_correlations %>%
  select(feature_id, feature, BodySite, method, estimate, p_value, p_FDR) %>%
  pivot_wider(names_from = method,
              values_from = c(estimate, p_value, p_FDR),
              names_sep = "_") %>%
  filter(passes_threshold(p_value_pearson, p_FDR_pearson),
         passes_threshold(p_value_spearman, p_FDR_spearman),
         sign(estimate_pearson) == sign(estimate_spearman)) %>%
  arrange(pmax(p_value_pearson, p_value_spearman), feature, BodySite)

correlation_plot_diagnostics <- list()
plotted_correlation_keys <- character()
for (ii in seq_len(nrow(correlation_hits))) {
  hit <- correlation_hits[ii, ]
  if (hit$feature_id %in% plotted_correlation_keys) next
  plotted_correlation_keys <- c(plotted_correlation_keys, hit$feature_id)

  feature_correlation_hits <- correlation_hits %>%
    filter(feature_id == hit$feature_id) %>%
    mutate(BodySite = factor(BodySite, levels = sites))

  plot_data <- feature_long %>%
    filter(feature_id == hit$feature_id) %>%
    filter(is.finite(SS10_total), is.finite(abundance))

  panel_y <- plot_data %>%
    group_by(BodySite) %>%
    summarise(y_base = max(abundance, na.rm = TRUE), .groups = "drop") %>%
    mutate(y_base = ifelse(is.finite(y_base) & y_base > 0, y_base, 0.05))
  panel_stats <- feature_correlation_hits %>%
    left_join(panel_y, by = "BodySite") %>%
    mutate(
      x = min(plot_data$SS10_total, na.rm = TRUE),
      y = y_base * 1.08,
      label = paste0(
        "Pearson r=", round(estimate_pearson, 3),
        ", p=", p_label(p_value_pearson),
        sig_mark(p_value_pearson, p_FDR_pearson),
        ", FDR=", p_label(p_FDR_pearson),
        "\nSpearman rho=", round(estimate_spearman, 3),
        ", p=", p_label(p_value_spearman),
        sig_mark(p_value_spearman, p_FDR_spearman),
        ", FDR=", p_label(p_FDR_spearman)
      )
    )

  subtitle_text <- paste0(
    "One plot per feature; panels are labeled only where Pearson and Spearman both meet ",
    threshold_label(), " and have the same direction",
    " (* p<0.10 ** p<0.05 *** FDR<0.10 **** FDR<0.05)"
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
         title = paste0(hit$feature, " vs SS-10"),
         subtitle = subtitle_text,
         caption = "Black regression line is fitted to Sensitive and NotSensitive samples together.") +
    theme_pub()

  file_name <- paste0(
    safe_feature_filename(hit$feature_id),
    ".SS10PearsonSpearmanCor.png"
  )
  file_path <- file.path(correlation_dir, file_name)
  ggsave(file_path, p, width = 10, height = 4.8, dpi = 300, bg = "white")
  record_plot("correlation", hit$feature_id, hit$feature, "all_body_sites",
              "Pearson and Spearman concordant SS10 correlation summary plot",
              min(pmax(feature_correlation_hits$p_value_pearson,
                       feature_correlation_hits$p_value_spearman), na.rm = TRUE),
              min(pmax(feature_correlation_hits$p_FDR_pearson,
                       feature_correlation_hits$p_FDR_spearman), na.rm = TRUE),
              file_path)

  for (jj in seq_len(nrow(feature_correlation_hits))) {
    dhit <- feature_correlation_hits[jj, ]
    diag_data <- plot_data %>% filter(BodySite == dhit$BodySite)
    loo_pearson <- loo_cor_diagnostics(diag_data$SS10_total, diag_data$abundance,
                                       "pearson", dhit$estimate_pearson)
    loo_spearman <- loo_cor_diagnostics(diag_data$SS10_total, diag_data$abundance,
                                        "spearman", dhit$estimate_spearman)
    correlation_plot_diagnostics[[length(correlation_plot_diagnostics) + 1]] <- tibble(
      feature_id = dhit$feature_id,
      feature = dhit$feature,
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
  tibble(plot_type = character(), feature_id = character(), feature = character(),
         scope = character(), test = character(), p_value = numeric(),
         fdr_value = numeric(), file = character())
} else {
  bind_rows(individual_plot_records) %>%
    arrange(plot_type, p_value, feature, scope, test)
}

correlation_plot_diagnostics <- if (length(correlation_plot_diagnostics) == 0) {
  tibble(feature_id = character(), feature = character(), BodySite = character(),
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
    arrange(pmax(pearson_p, spearman_p), feature, BodySite)
}

## ---------------------------------------------------------------------------
## 5. Export TSV and Excel results
## ---------------------------------------------------------------------------
write.table(feature_ranking, out(basename(data_file), ".feature_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(feature_group_summary, out(basename(data_file), ".group_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(feature_omnibus, out(basename(data_file), ".omnibus_stats.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(feature_contrasts, out(basename(data_file), ".contrast_stats.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(feature_correlations, out(basename(data_file), ".SS10_correlations.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(individual_plot_index, out(basename(data_file), ".individual_plot_index.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(correlation_plot_diagnostics,
            out(basename(data_file), ".correlation_plot_diagnostics.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(skin_correlations,
            out(basename(data_file), ".skin_measure_correlations.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(feature_master, out(basename(data_file), ".master_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
write.table(consistent_details, out(basename(data_file), ".consistent_SvsNS_SS10cor_details.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
write.table(consistent_summary, out(basename(data_file), ".consistent_SvsNS_SS10cor_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "")

wb <- createWorkbook()

## The workbook contains the master summary and the consistent-feature summary;
## detailed tables remain TSVs.
## ---- master_summary sheet with colour highlighting --------------------------
## Colour key:
##   yellow  significant p <= 0.05 (Kruskal-Wallis, contrast and correlation p)
##   red     up-regulated fold change (Sensitive > NotSensitive) / positive rho
##   blue    down-regulated fold change (Sensitive < NotSensitive) / negative rho
##   green   correlation flagged as solid (Pearson and Spearman agree, p <= 0.05)
ms_sheet <- "master_summary"
addWorksheet(wb, ms_sheet)
writeData(wb, ms_sheet, feature_master, headerStyle =
            createStyle(textDecoration = "bold", fgFill = "#DCE6F1",
                        halign = "center", wrapText = TRUE, border = "TopBottomLeftRight"))
freezePane(wb, ms_sheet, firstActiveRow = 2, firstActiveCol = 3)
setColWidths(wb, ms_sheet, cols = 1:ncol(feature_master), widths = "auto")

ms_rows <- 2:(nrow(feature_master) + 1)
ms_names <- colnames(feature_master)
style_sig   <- createStyle(bgFill = "#FFE699", fontColour = "#7F6000", textDecoration = "bold")
style_up    <- createStyle(bgFill = "#F8CBAD", fontColour = "#843C0C")
style_down  <- createStyle(bgFill = "#BDD7EE", fontColour = "#1F4E79")
style_solid <- createStyle(bgFill = "#C6EFCE", fontColour = "#006100", textDecoration = "bold")

cf <- function(cols, rule, style) {
  for (cl in cols) {
    conditionalFormatting(wb, ms_sheet, cols = cl, rows = ms_rows,
                          rule = rule, style = style, type = "expression")
  }
}

p_cols  <- which((grepl("^KW_", ms_names) | grepl("_p$", ms_names)) &
                 !grepl("_FDR$", ms_names))
fdr_cols <- which(grepl("_FDR$", ms_names))
fc_cols <- which(grepl("_fold_change$", ms_names))
r_cols  <- which(grepl("_pearson_r$|_spearman_rho$", ms_names))
solid_cols <- which(grepl("_solid$", ms_names))

cf(c(p_cols, fdr_cols), "<=0.05", style_sig)
cf(fc_cols, ">0", style_up)
cf(fc_cols, "<0", style_down)
cf(r_cols, ">0", style_up)
cf(r_cols, "<0", style_down)

## "solid" flags are logical, so highlight the TRUE cells directly instead of
## relying on a conditional-formatting expression.
for (cl in solid_cols) {
  hits <- which(vapply(feature_master[[ms_names[cl]]], isTRUE, logical(1)))
  if (length(hits) > 0) {
    addStyle(wb, ms_sheet, style_solid, rows = hits + 1, cols = cl,
             gridExpand = TRUE, stack = TRUE)
  }
}

## Legend for the colour coding, written below the table.
legend_row <- nrow(feature_master) + 3
legend_df <- data.frame(Colour_key = c(
  "Yellow: p <= 0.05 (raw p or FDR)",
  "Red: fold change > 0 (Sensitive > NotSensitive) or positive correlation",
  "Blue: fold change < 0 (Sensitive < NotSensitive) or negative correlation",
  "Green: solid correlation (Pearson and Spearman both p <= 0.05, same direction)"))
writeData(wb, ms_sheet, legend_df, startRow = legend_row, colNames = TRUE)
addStyle(wb, ms_sheet, createStyle(textDecoration = "bold"),
         rows = legend_row, cols = 1)
addStyle(wb, ms_sheet, style_sig,   rows = legend_row + 1, cols = 1)
addStyle(wb, ms_sheet, style_up,    rows = legend_row + 2, cols = 1)
addStyle(wb, ms_sheet, style_down,  rows = legend_row + 3, cols = 1)
addStyle(wb, ms_sheet, style_solid, rows = legend_row + 4, cols = 1)

## ---- consistent_summary sheet ----------------------------------------------
cs_sheet <- "consistent_summary"
addWorksheet(wb, cs_sheet)
cs_note <- data.frame(Note = c(
  paste0("Features with Sensitive vs NotSensitive Wilcoxon p<=", consistent_p,
         " in at least one group (All samples, Leg, Hand, Back)."),
  paste0("SvsNS_log2FC = log2(mean Sensitive / mean NotSensitive); blank if p>", consistent_p,
         ". Red = higher in Sensitive, blue = lower. * p<0.1  ** p<0.05  *** FDR<0.1  **** FDR<0.05"),
  paste0("Correlation direction columns use the master-sheet rule within each group: ",
         "Pearson and Spearman both p<=", solid_p, " with the same sign."),
  "Feature ID coloured (red/blue) when the shown fold change agrees with the SS10 direction in that group; purple = agreement in opposite directions across groups."))
writeData(wb, cs_sheet, cs_note, colNames = FALSE)

cs_hdr <- nrow(cs_note) + 2
cs_fixed <- c("feature_id", "feature", "SS10_agree_groups", "SS10_agree_direction",
              "n_SvsNS_groups")
short_block <- c("SvsNS_log2FC", names(cor_vars))
cs_style_hdr <- createStyle(textDecoration = "bold", fgFill = "#DCE6F1", halign = "center",
                            valign = "center", border = "TopBottomLeftRight", wrapText = TRUE)
writeData(wb, cs_sheet, as.data.frame(t(c(cs_fixed,
                                         rep(short_block, length(consistent_scopes))))),
          startRow = cs_hdr + 1, colNames = FALSE)
for (i in seq_along(consistent_scopes)) {
  c1 <- length(cs_fixed) + (i - 1) * length(short_block) + 1
  c2 <- c1 + length(short_block) - 1
  writeData(wb, cs_sheet, ifelse(consistent_scopes[i] == "All", "All samples",
                                 consistent_scopes[i]),
            startRow = cs_hdr, startCol = c1)
  mergeCells(wb, cs_sheet, cols = c1:c2, rows = cs_hdr)
}
addStyle(wb, cs_sheet, cs_style_hdr, rows = cs_hdr:(cs_hdr + 1),
         cols = seq_len(ncol(consistent_summary)), gridExpand = TRUE)
writeData(wb, cs_sheet, consistent_summary, startRow = cs_hdr + 2, colNames = FALSE)
freezePane(wb, cs_sheet, firstActiveRow = cs_hdr + 2, firstActiveCol = 3)
setColWidths(wb, cs_sheet, cols = 1:2, widths = 30)
setColWidths(wb, cs_sheet, cols = 3:ncol(consistent_summary), widths = 12)

if (nrow(consistent_summary) > 0) {
  cs_rows <- seq_len(nrow(consistent_summary)) + cs_hdr + 1
  style_red  <- createStyle(fgFill = "#F8CBAD", fontColour = "#843C0C", halign = "center")
  style_blue <- createStyle(fgFill = "#BDD7EE", fontColour = "#1F4E79", halign = "center")
  style_mix  <- createStyle(fgFill = "#D9C3E9")
  for (s in consistent_scopes) {
    cl <- match(paste0(s, "_SvsNS_log2FC"), colnames(consistent_summary))
    v <- consistent_summary[[cl]]
    up <- which(v != "" & !startsWith(v, "-"))
    dn <- which(startsWith(v, "-"))
    if (length(up)) addStyle(wb, cs_sheet, style_red, rows = cs_rows[up], cols = cl, stack = TRUE)
    if (length(dn)) addStyle(wb, cs_sheet, style_blue, rows = cs_rows[dn], cols = cl, stack = TRUE)
  }
  ad <- consistent_summary$SS10_agree_direction
  for (k in list(list("up", style_red), list("down", style_blue), list("mixed", style_mix))) {
    rr <- which(ad == k[[1]])
    if (length(rr)) addStyle(wb, cs_sheet, k[[2]], rows = cs_rows[rr], cols = 1, stack = TRUE)
  }
}
saveWorkbook(wb, out(basename(data_file), ".stats.xlsx"), overwrite = TRUE)

## ---------------------------------------------------------------------------
## 6. Console summary
## ---------------------------------------------------------------------------
cat("\n==================== DONE ====================\n")
cat("Abundance data file:  ", normalizePath(data_file), "\n", sep = "")
cat("Metadata file:        ", normalizePath(meta_file), "\n", sep = "")
cat("Samples: ", n_distinct(feature_long$sample_id),
    " | Features: ", n_distinct(feature_long$feature_id),
    " | Top features plotted: ", length(top_features), "\n", sep = "")
cat("Individual plot threshold: ", threshold_label(),
    " | Individual plots: ", nrow(individual_plot_index), "\n", sep = "")
cat("Output folder:        ", normalizePath(outdir), "\n", sep = "")
cat("Master summary:       ", nrow(feature_master), " features x ",
    ncol(feature_master), " columns\n", sep = "")
cat("Consistent summary:   ", nrow(consistent_summary),
    " features with S vs NS p<=", consistent_p, "; ",
    sum(consistent_summary$SS10_agree_direction != ""),
    " agree with SS-10 direction\n", sep = "")
cat("\nMost abundant features:\n")
print(head(feature_ranking %>% select(feature_id, feature, mean_abundance, prevalence), 10),
      row.names = FALSE)
