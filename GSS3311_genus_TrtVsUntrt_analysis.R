###############################################################################
## GSS3311 Falcon2 — taxon-level (genus, species, ...) Treated vs Untreated analysis
## -----------------------------------------------------------------------------
## Modelled on get_paired_wilcox_site_GSS3318_analysis15.R
##
## Design (from sample names / metadata)
##   Sample name prefix = CohortTrt: AA, BB, CC, DD (treated side) and
##   AN, BN, CN, DN (untreated side of the same subject).  Subjects are paired
##   by the metadata column "Identity".
##
## Per feature (row of the data, e.g. genus or species) the script computes
##   1. Paired Wilcoxon signed-rank test (by Identity) + fold change for
##      AA:AN, BB:BN, CC:CN, DD:DN  (FC = mean treated / mean untreated,
##      reported as log2FC and signed "true" FC; also median paired diff).
##   2. Kruskal-Wallis p values for treatment differences:
##        KW_8groups      : AA,AN,BB,BN,CC,CN,DD,DN
##        KW_Treatment5   : A,B,C,D (treated) and N (all untreated pooled)
##        KW_TreatedOnly  : AA,BB,CC,DD
##        KW_UntreatedOnly: AN,BN,CN,DN  (baseline / cohort check)
##        KW_TrtVsUntrt   : Treated vs Untreated (all cohorts pooled)
##   3. Pearson and Spearman correlation (r/rho, p, BH-FDR) with
##      PH.average, Irritation, TotalDNAYield_ng, MaxIrritationPerRandSide,
##      ss10_Total, ss10_Itch, ss10_Redness, and ss10_Pain for All samples,
##      Treated only, Untreated only.
##   4. Unpaired Wilcoxon responder vs nonresponder tests, separately for all
##      treated and all untreated samples, with group means, SE, counts, FC, p,
##      and BH-FDR.
##   BH-FDR is computed across features for every test column.
##   Significance marks: * p<0.1, ** p<0.05, *** FDR<0.1, **** FDR<0.05.
##   consistent_* : Pearson & Spearman both p<=0.1 with the same sign ->
##     "positive"/"negative" + stars from the weaker (larger) p/FDR of the two.
##
## Input data: tab-delimited, first column = feature (taxon) ID, other columns =
##   samples. Any taxonomic level works (genus_relab10_data, species_relab10_data, ...).
##
## Usage
##   Rscript GSS3311_genus_TrtVsUntrt_analysis.R [data_file] [meta_file] [outdir] [p_plot] [log_dir]
##   defaults: ../genus_relab10_data GSS3311_meta_PH101326.txt  ./  0.05  <outdir>/<base>_result/logs
##   All output goes to <outdir>/<basename(data_file)>_result/
##     e.g. ./genus_relab10_data_result/, ./species_relab10_data_result/
##     containing results/, figures/ (including separate responder treated/untreated plots).
##   The stdout/summary log is written to its own directory (log_dir), which
##   defaults to <base>_result/logs/ but can be redirected anywhere via the 5th arg.
###############################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(openxlsx)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
data_file <- if (length(args) >= 1) args[1] else "../genus_relab10_data"
meta_file <- if (length(args) >= 2) args[2] else "GSS3311_meta_PH101326.txt"
outdir    <- if (length(args) >= 3) args[3] else "."
p_plot    <- if (length(args) >= 4) as.numeric(args[4]) else 0.05
log_dir   <- if (length(args) >= 5) args[5] else NA   # where the stdout summary log goes
if (!file.exists(data_file)) stop("Data file not found: ", data_file)
if (!file.exists(meta_file)) stop("Metadata file not found: ", meta_file)

prefix  <- basename(data_file)
run_dir <- file.path(outdir, paste0(prefix, "_result"))
res_dir <- file.path(run_dir, "results")
fig_dir <- file.path(run_dir, "figures")
ind_dir <- file.path(fig_dir, "individual")
ind_trt_dir <- file.path(ind_dir, "treatment_comparison")
ind_cor_dir <- file.path(ind_dir, "correlation")
ind_bar_dir <- file.path(ind_dir, "group_barplot")
responder_taxon_dir <- file.path(fig_dir, "responder_analysis", "by_taxon")
## Standard-out summary log goes to its own directory, separate from results/figures.
if (is.na(log_dir)) log_dir <- file.path(run_dir, "logs")
## ind_trt_dir is intentionally not created while its plot block is disabled.
for (d in c(res_dir, fig_dir, ind_cor_dir, ind_bar_dir, log_dir))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
dir.create(responder_taxon_dir, recursive = TRUE, showWarnings = FALSE)
## remove stale individual plots so the subfolders reflect the current run only
unlink(list.files(ind_dir, pattern = "\\.png$", full.names = TRUE, recursive = TRUE))
unlink(list.files(file.path(fig_dir, "responder_analysis"),
                  pattern = "\\.png$", full.names = TRUE, recursive = TRUE))

log_file <- file.path(log_dir, paste0(prefix, ".run_log.txt"))
sink(log_file, split = TRUE)
cat("Run:", format(Sys.time()), "\nData:", normalizePath(data_file), "\nMeta:", normalizePath(meta_file),
    "\nOutput:", normalizePath(run_dir), "\nLog dir:", normalizePath(log_dir), "\np_plot:", p_plot, "\n\n")

rout <- function(...) file.path(res_dir, paste0(prefix, ".", ...))
fout <- function(...) file.path(fig_dir, paste0(prefix, ".", ...))
unlink(c(rout("PH_refined_correlations.tsv"),
         rout("PH_refined_removed_samples.tsv")))

pairs     <- c(A = "AA:AN", B = "BB:BN", C = "CC:CN", D = "DD:DN")
groups8   <- c("AA", "AN", "BB", "BN", "CC", "CN", "DD", "DN")
treated4  <- c("AA", "BB", "CC", "DD")                 # treated side, one per cohort
untrt4    <- c("AN", "BN", "CN", "DN")                 # untreated side
trt_pairs <- combn(treated4, 2, simplify = FALSE)      # 6 between-treatment pairs
untrt_pairs <- combn(untrt4, 2, simplify = FALSE)      # 6 between-cohort untreated pairs
cor_vars  <- c("PH.average", "Irritation", "TotalDNAYield_ng",
               "MaxIrritationPerRandSide", "ss10_Total",
               "ss10_Itch", "ss10_Redness", "ss10_Pain")
screened_cor_vars <- c("PH.average", "Irritation",
                       "ss10_Itch", "ss10_Redness", "ss10_Pain")
dedicated_cor_vars <- setdiff(screened_cor_vars, "PH.average")
cor_sets  <- c("All", "Treated", "Untreated")
pal_grp   <- c(AA = "#C44E52", AN = "#F2A6A8", BB = "#4C72B0", BN = "#A7BEE0",
               CC = "#55A868", CN = "#AAD8B4", DD = "#8172B3", DN = "#C9BFE6")

## --------------------------------------------------------------------------
## Helpers (same conventions as the analysis15 example)
## --------------------------------------------------------------------------
theme_pub <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = "grey92"),
          strip.background = element_rect(fill = "grey95", colour = NA),
          strip.text = element_text(face = "bold"),
          plot.title = element_text(face = "bold", hjust = 0),
          legend.position = "top")
}
safe_p <- function(expr) tryCatch(expr, error = function(e) NA_real_)

paired_wilcox_p <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]; y <- y[keep]
  if (length(x) < 2 || all(x - y == 0)) return(NA_real_)
  safe_p(suppressWarnings(wilcox.test(x, y, paired = TRUE, exact = FALSE)$p.value))
}

## Unpaired Wilcoxon / Mann-Whitney for comparing two independent groups.
mwu_p <- function(x, y) {
  x <- x[is.finite(x)]; y <- y[is.finite(y)]
  if (length(x) < 2 || length(y) < 2 || length(unique(c(x, y))) < 2) return(NA_real_)
  safe_p(suppressWarnings(wilcox.test(x, y, exact = FALSE)$p.value))
}

kw_p <- function(v, g) {
  keep <- is.finite(v) & !is.na(g)
  v <- v[keep]; g <- droplevels(as.factor(g[keep]))
  if (length(v) < 3 || nlevels(g) < 2 || length(unique(v)) < 2) return(NA_real_)
  safe_p(kruskal.test(v, g)$p.value)
}

cor_one <- function(x, y, method) {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]; y <- y[keep]
  if (length(x) < 4 || length(unique(x)) < 2 || length(unique(y)) < 2)
    return(c(n = length(x), r = NA_real_, p = NA_real_))
  ct <- tryCatch(suppressWarnings(
    if (method == "spearman") cor.test(x, y, method = method, exact = FALSE)
    else cor.test(x, y, method = method)), error = function(e) NULL)
  if (is.null(ct)) return(c(n = length(x), r = NA_real_, p = NA_real_))
  c(n = length(x), r = unname(ct$estimate), p = ct$p.value)
}

truefc <- function(ratio) {
  ## Signed fold change: 2x higher = +2; 2x lower = -2.
  ifelse(!is.finite(ratio) | ratio == 0, NA_real_, ifelse(ratio < 1, -1 / ratio, ratio))
}

sig_mark <- function(p, fdr) {
  case_when(!is.na(fdr) & fdr < 0.05 ~ "****",
            !is.na(fdr) & fdr < 0.10 ~ "***",
            !is.na(p) & p < 0.05 ~ "**",
            !is.na(p) & p < 0.10 ~ "*",
            TRUE ~ "")
}
sig_legend <- "* p<0.1  ** p<0.05  *** FDR<0.1  **** FDR<0.05"

## Consistent correlation: Pearson & Spearman both p<=0.1 with the same sign.
## Stars use the weaker (larger) p and FDR of the two methods.
consistent_mark <- function(pr, pp, pf, sr, sp, sf) {
  ok <- !is.na(pr) & !is.na(sr) & !is.na(pp) & !is.na(sp) &
    pp <= 0.1 & sp <= 0.1 & sign(pr) == sign(sr) & sign(sr) != 0
  stars <- sig_mark(pmax(pp, sp), pmax(pf, sf))
  ifelse(ok, paste0(ifelse(sr > 0, "positive", "negative"), stars), "")
}

clean_name <- function(x) sub("^[a-z]__", "", x)
safe_file  <- function(x) gsub("[^A-Za-z0-9_.-]+", "_", x)

## --------------------------------------------------------------------------
## Read data
## --------------------------------------------------------------------------
dat <- read.delim(data_file, check.names = FALSE, stringsAsFactors = FALSE)
feature_ids <- dat[[1]]
mat <- as.matrix(dat[, -1, drop = FALSE])
storage.mode(mat) <- "numeric"
rownames(mat) <- feature_ids

meta <- read.delim(meta_file, check.names = FALSE, stringsAsFactors = FALSE)
names(meta)[1] <- "Sample"
responder_col <- intersect(c("responder_Y2_N1", "responder_N1_Y2"), names(meta))
if (!length(responder_col))
  stop("Metadata must contain responder_Y2_N1 or responder_N1_Y2.")
responder_col <- responder_col[1]
meta <- meta %>%
  mutate(Group = substr(Sample, 1, 2)) %>%
  filter(Group %in% groups8, Sample %in% colnames(mat)) %>%
  mutate(Cohort = substr(Group, 1, 1),
         Status = ifelse(substr(Group, 2, 2) == "N", "Untreated", "Treated"),
         ResponderStatus = case_when(
           as.character(.data[[responder_col]]) == "1" ~ "Nonresponder",
           as.character(.data[[responder_col]]) == "2" ~ "Responder",
           TRUE ~ NA_character_
         ),
         Treatment5 = ifelse(Status == "Treated", Cohort, "N"),
         Group = factor(Group, levels = groups8))
## Responder status is a subject-level outcome recorded on treated samples;
## carry it to each matched untreated sample within cohort and Identity.
resp_map <- meta %>%
  filter(Status == "Treated", !is.na(ResponderStatus)) %>%
  distinct(Cohort, Identity, ResponderStatus) %>%
  rename(TreatedResponderStatus = ResponderStatus)
resp_conflicts <- resp_map %>% count(Cohort, Identity) %>% filter(n > 1)
if (nrow(resp_conflicts))
  stop("Conflicting responder labels among treated samples for cohort/Identity: ",
       paste(paste(resp_conflicts$Cohort, resp_conflicts$Identity, sep = "/"), collapse = ", "))
meta <- meta %>%
  left_join(resp_map, by = c("Cohort", "Identity")) %>%
  mutate(ResponderStatus = ifelse(Status == "Untreated",
                                  TreatedResponderStatus, ResponderStatus)) %>%
  select(-TreatedResponderStatus)

missing_meta <- setdiff(colnames(mat), meta$Sample)
if (length(missing_meta))
  warning(length(missing_meta), " data samples have no metadata and are dropped: ",
          paste(head(missing_meta, 5), collapse = ", "))
mat <- mat[, meta$Sample, drop = FALSE]
for (v in cor_vars) meta[[v]] <- suppressWarnings(as.numeric(meta[[v]]))
invalid_ph <- is.finite(meta$PH.average) & (meta$PH.average < 0 | meta$PH.average > 14)
ph_invalid_samples <- meta[invalid_ph, c("Sample", "Identity", "Group", "Status", "PH.average")]
write.table(ph_invalid_samples, rout("PH_invalid_removed_samples.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
if (any(invalid_ph)) {
  msg <- paste0(
    "PH.average values outside the valid [0,14] range were removed from all PH correlations: ",
    paste0(ph_invalid_samples$Sample, " [Identity=", ph_invalid_samples$Identity,
           ", Group=", ph_invalid_samples$Group, ", PH.average=",
           ph_invalid_samples$PH.average, "]", collapse = "; ")
  )
  cat("ERROR:", msg, "\n")
  message("ERROR: ", msg)
  meta$PH.average[invalid_ph] <- NA_real_
}

## Remove empty features (rows) and empty samples (columns): all NA or all 0.
## Removed items are reported to the log/console (stderr too) and saved to results/.
is_empty <- function(x) all(is.na(x) | x == 0)
removed <- data.frame(Type = character(), ID = character(), Reason = character())
empty_cols <- apply(mat, 2, is_empty)
if (any(empty_cols)) {
  ids <- colnames(mat)[empty_cols]
  rsn <- ifelse(apply(mat[, empty_cols, drop = FALSE], 2, function(x) all(is.na(x))), "all NA", "all 0 / NA")
  msg <- paste0("Removed ", sum(empty_cols), " sample column(s) that are all NA/0 before analysis: ",
                paste0(ids, " [", rsn, "]", collapse = ", "))
  cat("WARNING:", msg, "\n"); message("WARNING: ", msg)
  removed <- rbind(removed, data.frame(Type = "sample_column", ID = ids, Reason = rsn))
  mat  <- mat[, !empty_cols, drop = FALSE]
  meta <- meta[!empty_cols, , drop = FALSE]
}
empty_rows <- apply(mat, 1, is_empty)
if (any(empty_rows)) {
  ids <- rownames(mat)[empty_rows]
  msg <- paste0("Removed ", sum(empty_rows), " feature row(s) that are all NA/0: ",
                paste(head(ids, 10), collapse = ", "), if (length(ids) > 10) " ..." else "")
  cat("NOTE:", msg, "\n"); message("NOTE: ", msg)
  removed <- rbind(removed, data.frame(Type = "feature_row", ID = ids, Reason = "all NA/0"))
  mat <- mat[!empty_rows, , drop = FALSE]
}
write.table(removed, rout("removed_empty.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
## Human-readable data-problem report at the top of the output folder
qc_file <- file.path(run_dir, paste0(prefix, ".DATA_QC_removed.txt"))
writeLines(c(
  paste0("Data QC report for: ", normalizePath(data_file)),
  paste0("Generated: ", format(Sys.time())),
  "Rule: any sample column or feature row that is entirely NA or 0 is removed before analysis.",
  "",
  paste0("Sample columns removed: ", sum(removed$Type == "sample_column")),
  if (any(removed$Type == "sample_column"))
    paste0("  ", removed$ID[removed$Type == "sample_column"], "\t[", removed$Reason[removed$Type == "sample_column"], "]"),
  "",
  paste0("Feature rows removed: ", sum(removed$Type == "feature_row")),
  if (any(removed$Type == "feature_row")) paste0("  ", removed$ID[removed$Type == "feature_row"]),
  "",
  if (!nrow(removed)) "No data problems found: nothing removed." else
    "Note: a removed sample's paired partner is kept for unpaired tests (KW, correlation) but drops out of the paired Wilcoxon."
), qc_file)
if (anyNA(mat)) cat("Note:", sum(is.na(mat)), "remaining NA cells are ignored per test\n")
feature_ids <- rownames(mat)

## Treatment names from the 2nd "."-separated field of treated samples (AA, BB, CC, DD)
## Untreated samples (AN..DN) are all "NoTreatment" and are not used for naming.
trt_name <- sapply(names(pairs), function(co) {
  s <- meta$Sample[meta$Group == paste0(co, co)]
  nm <- unique(sapply(strsplit(s, ".", fixed = TRUE), function(z) if (length(z) >= 2) z[2] else NA))
  nm <- nm[!is.na(nm) & nzchar(nm)]
  if (length(nm) > 1) warning("Cohort ", co, " has multiple treatment names: ", paste(nm, collapse = ", "))
  if (length(nm)) nm[1] else co
})
cmp_lab <- setNames(paste0(pairs, " (", trt_name, ")"), pairs)
grp_lab <- setNames(c(rbind(paste0(names(pairs), names(pairs), "\n", trt_name),
                            paste0(names(pairs), "N\n", trt_name, " untrt"))),
                    c(rbind(paste0(names(pairs), names(pairs)), paste0(names(pairs), "N"))))
grp_lab1 <- gsub("\n", " ", grp_lab)
cat("Treatment names:", paste0(names(trt_name), "=", trt_name, collapse = ", "), "\n")

cat("Features:", nrow(mat), " Samples:", ncol(mat), "\n")
print(table(meta$Group))
cat("Responder status by treatment status (untreated labels inherited from matched treated subject):\n")
print(table(meta$Status, meta$ResponderStatus, useNA = "ifany"))
if (!any(!is.na(meta$ResponderStatus) & meta$Status == "Treated") ||
    !any(!is.na(meta$ResponderStatus) & meta$Status == "Untreated"))
  stop("Responder status must be coded 1/2 for at least one treated and one untreated sample.")

## Paired sample indices by Identity for each cohort
pair_index <- lapply(names(pairs), function(co) {
  tg <- paste0(co, co); ug <- paste0(co, "N")
  t <- meta %>% filter(Group == tg) %>% select(Identity, t_sample = Sample)
  u <- meta %>% filter(Group == ug) %>% select(Identity, u_sample = Sample)
  inner_join(t, u, by = "Identity")
})
names(pair_index) <- names(pairs)
cat("Paired subjects per cohort:", sapply(pair_index, nrow), "\n")

## --------------------------------------------------------------------------
## Per-feature statistics
## --------------------------------------------------------------------------
eps <- min(mat[mat > 0], na.rm = TRUE) / 2   # pseudocount for log2FC

one_feature <- function(i) {
  v <- mat[i, ]
  row <- list(FeatureID = feature_ids[i], Taxon = clean_name(feature_ids[i]),
              Mean_All = mean(v, na.rm = TRUE), Prevalence_All = mean(v > 0, na.rm = TRUE))
  ## group descriptives
  for (g in groups8) {
    vg <- v[meta$Group == g]; vg <- vg[!is.na(vg)]
    row[[paste0("n_", g)]]    <- length(vg)
    row[[paste0("Mean_", g)]] <- mean(vg)
    row[[paste0("SE_", g)]]   <- sd(vg) / sqrt(length(vg))
    row[[paste0("Prev_", g)]] <- mean(vg > 0)
  }
  ## paired Wilcoxon + fold change (complete pairs only)
  for (co in names(pairs)) {
    pi <- pair_index[[co]]
    x <- v[pi$t_sample]; y <- v[pi$u_sample]
    ok <- !is.na(x) & !is.na(y); x <- x[ok]; y <- y[ok]
    lab <- paste0(co, co, "vs", co, "N")
    mt <- mean(x); mu <- mean(y)
    row[[paste0("Pairs_", lab)]]       <- length(x)
    row[[paste0("MeanTrt_", lab)]]     <- mt
    row[[paste0("MeanUntrt_", lab)]]   <- mu
    row[[paste0("MedianDiff_", lab)]]  <- median(x - y)
    row[[paste0("nUp_", lab)]]         <- sum(x > y)
    row[[paste0("nDown_", lab)]]       <- sum(x < y)
    row[[paste0("log2FC_", lab)]]      <- log2((mt + eps) / (mu + eps))
    row[[paste0("FC_", lab)]]          <- truefc((mt + eps) / (mu + eps))
    row[[paste0("WilcoxP_", lab)]]     <- paired_wilcox_p(x, y)
  }
  ## Unpaired responder vs nonresponder comparisons, pooled over all cohorts
  for (s in c("Treated", "Untreated")) {
    sel <- meta$Status == s & !is.na(meta$ResponderStatus)
    x <- v[sel & meta$ResponderStatus == "Responder"]
    y <- v[sel & meta$ResponderStatus == "Nonresponder"]
    x <- x[is.finite(x)]; y <- y[is.finite(y)]
    mx <- if (length(x)) mean(x) else NA_real_
    my <- if (length(y)) mean(y) else NA_real_
    row[[paste0("Responder_N_", s)]] <- length(x)
    row[[paste0("Responder_Mean_", s)]] <- mx
    row[[paste0("Responder_SE_", s)]] <- if (length(x) > 1) sd(x) / sqrt(length(x)) else NA_real_
    row[[paste0("Nonresponder_N_", s)]] <- length(y)
    row[[paste0("Nonresponder_Mean_", s)]] <- my
    row[[paste0("Nonresponder_SE_", s)]] <- if (length(y) > 1) sd(y) / sqrt(length(y)) else NA_real_
    row[[paste0("ResponderFC_", s)]] <- if (is.finite(mx) && is.finite(my))
      (mx + eps) / (my + eps) else NA_real_
    row[[paste0("WilcoxRespP_", s)]] <- mwu_p(x, y)
  }
  ## Responder vs nonresponder within each individual sample group.
  ## Suppress the test unless the group has >=3 responders and >=2 nonresponders.
  for (g in groups8) {
    sel <- meta$Group == g & !is.na(meta$ResponderStatus)
    x <- v[sel & meta$ResponderStatus == "Responder"]
    y <- v[sel & meta$ResponderStatus == "Nonresponder"]
    x <- x[is.finite(x)]; y <- y[is.finite(y)]
    mx <- if (length(x)) mean(x) else NA_real_
    my <- if (length(y)) mean(y) else NA_real_
    row[[paste0("ResponderGroup_N_", g)]] <- length(x)
    row[[paste0("ResponderGroup_Mean_", g)]] <- mx
    row[[paste0("ResponderGroup_SE_", g)]] <- if (length(x) > 1) sd(x) / sqrt(length(x)) else NA_real_
    row[[paste0("NonresponderGroup_N_", g)]] <- length(y)
    row[[paste0("NonresponderGroup_Mean_", g)]] <- my
    row[[paste0("NonresponderGroup_SE_", g)]] <- if (length(y) > 1) sd(y) / sqrt(length(y)) else NA_real_
    row[[paste0("ResponderGroup_FC_", g)]] <- if (is.finite(mx) && is.finite(my))
      (mx + eps) / (my + eps) else NA_real_
    enough_responders <- sum(sel & meta$ResponderStatus == "Responder") >= 3
    enough_nonresponders <- sum(sel & meta$ResponderStatus == "Nonresponder") >= 2
    row[[paste0("WilcoxRespGroupP_", g)]] <- if (enough_responders && enough_nonresponders)
      mwu_p(x, y) else NA_real_
  }
  ## Kruskal-Wallis
  row$KW_8groups       <- kw_p(v, meta$Group)
  row$KW_Treatment5    <- kw_p(v, meta$Treatment5)
  tr <- meta$Status == "Treated"
  ## (1) treatment effect on the TREATED side: AA,BB,CC,DD
  row$KW_TreatedOnly   <- kw_p(v[tr], meta$Group[tr])
  ## (2) treatment effect on the UNTREATED side: AN,BN,CN,DN
  row$KW_UntreatedOnly <- kw_p(v[!tr], meta$Group[!tr])
  row$KW_TrtVsUntrt    <- kw_p(v, meta$Status)
  ## Pairwise unpaired Wilcoxon between treatments on the treated side (AA,BB,CC,DD)
  for (pr in trt_pairs) {
    lab <- paste0(pr[1], "vs", pr[2])
    row[[paste0("WilcoxTrt_", lab)]] <- mwu_p(v[meta$Group == pr[1]], v[meta$Group == pr[2]])
  }
  ## Pairwise unpaired Wilcoxon between untreated cohorts (AN,BN,CN,DN).
  for (pr in untrt_pairs) {
    lab <- paste0(pr[1], "vs", pr[2])
    row[[paste0("WilcoxUntreat_", lab)]] <- mwu_p(v[meta$Group == pr[1]], v[meta$Group == pr[2]])
  }
  ## correlations
  for (s in cor_sets) {
    sel <- switch(s, All = rep(TRUE, ncol(mat)), Treated = tr, Untreated = !tr)
    for (cv in cor_vars) {
      for (m in c("pearson", "spearman")) {
        cs <- cor_one(v[sel], meta[[cv]][sel], m)
        tag <- paste0(cv, "_", s)
        row[[paste0(m, "_r_", tag)]] <- cs[["r"]]
        row[[paste0(m, "_p_", tag)]] <- cs[["p"]]
      }
    }
  }
  as_tibble(row)
}

master <- bind_rows(lapply(seq_len(nrow(mat)), one_feature))

## BH-FDR across taxa for every p-value column
p_cols <- grep("^(WilcoxP_|WilcoxTrt_|WilcoxUntreat_|WilcoxRespP_|WilcoxRespGroupP_|KW_|pearson_p_|spearman_p_)", names(master), value = TRUE)
for (pc in p_cols) {
  fdr_name <- sub("WilcoxP_", "WilcoxFDR_", sub("WilcoxTrt_", "WilcoxTrtFDR_", sub("_p_", "_FDR_", pc)))
  if (grepl("^WilcoxUntreat_", pc)) fdr_name <- sub("^WilcoxUntreat_", "WilcoxUntreatFDR_", pc)
  if (grepl("^WilcoxRespP_", pc)) fdr_name <- sub("^WilcoxRespP_", "WilcoxRespFDR_", pc)
  if (grepl("^WilcoxRespGroupP_", pc))
    fdr_name <- sub("^WilcoxRespGroupP_", "WilcoxRespGroupFDR_", pc)
  if (grepl("^KW_", pc)) fdr_name <- paste0(pc, "_FDR")
  master[[fdr_name]] <- p.adjust(master[[pc]], method = "BH")
}
## re-order: put FDR right after its p column
ord <- names(master)[!names(master) %in% grep("FDR", names(master), value = TRUE)]
final_ord <- unlist(lapply(ord, function(n) {
  f <- if (grepl("^WilcoxP_", n)) sub("WilcoxP_", "WilcoxFDR_", n)
       else if (grepl("^WilcoxTrt_", n)) sub("WilcoxTrt_", "WilcoxTrtFDR_", n)
       else if (grepl("^WilcoxUntreat_", n)) sub("^WilcoxUntreat_", "WilcoxUntreatFDR_", n)
       else if (grepl("^WilcoxRespP_", n)) sub("^WilcoxRespP_", "WilcoxRespFDR_", n)
       else if (grepl("^WilcoxRespGroupP_", n)) sub("^WilcoxRespGroupP_", "WilcoxRespGroupFDR_", n)
       else if (grepl("^KW_", n)) paste0(n, "_FDR")
       else if (grepl("_p_", n)) sub("_p_", "_FDR_", n) else NULL
  c(n, f)
}))
master <- master[, final_ord]

## "solid" correlation flag: Pearson and Spearman agree in sign and both p<=0.05
for (s in cor_sets) for (cv in cor_vars) {
  tag <- paste0(cv, "_", s)
  pr <- master[[paste0("pearson_r_", tag)]];  pp <- master[[paste0("pearson_p_", tag)]]
  sr <- master[[paste0("spearman_r_", tag)]]; sp <- master[[paste0("spearman_p_", tag)]]
  master[[paste0("solid_", tag)]] <- ifelse(!is.na(pp) & !is.na(sp) & pp <= 0.05 & sp <= 0.05 &
                                              sign(pr) == sign(sr),
                                            ifelse(sr > 0, "pos", "neg"), "")
  master[[paste0("consistent_", tag)]] <- consistent_mark(
    pr, pp, master[[paste0("pearson_FDR_", tag)]],
    sr, sp, master[[paste0("spearman_FDR_", tag)]])
}
## Significance marks for paired Wilcoxon and KW tests
for (co in names(pairs)) {
  lab <- paste0(co, co, "vs", co, "N")
  master[[paste0("WilcoxSig_", lab)]] <- sig_mark(master[[paste0("WilcoxP_", lab)]],
                                                  master[[paste0("WilcoxFDR_", lab)]])
}
## Significance marks for the between-treatment pairwise tests (AA,BB,CC,DD)
for (pr in trt_pairs) {
  lab <- paste0(pr[1], "vs", pr[2])
  master[[paste0("WilcoxTrtSig_", lab)]] <- sig_mark(master[[paste0("WilcoxTrt_", lab)]],
                                                     master[[paste0("WilcoxTrtFDR_", lab)]])
}
for (pr in untrt_pairs) {
  lab <- paste0(pr[1], "vs", pr[2])
  master[[paste0("WilcoxUntreatSig_", lab)]] <- sig_mark(
    master[[paste0("WilcoxUntreat_", lab)]],
    master[[paste0("WilcoxUntreatFDR_", lab)]])
}
for (s in c("Treated", "Untreated"))
  master[[paste0("WilcoxRespSig_", s)]] <- sig_mark(master[[paste0("WilcoxRespP_", s)]],
                                                    master[[paste0("WilcoxRespFDR_", s)]])
for (g in groups8)
  master[[paste0("WilcoxRespGroupSig_", g)]] <- sig_mark(
    master[[paste0("WilcoxRespGroupP_", g)]],
    master[[paste0("WilcoxRespGroupFDR_", g)]])
for (k in grep("^KW_[A-Za-z0-9]+$", names(master), value = TRUE))
  master[[paste0(k, "_Sig")]] <- sig_mark(master[[k]], master[[paste0(k, "_FDR")]])
cons_cols <- grep("^consistent_", names(master), value = TRUE)
master$n_consistent_cor <- rowSums(master[, cons_cols] != "")
master$consistent_cor_list <- apply(master[, cons_cols], 1, function(x) {
  x <- x[x != ""]
  if (!length(x)) "" else paste0(sub("^consistent_", "", names(x)), ":", x, collapse = "; ")
})

write.table(master, rout("master_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

## --------------------------------------------------------------------------
## Long tables
## --------------------------------------------------------------------------
paired_long <- bind_rows(lapply(names(pairs), function(co) {
  lab <- paste0(co, co, "vs", co, "N")
  tibble(FeatureID = master$FeatureID, Taxon = master$Taxon, Comparison = pairs[[co]],
         Treatment = trt_name[[co]],
         Pairs = master[[paste0("Pairs_", lab)]],
         MeanTreated = master[[paste0("MeanTrt_", lab)]],
         MeanUntreated = master[[paste0("MeanUntrt_", lab)]],
         MedianPairedDiff = master[[paste0("MedianDiff_", lab)]],
         nUp = master[[paste0("nUp_", lab)]], nDown = master[[paste0("nDown_", lab)]],
         log2FC = master[[paste0("log2FC_", lab)]], FC = master[[paste0("FC_", lab)]],
         Wilcox_p = master[[paste0("WilcoxP_", lab)]],
         Wilcox_FDR = master[[paste0("WilcoxFDR_", lab)]])
})) %>% mutate(Sig = sig_mark(Wilcox_p, Wilcox_FDR)) %>% arrange(Comparison, Wilcox_p)
write.table(paired_long, rout("paired_wilcox_TrtVsUntrt.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

kw_tab <- master %>% select(FeatureID, Taxon, Mean_All, Prevalence_All, starts_with("Mean_"),
                            starts_with("KW_")) %>% arrange(KW_8groups)
write.table(kw_tab, rout("kruskal_wallis_treatment.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

## Between-treatment pairwise comparisons on the treated side (AA,BB,CC,DD)
trt_pairwise_long <- bind_rows(lapply(trt_pairs, function(pr) {
  lab <- paste0(pr[1], "vs", pr[2])
  tibble(FeatureID = master$FeatureID, Taxon = master$Taxon,
         Comparison = paste0(pr[1], " (", trt_name[[substr(pr[1], 1, 1)]], ") vs ",
                             pr[2], " (", trt_name[[substr(pr[2], 1, 1)]], ")"),
         Group1 = pr[1], Group2 = pr[2],
         Mean1 = master[[paste0("Mean_", pr[1])]], Mean2 = master[[paste0("Mean_", pr[2])]],
         Wilcox_p = master[[paste0("WilcoxTrt_", lab)]],
         Wilcox_FDR = master[[paste0("WilcoxTrtFDR_", lab)]])
})) %>% mutate(Sig = sig_mark(Wilcox_p, Wilcox_FDR)) %>% arrange(Comparison, Wilcox_p)
write.table(trt_pairwise_long, rout("treated_pairwise_wilcox.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

untreated_pairwise_long <- bind_rows(lapply(untrt_pairs, function(pr) {
  lab <- paste0(pr[1], "vs", pr[2])
  tibble(FeatureID = master$FeatureID, Taxon = master$Taxon,
         Comparison = paste(pr[1], "vs", pr[2]),
         Group1 = pr[1], Group2 = pr[2],
         Mean1 = master[[paste0("Mean_", pr[1])]],
         Mean2 = master[[paste0("Mean_", pr[2])]],
         Wilcox_p = master[[paste0("WilcoxUntreat_", lab)]],
         Wilcox_FDR = master[[paste0("WilcoxUntreatFDR_", lab)]])
})) %>% mutate(Sig = sig_mark(Wilcox_p, Wilcox_FDR)) %>% arrange(Comparison, Wilcox_p)
write.table(untreated_pairwise_long, rout("untreated_pairwise_wilcox.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

responder_long <- function(s) {
  tibble(
    FeatureID = master$FeatureID,
    Taxon = master$Taxon,
    Comparison = paste0(s, ": Responder vs Nonresponder (unpaired Wilcoxon)"),
    Responder_n = master[[paste0("Responder_N_", s)]],
    Responder_Mean = master[[paste0("Responder_Mean_", s)]],
    Responder_SE = master[[paste0("Responder_SE_", s)]],
    Nonresponder_n = master[[paste0("Nonresponder_N_", s)]],
    Nonresponder_Mean = master[[paste0("Nonresponder_Mean_", s)]],
    Nonresponder_SE = master[[paste0("Nonresponder_SE_", s)]],
    FoldChange_Responder_Nonresponder = master[[paste0("ResponderFC_", s)]],
    Wilcox_p = master[[paste0("WilcoxRespP_", s)]],
    Wilcox_FDR = master[[paste0("WilcoxRespFDR_", s)]]
  ) %>%
    mutate(Sig = sig_mark(Wilcox_p, Wilcox_FDR)) %>%
    arrange(Wilcox_p)
}
responder_treated_long <- responder_long("Treated")
responder_untreated_long <- responder_long("Untreated")
write.table(responder_treated_long, rout("responder_vs_nonresponder_treated.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(responder_untreated_long, rout("responder_vs_nonresponder_untreated.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
responder_group_long <- bind_rows(lapply(groups8, function(g) {
  tibble(
    FeatureID = master$FeatureID,
    Taxon = master$Taxon,
    Group = g,
    Responder_n = master[[paste0("ResponderGroup_N_", g)]],
    Responder_Mean = master[[paste0("ResponderGroup_Mean_", g)]],
    Responder_SE = master[[paste0("ResponderGroup_SE_", g)]],
    Nonresponder_n = master[[paste0("NonresponderGroup_N_", g)]],
    Nonresponder_Mean = master[[paste0("NonresponderGroup_Mean_", g)]],
    Nonresponder_SE = master[[paste0("NonresponderGroup_SE_", g)]],
    FoldChange_Responder_Nonresponder = master[[paste0("ResponderGroup_FC_", g)]],
    TestEligible = Responder_n >= 3 & Nonresponder_n >= 2,
    Wilcox_p = master[[paste0("WilcoxRespGroupP_", g)]],
    Wilcox_FDR = master[[paste0("WilcoxRespGroupFDR_", g)]]
  )
})) %>%
  mutate(Sig = sig_mark(Wilcox_p, Wilcox_FDR)) %>%
  arrange(factor(Group, levels = groups8), Wilcox_p)
write.table(responder_group_long, rout("responder_vs_nonresponder_by_group.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cor_long <- bind_rows(lapply(cor_sets, function(s) bind_rows(lapply(cor_vars, function(cv) {
  tag <- paste0(cv, "_", s)
  tibble(FeatureID = master$FeatureID, Taxon = master$Taxon, Subset = s, Variable = cv,
         pearson_r = master[[paste0("pearson_r_", tag)]],
         pearson_p = master[[paste0("pearson_p_", tag)]],
         pearson_FDR = master[[paste0("pearson_FDR_", tag)]],
         spearman_rho = master[[paste0("spearman_r_", tag)]],
         spearman_p = master[[paste0("spearman_p_", tag)]],
         spearman_FDR = master[[paste0("spearman_FDR_", tag)]],
         pearson_Sig = sig_mark(pearson_p, pearson_FDR),
         spearman_Sig = sig_mark(spearman_p, spearman_FDR),
         solid = master[[paste0("solid_", tag)]],
         consistent = master[[paste0("consistent_", tag)]])
})))) %>% arrange(Subset, Variable, spearman_p)
write.table(cor_long, rout("correlations.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
dedicated_cor_tables <- setNames(lapply(dedicated_cor_vars, function(cv) {
  out <- cor_long %>%
    filter(Variable == cv) %>%
    arrange(Subset, pmin(pearson_p, spearman_p, na.rm = TRUE))
  write.table(out, rout(paste0(cv, "_correlations.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  out
}), dedicated_cor_vars)
cor_consistent <- cor_long %>% filter(consistent != "") %>%
  mutate(Direction = sub("\\*+$", "", consistent)) %>%
  arrange(Subset, Variable, pmax(pearson_p, spearman_p))
write.table(cor_consistent, rout("correlations_consistent.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

## Sample-level merged table (meta + abundance) for reference
sample_tab <- bind_cols(meta %>% select(Sample, Identity, Group, Cohort, Status, ResponderStatus, all_of(cor_vars)),
                        as.data.frame(t(mat), check.names = FALSE))
write.table(sample_tab, rout("sample_level_merged.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

## --------------------------------------------------------------------------
## Excel workbook
## --------------------------------------------------------------------------
wb <- createWorkbook()
hs <- createStyle(textDecoration = "bold", fgFill = "#DDEBF7", border = "Bottom", wrapText = TRUE)
st_sig  <- createStyle(fgFill = "#FFEB9C")
st_up   <- createStyle(fontColour = "#C00000", textDecoration = "bold")
st_down <- createStyle(fontColour = "#1F4E79", textDecoration = "bold")

add_sheet <- function(name, df, p_like = NULL, signed = NULL, dir_text = NULL) {
  addWorksheet(wb, name)
  writeData(wb, name, df, headerStyle = hs)
  freezePane(wb, name, firstRow = TRUE, firstCol = TRUE)
  setColWidths(wb, name, cols = 1:ncol(df), widths = 13)
  setColWidths(wb, name, cols = 1:2, widths = 28)
  if (!nrow(df)) return(invisible())
  rows <- 2:(nrow(df) + 1)
  for (cl in which(names(df) %in% p_like))
    conditionalFormatting(wb, name, cols = cl, rows = rows, rule = "<=0.05", style = st_sig)
  for (cl in which(names(df) %in% signed)) {
    conditionalFormatting(wb, name, cols = cl, rows = rows, rule = ">0", style = st_up)
    conditionalFormatting(wb, name, cols = cl, rows = rows, rule = "<0", style = st_down)
  }
  for (cl in which(names(df) %in% dir_text)) {
    conditionalFormatting(wb, name, cols = cl, rows = rows, type = "contains", rule = "positive", style = st_up)
    conditionalFormatting(wb, name, cols = cl, rows = rows, type = "contains", rule = "negative", style = st_down)
  }
}
pm <- grep("(WilcoxP_|WilcoxFDR_|^KW_.*[^g]$|_p_|_FDR_)", names(master), value = TRUE)
sm <- grep("(log2FC_|^FC_|_r_)", names(master), value = TRUE)
add_sheet("master_summary", master, pm, sm, c(cons_cols, "consistent_cor_list"))
add_sheet("paired_wilcox", paired_long, c("Wilcox_p", "Wilcox_FDR"), c("log2FC", "FC"))
add_sheet("kruskal_wallis", kw_tab, grep("^KW_.*[^g]$", names(kw_tab), value = TRUE))
add_sheet("treated_pairwise", trt_pairwise_long, c("Wilcox_p", "Wilcox_FDR"))
add_sheet("untreated_pairwise", untreated_pairwise_long, c("Wilcox_p", "Wilcox_FDR"))
add_sheet("responder_treated", responder_treated_long, c("Wilcox_p", "Wilcox_FDR"),
          "FoldChange_Responder_Nonresponder")
add_sheet("responder_untreated", responder_untreated_long, c("Wilcox_p", "Wilcox_FDR"),
          "FoldChange_Responder_Nonresponder")
add_sheet("responder_by_group", responder_group_long, c("Wilcox_p", "Wilcox_FDR"),
          "FoldChange_Responder_Nonresponder")
add_sheet("correlations", cor_long, c("pearson_p", "pearson_FDR", "spearman_p", "spearman_FDR"),
          c("pearson_r", "spearman_rho"), "consistent")
for (cv in dedicated_cor_vars) {
  add_sheet(paste0(cv, "_cor"), dedicated_cor_tables[[cv]],
            c("pearson_p", "pearson_FDR", "spearman_p", "spearman_FDR"),
            c("pearson_r", "spearman_rho"), "consistent")
}
add_sheet("consistent_correlations", cor_consistent,
          c("pearson_p", "pearson_FDR", "spearman_p", "spearman_FDR"),
          c("pearson_r", "spearman_rho"), c("consistent", "Direction"))
add_sheet("PH_invalid_samples", ph_invalid_samples)
legend_df <- data.frame(
  Item = c("Yellow fill", "Red text", "Blue text", "solid_*", "consistent_* / consistent",
           "FC", "log2FC", "Sig marks", "Pairing", "FDR", "Treatments", "Responder analysis",
           "Group responder test eligibility"),
  Meaning = c("p or FDR <= 0.05", "positive FC / correlation", "negative FC / correlation",
              "Pearson & Spearman p<=0.05 with same sign",
              paste("Pearson & Spearman p<=0.1 with same sign (rho/r); 'positive' or 'negative' +",
                    "stars from the larger p and larger FDR of the two methods"),
              "signed fold change mean(Treated)/mean(Untreated); 2x lower = -2",
              paste0("log2((meanT+eps)/(meanU+eps)), eps=", signif(eps, 3)),
              sig_legend,
              "Treated vs untreated side of same subject (Identity); Wilcoxon signed-rank",
              "Benjamini-Hochberg across taxa per test column",
              paste0(paste0(names(trt_name), "=", trt_name, collapse = ", "),
                     "; xN = untreated (NoTreatment) side of the same subjects"),
              paste0("Unpaired Wilcoxon on all treated or all untreated samples; fold change is ",
                     "(Responder mean + eps)/(Nonresponder mean + eps), eps=", signif(eps, 3),
                     "; ", responder_col, ": 1=Nonresponder, 2=Responder"),
              "Within-group tests run only when that group has >=3 responders and >=2 nonresponders; BH-FDR per group across features."))
add_sheet("legend", legend_df)
saveWorkbook(wb, rout("results.xlsx"), overwrite = TRUE)

## --------------------------------------------------------------------------
## Figures
## --------------------------------------------------------------------------
## Heatmaps show at most max_heat features (smallest p first) so that large
## feature sets (EC, GO, pathways) stay readable and within the PNG size limit.
max_heat <- 60
row_minp <- function(df, pcols) {
  pm <- as.matrix(as.data.frame(df)[, pcols, drop = FALSE])
  pm[is.na(pm)] <- Inf
  out <- apply(pm, 1, min); out[!is.finite(out)] <- NA; out
}
top_ids <- function(df, pcols, cutoff) {
  p <- row_minp(df, pcols)
  keep <- !is.na(p) & p <= cutoff
  if (!any(keep)) return(character())
  best <- tapply(p[keep], df$FeatureID[keep], min)
  head(names(sort(best)), max_heat)
}
n_hits <- function(df, pcols, cutoff) {
  p <- row_minp(df, pcols); length(unique(df$FeatureID[!is.na(p) & p <= cutoff]))
}
trunc_lab <- function(x, n = 60) ifelse(nchar(x) > n, paste0(substr(x, 1, n - 3), "..."), x)
cap_note <- function(n_all) if (n_all > max_heat) paste0(" (top ", max_heat, " of ", n_all, " by p)") else ""

## 1) Paired FC heatmap for taxa with Wilcoxon p<=0.1 in any cohort
hit_ids <- top_ids(paired_long, "Wilcox_p", 0.1)
n_all <- n_hits(paired_long, "Wilcox_p", 0.1)
if (length(hit_ids)) {
  hd <- paired_long %>% filter(FeatureID %in% hit_ids) %>%
    mutate(log2FC_c = pmax(pmin(log2FC, 3), -3),
           Comparison = factor(cmp_lab[Comparison], levels = cmp_lab))
  ord_g <- hd %>% group_by(Taxon) %>% summarise(m = mean(log2FC, na.rm = TRUE)) %>% arrange(m) %>% pull(Taxon)
  hd$Taxon <- factor(hd$Taxon, levels = ord_g)
  p <- ggplot(hd, aes(Comparison, Taxon, fill = log2FC_c)) + geom_tile(colour = "white") +
    geom_text(aes(label = Sig), size = 3, vjust = 0.75) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                         name = "log2FC\n(Trt/Untrt)", limits = c(-3, 3)) +
    scale_y_discrete(labels = trunc_lab) +
    labs(title = paste0("Paired Wilcoxon: Treated vs Untreated (taxa with p<=0.1)", cap_note(n_all)),
         subtitle = sig_legend, x = NULL, y = NULL) +
    theme_pub(10) + theme(legend.position = "right", axis.text.x = element_text(angle = 30, hjust = 1))
  ggsave(fout("paired_FC_heatmap.png"), p, width = 7, height = max(3, 0.22 * length(hit_ids) + 2.5),
         dpi = 300, bg = "white", limitsize = FALSE)
}

## 2) Volcano plots per cohort
vd <- paired_long %>% filter(!is.na(Wilcox_p)) %>%
  mutate(Comparison = factor(cmp_lab[Comparison], levels = cmp_lab))
vd$lab <- ifelse(vd$Wilcox_p < p_plot, paste0(vd$Taxon, vd$Sig), "")
p <- ggplot(vd, aes(log2FC, -log10(Wilcox_p))) +
  geom_point(aes(colour = ifelse(Wilcox_p < p_plot, ifelse(log2FC > 0, "Up", "Down"), "NS")),
             size = 1.6, alpha = 0.8) +
  geom_hline(yintercept = -log10(p_plot), linetype = 2, colour = "grey50") +
  geom_vline(xintercept = 0, colour = "grey70") +
  geom_text(aes(label = lab), size = 2.3, vjust = -0.6, check_overlap = TRUE) +
  scale_colour_manual(values = c(Up = "#B2182B", Down = "#2166AC", NS = "grey70"), name = NULL) +
  facet_wrap(~Comparison, nrow = 1) +
  labs(title = "Paired Wilcoxon volcano (Treated vs Untreated)", subtitle = sig_legend,
       x = "log2 FC", y = "-log10 p") +
  theme_pub(10)
ggsave(fout("paired_volcano.png"), p, width = 13, height = 4.5, dpi = 300, bg = "white")

## 3) Correlation heatmaps (Spearman and Pearson) for taxa with p<0.05 anywhere
for (m in c("spearman", "pearson")) {
  rcol <- if (m == "spearman") "spearman_rho" else "pearson_r"
  pcol <- paste0(m, "_p"); fcol <- paste0(m, "_FDR")
  ids <- top_ids(cor_long, pcol, p_plot)
  if (!length(ids)) next
  n_all <- n_hits(cor_long, pcol, p_plot)
  cd <- cor_long %>% filter(FeatureID %in% ids) %>%
    mutate(r = .data[[rcol]], Sig = sig_mark(.data[[pcol]], .data[[fcol]]),
           Subset = factor(Subset, levels = cor_sets),
           Variable = factor(Variable, levels = cor_vars))
  ord_g <- cd %>% group_by(Taxon) %>% summarise(m = mean(r, na.rm = TRUE)) %>% arrange(m) %>% pull(Taxon)
  cd$Taxon <- factor(cd$Taxon, levels = ord_g)
  p <- ggplot(cd, aes(Variable, Taxon, fill = r)) + geom_tile(colour = "white") +
    geom_text(aes(label = Sig), size = 2.6, vjust = 0.75) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                         limits = c(-1, 1), name = ifelse(m == "spearman", "rho", "r")) +
    facet_wrap(~Subset, nrow = 1) +
    scale_y_discrete(labels = trunc_lab) +
    labs(title = paste0(tools::toTitleCase(m), " correlation (taxa with p<=", p_plot, " in any test)",
                        cap_note(n_all)),
         subtitle = sig_legend, x = NULL, y = NULL) +
    theme_pub(9) + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")
  ggsave(fout(m, "_correlation_heatmap.png"), p, width = 12,
         height = max(4, 0.18 * length(ids) + 2.5), dpi = 300, bg = "white", limitsize = FALSE)
}

## 3b) Consistent correlation heatmap (Pearson & Spearman p<=0.1, same sign)
cons_ids <- cor_consistent %>% mutate(pmaxc = pmax(pearson_p, spearman_p)) %>%
  top_ids("pmaxc", 1)
n_all <- length(unique(cor_consistent$FeatureID))
if (length(cons_ids)) {
  cd <- cor_long %>% filter(FeatureID %in% cons_ids) %>%
    mutate(Stars = sub("^(positive|negative)", "", consistent),
           Consistent = consistent != "",
           Subset = factor(Subset, levels = cor_sets),
           Variable = factor(Variable, levels = cor_vars))
  ord_g <- cd %>% group_by(Taxon) %>% summarise(m = mean(spearman_rho, na.rm = TRUE)) %>%
    arrange(m) %>% pull(Taxon)
  cd$Taxon <- factor(cd$Taxon, levels = ord_g)
  p <- ggplot(cd, aes(Variable, Taxon)) +
    geom_tile(aes(fill = ifelse(Consistent, spearman_rho, NA)), colour = "white") +
    geom_text(aes(label = Stars), size = 2.6, vjust = 0.75) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                         limits = c(-1, 1), na.value = "grey95", name = "Spearman\nrho") +
    facet_wrap(~Subset, nrow = 1) +
    scale_y_discrete(labels = trunc_lab) +
    labs(title = paste0("Consistent correlations (Pearson & Spearman p<=0.1, same direction)", cap_note(n_all)),
         subtitle = paste0(sig_legend, "  (stars from the weaker of the two methods)"),
         x = NULL, y = NULL) +
    theme_pub(9) + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")
  ggsave(fout("consistent_correlation_heatmap.png"), p, width = 12,
         height = max(4, 0.18 * length(cons_ids) + 2.5), dpi = 300, bg = "white", limitsize = FALSE)
}

## 4) KW boxplots of the top taxa (8 groups)
kw_top <- kw_tab %>% filter(!is.na(KW_8groups), KW_8groups < p_plot) %>% head(24) %>% pull(FeatureID)
if (length(kw_top)) {
  bd <- as.data.frame(t(mat[kw_top, , drop = FALSE]), check.names = FALSE) %>%
    mutate(Group = meta$Group) %>% pivot_longer(-Group, names_to = "FeatureID", values_to = "Abundance") %>%
    left_join(kw_tab %>% select(FeatureID, Taxon, KW_8groups, KW_8groups_Sig), by = "FeatureID") %>%
    mutate(Panel = paste0(Taxon, "\nKW p=", signif(KW_8groups, 2), " ", KW_8groups_Sig),
           Panel = factor(Panel, levels = unique(Panel[order(KW_8groups)])))
  p <- ggplot(bd, aes(Group, Abundance, fill = Group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.8) +
    geom_jitter(width = 0.15, size = 0.6, alpha = 0.6) +
    scale_fill_manual(values = pal_grp) + scale_y_continuous(trans = "sqrt") +
    scale_x_discrete(labels = grp_lab) +
    facet_wrap(~Panel, scales = "free_y", ncol = 6) +
    labs(title = paste0("Kruskal-Wallis across 8 groups: top taxa (p<", p_plot, ")"),
         subtitle = paste0(sig_legend, "   |   ", paste0(names(trt_name), "=", trt_name, collapse = ", ")),
         y = "Relative abundance (sqrt scale)", x = NULL) +
    theme_pub(8) + theme(legend.position = "none",
                         axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6))
  ggsave(fout("KW_top_boxplots.png"), p, width = 15, height = 3 * ceiling(length(kw_top) / 6) + 1,
         dpi = 300, bg = "white", limitsize = FALSE)
}

## Separate responder plots grouped by treated, untreated, or individual group
## are intentionally disabled. Only the combined by-taxon figures below are generated.
##
## 5) One composite responder figure per taxon with any p<=0.01 result.
##    Each panel shows a pooled treated/untreated or eligible individual-group
##    comparison for that taxon, regardless of the individual panel p-value.
responder_plot_tests <- bind_rows(
  responder_treated_long %>%
    transmute(FeatureID, Taxon, Analysis = "All treated", Responder_n, Responder_Mean,
              Responder_SE, Nonresponder_n, Nonresponder_Mean, Nonresponder_SE,
              Wilcox_p, Wilcox_FDR, Sig),
  responder_untreated_long %>%
    transmute(FeatureID, Taxon, Analysis = "All untreated", Responder_n, Responder_Mean,
              Responder_SE, Nonresponder_n, Nonresponder_Mean, Nonresponder_SE,
              Wilcox_p, Wilcox_FDR, Sig),
  responder_group_long %>%
    filter(TestEligible) %>%
    transmute(FeatureID, Taxon,
              Analysis = paste0("Group ", Group, " (",
                                trt_name[substr(Group, 1, 1)],
                                ifelse(grepl("N$", Group), ", untreated)", ")")),
              Responder_n,
              Responder_Mean, Responder_SE, Nonresponder_n, Nonresponder_Mean,
              Nonresponder_SE, Wilcox_p, Wilcox_FDR, Sig)
)
responder_trigger_ids <- responder_plot_tests %>%
  filter(!is.na(Wilcox_p), Wilcox_p <= 0.01) %>%
  distinct(FeatureID) %>%
  pull(FeatureID)
responder_test_order <- c(
  "All treated", "All untreated",
  paste0("Group ", groups8, " (",
         trt_name[substr(groups8, 1, 1)],
         ifelse(grepl("N$", groups8), ", untreated)", ")"))
)
for (id in responder_trigger_ids) {
  plot_tab <- responder_plot_tests %>%
    filter(FeatureID == id) %>%
    mutate(Analysis = factor(Analysis, levels = responder_test_order),
           Panel = paste0(as.character(Analysis), "\nResponder n=", Responder_n,
                          ", Nonresponder n=", Nonresponder_n),
           Panel = factor(Panel, levels = unique(Panel)))
  if (!nrow(plot_tab)) next
  bar_df <- bind_rows(
    plot_tab %>% transmute(FeatureID, Taxon, Analysis, Panel, Category = "Nonresponder",
                           Mean = Nonresponder_Mean, SE = Nonresponder_SE),
    plot_tab %>% transmute(FeatureID, Taxon, Analysis, Panel, Category = "Responder",
                           Mean = Responder_Mean, SE = Responder_SE)
  ) %>%
    mutate(Category = factor(Category, levels = c("Nonresponder", "Responder")),
           SE = ifelse(is.na(SE), 0, SE))
  stopifnot(nrow(bar_df) == 2 * nrow(plot_tab),
            all(table(bar_df$Panel) == 2))
  ann <- plot_tab %>%
    mutate(y = pmax(Responder_Mean + Responder_SE, Nonresponder_Mean + Nonresponder_SE,
                    na.rm = TRUE),
           y = ifelse(is.finite(y) & y > 0, y * 1.2, 1),
           label = paste0("p=", format.pval(Wilcox_p, digits = 2, eps = 0.001),
                          "\nFDR=", ifelse(is.na(Wilcox_FDR), "NA",
                                            format.pval(Wilcox_FDR, digits = 2, eps = 0.001)),
                          " ", Sig))
  p <- ggplot(bar_df, aes(as.numeric(Category), Mean, fill = Category)) +
    geom_col(width = 0.72, colour = "grey25", alpha = 0.9) +
    geom_errorbar(aes(ymin = pmax(Mean - SE, 0), ymax = Mean + SE), width = 0.22) +
    geom_text(data = ann, aes(x = 1.5, y = y, label = label),
              inherit.aes = FALSE, size = 2.8, lineheight = 0.95) +
    facet_wrap(~Panel, ncol = 3, scales = "free_y") +
    scale_x_continuous(breaks = 1:2, labels = levels(bar_df$Category)) +
    scale_fill_manual(values = c(Nonresponder = "#C44E52", Responder = "#4C72B0")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.32))) +
    labs(title = paste0(clean_name(unique(plot_tab$Taxon)), ": responder vs nonresponder"),
         subtitle = paste0("All eligible comparisons shown; figure trigger: any p<=0.01; ",
                           "unpaired Wilcoxon; bars = mean +/- SE; ", sig_legend),
         x = NULL, y = "Relative abundance") +
    theme_pub(9) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 25, hjust = 1, size = 8),
          strip.text = element_text(size = 8))
  ggsave(file.path(responder_taxon_dir,
                   paste0(safe_file(clean_name(unique(plot_tab$Taxon))),
                          ".responder_composite.png")),
         p, width = 12, height = max(5, 2.8 * ceiling(nrow(plot_tab) / 3)),
         dpi = 300, bg = "white", limitsize = FALSE)
}
cat("Composite responder figures (any test p<=0.01): ", length(responder_trigger_ids), "\n", sep = "")

## 6) Individual paired treatment-comparison plots are intentionally disabled.
## To restore them, change FALSE to TRUE.
if (FALSE) {
  paired_hits <- paired_long %>% filter(Wilcox_p < p_plot) %>% pull(FeatureID) %>% unique()
  for (id in paired_hits) {
    pdl <- bind_rows(lapply(names(pairs), function(co) {
      pi <- pair_index[[co]]
      tibble(Identity = rep(pi$Identity, 2),
             Side = rep(c("Treated", "Untreated"), each = nrow(pi)),
             Group = rep(c(paste0(co, co), paste0(co, "N")), each = nrow(pi)),
             Comparison = pairs[[co]],
             Abundance = c(mat[id, pi$t_sample], mat[id, pi$u_sample]))
    }))
    labs_df <- paired_long %>% filter(FeatureID == id) %>%
      mutate(lab = paste0("p=", signif(Wilcox_p, 2), " FDR=", signif(Wilcox_FDR, 2), " ", Sig,
                          "\nFC=", round(FC, 2)))
    pdl <- pdl %>% left_join(labs_df %>% select(Comparison, lab), by = "Comparison") %>%
      mutate(Facet = paste0(cmp_lab[Comparison], "\n", lab),
             Facet = factor(Facet, levels = unique(Facet)))
    p <- ggplot(pdl, aes(Side, Abundance)) +
      geom_boxplot(aes(fill = Group), outlier.shape = NA, alpha = 0.7) +
      geom_line(aes(group = Identity), colour = "grey55", alpha = 0.6) +
      geom_point(aes(colour = Group), size = 1.4) +
      scale_fill_manual(values = pal_grp) + scale_colour_manual(values = pal_grp) +
      facet_wrap(~Facet, nrow = 1) +
      labs(title = paste0(clean_name(id), ": paired Treated vs Untreated"), subtitle = sig_legend,
           x = NULL, y = "Relative abundance") +
      theme_pub(9) + theme(legend.position = "none")
    ggsave(file.path(ind_trt_dir, paste0(safe_file(clean_name(id)), ".paired_TrtVsUntrt.png")),
           p, width = 9, height = 4, dpi = 200, bg = "white")
  }
}

## 7) Individual scatter plots for consistent correlations. PH, Irritation,
##    ss10_Itch, ss10_Redness, and ss10_Pain plots are also generated when
##    either Pearson or Spearman p<=0.1.
screened_plot_hits <- cor_long %>%
  filter(Variable %in% screened_cor_vars) %>%
  filter((!is.na(pearson_p) & pearson_p <= 0.1) |
           (!is.na(spearman_p) & spearman_p <= 0.1)) %>%
  distinct(FeatureID, Variable)
solid_hits <- bind_rows(
  cor_long %>% filter(consistent != "") %>% distinct(FeatureID, Variable),
  screened_plot_hits
) %>% distinct(FeatureID, Variable)
for (k in seq_len(nrow(solid_hits))) {
  id <- solid_hits$FeatureID[k]; cv <- solid_hits$Variable[k]
  sd0 <- tibble(Sample = meta$Sample, Abundance = mat[id, ], Value = meta[[cv]],
                Group = meta$Group, Status = meta$Status) %>%
    filter(is.finite(Value), is.finite(Abundance))
  sdl <- bind_rows(sd0 %>% mutate(Subset = "All"),
                   sd0 %>% filter(Status == "Treated") %>% mutate(Subset = "Treated"),
                   sd0 %>% filter(Status == "Untreated") %>% mutate(Subset = "Untreated"))
  lab <- cor_long %>% filter(FeatureID == id, Variable == cv) %>%
    mutate(Facet = paste0(Subset, ifelse(consistent != "", paste0(" [", consistent, "]"), ""),
                          "\nr=", round(pearson_r, 2), " p=", signif(pearson_p, 2), " ", pearson_Sig,
                          "; rho=", round(spearman_rho, 2), " p=", signif(spearman_p, 2), " ",
                          spearman_Sig))
  sdl <- sdl %>% left_join(lab %>% select(Subset, Facet), by = "Subset") %>%
    mutate(Facet = factor(Facet, levels = lab$Facet[match(cor_sets, lab$Subset)]))
  p <- ggplot(sdl, aes(Value, Abundance)) +
    geom_jitter(aes(colour = Group), width = 0.03 * diff(range(sd0$Value, na.rm = TRUE)),
                height = 0, size = 1.6, alpha = 0.8) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = "black", linewidth = 0.6) +
    scale_colour_manual(values = pal_grp, labels = grp_lab1, name = NULL) +
    facet_wrap(~Facet, nrow = 1) +
    labs(title = paste0(clean_name(id), " vs ", cv),
         subtitle = if (cv == "PH.average")
           paste0("PH values outside [0,14] excluded (n=", nrow(ph_invalid_samples), "); ",
                  sig_legend)
         else sig_legend,
         x = cv, y = "Relative abundance") +
    theme_pub(9)
  ggsave(file.path(ind_cor_dir, paste0(safe_file(clean_name(id)), ".cor_", safe_file(cv), ".png")),
         p, width = 11, height = 4.8, dpi = 200, bg = "white")
}

## 8) Grouped bar graphs (mean +/- SE) of all 8 groups with pairwise p-value
##    brackets: the 4 treated-vs-untreated pairs (AA:AN, BB:BN, CC:CN, DD:DN)
##    and the 6 between-treatment pairs (AA,BB,CC,DD). One plot per feature that
##    is significant (p<p_plot) in any of those comparisons or in KW_TreatedOnly
##    / KW_UntreatedOnly.
xpos <- setNames(seq_along(groups8), groups8)
paired_bracket <- lapply(names(pairs), function(co)
  list(g1 = paste0(co, co), g2 = paste0(co, "N"), col = paste0("WilcoxP_", co, co, "vs", co, "N")))
treat_bracket <- lapply(trt_pairs, function(pr)
  list(g1 = pr[1], g2 = pr[2], col = paste0("WilcoxTrt_", pr[1], "vs", pr[2])))
untreat_bracket <- lapply(untrt_pairs, function(pr)
  list(g1 = pr[1], g2 = pr[2], col = paste0("WilcoxUntreat_", pr[1], "vs", pr[2])))
all_brackets <- c(paired_bracket, treat_bracket, untreat_bracket)
bar_hit_cols <- c(unlist(lapply(all_brackets, `[[`, "col")), "KW_TreatedOnly", "KW_UntreatedOnly")
bar_hits <- master$FeatureID[apply(
  sapply(bar_hit_cols, function(cc) master[[cc]] < p_plot), 1, any, na.rm = TRUE)]
bar_hits <- head(bar_hits[!is.na(bar_hits)], 60)
for (id in bar_hits) {
  mr <- master[master$FeatureID == id, ]
  bd <- data.frame(Group = factor(groups8, levels = groups8),
                   x = xpos[groups8],
                   Mean = sapply(groups8, function(g) mr[[paste0("Mean_", g)]]),
                   SE   = sapply(groups8, function(g) mr[[paste0("SE_", g)]]),
                   Side = ifelse(grepl("N$", groups8), "Untreated", "Treated"))
  grp_plot_labels <- setNames(vapply(groups8, function(g) {
    in_group <- meta$Group == g & !is.na(meta$ResponderStatus)
    paste0(grp_lab[g], "\nR n=", sum(in_group & meta$ResponderStatus == "Responder"),
           ", NR n=", sum(in_group & meta$ResponderStatus == "Nonresponder"))
  }, character(1)), groups8)
  bd$SE[is.na(bd$SE)] <- 0
  ytop <- max(bd$Mean + bd$SE, na.rm = TRUE); if (!is.finite(ytop) || ytop <= 0) ytop <- 1
  ## build bracket data, shortest span lowest, increasing heights
  brk <- lapply(all_brackets, function(b) {
    p <- mr[[b$col]]
    if (is.na(p) || p > 0.1) return(NULL)
    data.frame(x1 = xpos[b$g1], x2 = xpos[b$g2], p = p,
               lab = paste0("p=", signif(p, 2), " ", sig_mark(p, NA)))
  })
  brk <- do.call(rbind, brk[!sapply(brk, is.null)])
  p <- ggplot(bd, aes(x, Mean, fill = Group)) +
    geom_col(width = 0.72, colour = "grey25", alpha = 0.9) +
    geom_errorbar(aes(ymin = pmax(Mean - SE, 0), ymax = Mean + SE), width = 0.25) +
    scale_fill_manual(values = pal_grp) +
    scale_x_continuous(breaks = xpos, labels = grp_plot_labels[groups8]) +
    labs(title = paste0(clean_name(id), ": group means +/- SE"),
         subtitle = paste0("KW treated(AA,BB,CC,DD) p=", signif(mr$KW_TreatedOnly, 2), " ",
                           sig_mark(mr$KW_TreatedOnly, NA),
                           "  |  KW untreated(AN,BN,CN,DN) p=", signif(mr$KW_UntreatedOnly, 2), " ",
                           sig_mark(mr$KW_UntreatedOnly, NA), "\n", sig_legend, "   |   ",
                           paste0(names(trt_name), "=", trt_name, collapse = ", ")),
         x = NULL, y = "Relative abundance (mean +/- SE)") +
    theme_pub(9) + theme(legend.position = "none",
                         axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
  if (!is.null(brk) && nrow(brk)) {
    brk <- brk[order(abs(brk$x2 - brk$x1), brk$x1), ]
    step <- 0.09 * ytop
    brk$y <- ytop + step * seq_len(nrow(brk))
    p <- p +
      geom_segment(data = brk, aes(x = x1, xend = x2, y = y, yend = y),
                   inherit.aes = FALSE, linewidth = 0.4) +
      geom_text(data = brk, aes(x = (x1 + x2) / 2, y = y + step * 0.25, label = lab),
                inherit.aes = FALSE, size = 2.4) +
      coord_cartesian(ylim = c(0, max(brk$y) + step), clip = "off")
  }
  ggsave(file.path(ind_bar_dir, paste0(safe_file(clean_name(id)), ".group_barplot.png")),
         p, width = 8, height = 6, dpi = 200, bg = "white")
}

## --------------------------------------------------------------------------
## Console summary (captured into summary)
## --------------------------------------------------------------------------
cat("\n== Paired Wilcoxon hits (p<0.05 / FDR<0.1) ==\n")
print(paired_long %>% group_by(Comparison) %>%
        summarise(p05 = sum(Wilcox_p < 0.05, na.rm = TRUE), FDR10 = sum(Wilcox_FDR < 0.1, na.rm = TRUE),
                  up_p05 = sum(Wilcox_p < 0.05 & log2FC > 0, na.rm = TRUE),
                  down_p05 = sum(Wilcox_p < 0.05 & log2FC < 0, na.rm = TRUE)))
cat("\n== Kruskal-Wallis hits (p<0.05 / FDR<0.1) ==\n")
print(sapply(grep("^KW_[A-Za-z0-9]+$", names(kw_tab), value = TRUE), function(k)
  c(p05 = sum(kw_tab[[k]] < 0.05, na.rm = TRUE),
    FDR10 = sum(kw_tab[[paste0(k, "_FDR")]] < 0.1, na.rm = TRUE))))
cat("\n== Correlation hits: Spearman p<0.05 / FDR<0.1 / solid / consistent ==\n")
print(as.data.frame(cor_long %>% group_by(Subset, Variable) %>%
        summarise(sp_p05 = sum(spearman_p < 0.05, na.rm = TRUE),
                  sp_FDR10 = sum(spearman_FDR < 0.1, na.rm = TRUE),
                  pe_p05 = sum(pearson_p < 0.05, na.rm = TRUE),
                  solid = sum(solid != ""),
                  consistent_pos = sum(grepl("^positive", consistent)),
                  consistent_neg = sum(grepl("^negative", consistent)), .groups = "drop")))
cat("\n== Responder vs nonresponder unpaired Wilcoxon hits ==\n")
print(bind_rows(
  responder_treated_long %>% summarise(Comparison = "Treated",
    p_lt_0_05 = sum(Wilcox_p < 0.05, na.rm = TRUE),
    FDR_lt_0_1 = sum(Wilcox_FDR < 0.1, na.rm = TRUE),
    p_le_0_1_plot = sum(Wilcox_p <= 0.1, na.rm = TRUE)),
  responder_untreated_long %>% summarise(Comparison = "Untreated",
    p_lt_0_05 = sum(Wilcox_p < 0.05, na.rm = TRUE),
    FDR_lt_0_1 = sum(Wilcox_FDR < 0.1, na.rm = TRUE),
    p_le_0_1_plot = sum(Wilcox_p <= 0.1, na.rm = TRUE))
))
cat("\n== Responder vs nonresponder within-group tests ==\n")
print(responder_group_long %>%
        group_by(Group) %>%
        summarise(Responder_n = first(Responder_n),
                  Nonresponder_n = first(Nonresponder_n),
                  Eligible = first(TestEligible),
                  p_lt_0_05 = sum(Wilcox_p < 0.05, na.rm = TRUE),
                  FDR_lt_0_1 = sum(Wilcox_FDR < 0.1, na.rm = TRUE),
                  p_le_0_1_plot = sum(Wilcox_p <= 0.1, na.rm = TRUE),
                  .groups = "drop"))
cat("\nDone. All output in", normalizePath(run_dir), "\nSummary log:", normalizePath(log_file), "\n")
sink()
