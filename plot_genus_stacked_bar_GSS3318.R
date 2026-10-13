###############################################################################
## GSS3318 SensiScope — stacked bar plot of taxon relative abundance
## -----------------------------------------------------------------------------
## Taxa whose maximum abundance in any sample is > min_max_pct (default 1%)
## are shown individually; all others are pooled into "Other".
## The taxonomic rank used in the title and legend is taken from the input file
## name: genus_cleandata -> "Genus"/"genera", species_cleandata -> "Species".
## Bars show the mean relative abundance (%) per group:
##   All, Sensitive, NotSensitive, Leg, Back, Hand,
##   Leg/Back/Hand x Sensitive/NotSensitive
##
## Outputs (in outdir):
##   *.stacked_bar.png   — static, publication-quality (300 dpi)
##   *.stacked_bar.pdf   — vector version
##   *.stacked_bar.html  — interactive (hover, zoom, legend toggling);
##                         camera icon in toolbar downloads a hi-res PNG
##   *.stacked_bar_data.tsv — plotted values
##
## Usage
##   Rscript plot_genus_stacked_bar_GSS3318.R [genus_file] [min_max_pct] [outdir]
## Example
##   Rscript plot_genus_stacked_bar_GSS3318.R genus_cleandata 1 stacked_bar_result
###############################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(plotly)
  library(htmlwidgets)
  library(RColorBrewer)
})

args <- commandArgs(trailingOnly = TRUE)

find_file <- function(path) {
  candidates <- unique(c(path, file.path("..", path)))
  hit <- candidates[file.exists(candidates)]
  if (length(hit) == 0) stop("Cannot find file: ", path, call. = FALSE)
  hit[1]
}

genus_file  <- find_file(if (length(args) >= 1) args[1] else "genus_cleandata")
min_max_pct <- if (length(args) >= 2) as.numeric(args[2]) else 1
outdir      <- if (length(args) >= 3) args[3] else
  paste0(basename(genus_file), "_stacked_bar_", format(Sys.time(), "%Y%m%d_%H%M%S"))
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
prefix <- file.path(outdir, paste0(basename(genus_file), ".stacked_bar"))

## ---------------------------------------------------------------------------
## Taxonomic rank taken from the input file name, e.g. genus_cleandata -> genus,
## species_cleandata -> species.  Used for all figure titles and legend labels.
## ---------------------------------------------------------------------------
rank_raw <- tolower(gsub("[_.-]?clean.*$", "",
                         tools::file_path_sans_ext(basename(genus_file))))
rank_raw <- gsub("[^a-z]", "", rank_raw)
if (rank_raw == "") rank_raw <- "taxon"
plural_map <- c(genus = "genera", species = "species", phylum = "phyla",
                family = "families", order = "orders", class = "classes",
                kingdom = "kingdoms", domain = "domains", taxon = "taxa",
                otu = "OTUs", asv = "ASVs")
rank_plural <- if (rank_raw %in% names(plural_map)) plural_map[[rank_raw]] else paste0(rank_raw, "s")
title_case <- function(x) paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
## "OTU"/"ASV" are already correctly cased; only capitalize normal words.
rank_sing_lab <- if (rank_raw %in% c("otu", "asv")) toupper(rank_raw) else title_case(rank_raw)
rank_plur_lab <- rank_plural

## ---------------------------------------------------------------------------
## Load data (rows = genera, columns = BodySite.Sensitivity.SID, values in %)
## ---------------------------------------------------------------------------
gw <- read.table(genus_file, header = TRUE, row.names = 1, sep = "\t",
                 quote = "", comment.char = "", check.names = FALSE)
gw[] <- lapply(gw, function(x) { z <- suppressWarnings(as.numeric(x)); z[is.na(z)] <- 0; z })

## Rescale to percent if table is in proportions (column sums ~ 1).
if (median(colSums(gw)) <= 1.5) gw <- gw * 100

clean_name <- function(x) gsub("_", " ", gsub("^[a-z]__", "", x))

keep <- apply(gw, 1, max) > min_max_pct
if (!any(keep)) stop("No genus has max abundance > ", min_max_pct, "%", call. = FALSE)

long <- gw %>%
  rownames_to_column("genus_id") %>%
  mutate(genus = ifelse(keep, clean_name(genus_id), "Other")) %>%
  select(-genus_id) %>%
  pivot_longer(-genus, names_to = "sample_id", values_to = "abundance") %>%
  group_by(genus, sample_id) %>%
  summarise(abundance = sum(abundance), .groups = "drop") %>%
  tidyr::separate(sample_id, c("BodySite", "Sensitivity", "SID"),
                  sep = "[.]", remove = FALSE, extra = "merge")

## ---------------------------------------------------------------------------
## Groups
## ---------------------------------------------------------------------------
sites <- c("Leg", "Back", "Hand")
sens  <- c("Sensitive", "NotSensitive")
sens_lab <- c(Sensitive = "Sensitive", NotSensitive = "Not sensitive")

group_defs <- c(
  list(list(label = "All samples", f = function(d) rep(TRUE, nrow(d)))),
  lapply(sens, function(s) list(label = paste("All", sens_lab[[s]]),
                                f = function(d) d$Sensitivity == s)),
  lapply(sites, function(b) list(label = paste("All", b),
                                 f = function(d) d$BodySite == b)),
  unlist(lapply(sites, function(b) lapply(sens, function(s)
    list(label = paste(b, sens_lab[[s]]),
         f = function(d) d$BodySite == b & d$Sensitivity == s))),
    recursive = FALSE)
)
group_levels <- vapply(group_defs, `[[`, character(1), "label")

n_samples <- n_distinct(long$sample_id)
bar_df <- lapply(group_defs, function(g) {
  d <- long[g$f(long), ]
  ns <- n_distinct(d$sample_id)
  d %>% group_by(genus) %>%
    summarise(mean_pct = sum(abundance) / ns, .groups = "drop") %>%
    mutate(group = g$label, n = ns)
}) %>% bind_rows()

## Genus order: by overall mean abundance; "Other" on top of the stack.
genus_order <- bar_df %>% filter(group == "All samples", genus != "Other") %>%
  arrange(desc(mean_pct)) %>% pull(genus)
if ("Other" %in% bar_df$genus) genus_order <- c(genus_order, "Other")

n_lab <- bar_df %>% distinct(group, n) %>%
  mutate(group = factor(group, levels = group_levels)) %>% arrange(group)
x_labels <- setNames(paste0(n_lab$group, "\n(n=", n_lab$n, ")"), as.character(n_lab$group))

bar_df <- bar_df %>%
  mutate(genus = factor(genus, levels = genus_order),
         group = factor(group, levels = group_levels)) %>%
  arrange(group, genus)

write.table(bar_df %>% select(group, n, genus, mean_pct),
            paste0(prefix, "_data.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

## ---------------------------------------------------------------------------
## Palette: distinct qualitative colours, grey for "Other"
## ---------------------------------------------------------------------------
n_gen <- sum(genus_order != "Other")
base_pal <- c(brewer.pal(12, "Paired"), brewer.pal(8, "Dark2"), brewer.pal(9, "Set1"),
              brewer.pal(8, "Set2"), brewer.pal(12, "Set3"), brewer.pal(8, "Accent"),
              brewer.pal(9, "Pastel1"), brewer.pal(8, "Pastel2"))
base_pal <- unique(base_pal[!toupper(base_pal) %in% c("#999999", "#B3B3B3", "#CCCCCC", "#D9D9D9", "#F2F2F2", "#666666")])
if (n_gen > length(base_pal)) {
  base_pal <- c(base_pal, hcl.colors(n_gen - length(base_pal), "Dynamic"))
}
pal <- setNames(base_pal[seq_len(n_gen)], genus_order[genus_order != "Other"])
if ("Other" %in% genus_order) pal <- c(pal, Other = "grey75")

## ---------------------------------------------------------------------------
## Static figure (PNG + PDF)
## ---------------------------------------------------------------------------
legend_cols <- ceiling(length(genus_order) / 30)
p <- ggplot(bar_df, aes(x = group, y = mean_pct, fill = genus)) +
  geom_col(width = 0.8, colour = "white", linewidth = 0.1,
           position = position_stack(reverse = TRUE)) +
  geom_vline(xintercept = c(1.5, 3.5, 6.5), linetype = "dashed",
             colour = "grey50", linewidth = 0.3) +
  scale_fill_manual(values = pal, name = rank_sing_lab) +
  scale_x_discrete(labels = x_labels) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.01)),
                     labels = function(x) paste0(x, "%")) +
  labs(x = NULL, y = "Mean relative abundance",
       title = paste0(rank_sing_lab, " composition by group"),
       subtitle = paste0(title_case(rank_plur_lab), " with max abundance > ", min_max_pct,
                         "% in any sample (", n_gen, " ", rank_plur_lab, "); others pooled as 'Other'")) +
  guides(fill = guide_legend(ncol = legend_cols, keywidth = unit(0.35, "cm"),
                             keyheight = unit(0.35, "cm"))) +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
        axis.text.y = element_text(colour = "black"),
        legend.text = element_text(size = 7, face = "italic"),
        legend.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold"),
        legend.position = "right")

fig_w <- 7 + 1.9 * legend_cols
ggsave(paste0(prefix, ".png"), p, width = fig_w, height = 8, dpi = 300, bg = "white")
ggsave(paste0(prefix, ".pdf"), p, width = fig_w, height = 8, bg = "white")

## ---------------------------------------------------------------------------
## Interactive figure (HTML, with PNG download button)
## ---------------------------------------------------------------------------
x_html <- gsub("\n", "<br>", x_labels)
ip <- plot_ly()
for (g in genus_order) {
  d <- bar_df %>% filter(genus == g)
  ip <- ip %>% add_trace(
    data = d, x = ~x_html[as.character(group)], y = ~mean_pct, type = "bar",
    name = g, marker = list(color = pal[[g]], line = list(color = "white", width = 0.3)),
    text = ~sprintf("<b><i>%s</i></b><br>%s (n=%d)<br>Mean: %.2f%%", genus, group, n, mean_pct),
    hoverinfo = "text", textposition = "none")
}
ip <- ip %>% layout(
  barmode = "stack",
  title = list(text = paste0("<b>", rank_sing_lab, " composition by group</b><br><sup>",
                             title_case(rank_plur_lab), " with max abundance > ",
                             min_max_pct, "% in any sample (", n_gen, " ",
                             rank_plur_lab, "); others pooled as 'Other'</sup>")),
  xaxis = list(title = "", categoryorder = "array",
               categoryarray = unname(x_html[group_levels]), tickangle = -45),
  yaxis = list(title = "Mean relative abundance (%)", range = c(0, 100.5),
               ticksuffix = "%"),
  legend = list(title = list(text = paste0("<b>", rank_sing_lab, "</b>")), traceorder = "normal",
                font = list(size = 10)),
  font = list(family = "Arial", size = 12),
  margin = list(t = 90, b = 140),
  plot_bgcolor = "white", paper_bgcolor = "white"
) %>% config(
  displaylogo = FALSE,
  toImageButtonOptions = list(format = "png", filename = basename(prefix),
                              width = 1600, height = 1000, scale = 3)
)
saveWidget(ip, paste0(prefix, ".html"), selfcontained = TRUE,
           title = paste0("GSS3318 ", rank_raw, " stacked bar"))
unlink(paste0(prefix, "_files"), recursive = TRUE)

cat(title_case(rank_plur_lab), "shown individually (max >", min_max_pct, "%):", n_gen, "\n")
cat("Output folder:", normalizePath(outdir), "\n")
