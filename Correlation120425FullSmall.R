###################################################
#library(proxy)
#myDist <- dist(W_matrix[c(i, row_index), ], method = "euclidean")
#system("git config --global user.name PingPnG")
#system("git config --global user.email 'hu.p@pg.com'")
###First try with correlation network layout show significant
####################################################
rm(list=ls())
library(ggplot2)
library(ggrepel)
library(corrplot)
library(officer)
library(rvg)
library(dplyr)
PCUTOFF=0.05
# Create PowerPoint to fill in important result


args <- commandArgs(trailingOnly = TRUE)
print(args)
file <- args[1]

rm(args)

#file="SmallSet_ArmCheek67.txt"
A <- read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
d <- dim(A)

heatmap_data <- A[, 2:d[2]]
row_names <- A[, 1]
row.names(heatmap_data) <- A[, 1]
data_matrix <- as.matrix(t(heatmap_data))

cor_matrix_spearman <- cor(data_matrix, method = "spearman", use = "pairwise.complete.obs")
testRes <- cor.mtest(as.matrix(data_matrix), method = "spearman", conf.level = 0.95)
p.mat = testRes$p

library(ggplot2)
library(ggpubr)
library(rstatix)
library(dplyr)
library(tidyverse)
library(openxlsx) 
library(circlize)
library(ComplexHeatmap)
library(lattice)   
library(reshape2)
create_pk_mat <- function(pk, adjpk = NULL) {
  # Coerce to numeric matrices
  pk    <- as.matrix(pk)
  storage.mode(pk) <- "double"
  
  if (!is.null(adjpk)) {
    adjpk <- as.matrix(adjpk)
    storage.mode(adjpk) <- "double"
  }
  
  # Ensure pk has row/col names (required for alignment)
  if (is.null(rownames(pk))) rownames(pk) <- seq_len(nrow(pk))
  if (is.null(colnames(pk))) colnames(pk) <- seq_len(ncol(pk))
  
  # Align adjpk to pk's shape and ordering (by names)
  if (!is.null(adjpk)) {
    # If adjpk lacks names, fall back to positional but warn
    if (is.null(rownames(adjpk)) || is.null(colnames(adjpk))) {
      warning("adjpk has no row and/or column names; aligning by position.")
      # recycle/trim to pk dims by position
      adjpk <- adjpk[seq_len(min(nrow(adjpk), nrow(pk))),
                     seq_len(min(ncol(adjpk), ncol(pk))), drop = FALSE]
      # Expand to pk dims with NA if needed
      tmp <- matrix(NA_real_, nrow = nrow(pk), ncol = ncol(pk),
                    dimnames = list(rownames(pk), colnames(pk)))
      tmp[seq_len(nrow(adjpk)), seq_len(ncol(adjpk))] <- adjpk
      adjpk <- tmp
    } else {
      # Proper name-based realignment (introduces NA if missing)
      adjpk <- adjpk[rownames(pk), colnames(pk), drop = FALSE]
    }
  }
  
  # Initialize output with pk's dims and names
  pk_mat <- matrix("",
                   nrow = nrow(pk), ncol = ncol(pk),
                   dimnames = list(rownames(pk), colnames(pk)))
  
  # Apply thresholds with clear precedence:
  # pk p<=0.10 -> "*", pk p<=0.05 -> "**",
  # adjpk q<=0.10 -> "***", adjpk q<=0.05 -> "****" (highest precedence)
  idx <- !is.na(pk) & pk <= 0.05
  pk_mat[idx] <- "*"
  
  #idx <- !is.na(pk) & pk <= 0.05
  #pk_mat[idx] <- "**"
  
  if (!is.null(adjpk)) {
    idx <- !is.na(adjpk) & adjpk <= 0.10
    pk_mat[idx] <- "**"
    
    idx <- !is.na(adjpk) & adjpk <= 0.05
    pk_mat[idx] <- "***"
  }
  
  pk_mat
}

pmMat<-create_pk_mat(p.mat)
# Example: cor_matrix_spearman = numeric correlation matrix
#          pmMat = "*" where p <= 0.05, "" otherwise

# Round correlations to 2 digits
cor_2d <- round(cor_matrix_spearman, 2)

# Convert to character
cor_char <- formatC(cor_2d, format = "f", digits = 2)

# Paste correlation with significance mark
combinedMat <- matrix(
  paste0(cor_char, pmMat),
  nrow = nrow(cor_matrix_spearman),
  ncol = ncol(cor_matrix_spearman),
  dimnames = dimnames(cor_matrix_spearman)
)

# Check result
print(combinedMat)



# row_groups <- c(rep("Rash", 2),   # rows 1–2
#                 rep("PH", 3),   # rows 3–5
#                 rep("DNA", 3),
#                 rep("HuDNA", 2),
#                 rep("MicDNA", 2),
#                 rep("Shannon", 2),
#                 rep("Species", 2)) 

library(ComplexHeatmap)
library(circlize)
library(grid)

row_groups <- c(rep("G1", 2),   # rows 1–2
                rep("G2", 3),   # rows 3–5
                rep("G3", 3),
                rep("G4", 2),
                rep("G5", 2),
                rep("G6", 2),
                rep("G7", 2))  

# col_fun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
# 
# ht <- Heatmap(
#   cor_matrix_spearman, name = "Cor",
#   col = col_fun,
#   na_col = "white",
#   column_names_side = "top",
#   border_gp = gpar(col = "black", lty = 1),
#   cluster_rows = FALSE,
#   cluster_columns = FALSE,
#   row_names_side = "left", 
#   row_split = row_groups,
#   show_row_dend = FALSE, 
#   row_gap = unit(2, "mm"),
#   rect_gp = gpar(col = "white", lwd = 2),
#   cell_fun = function(j, i, x, y, w, h, fill) {
#     lab <- combinedMat[i, j]
#     if (!is.na(lab) && nzchar(lab)) {
#       grid.text(lab, x = x, y = y, gp = gpar(fontsize = 7))
#     }
#   }
# )
# 
# png("meta.spearman.png", width = 2600, height = 2300, res = 300) 
# draw(ht, 
#      heatmap_legend_side = "left", 
#      row_title = NULL)   # 🚨 removes the "G1 G2..." labels this is not right
# dev.off()
# 
#Increasing font weight and size, and
#Changing the text color dynamically depending on the background brightness (so white text on dark colors, black text on light colors).
library(ComplexHeatmap)
library(circlize)
library(grid)

col_fun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

ht <- Heatmap(
  cor_matrix_spearman, name = "Cor",
  col = col_fun,
  na_col = "white",
  column_names_side = "top",
  border_gp = gpar(col = "black", lty = 1),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_side = "left",
  #row_split = row_groups,
  show_row_dend = FALSE,
  row_gap = unit(2, "mm"),
  rect_gp = gpar(col = "white", lwd = 2),
  cell_fun = function(j, i, x, y, w, h, fill) {
    lab <- combinedMat[i, j]
    if (!is.na(lab) && nzchar(lab)) {
      # Determine brightness of fill color to pick text color
      rgb_vals <- col2rgb(fill) / 255
      brightness <- 0.299 * rgb_vals[1] + 0.587 * rgb_vals[2] + 0.114 * rgb_vals[3]
      text_col <- if (brightness < 0.5) "white" else "black"
      
      grid.text(
        lab,
        x = x, y = y,
        gp = gpar(fontsize = 9, fontface = "bold", col = text_col)
      )
    }
  }
)

png(paste0(file,".spearman.png"), width = 8800, height = 8800, res = 300)
draw(ht, heatmap_legend_side = "left", row_title = NULL)
dev.off()







##################################
 combinedMat[p.mat > 0.05] <- ""
# 
# ht <- Heatmap(
#   cor_matrix_spearman, name = "Cor",
#   col = col_fun,
#   na_col = "white",
#   column_names_side = "top",
#   border_gp = gpar(col = "black", lty = 1),
#   cluster_rows = FALSE,
#   cluster_columns = FALSE,
#   row_names_side = "left", 
#   row_split = row_groups,
#   show_row_dend = FALSE, 
#   row_gap = unit(2, "mm"),
#   rect_gp = gpar(col = "white", lwd = 2),
#   cell_fun = function(j, i, x, y, w, h, fill) {
#     lab <- combinedMat[i, j]
#     if (!is.na(lab) && nzchar(lab)) {
#       grid.text(lab, x = x, y = y, gp = gpar(fontsize = 7))
#     }
#   }
# )
# 
# png("meta.spearmanSig.png", width = 2600, height = 2300, res = 300) 
# draw(ht, 
#      heatmap_legend_side = "left", 
#      row_title = NULL)   # 🚨 removes the "G1 G2..." labels this is not right
# dev.off()


library(ComplexHeatmap)
library(circlize)
library(grid)

col_fun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

ht <- Heatmap(
  cor_matrix_spearman, name = "Cor",
  col = col_fun,
  na_col = "white",
  column_names_side = "top",
  border_gp = gpar(col = "black", lty = 1),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_side = "left", 
  #row_split = row_groups,
  show_row_dend = FALSE, 
  row_gap = unit(2, "mm"),
  rect_gp = gpar(col = "white", lwd = 2),
  cell_fun = function(j, i, x, y, w, h, fill) {
    lab <- combinedMat[i, j]
    if (!is.na(lab) && nzchar(lab)) {
      # Compute brightness of the cell background
      rgb_vals <- col2rgb(fill) / 255
      brightness <- 0.299 * rgb_vals[1] + 0.587 * rgb_vals[2] + 0.114 * rgb_vals[3]
      text_col <- if (brightness < 0.5) "white" else "black"
      
      grid.text(
        lab,
        x = x, y = y,
        gp = gpar(fontsize = 9, fontface = "bold", col = text_col)
      )
    }
  }
)

png(paste0(file, ".spearmanSig.png"), width = 3800, height = 3800, res = 300)
draw(ht,
     heatmap_legend_side = "left",
     row_title_gp = gpar(fontsize = 10, fontface = "bold"))  # keeps row titles visible and bold
dev.off()

myCorChar=cor_char
myCorChar <- gsub("-", "\u2212", cor_char)
combinedMat2 <- matrix(
  paste0(myCorChar, pmMat),
  nrow = nrow(cor_matrix_spearman),
  ncol = ncol(cor_matrix_spearman),
  dimnames = dimnames(cor_matrix_spearman)
)
combinedMat2[p.mat > 0.05] <- ""

library(ComplexHeatmap)
library(circlize)
library(grid)

col_fun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

ht <- Heatmap(
  cor_matrix_spearman, name = "Cor",
  col = col_fun,
  na_col = "white",
  column_names_side = "top",
  border_gp = gpar(col = "black", lty = 1),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_side = "left", 
  #row_split = row_groups,
  show_row_dend = FALSE, 
  row_gap = unit(2, "mm"),
  rect_gp = gpar(col = "white", lwd = 2),
  cell_fun = function(j, i, x, y, w, h, fill) {
    lab <- combinedMat2[i, j]
    if (!is.na(lab) && nzchar(lab)) {
      # Compute brightness of the cell background
      rgb_vals <- col2rgb(fill) / 255
      brightness <- 0.299 * rgb_vals[1] + 0.587 * rgb_vals[2] + 0.114 * rgb_vals[3]
      text_col <- if (brightness < 0.5) "white" else "black"
      
      grid.text(
        lab,
        x = x, y = y,
        gp = gpar(fontsize = 9, fontface = "bold", col = text_col)
      )
    }
  }
)

png(paste0(file,".spearmanSig2.png"), width = 3800, height = 3800, res = 300)
draw(ht,
     heatmap_legend_side = "left",
     row_title_gp = gpar(fontsize = 10, fontface = "bold"))  # keeps row titles visible and bold
dev.off()




