##########################################
#Ping Hu April 16, 2025
### Would like to add into the PPT summary
#### Need to filter out the summary graph with the %zero to be less than 75%
##########################################
rm(list=ls())
library(dplyr)
library(ggplot2)
library(ggpubr)

args <- commandArgs(trailingOnly = TRUE)
#filename1="GSS3215Pam_meta.txt"
#X<-read.table(filename1, sep="\t", header=TRUE)
filename=args[1]
filename="full_Paired_Data"

A0<-read.table(filename, sep="\t", header=TRUE)
d <- dim(A0);
B=A0[1:d[1], 2:d[2]]
rownames(B)=A0[,1]

stopifnot(is.matrix(B) || is.data.frame(B))
B <- as.matrix(B)
d<-dim(B)
stopifnot(nrow(B) == d[1], ncol(B) == d[2])
p <- d[2]/2
baseline <- B[, 1:p, drop = FALSE]
week8    <- B[, (p + 1):(2 * p), drop = FALSE]
## delta (week8 - baseline), same size 37951 x 175
delta <- week8 - baseline

## Make sure delta is a matrix
delta <- as.matrix(delta)
## Make sure delta has column names
colnames(delta) <- gsub("W8", "W8-BL", colnames(delta))
## Write to tab-delimited text file
write.table(delta,
            file = "delta_matrix.txt",
            sep = "\t",
            quote = FALSE,
            row.names = TRUE,
            col.names = NA)

write.table(week8,
            file = "week8_matrix.txt",
            sep = "\t",
            quote = FALSE,
            row.names = TRUE,
            col.names = NA)
write.table(baseline,
            file = "baseline_matrix.txt",
            sep = "\t",
            quote = FALSE,
            row.names = TRUE,
            col.names = NA)

