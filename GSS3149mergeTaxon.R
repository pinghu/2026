rm(list=ls())

library(ggpubr)

library(dplyr)
library(tidyverse)
args <- commandArgs(trailingOnly = TRUE)
filename <- args[1]


#filename="SalmonOral83.clean.ReadsNum.ann"
cnts<-read.table(filename, sep="\t", header=TRUE)
d=dim(cnts)
geneID=cnts[,1]
samplenames <- colnames(cnts)[2:d[2]]
Clen=length(samplenames)
df <-data.frame(cnts)



# Sum counts per species
keggID_counts <- df %>%
  group_by(KeggID) %>%
  summarise(across(starts_with("GSS"), sum), .groups = "drop")
library(openxlsx)

# Save to Excel
write.xlsx(keggID_counts, "KO_counts.xlsx", rowNames = FALSE)
# Save to tab-delimited file
write.table(
  keggID_counts,
  file = "keggID_counts.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

