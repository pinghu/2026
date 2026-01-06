rm(list=ls())

library(ggpubr)

library(dplyr)
library(tidyverse)
args <- commandArgs(trailingOnly = TRUE)
filename <- args[1]


filename="SalmonOral83.clean.ReadsNum.10filter.ann"
cnts<-read.table(filename, sep="\t", header=TRUE)
d=dim(cnts)
geneID=cnts[,1]
samplenames <- colnames(cnts)[2:d[2]]
Clen=length(samplenames)
df <-data.frame(cnts)




# Extract species from Gene
df <- df %>%
  mutate(Species = sub(":.*", "", gene))  # take part before ":"

# Sum counts per species
species_counts <- df %>%
  group_by(Species) %>%
  summarise(across(starts_with("GSS"), sum), .groups = "drop")
library(openxlsx)

# Save to Excel
write.xlsx(species_counts, "species_counts.xlsx", rowNames = FALSE)
# Save to tab-delimited file
write.table(
  species_counts,
  file = "species_counts.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

