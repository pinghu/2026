rm(list=ls())
library(ggplot2)
library(ggpubr)
library(rstatix)
library(dplyr)
library(tidyverse)
library("plotrix")
library(vegan)
library(openxlsx)

library(officer)
library(rvg)
library(flextable)

wb <- createWorkbook()
filename1="GSS3215Pam_meta.txt"

args <- commandArgs(trailingOnly = TRUE)
filename<- args[1]
#filename="salmon_regular.clean.NumReads"
A<-read.table(filename1, sep="\t", header=TRUE)
d <- dim(A);

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE),
      se= std.error(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

my.t.test.p.value <- function(...) {
  obj<-try(t.test(...), silent=TRUE)
  if (inherits(obj, "try-error")) return(NA) else return(obj$p.value)
}

my.wilcox.p.value <- function(...) {
  obj<-try(wilcox.test(...), silent=TRUE)
  if (inherits(obj, "try-error")) return(NA) else return(obj$p.value)
}
my.kruskal.p.value <- function(...) {
  obj<-try(kruskal.test(...), silent=TRUE)
  if (inherits(obj, "try-error")) return(NA) else return(obj$p.value)
}

truefc<-function(VVV){
  #print(VVV)
  if (is.finite(VVV )){
    XXX=VVV
    if(VVV==0){
      XXX=NA
    }else if(VVV<1){
      XXX=-1/VVV
    }
    return(XXX)
  }else{
    return("NA")
  }
}

test <- read.table(filename, header = TRUE, row.names = 1, sep="\t")
test[is.na(test)]<-0
d=dim(test)  #103133
d
X=apply(test, 1, sum)

test_filter2=test[(X/d[2] >=10),] #23419
idx <- match(A$ID, colnames(test_filter2))

# Optional safety check
if (any(is.na(idx))) {
  stop("Some A$cleanMeta values are not found in colnames(test_filter2).")
}

test_filter2_mat <- as.matrix(test_filter2[, idx])
#new_names <- as.character(A$NewID)  # use the real column, capital D
#stopifnot(length(new_names) == ncol(test_filter2_mat))

#colnames(test_filter2_mat) <- new_names

#colnames(test_filter2_mat)  # should now show your No.Arm..., Yes.Cheek..., etc.

relative_fraction_mat <- sweep(
  test_filter2_mat,
  MARGIN = 2,                    # operate over columns
  STATS  = A$original.total,
  FUN    = "/"
)
## 4. If you want *percentages* instead of fractions:
relative_percent_mat <- relative_fraction_mat * 100
# ---- Save as tab-delimited file ----
write.table(
  relative_percent_mat,
  file = paste0(filename, ".relative_percent_mat.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = TRUE,
  col.names = TRUE
)
write.table(
  test_filter2_mat,
  file = paste0(filename, ".test_filter2_mat.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = TRUE,
  col.names = TRUE
)
write.table(
  test_filter2,
  file = paste0(filename, ".test_filter2.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = TRUE,
  col.names = TRUE
)

# ---- Save as Excel file ----
library(openxlsx)

write.xlsx(
  relative_percent_mat,
  file = paste0(filename, ".relative_percent_mat.xlsx"),
  rowNames = TRUE
)
write.xlsx(
  test_filter2_mat,
  file = paste0(filename, ".test_filter2_mat.xlsx"),
  rowNames = TRUE
)
write.xlsx(
  test_filter2,
  file = paste0(filename, ".test_filter2.xlsx"),
  rowNames = TRUE
)

