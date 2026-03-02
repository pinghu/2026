rm(list=ls())
#remotes::install_github('YuLab-SMU/ggtree') 
####Later to check Mesh and cell marker data
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE) 
library(ReactomePA)
#library(MeSHDbi)
#library(AnnotationHub)
#library(meshes)
##check their visualization: https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html
######Given a list tell you GO and KEGG matching U133 entrenz gene id
perform_gene_enrichment_analysis <- function(gene_symbols, outname) {
  # Map gene symbols to Entrez Gene IDs using org.Hs.eg.db
  entrez_ids <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  # Extract unique Entrez Gene IDs
  entrez_ids <- unique(entrez_ids$ENTREZID)
  gene_symbols <- unique(gene_symbols)
  
  # Run GO enrichment analysis using clusterProfiler for Biological Process (BP), Molecular Function (MF), and Cellular Component (CC)
  go_BP <- enrichGO(gene = entrez_ids, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP", pvalueCutoff  = 0.05)
  go_MF <- enrichGO(gene = entrez_ids, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "MF", pvalueCutoff  = 0.05)
  go_CC <- enrichGO(gene = entrez_ids, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "CC", pvalueCutoff  = 0.05)
  kk <- enrichKEGG(gene= entrez_ids, organism     = 'hsa', pvalueCutoff = 0.05)
  wp<-enrichWP(entrez_ids, organism = "Homo sapiens", pvalueCutoff = 0.05) 
  react <-enrichPathway(gene=entrez_ids, pvalueCutoff = 0.05, readable=TRUE)
  dgn <- enrichDGN(entrez_ids) 
  x <- enrichDO(gene= entrez_ids, pvalueCutoff  = 0.05)
  bp_df <- as.data.frame(go_BP)
  if (nrow(bp_df) > 0) {
    bp_df$Category <- "goBP"
  } else {
    bp_df <- data.frame(Category = character(0))
  }
  if (nrow(dgn) > 0) {
    dgn$Category <- "DiseaseGeneNetwork"
  } else {
    dgn <- data.frame(Category = character(0))
  }
  if (nrow(x) > 0) {
    x$Category <- "DiseaseOntology"
  } else {
    x <- data.frame(Category = character(0))
  }
  if (nrow(kk) > 0) {
    kk$Category <- "kegg"
  } else {
    kk <- data.frame(Category = character(0))
  }
  
  if (nrow(wp) > 0) {
    wp$Category <- "wikiPathway"
  } else {
    wp <- data.frame(Category = character(0))
  }
  
  if (nrow(react) > 0) {
    react$Category <- "Reactone"
  } else {
    react <- data.frame(Category = character(0))
  }
  # Convert and add category for MF
  mf_df <- as.data.frame(go_MF)
  if (nrow(mf_df) > 0) {
    mf_df$Category <- "goMF"
  } else {
    mf_df <- data.frame(Category = character(0))
  }
  
  # Convert and add category for CC
  cc_df <- as.data.frame(go_CC)
  if (nrow(cc_df) > 0) {
    cc_df$Category <- "goCC"
  } else {
    cc_df <- data.frame(Category = character(0))
  }
  
  
  # Combine all into one table
  go_combined <- rbind(mf_df, bp_df, cc_df, kk, wp, react, dgn, x)
  
  # Write to file
  # write.table(go_combined,
  #             file = paste0(outname, ".GO_clusterprofiler.tsv"),
  #             sep = "\t",
  #             row.names = TRUE,
  #             col.names = NA,
  #             eol = "\n",
  #             na = "NA",
  #             quote = FALSE)
  
  return(go_combined)
}

args <- commandArgs(trailingOnly = TRUE)
print(args)
filename <- args[1]
outname <- args[2]
#filename="GSS3110.10filter.stat.xls.GER_vs_V1.stat"
rm(args)
#filename="agree" ###remove long gene name, there are issure with it
#filename ="GSS3049.T_V2.TNFalpha.xls"
A<-read.table(filename, sep="\t", header=TRUE)
A$gene <- gsub("_[A-Za-z0-9]+$", "", A[,1])
Asig=A[as.numeric(A[,2])<=0.05, ]
AsigU = Asig[as.numeric(Asig[,3])>0,]
AsigD = Asig[as.numeric(Asig[,3])<0,]
#5.2 GO classification "MF", "BP", and "CC" subontologies.
#myEZ=AsigU[,3]
ResultU=perform_gene_enrichment_analysis(AsigU$gene, paste0(filename, ".sigUp"))
ResultD=perform_gene_enrichment_analysis(AsigD$gene, paste0(filename, ".sigDown"))
ResultA=perform_gene_enrichment_analysis(Asig$gene, paste0(filename, ".sigAll"))

# Rename columns (keeping the first two columns untouched)
colnames(ResultA)[-c(1, 2)] <- paste0(colnames(ResultA)[-c(1, 2)], ".All")
colnames(ResultD)[-c(1, 2)] <- paste0(colnames(ResultD)[-c(1, 2)], ".Down")
colnames(ResultU)[-c(1, 2)] <- paste0(colnames(ResultU)[-c(1, 2)], ".Up")

# Get the names of the first two columns (e.g., "ID" and "Category")
merge_keys <- colnames(ResultA)[1:2]

# Merge all three tables by the first two columns
merged_AD <- merge(ResultA, ResultD, by = merge_keys, all = TRUE)
merged_all <- merge(merged_AD, ResultU, by = merge_keys, all = TRUE)


# Create a new numeric column with the extracted numerator
merged_all$SigAllNum <- ifelse(
  is.na(merged_all$GeneRatio.All),
  0,
  as.numeric(sub("/.*", "", merged_all$GeneRatio.All))
)

merged_all$SigUpNum <- ifelse(
  is.na(merged_all$GeneRatio.Up),
  0,
  as.numeric(sub("/.*", "", merged_all$GeneRatio.Up))
)

merged_all$SigDownNum <- ifelse(
  is.na(merged_all$GeneRatio.Down),
  0,
  as.numeric(sub("/.*", "", merged_all$GeneRatio.Down))
)

merged_all$direction <- ifelse(
  merged_all$SigUpNum > merged_all$SigDownNum, "Up",
  ifelse(merged_all$SigDownNum > merged_all$SigUpNum, "Down", "Sig")
)

write.table(merged_all, file = paste0(outname, ".clusterprofiler.txt"), sep = "\t", row.names = TRUE, col.names = TRUE, eol = "\n", na = "NA")
write.xlsx(merged_all, file = paste0(outname, ".clusterprofiler.xlsx"), rowNames = FALSE)

