##########################
#Ping Hu April 3rd, 2025
##############################
rm(list=ls())
args <- commandArgs(trailingOnly = TRUE)

filename <- args[1]

#filename="metaphlan.relab10.7"
#output_file <- paste0(filename, ".stat.csv")  # Define output file name
#results_df <- data.frame()
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
filename1="~/Desktop/Disk2/project/GSS3149OralBiofilm/GSS3149OralBiofilm.meta.xls"
A<-read.table(filename1, sep="\t", header=TRUE)
d <- dim(A);
order_levels <- c("NoTreat","NCColgCP","Formula1","Formula2","Formula4", "Formula5")
A$Treat <- factor(A$Treat, levels = order_levels, ordered = TRUE)

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

A0<-read.table(filename, sep="\t", header=TRUE)
d <- dim(A0);
B=A0[1:d[1], 2:d[2]]
ZZ=as.numeric(min(B[B>0&!is.na(B)]))/100
C=B+ZZ

Cname=colnames(A0)[2:d[2]]
Clen=length(Cname) ##there are 3 annotation columns
splitname<-strsplit(Cname, "[_]")
SID=rep("NA", Clen)
for(mm in  1:Clen ){
  SID[mm]=splitname[[mm]][3]
}
SID=as.numeric(SID)
matching_dataset <- A %>%
  filter(SID %in% SID) %>%
  arrange(match(SID, SID))


 trt_groups=unique(A$Treat)
# # Generate all possible pairwise comparisons
 trt_pairs <- combn(trt_groups, 2, simplify = FALSE)
big_results_list <- list()
for (i in 1:d[1]){
  genename=A0[i,1]
  gene<-as.numeric(C[i,])
  meanAll <-mean(gene)
  splitG<-strsplit(as.character(genename), "[|]")
  LLL=length(splitG[[1]])
  genus=splitG[[1]][LLL]
  cleaned_genus <- gsub("^\\w__", "", genus)  # Removes the prefix
  cleaned_genus <- gsub("_", " ", cleaned_genus)  # Replaces underscores with spaces
  percentZeroNA=(sum(B[i,]==0)+sum(is.na(B[i,])))/(d[2]-1)
  mydata0=data.frame(gene, SID)
  mydata <- inner_join(mydata0, matching_dataset, by = c("SID" = "SID"))

  KP_Treat_AR=kruskal.test(gene ~ Treat, data = mydata)$p.value
  
  summary_testdata2 <- mydata %>%
    group_by(Treat) %>%
    summarize(
      mean_RA = mean(gene, na.rm = TRUE),
      se_RA = sd(gene, na.rm = TRUE) / sqrt(n()),
      count=n()
    )
  
  # Initialize results list for pairwise tests
  results_list <- list()
  results_list[["TaxonName"]] <-cleaned_genus
  results_list[["TaxonFullName"]] <-genename
  results_list[["PercentZero"]] <-percentZeroNA
  results_list[["MeanAllAbsolute"]] <-meanAll
  
  for ( test in trt_groups) { 
    results_list[[paste0("mean_AR.", test)]] <- summary_testdata2$mean_RA[summary_testdata2$Treat == test]
    results_list[[paste0("se_AR.", test)]] <- summary_testdata2$se_RA[summary_testdata2$Treat == test]
    results_list[[paste0("sample_number.", test)]] <- summary_testdata2$count[summary_testdata2$Treat == test]
    
  }
  results_list[["Kruskal_P.AR"]] <-KP_Treat_AR
 
  for (pair in trt_pairs) {
    
    trt1 <- pair[1]
    trt2 <- pair[2]
    
    # Extract the gene values for each treatment group
    data1 <- mydata[mydata$Treat == trt1, "gene"]
    data2 <- mydata[mydata$Treat == trt2, "gene"]
    
    if (length(data1) > 1 && length(data2) > 1) {
      
      m1 <- mean(data1, na.rm = TRUE)
      m2 <- mean(data2, na.rm = TRUE)
      fc <- m1 / m2
      
      TFC_AR <- truefc(fc)
      
      # Check if there's variation before running Wilcoxon
      if (length(unique(data1)) > 1 && length(unique(data2)) > 1) {
        wilcox_AR <- wilcox.test(data1, data2, paired = FALSE)$p.value
      } else {
        wilcox_AR <- NA
      }
      
      # Store results with named keys
      label <- paste0(trt1, "_vs_", trt2)
      results_list[[paste0("wilcox_AR.", label)]] <- wilcox_AR
      results_list[[paste0("truefold_AR.", label)]] <- TFC_AR
    }
  }
  
  # Store the transformed results for this gene
  big_results_list[[genename]] <- results_list
  
  ### Now do the graph
      png(filename = paste0(genus, ".relative.png"), width = 1000, height = 1000, res = 300)
      bxp <- ggboxplot(mydata, x = "Treat", y = "gene", color = "Treat") +
        ggtitle(paste0(cleaned_genus,"\nKP=", sprintf("%.02f", KP_Treat_AR), " 0%=", sprintf("%.02f", 100*percentZeroNA))) + ####This is for pathway
        geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.7) +
        theme(
          legend.position = "none",
          axis.text.x = element_text(angle = 60, hjust = 1),
          plot.title = element_text(face = "italic") # Italicizes the entire title
        )
      print(bxp)
      dev.off()

      
}
# Convert big results list to a single data frame
big_results_table <- bind_rows(big_results_list)
pData <- big_results_table %>% select("Kruskal_P.AR", starts_with("wilcox_AR"))
pData[is.na(pData)] <- 1
d=dim(pData)
df<-data.frame(matrix(ncol=d[2], nrow=d[1]))
colnames(df)=paste0("fdr.",colnames(pData))

for (i in 1:d[2]) {
  p <- as.numeric(pData[[i]])       # Use [[i]] or as.numeric() to avoid list issues
  df[, i] <- p.adjust(p, method = "fdr")
}
final_result_table <- cbind(big_results_table, df)
write.table(final_result_table, file = paste0(filename, ".relative_stat"), sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)     
#library(openxlsx)  # Load the package
output_file <- paste0(filename, ".relative_stat.xlsx")
write.xlsx(final_result_table, file = output_file, rowNames = FALSE, colNames = TRUE)

print(paste("Results successfully written to:", output_file))


# Filter the big_results_table for rows where KP_Treat_AR <= 0.05
filtered_results1 <- final_result_table[final_result_table$Kruskal_P.AR <= 0.05,]
output_file2 <- paste0(filename, ".", dim(filtered_results1)[1], ".relative_sig.xlsx")
write.xlsx(filtered_results1, file = output_file2, rowNames = FALSE, colNames = TRUE)


# Filter the big_results_table for rows where KP_Treat_AR <= 0.05
filtered_results <- final_result_table[final_result_table$fdr.Kruskal_P.AR <= 0.1,]
output_file3 <- paste0(filename, ".", dim(filtered_results)[1], ".relative_fdr.1.xlsx")
write.xlsx(filtered_results, file = output_file3, rowNames = FALSE, colNames = TRUE)

meanData <- filtered_results %>% select( starts_with("mean_AR"))
rownames(meanData) =filtered_results$TaxonName
mean_zscore_mat <- t(apply(meanData, 1, function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)))
rownames(meanData) =filtered_results$TaxonName
fcData <- filtered_results %>% select(starts_with("truefold_AR"))
rownames(fcData) =filtered_results$TaxonName

# Convert fcData to a numeric matrix if it's not already
# fcData2 <- as.matrix(fcData)
# mode(fcData2) <- "numeric"  # force all elements to numeric if needed
# # Apply thresholding
# fcData2 <- pmin(pmax(fcData2, -5), 5)
pData<-as.matrix(pData)
mode(pData) <- "numeric" 
pData[is.na(pData)] <- 0
#if not switch to the matrix and change mode to numeric it will lead mistake for all the formular1 to notreat turn into -5, not sure why
fcData <- as.matrix(fcData)
mode(fcData) <- "numeric" 
fcData[fcData>5] <- 5
fcData[fcData <= -5] = -5
#pData <- filtered_results %>% select("TaxonFullName", starts_with("wilcox_AR"))
#adjData <-filtered_results %>% select("TaxonFullName", starts_with("adj.wilcox_AR"))
pData <- filtered_results %>% select(starts_with("wilcox_AR"))
adjData <-filtered_results %>% select(starts_with("adj.wilcox_AR"))
create_pk_mat <- function(pk, adjpk) {
  nrow_pk <- nrow(pk)
  ncol_pk <- ncol(pk)
  
  pk_mat <- matrix("", nrow = nrow_pk, ncol = ncol_pk)
  
  # Set the row and column names
  rownames(pk_mat) <- rownames(pk)
  colnames(pk_mat) <- colnames(pk)
  
  pk_mat[!is.na(pk) & pk <= 0.1] <- "*"
  pk_mat[!is.na(pk) & pk <= 0.05] <- "**"
  pk_mat[!is.na(adjpk) & adjpk <= 0.1] <- "***"
  pk_mat[!is.na(adjpk) & adjpk <= 0.05] <- "****"
  
  return(pk_mat)
}
pmat <-create_pk_mat(pData, adjData)

ht4 <- Heatmap(fcData,name="Fold Change", 
               column_names_gp = gpar(fontsize = 9), 
               cluster_columns = FALSE, 
               row_names_side = "left", 
               rect_gp = gpar(col = "white", lwd = 2),
               row_names_gp = gpar(fontsize = 8, fontface = "italic"), 
               cluster_rows = FALSE, 
               cell_fun = function(j, i, x, y, w, h, col) {
                 grid.text(pmat[i, j], x, y)
               },
               col = colorRamp2(c(-5, 0, 5), c("blue", "white", "red")))

png(paste0(filename, ".relativefoldheat.png"), width = 1200, height = 1200, res = 150)  # Adjust size and resolution as needed
#print(htcor)
draw(ht4)
# Close the PNG device
dev.off()

htcor<- Heatmap(mean_zscore_mat, name = "Zscore Relative Abundance", rect_gp = gpar(col = "white", lwd = 2),
                row_names_gp = gpar(fontsize = 8, fontface = "italic"),
                row_names_side = "left",
                #column_names_rot = 0,  # Keep column names horizontal
                column_names_gp = gpar(fontsize = 8),  # Reduce column name size
                column_title = "Relative Percentage",
                #,  cluster_rows = FALSE,cluster_columns = FALSE
                show_row_dend = FALSE,  # Hide row dendrogram
                show_column_dend = FALSE, 
                cluster_columns = FALSE
)
#height2300 for species
png(paste0(filename, ".relativeMeanHeat.png"), width = 800, height = 1200, res = 150)  # Adjust size and resolution as needed
#print(htcor)
draw(htcor)
# Close the PNG device
dev.off()
ht_list=htcor+ht4
jpeg(paste0(filename, ".relative.jpg"), width = 2500, height = 1200, res=150 )
draw(ht_list, heatmap_legend_side = "right", 
     column_title = paste( filename, "\n relative abundance"))
dev.off()

