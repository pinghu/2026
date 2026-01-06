rm(list=ls())
library(ggplot2)
library(ggpubr)
library(rstatix)
library(dplyr)
library(tidyverse)
library(ggtext)
library(xfun)

library(openxlsx) 
library(circlize)
library(ComplexHeatmap)
library(lattice)   
library(reshape2)
library(gridExtra)

library(ggplot2)
library(stringr)
library(officer)
library(rvg)

doc <- read_pptx()
drawSortedBarPlot <- function(testData, filename) {
  # Ensure testData is a data frame
  if (!is.data.frame(testData)) {
    testData <- as.data.frame(testData)
  }
  
  # Make sure testData has the expected columns
  if(!("Sample" %in% names(testData)) || !("similarityScore" %in% names(testData))) {
    stop("testData must have 'Sample' and 'similarityScore' columns")
  }
  
  # Convert similarityScore to numeric if it's not already
  testData$similarityScore <- as.numeric(as.character(testData$similarityScore))
  
  # In case of any NA introduced by conversion (e.g., non-numeric values present), handle or warn
  if(any(is.na(testData$similarityScore))) {
    warning("NAs introduced by coercion to numeric in similarityScore")
  }
  
  color_vector<- ifelse(grepl("etro",testData$Sample), "blue",
                        ifelse(grepl("SLS|TNF|Geraniol",testData$Sample), "red", "steelblue"))
  testData$Col =color_vector
  
  sorted_df <- testData[order(testData$similarityScore), ]
  p <- ggplot(sorted_df, aes(y = reorder(Sample, similarityScore), x = similarityScore)) +
    geom_bar(stat = "identity", fill = sorted_df$Col) +
    ylab(filename) +
    xlab("Similarity Score") +
    ggtitle("Rank") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size = 12))
  
  # Save the plot as a JPEG image file
  filename_with_extension <- paste0(filename, ".jpg")
  jpeg(filename_with_extension, width = 285, height = 380)
  print(p)
  dev.off()
 # caption_text = paste(outname, ": unique genes ", length(gene_symbols), "; database ", db)
    
  # Use rank with ties.method = "min" to handle scores that are the same
  sorted_df$order <- rank(sorted_df$similarityScore, ties.method = "min")
  
  # Specify your desired output file name
  output_tab_delimited_file <- paste0(filename, ".rank")
  
  # Write the specific columns to a tab-delimited text file
  write.table(sorted_df[, c("Sample", "order", "similarityScore")], 
              file = output_tab_delimited_file, 
              sep = "\t", 
              row.names = FALSE, 
              col.names = TRUE)
  doc <- doc %>%
    add_slide(layout = "Title and Content", master = "Office Theme") %>%
    ph_with(
      fpar(ftext(filename, prop = fp_text(font.size = 20))),  # Adjust size as needed
      location = ph_location_type(type = "title")
    )%>%
    ph_with(dml(ggobj = p), location = ph_location(left = 1, top = 1.5, width = 5, height = 5)) %>%
    ph_with(sorted_df[, c("Sample", "order", "similarityScore")], 
            location = ph_location(left = 6.2, top = 2.0, width = 3.5, height = 5))
  
  
}

calculateSimilarityScores <- function(distance_dict) {
  # Flatten the distance dictionary to a numeric vector
  myDist <- unlist(distance_dict)
  
  # Find minimum and maximum distances
  min_dist <- min(myDist)
  max_dist <- max(myDist)
  
  # Calculate similarity scores based on distances from -1 to 1
  #similarityScore <- 2 * (myDist - min_dist) / (max_dist - min_dist) - 1
  # Calculate similarity scores based on distances from 0 to 1
  similarityScore <- (myDist - min_dist) / (max_dist - min_dist)
  return(similarityScore)
}

#args <- commandArgs(trailingOnly = TRUE)
#filename <- args[1]
#outname <-arg[2]

filename="SalmonOral83.clean.ReadsNum.10filter"
metafile="~/Desktop/Disk2/project/GSS3149OralBiofilm/GSS3149OralBiofilm.meta.xls"
A2<-read.table(metafile, sep="\t", header=TRUE)
cnts<-read.table(filename, sep="\t", header=TRUE)
outname="salmon83"

d=dim(cnts)
geneID=cnts[,1]
samplenames <- colnames(cnts)[2:d[2]]
Clen=length(samplenames)

splitname<-strsplit(samplenames, "[_]")
SID=rep("NA", Clen)
for(mm in  1:Clen ){
  SID[mm]=splitname[[mm]][3]
}
SID=as.numeric(SID)
matching_dataset <- A2 %>%
  filter(SID %in% SID) %>%
  arrange(match(SID, SID))

trt_groups=unique(A2$Treat)

cnts = cnts[,2:d[2]]
colnames(cnts)=matching_dataset$NewID
rownames(cnts)=geneID
matching_dataset$category= matching_dataset$Treat


###########################################
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# # The following initializes usage of Bioc devel
# BiocManager::install(version='devel')
# 
# BiocManager::install("DESeq2")
library("DESeq2")
y=round(cnts)
ddsMat <- DESeqDataSetFromMatrix(countData = y,colData = matching_dataset, design = ~ Treat)

#nrow(ddsMat)###9381
#keep <- rowSums(counts(dds)) > 1
# at least 3 samples with a count of 10 or higher

keep <- rowSums(counts(ddsMat) >= 10) >= 3
dds <- ddsMat[keep,]

###nrow(dds)###9381
#install.packages("hexbin")
###################################################################
#DESeq2 offers two transformations for count data that stabilize the variance across the mean: the variance stabilizing transformation (VST) for negative binomial data with a dispersion-mean trend (Anders and Huber 2010), implemented in the vst function, and the regularized-logarithm transformation or rlog (Love, Huber, and Anders 2014).For genes with high counts, both the VST and the rlog will give similar result to the ordinary log2 transformation of normalized counts. For genes with lower counts, however, the values are shrunken towards a middle value. The VST or rlog-transformed data then become approximately homoskedastic (more flat trend in the meanSdPlot), and can be used directly for computing distances between samples, making PCA plots, or as input to downstream methods which perform best with homoskedastic data.The VST is much faster to compute and is less sensitive to high count outliers than the rlog. The rlog tends to work well on small datasets (n < 30), potentially outperforming the VST when there is a wide range of sequencing depth across samples (an order of magnitude difference). We therefore recommend the VST for medium-to-large datasets (n > 30). You can perform both transformations and compare the meanSdPlot or PCA plots generated, as described below.
##########################################################
## vst: variance stablizing transformation----------------------------
###In the above function calls, we specified blind = FALSE, which means that differences between cell lines and treatment (the variables in the design) will not contribute to the expected variance-mean trend of the experiment. The experimental design is not used directly in the transformation, only in estimating the global amount of variability in the counts. For a fully unsupervised transformation, one can set blind = TRUE (which is the default).
##############################################
library("vsn")
#vsd <- vst(dds, blind = FALSE)
vsd <- vst(dds, blind = FALSE, fitType = "local")

#head(assay(vsd), 3)
#colData(vsd)
############################################################################
## ----rlog is extremely time consuming, consider to get away from it. 
##########################################################################
#rld <- rlog(dds, blind = FALSE)
#head(assay(rld), 3)
######################################################################
## ----transformplot, fig.width = 6, fig.height = 2.5------------------------
############################################################
library("dplyr")
library("ggplot2")

dds <- estimateSizeFactors(dds)
new_cnt<-counts(dds, normalized=TRUE)
###################################
#These plot is not as good as the limma data
#I do not like the new count glimma plot, did not seperate the sample well. 
#glMDSPlot(new_cnt, labels=group, 
#          groups=trt, launch=TRUE)

library(Glimma)
glMDSPlot(log2(cnts+1), labels=matching_dataset$NewID, 
          groups=matching_dataset$Treat, launch=TRUE)
##############################################
write.csv(new_cnt, file = "Deseq.normalized.count.csv")
#write.csv(vsd, file = "Deseq.vsd.csv")
#write.csv(rld, file = "Deseq.rld.csv")
write.csv(cnts, file = "Desq.cleaned.count.csv")


###########################################################################
## https://davetang.org/muse/2018/05/15/making-a-heatmap-in-r-with-the-pheatmap-package/
#install.packages("pheatmap")
########################################################
sampleDists <- dist(t(assay(vsd)))
sampleDists
library("pheatmap")
library("RColorBrewer")

## ----distheatmap, fig.width = 6.1, fig.height = 4.5------------------------
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- colnames(vsd)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
my_sample_col <- data.frame(sample = matching_dataset$Treat)
row.names(my_sample_col) <- colnames(vsd)
jpeg('sampleDistance.jpg', width = 6000, height = 6000,  res = 100)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
	 #annotation_row = my_sample_col,
         col = colors)
dev.off()
doc <- doc %>%
  add_slide(layout = "Title and Content", master = "Office Theme") %>%
  ph_with(value = paste0("SampleDistance ", filename),
          location = ph_location_type(type = "title")) %>%
  ph_with(external_img("sampleDistance.jpg", width = 9, height = 5), 
          location = ph_location_type(type = "body"))

#######################################################
## Another option for calculating sample distances is to use the Poisson Distance (Witten 2011), implemented in the PoiClaClu package. This measure of dissimilarity between counts also takes the inherent variance structure of counts into consideration when calculating the distances between samples. The PoissonDistance function takes the original count matrix (not normalized) with samples as rows instead of columns, so we need to transpose the counts in dds.

#install.packages("PoiClaClu")
###################################################
library("PoiClaClu")
poisd <- PoissonDistance(t(counts(dds)))
## ----poisdistheatmap, fig.width = 6.1, fig.height = 4.5--------------------
samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- colnames(dds)
colnames(samplePoisDistMatrix) <- NULL

jpeg('samplePoisonDistance.jpg', width = 3000, height = 3000,  res = 300)
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,	 
         col = colors)
dev.off()

doc <- doc %>%
  add_slide(layout = "Title and Content", master = "Office Theme") %>%
  ph_with(value = paste0("SamplePoisonDistance ", filename),
          location = ph_location_type(type = "title")) %>%
  ph_with(external_img("samplePoisonDistance.jpg", width = 9, height = 5), 
          location = ph_location_type(type = "body"))
###################################################################
## ----plotpca, fig.width=6, fig.height=4.5----------------------------------
###############################################################
jpeg('PCA.jpg')
plotPCA(vsd, intgroup = "Treat")
dev.off()

doc <- doc %>%
  add_slide(layout = "Title and Content", master = "Office Theme") %>%
  ph_with(value = paste0("PCA ", filename),
          location = ph_location_type(type = "title")) %>%
  ph_with(external_img("PCA.jpg", width = 9, height = 5), 
          location = ph_location_type(type = "body"))


## --------------------------------------------------------------------------
pcaData <- plotPCA(vsd, intgroup = "Treat", returnData = TRUE)
pcaData
percentVar <- round(100 * attr(pcaData, "percentVar"))

## ----ggplotpca, fig.width=6, fig.height=4.5--------------------------------


# Calculate the reference point
reference_point <- pcaData %>%
  filter(Treat == "NCColgCP") %>%  # Filter for control group
  summarise(
    avg_PC1 = mean(PC1),  # Calculate average of PC1
    avg_PC2 = mean(PC2)   # Calculate average of PC2
  )

# Calculate the average PC1 and PC2 for each treatment group
avePCAData <- pcaData %>%
  group_by(Treat) %>%
  summarise(
    avePC1 = mean(PC1, na.rm = TRUE), # Calculate the average of PC1, remove NA values
    avePC2 = mean(PC2, na.rm = TRUE)  # Calculate the average of PC2, remove NA values
  )
library(rstatix)
library(ggplot2)
library(ggpubr)
custom_colors <- c("NCColgCP" = "black", 
                   "Formula1" = "cyan", 
                   "Formula2" = "green",
                   "Formula4" = "hotpink",   
                   "Formula5" = "purple",
                   "NoTreat"="orange"
)

#avePCAData_filtered <- avePCAData[!(avePCAData$Treat %in% c("GER", "SLS", "TNF")), ]

jpeg(paste0(outname, '.PCA2.jpg'), width=800, height=880)
ggplot(pcaData, aes(x = PC1, y = PC2, color = Treat)) +
geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed()+theme_bw()+scale_color_manual(values = custom_colors) + 
  stat_chull(aes(color=Treat, fill=Treat), alpha=0.1, geom="polygon") +
  geom_text(data = avePCAData, aes(label = Treat, x = avePC1, y = avePC2), nudge_x =0.2, nudge_y = 0.2) + # Label first points
 # geom_text(data = first_points, aes(label = trt, x = PC1, y = PC2), nudge_x =2, nudge_y = 2) + # Label first points
  #geom_text(data = pcaData, aes(label = trt, x = PC1, y = PC2), nudge_x = 0.05, nudge_y = 0.05) # Adjust nudge_x and nudge_y as needed
  ggtitle("PCA with VST data")
dev.off()

doc <- doc %>%
  add_slide(layout = "Title and Content", master = "Office Theme") %>%
  ph_with(value = paste0("PCA2 ", filename),
          location = ph_location_type(type = "title")) %>%
  ph_with(external_img(paste0(outname, '.PCA2.jpg'), width = 9, height = 5), 
          location = ph_location_type(type = "body"))


#pcaData_filtered <- pcaData[!(pcaData$group %in% c("GER", "SLS", "TNF")), ]


# jpeg(paste0(outname, '.PCA2-filter.jpg'), width=800, height=880)
# ggplot(pcaData_filtered, aes(x = PC1, y = PC2, color = Treat)) +
#   geom_point(size = 5) +
#   scale_color_manual(values = custom_colors) + 
#   xlab(paste0("PC1: ", percentVar[1], "% variance")) +
#   ylab(paste0("PC2: ", percentVar[2], "% variance")) +
#   coord_fixed() +
#   theme_bw() +stat_chull(aes(color=Treat, fill=Treat), alpha=0.05, geom="polygon") +
#   #stat_chull(aes(color=trt, fill=trt), alpha=0.1, geom="polygon") +
#   geom_text(data = avePCAData_filtered, aes(label = Treat, x = avePC1, y = avePC2), 
#             nudge_x = 0, nudge_y = 0, size = 8) +  # Increased font size here
#   ggtitle("PCA with VST data-no irritant")
# dev.off()
# doc <- doc %>%
#   add_slide(layout = "Title and Content", master = "Office Theme") %>%
#   ph_with(value = paste0("PCA2-filter ", filename),
#           location = ph_location_type(type = "title")) %>%
#   ph_with(external_img(paste0(outname, '.PCA2-filter.jpg'), width = 9, height = 5), 
#           location = ph_location_type(type = "body"))


jpeg(paste0(outname, '.avePCA2.jpg'), width=1800, height=1800, res=300)
nudge_value_x <- (max(avePCAData$avePC1) - min(avePCAData$avePC1)) / 25
nudge_value_y <- (max(avePCAData$avePC2) - min(avePCAData$avePC2)) / 25 
ggplot(avePCAData, aes(x = avePC1, y = avePC2, color = Treat)) +
  geom_point(size = 3) +
  scale_color_manual(values = custom_colors) +  # Apply custom colors
  #xlab(pc1_label) +
  #ylab(pc2_label) +
  coord_fixed() +
  theme_bw() +
  theme(legend.position = "top") +  # Position legend at the top
  stat_chull(aes(color = Treat), alpha = 0.05, geom = "polygon") +  # Convex hull for groups
  geom_line(aes(group = Treat), linetype = "solid") +  # Add a line based on 'trt'
  geom_text(data = avePCAData, aes(label = Treat, x = avePC1, y = avePC2), size=2,
            nudge_x = nudge_value_x, nudge_y = - nudge_value_y ) +  # Further nudge the labels
  ggtitle(paste(outname, "PCA"))
dev.off()

doc <- doc %>%
  add_slide(layout = "Title and Content", master = "Office Theme") %>%
  ph_with(value = paste0("AvePCA2 ", filename),
          location = ph_location_type(type = "title")) %>%
  ph_with(external_img(paste0(outname, '.avePCA2.jpg'), width = 9, height = 5), 
          location = ph_location_type(type = "body"))


# jpeg(paste0(outname, '.avePCA2_filtered.jpg'), width=1800, height=1800, res=300)
# #nudge_value_x <- (max(avePCAData$avePC1) - min(avePCAData$avePC1)) / 25
# #nudge_value_y <- (max(avePCAData$avePC2) - min(avePCAData$avePC2)) / 25 
# ggplot(avePCAData_filtered, aes(x = avePC1, y = avePC2, color = Treat)) +
#   geom_point(size = 3) +
#   scale_color_manual(values = custom_colors) +  # Apply custom colors
#   #xlab(pc1_label) +
#   #ylab(pc2_label) +
#   coord_fixed() +
#   theme_bw() +
#   theme(legend.position = "top") +  # Position legend at the top
#   stat_chull(aes(color = Treat), alpha = 0.05, geom = "polygon") +  # Convex hull for groups
#   geom_line(aes(group = Treat), linetype = "solid") +  # Add a line based on 'trt'
#   geom_text(data = avePCAData_filtered, aes(label = Treat, x = avePC1, y = avePC2), size=2,
#             nudge_x = 0, nudge_y = - 0.05 ) +  # Further nudge the labels
#   ggtitle(paste(outname, "PCA"))
# dev.off()
# 
# doc <- doc %>%
#   add_slide(layout = "Title and Content", master = "Office Theme") %>%
#   ph_with(value = paste0("AvePCA2 filtered ", filename),
#           location = ph_location_type(type = "title")) %>%
#   ph_with(external_img(paste0(outname, '.avePCA2_filtered.jpg'), width = 9, height = 5), 
#           location = ph_location_type(type = "body"))
# 


library(geosphere)
distance_dict <- list()
for (i in 1:dim(avePCAData)[1]) {
  myDist<-distm (avePCAData[i,2:3], reference_point, fun = distGeo)
  #myDist <- rdist(as.matrix(result$rotation[i,1:2, drop = FALSE]), as.matrix(reference_point))
  myName=as.character(avePCAData$trt[i])
  distance_dict[myName] <- myDist
}
pca_similarityScore=calculateSimilarityScores(distance_dict)
myResult=cbind(names(distance_dict),unlist(distance_dict), pca_similarityScore )
colnames(myResult)=c("Sample", "Distance", "similarityScore")
write.table(myResult, file = outname, sep = "\t", row.names = FALSE)
drawSortedBarPlot(myResult, outname)

################################################
#devtools::install_github("willtownes/glmpca") ===>not work
####################################################
# library("glmpca")
# gpca <- glmpca(counts(dds), L=2)
# gpca.dat <- gpca$factors
# jpeg('PCA3.jpg')
# 
# ggplot(gpca.dat, aes(x = dim1, y = dim2, color = Treat)) +
#   geom_point(size =3) + coord_fixed()  +  ggtitle("glmpca - Generalized PCA")
# dev.off()


################################################
##MDS plot
############################################
mds <- as.data.frame(colData(vsd))  %>%
         cbind(cmdscale(sampleDistMatrix))
jpeg('MDS.jpg')
ggplot(mds, aes(x = `1`, y = `2`, color = Treat)) +
geom_point(size = 3) + coord_fixed()+ ggtitle("MDS with VST data")
dev.off()
doc <- doc %>%
  add_slide(layout = "Title and Content", master = "Office Theme") %>%
  ph_with(value = paste0("MDS ", filename),
          location = ph_location_type(type = "title")) %>%
  ph_with(external_img('MDS.jpg', width = 9, height = 5), 
          location = ph_location_type(type = "body"))



####################################################
mdsPois <- as.data.frame(colData(dds)) %>%
   cbind(cmdscale(samplePoisDistMatrix))
jpeg('MDS2.jpg')
ggplot(mdsPois, aes(x = `1`, y = `2`, color = Treat)) +
geom_point(size = 3) + coord_fixed()+ ggtitle("MDS with PoissonDistances")
dev.off()

doc <- doc %>%
  add_slide(layout = "Title and Content", master = "Office Theme") %>%
  ph_with(value = paste0("MDS2 ", filename),
          location = ph_location_type(type = "title")) %>%
  ph_with(external_img('MDS2.jpg', width = 9, height = 5), 
          location = ph_location_type(type = "body"))

#############################
results <- DESeq(dds)
SigText="";
for (i in trt_groups) {
  if(i != "NCColgCP"){
    comparison1 <- results(results, pAdjustMethod = "BH", contrast = c("Treat",i,"NCColgCP"))
    comparison1a <- cbind(rownames(dds),comparison1)
    colnames(comparison1a)[1] <- "Genes"
    XX=2**comparison1a$log2FoldChange
    comparison1a$TrueFC <- sapply(XX,function(x){ifelse(x>=1, x, -1/x)})
    YY=cbind(rownames(dds), comparison1a$pvalue, comparison1a$TrueFC)
    write.table(YY, file = paste0(i,"_vs_NCColgCP.stat"), sep="\t")
    write.table(comparison1a, file = paste0(i, "_vs_NCColgCP.txt"), sep="\t")
    SigNumInfo= paste0(i, "_vs_NCColgCP", " p<=0.05 ",  sum(comparison1a$pvalue<=0.05, na.rm=TRUE), " padj<=0.05 ", sum(comparison1a$padj<=0.05, na.rm=TRUE))
    SigText=paste(SigText, SigNumInfo, "\n")
    print (SigNumInfo)
  }
}
doc <- doc %>%
  add_slide(layout = "Title and Content", master = "Office Theme") %>%
  ph_with(value = paste0("Significant Count ", filename),
          location = ph_location_type(type = "title")) %>%
  # Add the caption text under the plot
  ph_with(fpar(ftext(SigText, prop = fp_text(font.size = 10))), 
          location = ph_location(left = 0.5, top = 6.3, width = 5.5, height = 0.5))

SigText="";
for (i in trt_groups) {
  if(i != "NoTreat"){
    comparison1 <- results(results, pAdjustMethod = "BH", contrast = c("Treat",i,"NoTreat"))
    comparison1a <- cbind(rownames(dds),comparison1)
    colnames(comparison1a)[1] <- "Genes"
    XX=2**comparison1a$log2FoldChange
    comparison1a$TrueFC <- sapply(XX,function(x){ifelse(x>=1, x, -1/x)})
    YY=cbind(rownames(dds), comparison1a$pvalue, comparison1a$TrueFC)
    write.table(YY, file = paste0(i,"_vs_NoTreat.stat"), sep="\t")
    write.table(comparison1a, file = paste0(i, "_vs_NoTreat.txt"), sep="\t")
    SigNumInfo= paste0(i, "_vs_NoTreat", " p<=0.05 ",  sum(comparison1a$pvalue<=0.05, na.rm=TRUE), " padj<=0.05 ", sum(comparison1a$padj<=0.05, na.rm=TRUE))
    SigText=paste(SigText, SigNumInfo, "\n")
    print (SigNumInfo)
  }
}


    comparison1 <- results(results, pAdjustMethod = "BH", contrast = c("Treat","Formula2","Formula1"))
    comparison1a <- cbind(rownames(dds),comparison1)
    colnames(comparison1a)[1] <- "Genes"
    XX=2**comparison1a$log2FoldChange
    comparison1a$TrueFC <- sapply(XX,function(x){ifelse(x>=1, x, -1/x)})
    YY=cbind(rownames(dds), comparison1a$pvalue, comparison1a$TrueFC)
    write.table(YY, file = paste0(i,"_vs_NCColgCP.stat"), sep="\t")
    write.table(comparison1a, file = paste0("Formula2", "_vs_Formula1.txt"), sep="\t")
    SigNumInfo= paste0("Formula2_vs_Formula1", " p<=0.05 ",  sum(comparison1a$pvalue<=0.05, na.rm=TRUE), " padj<=0.05 ", sum(comparison1a$padj<=0.05, na.rm=TRUE))
    SigText=paste(SigText, SigNumInfo, "\n")
    print (SigNumInfo)

doc <- doc %>%
  add_slide(layout = "Title and Content", master = "Office Theme") %>%
  ph_with(value = paste0("Significant Count ", filename),
          location = ph_location_type(type = "title")) %>%
  # Add the caption text under the plot
  ph_with(fpar(ftext(SigText, prop = fp_text(font.size = 10))), 
          location = ph_location(left = 0.5, top = 6.3, width = 5.5, height = 0.5))


print(doc, target = paste0(filename, ".pptx"))
sessionInfo()