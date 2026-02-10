#############################################
#Need to add in the statistical analysis
#Need to add the ranking
#Need to add statistical figure into pptx
#############################################
rm(list=ls())

library(dplyr)
library(ggplot2)
library(ggpubr)
library(openxlsx)
#library(cowplot)

args <- commandArgs(trailingOnly = TRUE)
print(args)
file <- args[1]

#rm(args)

#file <-"GSS3225_compass_ddct.txt"
file <-"GSS3240_ddCt"
A <- read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
d <- dim(A)
df <- A[,2:d[2]]
rownames(df) <-A[,1]
housekeeping <- c("ACTB", "GAPDH", "PPIA", "B2M")
##########Now start to check different sets######################
inflammation <- c("DDIT4","dnajb9","DSG3","HBEGF","PLAUR","pmaip1","SESN2","slc30a1","stc2","TRIB3")
TissueDamage <- c("CLDN1","cxcl8","mmp10","MMP3","mt1g", "ddit3", "il13ra2", "MMP1")
Redox <- c("AKR1C1","akr1c2","HMOX1","NQO1", "osgin1", "slc7a11", "SQSTM1")
USS<- c( "DDIT4","dnajb9","DSG3","HBEGF","PLAUR","pmaip1","SESN2","slc30a1",
         "stc2","TRIB3","akr1c2","HMOX1","NQO1","osgin1","slc7a11","SQSTM1","ddit3","il13ra2",
         "MMP1","ARRDC3","asns","BCL3","btg1","cbs","chac1","fth1","GCLC","gclm",
         "gpt2","HERPUD1","KRT16","mt1x","NDRG1","pck2","pnrc1","psat1","psph","slc1a4",
         "smim14","srxn1","trim16","txnrd1","UPP1","VEGFA"
)
All <- c(
  "abca12",  "AKR1C1", "angptl4", "ARNTL2", "ARRDC3", "asns", "ass1",
  "atf3",  "BCL3", "btg1", "CCL20", "cebpb", "cers3", "chac1", "CLDN1",
  "csta", "CXCL1", "cxcl14", "cxcl8", "ddit3", "DDIT4", "dnajb9", "DSG3", "dst",
  "dusp1", "dusp10", "eppk1", "FAT4", "fgfr2", "FRAT1", "fth1", "gadd45a",
  "GCLC", "gclm", "gdf15", "gjb2", "gpt2", "HBEGF", "HERPUD1",
  "hist1h2ac", "HMOX1", "hspa1a", "hspa1b", "IL1a", "il1rl1", "il23a", "il6r",
  "irak2", "IRF7", "KRT16", "KRT6B", "KRTAP2-3", "lamp3", "MMP1", "mmp10",
  "MMP3", "mt1g", "mt1x", "NCF2", "NDRG1", "NQO1", "osgin1", "pck2", "PLAU",
  "PLAUR", "pmaip1", "pnrc1",  "PPP1R15A", "psat1", "psph", "PTGS2",
  "scnn1a", "serpine1", "serpine2", "slc1a4", "slc30a1", "slc7a11", "smim14",
  "sod2", "sprr1b", "SQSTM1", "srxn1", "stc2", "TGM1", "TIMP3", "tnc", "TNIP2",
  "TPBG", "TRIB3", "trim16", "UPP1", "VEGFA"
)
# Combine all other gene sets
other_genes <- unique(c(inflammation, TissueDamage, Redox, USS))
Rest <- setdiff(All, other_genes)
NewDesign<-c(
  "serpine1","ACTB","TGM1","sod2","VEGFA","MMP3","PPIA","CLDN1","MMP1","ARRDC3",
  "HERPUD1","UPP1","DSG3","KRT16","NDRG1","AKR1C1","sprr1b","TRIB3","SQSTM1",
  "DDIT4","angptl4","PLAUR","CCL20","B2M","cers3","HMOX1","HBEGF","TIMP3",
  "PPP1R15A","PTGS2","PLAU","CXCL1","gadd45a","serpine2","GAPDH","ARNTL2",
  "abca12","gjb2","NQO1","GCLC","asns","ass1","atf3","BCL3","cebpb","chac1",
  "csta","cxcl14","cxcl8","ddit3","dnajb9","dusp1","dusp10","eppk1","fth1",
  "gclm","gdf15","gpt2","hist1h2ac","hspa1a","hspa1b","IL1a","il1rl1","il23a",
  "irak2","KRT6B","lamp3","mmp10","mt1x","NCF2","osgin1","pck2","pnrc1","psat1",
  "scnn1a","slc30a1","slc7a11","srxn1","tnc","TPBG","trim16","AREG","C1orf68",
  "CDC42BPG","FABP5","GRHL1","IL17RE","LCE1A","LCE1E","NFE2","RNA28S5",
  "SLC7A5","SLURP1","TEAD1","TEAD4","UGT3A2"
)
NewAdd <- c(
  "ACTB","PPIA","B2M","GAPDH","AREG","C1orf68","CDC42BPG","FABP5","GRHL1",
  "IL17RE","LCE1A","LCE1E","NFE2","RNA28S5","SLC7A5","SLURP1","TEAD1","TEAD4","UGT3A2"
)

gene_sets <- list(Inflammation = inflammation, TissueDamage = TissueDamage, Redox = Redox, USS = USS, All=All, Rest=Rest, NewDesign=NewDesign, NewAdd=NewAdd)


# Remove housekeeping genes by filtering out their row names
ddct_matrix_all <- df[!(rownames(df) %in% housekeeping),]
nrow(df)           # Before
nrow(ddct_matrix_all)  # After
orgID=colnames(A)[2:d[2]]

splitname<-strsplit(orgID, "[._]")
Clen=length(splitname)

trt=rep("NA", Clen)
dup=rep("NA", Clen)
Time=rep("NA", Clen)
tType=rep("NA", Clen)
for(mm in  1:Clen ){
 
  trt[mm]=splitname[[mm]][1]
  dup[mm]=splitname[[mm]][4]
  Time[mm]=splitname[[mm]][2]
  tType[mm]=splitname[[mm]][3]
}

SID=paste0(tType, ".", trt, ".", Time,  ".", dup)
TypeTrtTime=paste0(tType, ".", trt, ".", Time)
colnames(ddct_matrix_all)=SID
meta<- data.frame(orgID,Time,  dup, trt, SID, tType, TypeTrtTime)
flatten_dist_mat <- function(distance_mat) {
  myDist <- unlist(distance_mat)
  # Find minimum and maximum distances
  min_dist <- min(myDist)
  max_dist <- max(myDist)
  # Calculate similarity scores based on distances
  similarityScore <-  1 + 2 * (min_dist - myDist) / (max_dist - min_dist)
  return(similarityScore)
}

TimePoints <- sort(unique(as.numeric(meta$Time)))
RefTime <- TimePoints[1]

library(officer)
library(rvg)
# Create PowerPoint

####################Now Let's Work on the PCA and Ranking of the material#######################
###############################################
for (nm in names(gene_sets)) {
  #nm="All"
  doc <- read_pptx()
  ddct_matrix <- ddct_matrix_all[(rownames(ddct_matrix_all) %in% gene_sets[[nm]]), ]
  print(nm)
  print(nrow(df))           # Before
  print(nrow(ddct_matrix))  # After
  pca_result <- prcomp(t(ddct_matrix), center = TRUE, scale. = TRUE)
  loadings <- pca_result$scale
  write.csv(loadings, file = paste0(file,".",nm, ".ddct_PCAGeneImportance.csv") , row.names = TRUE)
  pca_df <- as.data.frame(pca_result$x)
  pca_df$SID <- rownames(pca_df)
  combined_df <- merge(pca_df, meta, by = "SID", all.x = TRUE)
  # Create a data frame with the PCA results
  write.csv(combined_df, file = paste0(file,".", nm, ".ddct_PCA.csv") , row.names = TRUE)
  variance_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)
  pc1_label <- paste0("PC1 (", round(variance_explained[1] * 100, 2), "%)")
  pc2_label <- paste0("PC2 (", round(variance_explained[2] * 100, 2), "%)")
  
  avePCAData <- combined_df %>%
    group_by(trt, Time, tType, TypeTrtTime) %>%
    summarise(
      avePC1 = mean(PC1, na.rm = TRUE),  # Calculate the average of PC1
      avePC2 = mean(PC2, na.rm = TRUE),   # Calculate the average of PC2
      se_PC1 = sd(PC1, na.rm = TRUE) / sqrt(n()),
      se_PC2 = sd(PC2, na.rm = TRUE) / sqrt(n())
    ) %>%
    ungroup() 
  
  avePCAData$Time <- as.numeric(avePCAData$Time)
 
  #####
  p<-ggplot(avePCAData, aes(x = avePC1, y = avePC2, color = trt, shape = trt)) +
    geom_point(size = 3) +
    xlab(pc1_label) +
    ylab(pc2_label) +
    coord_fixed() +
    theme_bw()  + 
    theme(legend.position = "top") +  
    geom_line(aes(group = trt), linetype = "solid") + 
    geom_text(data = avePCAData, aes(label = TypeTrtTime, x = avePC1, y = avePC2), size=2)+#,
    ggtitle(paste(file, nm, " average PCA"))
  
  img_path <- paste0(file, ".", nm, ".avePCA2.jpg")
  
  jpeg(img_path, width = 2000, height = 2600, res = 300)
  print(p)
  dev.off()
  
  
  doc <- doc %>%
    add_slide(layout = "Title and Content", master = "Office Theme") %>%
    ph_with(value = paste0(file, ".", nm," Average PCA"),
            location = ph_location_type(type = "title")) %>%
    ph_with(
      external_img(img_path, width = 8, height = 4.5),  # size in inches
      location = ph_location(left = 0.5, top = 1.8)
    )
  #############################################################################
  
  plot5 <- ggplot(avePCAData, aes(x = Time, y = avePC1, color = trt, group = trt)) +
    geom_line(linewidth = 1) +
    geom_point(size = 6) +
    geom_errorbar(aes(ymin = avePC1 - se_PC1, ymax = avePC1 + se_PC1), width = 0.2) +
    labs(x = "Time", y = "Mean PCA PC1", title = paste(nm, "PC1:", pc1_label)) +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 7),        # shrink legend text
      legend.title = element_text(size = 6),       # shrink legend title (optional)
      panel.background = element_blank(),  # Remove panel background
      plot.background = element_blank(),  # Remove plot background
      panel.grid.major = element_blank(),  # Remove major grid lines
      panel.grid.minor = element_blank()   # Remove minor grid lines
      
    ) 
  
  
  jpeg(paste0(file, ".", nm, '.PC1.jpg'), width=1800, height=1800, res=300)
  print(plot5)
  dev.off()
  doc <- doc %>%
    add_slide(layout = "Title and Content", master = "Office Theme") %>%
    ph_with(value = paste0(file," ", nm, " PC1"),
            location = ph_location_type(type = "title")) %>%
    ph_with(
      external_img(paste0(file, ".", nm, '.PC1.jpg'), width = 8, height = 4.5),  # size in inches
      location = ph_location(left = 0.5, top = 1.8)
    )
  
  plot6 <- ggplot(avePCAData, aes(x = Time, y = avePC2, color = trt, group = trt)) +
    geom_line(linewidth = 1) +  # Draw line for each trt
    geom_point(size=6) +  # Draw points for each data point
    geom_errorbar(aes(ymin = avePC2 - se_PC2, ymax = avePC2 + se_PC2), width = 0.2) +  # Error bars
    labs(x = "Time", y = "Mean PCA PC2 ", title = paste(nm, "PC2:", pc2_label)) +
    theme(    legend.position = "bottom",
              legend.text = element_text(size = 7),        # shrink legend text
              legend.title = element_text(size = 6),  
              panel.background = element_blank(),  # Remove panel background
              plot.background = element_blank(),  # Remove plot background
              panel.grid.minor = element_blank()   # Remove minor grid lines
    ) 
  jpeg(paste0(file, ".", nm, '.PC2.jpg'), width=1800, height=1800, res=300)
  print(plot6)
  dev.off()
  
  doc <- doc %>%
    add_slide(layout = "Title and Content", master = "Office Theme") %>%
    ph_with(value = paste0(file, " ",nm,  " PC2"),
            location = ph_location_type(type = "title")) %>%
    ph_with(
      external_img(paste0(file, ".", nm,'.PC2.jpg'), width = 8, height = 4.5),  # size in inches
      location = ph_location(left = 0.5, top = 1.8)
    )
  
  p1<-ggplot(pca_df, aes(x = PC1, y = PC2, color = trt, shape=trt)) +
    geom_point(size = 3) +
    xlab(pc1_label) +
    ylab(pc2_label) +
    coord_fixed() +
    theme_bw()  +
    theme(legend.position = "bottom") +  # Position legend at the bottom
    stat_chull(aes(color = trt, fill = trt), alpha = 0.05, geom = "polygon") +
    geom_text(data = combined_df, aes(label = TypeTrtTime, x = PC1, y = PC2), size=2 ) + ggtitle(paste(file, nm, " PCA"))
  
  img_path <- paste0(file,".", nm,  ".PCA2.jpg")
  
  jpeg(img_path, width = 2000, height = 2600, res = 300)
  print(p1)
  dev.off()
  doc <- doc %>%
    add_slide(layout = "Title and Content", master = "Office Theme") %>%
    ph_with(value = paste(file,".", nm, " PCA2"),
            location = ph_location_type(type = "title")) %>%
    ph_with(
      external_img(img_path, width = 8, height = 4.5),  # size in inches
      location = ph_location(left = 0.5, top = 1.8)
    )
  
  ##########################################################
  for (myTypeTrtTime in unique(combined_df$TypeTrtTime)){
    
    parts <- strsplit(myTypeTrtTime, "\\.")[[1]]
    
    mytType <- parts[1]
    myTrt   <- parts[2]
    myTime  <- parts[3]
    reference_point <- combined_df %>%
      filter(TypeTrtTime == myTypeTrtTime) %>%  # Filter for control group
      summarise(
        avePC1 = mean(PC1),  # Calculate average of PC1
        avePC2 = mean(PC2)   # Calculate average of PC2
      )
    for (i in 1:dim(avePCAData)[1]) {
      point1<-avePCAData[i,c("avePC1","avePC2")]
      myDist <- sqrt(sum((point1 - reference_point)^2)) ###Eucledian distance
      avePCAData$EucDistance[i] <-myDist
    }
    avePCAData$similarityScore <- flatten_dist_mat(avePCAData$EucDistance)
    avePCAData <- avePCAData[order(avePCAData$similarityScore, decreasing = TRUE), ]
    
    for (i in 1:dim(combined_df)[1]) {
      point1<-pca_df[i,c("PC1","PC2")]
      myDist <- sqrt(sum((point1 - reference_point)^2)) ###Eucledian distance
      combined_df$EucDistance[i] <-myDist
    }
    combined_df$similarityScore <- flatten_dist_mat(combined_df$EucDistance)
    if(myTrt == "Untreated"){
    p <- ggplot(avePCAData, aes(y = reorder(TypeTrtTime, similarityScore), x = similarityScore)) +
      geom_bar(stat = "identity", aes(fill = Time)#, color = "black"
      ) +
      labs(y = "", x = "SimilarityScore", title = paste0(nm, ": Rank to ", myTypeTrtTime) )+
      theme(
        legend.position = "none",  # Remove the legend if not needed
        panel.background = element_rect(fill = "white", color = NA),  # White background
        plot.background = element_rect(fill = "white", color = NA),  # White plot area background
        panel.grid = element_blank(),  # Remove grid lines
        axis.line = element_line(color = "black")  # Add black axis lines
      )+
      ggtitle(paste(file, "\n",nm,  myTypeTrtTime, " Rank"))
    # 
    # Save the plot to a file
    file_name <- paste0(file,".", nm, ".", myTypeTrtTime, ".jpg")
    #ggsave(file_name, plot = p, width = 1000, height = 800, res=300)
    jpeg(file_name, width=2800, height=1600, res=300)
    print(p)
    dev.off()
    
    doc <- doc %>%
      add_slide(layout = "Title and Content", master = "Office Theme") %>%
      ph_with(value = paste(file,nm,  myTypeTrtTime, "Rank"), location = ph_location_type(type = "title")) %>%
      ph_with(dml(ggobj = p), location = ph_location(left = 0.5, top = 1.8, width = 8, height = 4.5))
    }
    names(avePCAData)[names(avePCAData) == "similarityScore"] <- paste0("SimilarityScore.",nm, ".",  myTypeTrtTime)
    names(avePCAData)[names(avePCAData) == "EucDistance"] <- paste0("EucDistance.", nm, ".",myTypeTrtTime)
    names(pca_df)[names(pca_df) == "similarityScore"] <- paste0("SimilarityScore.",nm, ".", myTypeTrtTime)
    names(pca_df)[names(pca_df) == "EucDistance"] <- paste0("EucDistance.", nm, ".",myTypeTrtTime)
  }
  write.xlsx(avePCAData,paste0(file,".", nm, ".rank.xlsx"))
  write.xlsx(pca_df,paste0(file,".", nm,".individual.rank.xlsx"))
  print(doc, target = paste0(file, ".", nm,".pptx"))
  

}

