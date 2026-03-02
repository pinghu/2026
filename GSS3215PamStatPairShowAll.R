##########################################
#Ping Hu April 16, 2025
### Would like to add into the PPT summary
#### Need to filter out the summary graph with the %zero to be less than 75%
##########################################
rm(list=ls())
library(dplyr)
library(ggplot2)
library(ggpubr)

my.wilcox.p.value <- function(...) {
  obj<-try(wilcox.test(...), silent=TRUE)
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

filename0="GSS3215Pam_meta_405Trans.txt.in.paired"
m0<-read.table(filename0, sep="\t", header=TRUE)
d <- dim(m0);
mB=m0[1:d[1], 2:d[2]]
rownames(mB)=m0[,1]
Ref_BrownSpotMelanin <- as.numeric(
  mB["BaselineSpotTracking.brownspot.SpotMelaMean.Sum", ]
)


args <- commandArgs(trailingOnly = TRUE)
#filename1="GSS3215Pam_meta.txt"
#X<-read.table(filename1, sep="\t", header=TRUE)
filename=args[1]
#filename="relab7.paired.data"

A0<-read.table(filename, sep="\t", header=TRUE)
d <- dim(A0);
B=A0[1:d[1], 2:d[2]]
rownames(B)=A0[,1]

ZZ=as.numeric(min(B[B>0&!is.na(B)]))/100
C=B+ZZ
rownames(C)=A0[,1]

Cname=colnames(A0)[2:d[2]]
Clen=length(Cname) ##there are 3 annotation columns
Visit=rep("NA", Clen)
Treat=rep("NA", Clen)
SID=rep("NA", Clen)
Side=rep("NA", Clen)
V0spot=rep("NA", Clen)
Fitz=rep("NA", Clen)
Age=rep("NA", Clen)
Race=rep("NA", Clen)
OldID=rep("NA", Clen)
splitname<-strsplit(Cname, "[.]")

for(mm in  1:Clen ){
  Visit[mm]=splitname[[mm]][1]
  Treat[mm]=splitname[[mm]][2]
  SID[mm]=splitname[[mm]][3]
  Side[mm]=splitname[[mm]][4]
  V0spot[mm]=splitname[[mm]][5]
  Fitz[mm]=splitname[[mm]][6]
  Age[mm]=splitname[[mm]][7]
  Race[mm]=splitname[[mm]][8]
  OldID[mm]=splitname[[mm]][9]
}

shortID=paste0(Visit,Treat,SID)
mydata <- data.frame(
  Visit, SID, Treat, Side, Fitz,
  Age, Race, OldID, V0spot,Cname,shortID
)
mydata$Age    <- as.numeric(mydata$Age)
mydata$V0spot <- as.numeric(mydata$V0spot)
mydata$VisitTreat =paste0(mydata$Treat, mydata$Visit)
order_levels <- c("ABL", "BBL","CBL","DBL", "EBL","AW8","BW8","CW8", "DW8", "EW8") 
mydata$VisitTreat <- factor(mydata$VisitTreat, levels = order_levels, ordered = TRUE)

library(officer)
library(rvg)
# Create PowerPoint
doc <- read_pptx()

sanitize_filename <- function(x) {
  x <- gsub("%", "pct", x, fixed = TRUE)        # replace % with 'pct'
  x <- gsub("[\\/:*?\"<>|]", "_", x)            # remove illegal path chars
  x <- gsub("\\s+", "_", x)                     # optional: replace spaces
  x
}

names_vec <- c(
  "Species", "Measurement","percentZeroNA",
  "spearman_rho.brownspotMelanin","pearson_rho.brownspotMelanin","spearman_pval.brownspotMelanin","pearson_pval.brownspotMelanin",
  "spearman_rho.V0Spot","pearson_rho.V0Spot","spearman_pval.V0Spot","pearson_pval.V0Spot",
  "KP_V0spot_BL","KP_SID","KP_SID_BL","KP_Side","KP_Side_BL",
  "KP_Treat","KP_Treat_BL","KP_Treat_W8",
  "KP_Visit","KP_VisitTreat",
  "KP_Age_BL","KP_Fitz_BL","KP_Race_BL",
  "pW8A_BLA","pW8B_BLB","pW8B_W8A","pBLB_BLA", "pW8C_BLC","pW8C_W8A", "pBLC_BLA",
  "pW8D_BLD","pW8D_W8A", "pBLD_BLA", "pW8E_BLE","pW8E_W8A", "pBLE_BLA",
  "meanAllOriginal","mBLA","mBLB","mBLC","mBLD","mBLE",
   "mW8A","mW8B","mW8C","mW8D","mW8E",
  "tfc_A_W8vBL","tfc_B_W8vBL","tfc_BvA_W8", "tfc_BvA_BL",
  "tfc_C_W8vBL","tfc_CvA_W8","tfc_CvA_BL", "tfc_D_W8vBL","tfc_DvA_W8","tfc_DvA_BL", "tfc_E_W8vBL","tfc_EvA_W8", "tfc_EvA_BL",
  "W8-BL differences", "pdB_dA","pdC_dA","pdD_dA","pdE_dA","mddB","mddC","mddD","mddE","mdA","mdB","mdC","mdD","mdE"
  #,
 # "Baseline_adjust_all", "p_BA_all", "p_CA_all", "p_DA_all", "p_EA_all",
  #"Baseline_adjust_individual", "p_BA", "p_CA", "p_DA", "p_EA"
)
print(paste(names_vec, collapse = ", "))
p_cutoff = 0.1
for (i in 1:d[1]){
  genename=A0[i,1]
  species <- sub(".*\\.(s__[^.]+)$", "\\1", genename)
  gene_original <-as.numeric(B[i,])
  gene_relative <-as.numeric(C[i,])
  if (length(unique(gene_original)) <= 1 || all(is.na(gene_relative))) {
    message(paste("Skipping", genename, "- no variation in relative_gene"))
    next  # skip to next iteration
  }
  
  meanAllRelative <-mean(gene_relative)
  meanAllOriginal <-mean(gene_original)
  
  percentZeroNA=(sum(B[i,]==0)+sum(is.na(B[i,])))/(d[2]-1)
  mydata0=data.frame(gene_relative,gene_original,  Cname, Ref_BrownSpotMelanin)
  mydata2 <- inner_join(mydata0, mydata, by = c("Cname" = "Cname"))
  #mydata2 <-inner_join(mydata1, X, by = c("Cname" = "ID"))
  #mydata2 <- mydata2[mydata2$gene_original != "NA", ]
  
  mydataBL <- mydata2 %>% filter(Visit == "BL")
  mydataW8 <- mydata2 %>% filter(Visit == "W8")
  yBL <- mydataBL$gene_relative
  # Correlations at BL: V0spot vs metric
  ct_s <- cor.test(mydataBL$Ref_BrownSpotMelanin, yBL, method = "spearman", exact = FALSE)
  ct_p <- cor.test(mydataBL$Ref_BrownSpotMelanin, yBL, method = "pearson",  exact = FALSE)
  #ct <- cor.test(mydataBL$Ref_BrownSpotMelanin, mydataBL$V0spot, method = "pearson",  exact = FALSE)
  ct0_s <- cor.test(mydataBL$V0spot, yBL, method = "spearman", exact = FALSE)
  ct0_p <- cor.test(mydataBL$V0spot, yBL, method = "pearson",  exact = FALSE)
  spearman_rho  <- sprintf("%.02f",unname(ct_s$estimate))
  spearman_pval <- sprintf("%.02f",ct_s$p.value)
  pearson_rho   <- sprintf("%.02f",unname(ct_p$estimate))
  pearson_pval  <- sprintf("%.02f",ct_p$p.value)
  
  spearman0_rho  <- sprintf("%.02f",unname(ct0_s$estimate))
  spearman0_pval <- sprintf("%.02f",ct0_s$p.value)
  pearson0_rho   <- sprintf("%.02f",unname(ct0_p$estimate))
  pearson0_pval  <- sprintf("%.02f",ct0_p$p.value)
  
  # Kruskal tests (use reformulate)
  KP_Visit       <- kruskal.test(gene_relative ~ Visit, data = mydata2)$p.value
  KP_Treat       <- kruskal.test(gene_relative ~ Treat, data = mydata2)$p.value
  KP_Side        <- kruskal.test(gene_relative ~ Side, data = mydata2)$p.value
  KP_SID         <- kruskal.test(gene_relative ~ SID, data = mydata2)$p.value
  KP_VisitTreat  <- kruskal.test(gene_relative ~ VisitTreat, data = mydata2)$p.value
  
  KP_Treat_BL    <- kruskal.test(gene_relative ~ Treat, data = mydataBL)$p.value
  KP_Treat_W8    <- kruskal.test(gene_relative ~ Treat, data = mydataW8)$p.value
  
  KP_V0spot_BL   <- kruskal.test(gene_relative ~ V0spot, data = mydataBL)$p.value
  KP_Fitz_BL     <- kruskal.test(gene_relative ~ Fitz, data = mydataBL)$p.value
  KP_Age_BL      <- kruskal.test(gene_relative ~ Age, data = mydataBL)$p.value
  KP_Race_BL     <- kruskal.test(gene_relative ~ Race, data = mydataBL)$p.value
  KP_SID_BL      <- kruskal.test(gene_relative ~ SID, data = mydataBL)$p.value
  KP_Side_BL     <- kruskal.test(gene_relative ~ Side, data = mydataBL)$p.value
  # Convenience getter
  get_group <- function(label) mydata2$gene_relative[mydata2$VisitTreat == label]
  
  gBLA <- get_group("ABL"); gBLB <- get_group("BBL"); gBLC <- get_group("CBL"); gBLD <- get_group("DBL"); gBLE <- get_group("EBL")
  gW8A <- get_group("AW8"); gW8B <- get_group("BW8"); gW8C <- get_group("CW8"); gW8D <- get_group("DW8"); gW8E <- get_group("EW8")
  ##########################Delta anlysis for baseline adjusted treatment effect
  gdA<-gW8A-gBLA;gdB<-gW8B-gBLB;gdC<-gW8C-gBLC;gdD<-gW8D-gBLD; gdE<-gW8E-gBLE;
  mdA=mean(gdA);mdB=mean(gdB); mdC=mean(gdC); mdE=mean(gdE);mdD=mean(gdD);
  pdB_dA <- my.wilcox.p.value(gdA, gdB, paired = FALSE)
  pdC_dA <- my.wilcox.p.value(gdA, gdC, paired = FALSE)
  pdD_dA <- my.wilcox.p.value(gdA, gdD, paired = FALSE)
  pdE_dA <- my.wilcox.p.value(gdA, gdE, paired = FALSE)
  mddB=mdB-mdA; mddC=mdC-mdA; mddD=mdD-mdA; mddE=mdE-mdA;
  mindP=min(pdB_dA,pdC_dA,pdC_dA,pdE_dA,na.rm=TRUE)
  

  ###############################################################
  pW8A_BLA <- my.wilcox.p.value(gBLA, gW8A, paired = TRUE)
  pW8B_BLB <- my.wilcox.p.value(gBLB, gW8B, paired = TRUE)
  pW8C_BLC <- my.wilcox.p.value(gBLC, gW8C, paired = TRUE)
  pW8D_BLD <- my.wilcox.p.value(gBLD, gW8D, paired = TRUE)
  pW8E_BLE <- my.wilcox.p.value(gBLE, gW8E, paired = TRUE)
  pW8B_W8A <- my.wilcox.p.value(gW8A, gW8B, paired = FALSE)
  pW8C_W8A <- my.wilcox.p.value(gW8A, gW8C, paired = FALSE)
  pW8D_W8A <- my.wilcox.p.value(gW8A, gW8D, paired = FALSE)
  pW8E_W8A <- my.wilcox.p.value(gW8A, gW8E, paired = FALSE)
  pBLB_BLA <- my.wilcox.p.value(gBLA, gBLB, paired = FALSE)
  pBLC_BLA <- my.wilcox.p.value(gBLA, gBLC, paired = FALSE)
  pBLD_BLA <- my.wilcox.p.value(gBLA, gBLD, paired = FALSE)
  pBLE_BLA <- my.wilcox.p.value(gBLA, gBLE, paired = FALSE)
  
  # Means
  mBLA <- mean(gBLA, na.rm = TRUE); mBLB <- mean(gBLB, na.rm = TRUE); mBLC <- mean(gBLC, na.rm = TRUE)
  mBLD <- mean(gBLD, na.rm = TRUE); mBLE <- mean(gBLE, na.rm = TRUE)
  mW8A <- mean(gW8A, na.rm = TRUE); mW8B <- mean(gW8B, na.rm = TRUE); mW8C <- mean(gW8C, na.rm = TRUE)
  mW8D <- mean(gW8D, na.rm = TRUE); mW8E <- mean(gW8E, na.rm = TRUE)
  
  # Fold-changes (your truefc)
  tfc_A_W8vBL <- truefc(mW8A / mBLA)
  tfc_B_W8vBL <- truefc(mW8B / mBLB)
  tfc_C_W8vBL <- truefc(mW8C / mBLC)
  tfc_D_W8vBL <- truefc(mW8D / mBLD)
  tfc_E_W8vBL <- truefc(mW8E / mBLE)
  
  tfc_BvA_W8  <- truefc(mW8B / mW8A)
  tfc_CvA_W8  <- truefc(mW8C / mW8A)
  tfc_DvA_W8  <- truefc(mW8D / mW8A)
  tfc_EvA_W8  <- truefc(mW8E / mW8A)
  
  tfc_BvA_BL  <- truefc(mBLB / mBLA)
  tfc_CvA_BL  <- truefc(mBLC / mBLA)
  tfc_DvA_BL  <- truefc(mBLD / mBLA)
  tfc_EvA_BL  <- truefc(mBLE / mBLA)
  values_vec <- c(
    species, genename, percentZeroNA, 
    spearman_rho, pearson_rho, spearman_pval, pearson_pval,
    spearman0_rho, pearson0_rho, spearman0_pval, pearson0_pval,
    KP_V0spot_BL, KP_SID, KP_SID_BL, KP_Side, KP_Side_BL,
    KP_Treat, KP_Treat_BL, KP_Treat_W8,
    KP_Visit, KP_VisitTreat,
    KP_Age_BL, KP_Fitz_BL, KP_Race_BL,
    pW8A_BLA, pW8B_BLB, pW8B_W8A, pBLB_BLA, pW8C_BLC, pW8C_W8A, pBLC_BLA,
    pW8D_BLD, pW8D_W8A, pBLD_BLA, pW8E_BLE, pW8E_W8A, pBLE_BLA,
    meanAllOriginal,mBLA, mBLB, mBLC, mBLD, mBLE,
    mW8A, mW8B, mW8C, mW8D, mW8E,
    tfc_A_W8vBL, tfc_B_W8vBL, tfc_BvA_W8, tfc_BvA_BL,
    tfc_C_W8vBL, tfc_CvA_W8, tfc_CvA_BL, tfc_D_W8vBL, tfc_DvA_W8, tfc_DvA_BL, tfc_E_W8vBL, tfc_EvA_W8, tfc_EvA_BL,
    "W8-BL delta to Vehicle", pdB_dA,pdC_dA,pdD_dA,pdE_dA,mddB,mddC,mddD,mddE,mdA,mdB,mdC,mdD,mdE
    #,
   # "Baseline_adjust_all", p_BA_all, p_CA_all, p_DA_all, p_EA_all,
   # "Baseline_adjust", p_BA, p_CA, p_DA, p_EA
  )
  print(paste(values_vec, collapse = ", "))
  word_metric <- paste(paste0(names_vec, "=", values_vec), collapse = ", ")
  minP=min(pW8A_BLA, pW8B_BLB, pW8B_W8A, pW8C_BLC, pW8C_W8A,
           pW8D_BLD, pW8D_W8A, pW8E_BLE, pW8E_W8A,
           spearman_pval, pearson_pval,
           spearman0_pval, pearson0_pval,
           KP_V0spot_BL,pdB_dA,pdC_dA,pdC_dA,pdE_dA, na.rm=TRUE)
  
  #print(word_metric)
  # stat tibble for annotation
  
  if(minP>p_cutoff){next}
  round_p = 2
  stat_metric <- tribble(
    ~.y.,     ~group1, ~group2, ~p,
    gene_relative,   "ABL",   "AW8",   pW8A_BLA,
    gene_relative,   "BBL",   "BW8",   pW8B_BLB,
    gene_relative,   "CBL",   "CW8",   pW8C_BLC,
    gene_relative,   "DBL",   "DW8",   pW8D_BLD,
    gene_relative,   "EBL",   "EW8",   pW8E_BLE,
    gene_relative,   "AW8",   "BW8",   pW8B_W8A,
    gene_relative,   "AW8",   "CW8",   pW8C_W8A,
    gene_relative,   "AW8",   "DW8",   pW8D_W8A,
    gene_relative,   "AW8",   "EW8",   pW8E_W8A,
    gene_relative,   "ABL",   "BBL",   pBLB_BLA,
    gene_relative,   "ABL",   "CBL",   pBLC_BLA,
    gene_relative,   "ABL",   "DBL",   pBLD_BLA,
    gene_relative,   "ABL",   "EBL",   pBLE_BLA
  ) %>%
    mutate(p = round(as.numeric(p), round_p)) %>%
    filter(p <= p_cutoff) %>%
    mutate(y.position = 1.05 * max(mydata2$gene_relative, na.rm = TRUE))
  
  # plot
  species_clean <- species |>
    sub("^s__", "", x = _) |>
    gsub("_", " ", x = _)
  
  # p <- ggboxplot(mydata2, x = "VisitTreat", y = "gene_relative", color = "Treat",
  #                add = "jitter", add.params = list(alpha = 0.6, width = 0.2)) +
  #   stat_pvalue_manual(stat_metric, label = "p", tip.length = 0.01, step.increase = 0.1) +
  #   labs(x = NULL, y=genename, title = genename) +
  #   theme(legend.position = "none",
  #     plot.title = element_text(size = 14, hjust = 0.5),
  #     axis.text.x = element_text(angle = 45, hjust = 1)
  #   )
  # 
  # 1) Summary (mean + sd) for bars/errorbars
  
  sum_df <- mydata2 %>%
    group_by(VisitTreat, Visit) %>%
    summarise(
      n    = sum(!is.na(gene_relative)),
      mean = mean(gene_relative, na.rm = TRUE),
      se   = sd(gene_relative, na.rm = TRUE) / sqrt(n),
      .groups = "drop"
    )
  stat_metric <- stat_metric %>%
    mutate(y.position = 1.05 * (max(sum_df$mean, na.rm = TRUE)+max(sum_df$se, na.rm=TRUE)))
  # 2) Barplot + SD + jittered dots + p-values
  p2 <- ggplot() +
    geom_col(
      data = sum_df,
      aes(x = VisitTreat, y = mean, fill = Visit),
      width = 0.7
    ) +
    geom_errorbar(
      data = sum_df,
      aes(x = VisitTreat, ymin = mean - se, ymax = mean + se),
      width = 0.2
    ) +
    stat_pvalue_manual(
      stat_metric,
      label = "p",
      tip.length = 0.01,
      step.increase = 0.1
    ) +
    labs(x = NULL, y=genename, title = genename) +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  #if(mindP <= p_cutoff){
    df_delta <- tibble(
      Arm = rep(c("B","C","D","E"), times = c(length(gdB), length(gdC), length(gdD), length(gdE))),
      delta_adj = c(gdB - mdA, gdC - mdA, gdD - mdA, gdE - mdA)
    )
    
    sum_df<- df_delta %>%
      group_by(Arm) %>%
      summarise(
        n = sum(!is.na(delta_adj)),
        mean = mean(delta_adj, na.rm = TRUE),
        se = sd(delta_adj, na.rm = TRUE) / sqrt(n),
        .groups = "drop"
      )
    
    stat_delta <- tibble(
      group = c("B","C","D","E"),
      p = c(pdB_dA, pdC_dA, pdD_dA, pdE_dA)
    ) %>%
      mutate(
        p = as.numeric(p),
        p.label = sprintf("%.2f", p)
      )
    
    # choose label height above bars
    ymax <- max(sum_df$mean + sum_df$se, na.rm = TRUE)
    stat_delta <- stat_delta %>%
      mutate(y.position = ymax * 1.10)
    
    p3 <- ggplot() +
      geom_col(data = sum_df,aes(x = Arm, y = mean, fill = Arm),width = 0.7) +
      geom_errorbar(data = sum_df,aes(x = Arm, ymin = mean - se, ymax = mean + se),width = 0.2) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_text(data = stat_delta,aes(x = group, y = y.position, label = paste0("p=", p.label)),vjust = 0) +
      labs(x = NULL,y = "W8-BL:(TreatΔ-VehicleAΔ)",title = genename) +
      theme_bw()+theme(legend.position = "none",axis.text.x = element_text(size = 12))
    #png(filename=paste0(species, ".BL_adjust.png"), width=1600, height=1600, res=300) 
    #print(p3)
    #dev.off()

    combined_plot <- ggarrange(p2, p3,ncol = 2,nrow = 1,labels = c("A", "B"))
    imgfile <- paste0(species, ".png")
    ggsave(imgfile,combined_plot,width = 8,height = 4)
    
    doc <- doc %>%
      add_slide(layout = "Title and Content", master = "Office Theme") %>%
      ph_with(
        fpar(ftext(paste0(species, ": baseline adjusted treatment effect"),
                   prop = fp_text(font.size = 20, bold = TRUE))),
        location = ph_location_type(type = "title")
      ) %>%
      ph_with(
        external_img(src = imgfile),
        location = ph_location(left = 0.5, top = 1.5, width = 5.1, height = 6)  # spans where p3+p2 used to be
      ) %>%
      ph_with(
        fpar(ftext(word_metric, prop = fp_text(font.size = 10))),
        location = ph_location(left = 5.5, top = 0.5, width = 4.6, height = 6)
      )
    
 # }
  
    message(word_metric)
    
  
}

print(doc, target = paste0(filename, ".paired.stat.pptx"))

