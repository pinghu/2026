##########################################
#Ping Hu April 16, 2025
### Would like to add into the PPT summary
#### Need to filter out the summary graph with the %zero to be less than 75%
##########################################
rm(list=ls())
args <- commandArgs(trailingOnly = TRUE)
filename="relative_percent_mat.txt"
filename1="GSS3201_P4.meta-12-5-2025-ArmCheek67.txt"
X<-read.table(filename1, sep="\t", header=TRUE)
d <- dim(X);
library(ggplot2)
library(ggpubr)
library(rstatix)
library(dplyr)
library(tidyverse)
library(ggtext)
library(xfun)

library(openxlsx) 
library(lattice)   
library(reshape2)
library(gridExtra)
library(stringr)
library(rlang)
library(officer)
library(magrittr)
library(rvg)
plot_bar_with_p_ttest <- function(data, var, group, filename, titlename, p_cutoff = 0.10,
                            comparisons = NULL, adjust = "BH") {
  var_sym   <- rlang::ensym(var)
  group_sym <- rlang::ensym(group)
  group_name <- rlang::as_name(group_sym)
  
  # --- summary stats for barplot ---
  sumdf <- data %>%
    group_by(!!group_sym) %>%
    summarise(
      mean  = mean(!!var_sym, na.rm = TRUE),
      se    = sd(!!var_sym,   na.rm = TRUE) / sqrt(sum(!is.na(!!var_sym))),
      n     = sum(!is.na(!!var_sym)),
      .groups = "drop"
    ) %>%
    mutate(ymax = mean + se)
  
  # --- build comparisons (pairs of groups) ---
  present_groups <- unique(as.character(dplyr::pull(sumdf, !!group_sym)))
  if (is.null(comparisons)) {
    # all pairwise among present groups
    cmp <- combn(present_groups, 2, simplify = FALSE)
  } else {
    # use provided comparisons but keep only those present in data
    cmp <- lapply(comparisons, function(x) x[ x %in% present_groups ])
    cmp <- Filter(function(x) length(x) == 2, cmp)
  }
  if (length(cmp) == 0) {
    message("No valid comparisons; will plot bars without p-values.")
  }
  
  # --- safe pairwise t-tests ---
  safe_p <- function(x1, x2) {
    n1 <- sum(!is.na(x1)); n2 <- sum(!is.na(x2))
    sd1 <- stats::sd(x1, na.rm = TRUE); sd2 <- stats::sd(x2, na.rm = TRUE)
    constant_or_small <- (n1 < 2) || (n2 < 2) || (isTRUE(sd1 == 0)) || (isTRUE(sd2 == 0))
    if (constant_or_small) return(1)
    out <- try(stats::t.test(x1, x2), silent = TRUE)
    if (inherits(out, "try-error")) 1 else out$p.value
  }
  
  stat.df <- if (length(cmp) > 0) {
    do.call(rbind, lapply(cmp, function(pr) {
      g1 <- pr[1]; g2 <- pr[2]
      x1 <- data %>% filter(!!group_sym == g1) %>% pull(!!var_sym)
      x2 <- data %>% filter(!!group_sym == g2) %>% pull(!!var_sym)
      data.frame(group1 = g1, group2 = g2, p = safe_p(x1, x2), stringsAsFactors = FALSE)
    })) %>% as_tibble()
  } else {
    tibble(group1 = character(), group2 = character(), p = numeric())
  }
  
  if (nrow(stat.df) > 0) {
    stat.df <- stat.df %>%
      mutate(p.adj = p.adjust(p, method = adjust)) %>%
      mutate(
        label = dplyr::if_else(p.adj < 0.001, "p < 0.001",
                               paste0("p = ", rstatix::p_format(p.adj, digits = 2)))
      ) %>%
      filter(p.adj <= p_cutoff)
  }
  
  # --- y position for p-value labels ---
  y_top <- 1.05 * max(sumdf$ymax, na.rm = TRUE)
  if (nrow(stat.df) > 0) stat.df$y.position <- y_top
  
  # --- plotting ---
  png(filename = paste0(filename, ".ttest.png"), width = 1600, height = 1600, res = 300)
  p <- ggplot(sumdf, aes(x = !!group_sym, y = mean, fill = !!group_sym)) +
    geom_col(width = 0.7, color = "black", alpha = 0.9) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.15) +
    labs(
      x = NULL,
      y = rlang::as_name(var_sym),
      title = titlename
    ) +
    theme_classic(base_size = 12) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  if (nrow(stat.df) > 0) {
    p <- p + stat_pvalue_manual(
      stat.df,
      label = "label",
      tip.length = 0.01,
      step.increase = 0.10
    )
  }
  
  print(p)
  dev.off()
  invisible(p)
}

plot_bar_with_p_wilcox <- function(data, var, group, filename, titlename,
                                   p_cutoff = 0.10, comparisons = NULL, adjust = "BH") {
  var_sym   <- ensym(var)
  group_sym <- ensym(group)
  group_name <- as_name(group_sym)
  
  # Summary stats (mean, SE, n)
  sumdf <- data %>%
    group_by(!!group_sym) %>%
    summarise(
      mean  = mean(!!var_sym, na.rm = TRUE),
      se    = sd(!!var_sym,   na.rm = TRUE) / sqrt(sum(!is.na(!!var_sym))),
      n     = sum(!is.na(!!var_sym)),
      .groups = "drop"
    ) %>%
    mutate(ymax = mean + se)
  
  # Ensure x order is stable
  sumdf[[group_name]] <- factor(sumdf[[group_name]], levels = unique(sumdf[[group_name]]))
  
  # Build comparisons (default: all pairwise among present groups)
  present <- levels(sumdf[[group_name]])
  if (is.null(comparisons)) {
    cmp <- combn(present, 2, simplify = FALSE)
  } else {
    # keep only valid pairs that are present
    cmp <- Filter(function(x) length(x) == 2 && all(x %in% present), comparisons)
  }
  
  # Safe Wilcoxon rank-sum p-value (fallback p=1)
  safe_wilcox_p <- function(x1, x2) {
    n1 <- sum(!is.na(x1)); n2 <- sum(!is.na(x2))
    if (n1 < 1 || n2 < 1) return(1)
    out <- try(wilcox.test(x1, x2, exact = FALSE, correct = TRUE), silent = TRUE)
    if (inherits(out, "try-error") || is.na(out$p.value)) 1 else out$p.value
  }
  
  # Pairwise tests
  stat.df <- if (length(cmp)) {
    do.call(rbind, lapply(cmp, function(pr) {
      g1 <- pr[1]; g2 <- pr[2]
      x1 <- data %>% filter(!!group_sym == g1) %>% pull(!!var_sym)
      x2 <- data %>% filter(!!group_sym == g2) %>% pull(!!var_sym)
      data.frame(group1 = g1, group2 = g2, p = safe_wilcox_p(x1, x2))
    })) %>% as_tibble()
  } else {
    tibble(group1 = character(), group2 = character(), p = numeric())
  }
  
  if (nrow(stat.df)) {
    stat.df <- stat.df %>%
      mutate(p.adj = p.adjust(p, method = adjust),
             label = if_else(p.adj < 0.001, "p < 0.001",
                             paste0("p = ", rstatix::p_format(p.adj, digits = 2)))) %>%
      filter(p.adj <= p_cutoff)
    # y position for brackets
    y_top <- 1.05 * max(sumdf$ymax, na.rm = TRUE)
    stat.df$y.position <- y_top
  }
  
  # Plot
  p <- ggplot(sumdf, aes(x = !!group_sym, y = mean, fill = !!group_sym)) +
    geom_col(width = 0.7, color = "black") +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.15) +
    labs(x = NULL, y = as_name(var_sym), title = titlename) +
    theme_classic(base_size = 12) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  if (nrow(stat.df)) {
    p <- p + stat_pvalue_manual(stat.df, label = "label", tip.length = 0.01, step.increase = 0.10)
  }
  
  # Save safely
  ggsave(filename = paste0(filename, ".wilcox.png"), plot = p, width = 6, height = 6, dpi = 300)
  invisible(p)
}

plot_bar_with_p_wilcox_withDots <- function(data, var, group, filename, titlename,
                                   p_cutoff = 0.10, comparisons = NULL, adjust = "BH") {
  var_sym   <- ensym(var)
  group_sym <- ensym(group)
  group_name <- as_name(group_sym)
  
  # summary stats
  sumdf <- data %>%
    group_by(!!group_sym) %>%
    summarise(
      mean = mean(!!var_sym, na.rm = TRUE),
      se   = sd(!!var_sym,   na.rm = TRUE) / sqrt(sum(!is.na(!!var_sym))),
      n    = sum(!is.na(!!var_sym)),
      .groups = "drop"
    ) %>%
    mutate(ymax = mean + se)
  
  # stable x order
  sumdf[[group_name]] <- factor(sumdf[[group_name]], levels = unique(sumdf[[group_name]]))
  
  # default: all pairwise among present groups
  present <- levels(sumdf[[group_name]])
  if (is.null(comparisons)) {
    cmp <- combn(present, 2, simplify = FALSE)
  } else {
    cmp <- Filter(function(x) length(x) == 2 && all(x %in% present), comparisons)
  }
  
  # safe wilcoxon
  safe_wilcox_p <- function(x1, x2) {
    n1 <- sum(!is.na(x1)); n2 <- sum(!is.na(x2))
    if (n1 < 1 || n2 < 1) return(1)
    out <- try(wilcox.test(x1, x2, exact = FALSE, correct = TRUE), silent = TRUE)
    if (inherits(out, "try-error") || is.na(out$p.value)) 1 else out$p.value
  }
  
  # tests
  stat.df <- if (length(cmp)) {
    do.call(rbind, lapply(cmp, function(pr) {
      g1 <- pr[1]; g2 <- pr[2]
      x1 <- data %>% filter(!!group_sym == g1) %>% pull(!!var_sym)
      x2 <- data %>% filter(!!group_sym == g2) %>% pull(!!var_sym)
      data.frame(group1 = g1, group2 = g2, p = safe_wilcox_p(x1, x2))
    })) %>% tibble::as_tibble()
  } else {
    tibble(group1 = character(), group2 = character(), p = numeric())
  }
  
  if (nrow(stat.df)) {
    stat.df <- stat.df %>%
      mutate(p.adj = p.adjust(p, method = adjust),
             label = if_else(p.adj < 0.001, "p < 0.001",
                             paste0("p = ", rstatix::p_format(p.adj, digits = 2)))) %>%
      filter(p.adj <= p_cutoff)
  }
  
  # raw points for dots
  pts <- data %>%
    select(!!group_sym, value = !!var_sym) %>%
    filter(!is.na(value)) %>%
    mutate(!!group_name := factor(!!group_sym, levels = present))
  
  # bracket height considers bars and dots
  y_top <- 1.05 * max(c(sumdf$ymax, pts$value), na.rm = TRUE)
  if (nrow(stat.df)) stat.df$y.position <- y_top
  
  # plot: bars + dots + errorbars (+ p-values)
  p <- ggplot(sumdf, aes(x = !!group_sym, y = mean, fill = !!group_sym)) +
    geom_col(width = 0.7, color = "black") +
    geom_point(data = pts,
               aes(x = !!group_sym, y = value, group = !!group_sym),
               position = position_jitter(width = 0.08, height = 0),
               size = 1.6, alpha = 0.7, inherit.aes = FALSE) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.15) +
    labs(x = NULL, y = as_name(var_sym), title = titlename) +
    theme_classic(base_size = 12) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  if (nrow(stat.df)) {
    p <- p + stat_pvalue_manual(stat.df, label = "label", tip.length = 0.01, step.increase = 0.10)
  }
  
  ggsave(filename = paste0(filename, ".wilcox.png"), plot = p, width = 6, height = 6, dpi = 300)
  invisible(p)
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

A0<-read.table(filename, sep="\t", header=TRUE)
d <- dim(A0);
B=A0[1:d[1], 2:d[2]]
rownames(B)=A0[,1]
ZZ=as.numeric(min(B[B>0&!is.na(B)]))/100
C=B+ZZ
rownames(C)=A0[,1]

Cname=colnames(A0)[2:d[2]]
colnames(A0)[2:d[2]] <-Cname
colnames(C) <- Cname

Clen=length(Cname) ##there are 3 annotation columns
Disease=rep("NA", Clen)
Lesion=rep("NA", Clen)
Site=rep("NA", Clen)
SID=rep("NA", Clen)
splitname<-strsplit(Cname, "[.]")

for(mm in  1:Clen ){
  Disease[mm]=splitname[[mm]][1]
  Site[mm]=splitname[[mm]][2]
  SID[mm]=splitname[[mm]][4]
  Lesion[mm]=splitname[[mm]][3]
}
Disease[Disease=="No"]="HL"
Disease[Disease=="Yes"]="AD"
SiteLesion=paste0(Site, Lesion)
SiteDisease=paste0(Site, Disease)
SiteDiseaseLesion=paste0(SiteDisease, Lesion)
DiseaseLesion=paste0(Disease, Lesion)
SID[SID=="R28"]=28
ID2=paste0(Disease, ".",Site, ".", Lesion, ".", SID)
#IGA <- as.numeric(B["IGA_score", ])
#SampleDNA<-as.numeric(B["DNAYield_ng", ])

mydata=data.frame(Cname, SID,Disease, Site, Lesion, SiteLesion, SiteDisease, SiteDiseaseLesion,  DiseaseLesion, ID2)

order_levels <- c("ArmADLesion", "ArmADNonlesion", "ArmHLNonlesion",   "CheekADLesion",  "CheekADNonlesion",   "CheekHLNonlesion") 
mydata$SiteDiseaseLesion <- factor(mydata$SiteDiseaseLesion, levels = order_levels, ordered = TRUE)


# Create PowerPoint
doc <- read_pptx()
print(paste("genename", 
      "IGASpearmanCorrelation_rho","IGAPearsonCorrelation_rho", 
      "Arm-IGASpearmanCorrelation_rho",
      "Arm-IGAPearsonCorrelation_rho", 
      "Cheek-IGASpearmanCorrelation_rho", 
      "Cheek-IGAPearsonCorrelation_rho", 
      "IGASpearmanCorrelation_pval","IGAPearsonCorrelation_pval",
      "Arm-IGASpearmanCorrelation_pval",
      "Arm-PearsonCorrelation_pval", 
      "Cheek-IGASpearmanCorrelation_pval", 
      "Cheek-PearsonCorrelation_pval", 
      "ZeroNotDetectable%",
      "KP.SiteDiseaseLesion", 
      "KP.site", 
      "KP.disease", 
      "KP.lesion", 
      "KP.SID", 
      "KP.SiteDisease", 
      "KP.SiteLesion", 
      "KP.DiseaseLesion", 
      "Wilcox.ArmAD.Lesion_vs_Nonlesion", 
      "Wilcox.Arm.ADLesion_vs_HL", 
      "Wilcox.Arm.ADNonlesion_vs_HL", 
      "Wilcox.Arm.AD_vs_HL", 
      "Wilcox.Arm.Lesion_vs_Nonlesion", 
      "Wilcox.CheekAD.Lesion_vs_Nonlesion", 
      "Wilcox.Cheek.ADLesion_vs_HL", 
      "Wilcox.Cheek.ADNonlesion_vs_HL", 
      "Wilcox.Cheek.AD_vs_HL", 
      "Wilcox.Cheek.Lesion_vs_Nonlesion", 
      "Wilcox.AllAD.Lesion_vs_Nonlesion",
      "Wilcox.All.ADLesion_vs_HL", 
      "Wilcox.All.ADNonlesion_vs_HL", 
      "Wilcox.Cheek_vs_Arm", 
      "Wilcox.Lesion.Cheek_vs_Arm", 
      "Wilcox.ADNonlesion.Cheek_vs_Arm", 
      "Wilcox.Healthy.Cheek_vs_Arm",
      "pFisher.AllDisease", 
      "pFisher.ArmDisease", 
      "pFisher.CheekDisease",
      "pFisher.AllLesion",
      "pFisher.ArmLesion", 
      "pFisher.CheekLesion", 
      "pFisher.AllSite", 
      "meanAlloriginal",
      "meanAll",
      "meanArm",  "meanArmAD", "meanArmHealthy",
      "meanArmADNonlesion",
      "meanArmADLesion",
      "meanArmNonlesion",
      "meanCheek",  "meanCheekAD", "meanCheekHealthy",
      "meanCheekADNonlesion",
      "meanCheekADLesion",
      "meanCheekNonlesion",
      "meanADLesion", 
      "meanADNonlesion", 
      "meanHealthy",
      "Present_ArmADLesion", 
      "Present_ArmADNonlesion",
      "Present_CheekADLesion",
      "Present_CheekADNonlesion",
      "Present_ArmHLNonlesion",
      "Present_CheekHLNonlesion",
      collapse = "\t"
))


for (i in 1:d[1]){
  genename=A0[i,1]
  
  sanitize_filename <- function(x) {
    x <- gsub("%", "pct", x, fixed = TRUE)        # replace % with 'pct'
    x <- gsub("[\\/:*?\"<>|]", "_", x)            # remove illegal path chars
    x <- gsub("\\s+", "_", x)                     # optional: replace spaces
    x
  }
  
  safe_genename <- sanitize_filename(genename)
  gene_original <-as.numeric(B[i,])
  gene_relative <-as.numeric(C[i,])
  if (length(unique(gene_original)) <= 1 || all(is.na(gene_relative))) {
    message(paste("Skipping", genename, "- no variation in relative_gene"))
    next  # skip to next iteration
  }
  
  meanAllRelative <-mean(gene_relative)
  meanAllOriginal <-mean(gene_original)
 
  percentZeroNA=(sum(B[i,]==0)+sum(is.na(B[i,])))/(d[2]-1)
  mydata1=data.frame(gene_relative,gene_original,  Cname)
  
  mydata11 <- inner_join(mydata1, X, by = c("Cname" = "NewID"))
  mydata2 <- inner_join(mydata, mydata11, by = c("Cname" = "Cname"))
  mydata2 <- mydata2[mydata2$gene_original != "NA", ]
  result <- cor.test(mydata2$IGA, mydata2$gene_relative, method = "spearman", exact = FALSE)
  spearman_rho <- result$estimate
  spearman_pval <- result$p.value
  result2 <- cor.test(mydata2$IGA, mydata2$gene_relative, method = "pearson", exact = FALSE)
  pearson_rho <- result2$estimate
  pearson_pval <- result2$p.value
  
  KP_relative_SiteDiseaseLesion=kruskal.test(gene_relative ~ SiteDiseaseLesion, data = mydata2)$p.value
  KP_relative_Site=kruskal.test(gene_relative ~ Site, data = mydata2)$p.value
  KP_relative_SID=kruskal.test(gene_relative ~ SID, data = mydata2)$p.value
  KP_relative_Lesion=kruskal.test(gene_relative ~ Lesion, data = mydata2)$p.value
  KP_relative_Disease=kruskal.test(gene_relative ~ Disease, data = mydata2)$p.value
  KP_relative_SiteLesion=kruskal.test(gene_relative ~ SiteLesion, data = mydata2)$p.value
  KP_relative_SiteDisease=kruskal.test(gene_relative ~ SiteDisease, data = mydata2)$p.value
  KP_relative_DiseaseLesion=kruskal.test(gene_relative ~ DiseaseLesion, data = mydata2)$p.value
  
  gArmADLesion=mydata2$gene_relative[mydata2$SiteDiseaseLesion =="ArmADLesion"]
  
  gArmADNonlesion=mydata2$gene_relative[mydata2$SiteDiseaseLesion =="ArmADNonlesion"]
  pArm=my.wilcox.p.value(gArmADLesion, gArmADNonlesion, paired=FALSE)
  
  gCheekADLesion=mydata2$gene_relative[mydata2$SiteDiseaseLesion =="CheekADLesion"]
  gCheekADNonlesion=mydata2$gene_relative[mydata2$SiteDiseaseLesion =="CheekADNonlesion"]
  pCheek=my.wilcox.p.value(gCheekADLesion, gCheekADNonlesion, paired=FALSE)
  
  gCheekHL =mydata2$gene_relative[mydata2$SiteDiseaseLesion =="CheekHLNonlesion"]
  gArmHL =mydata2$gene_relative[mydata2$SiteDiseaseLesion =="ArmHLNonlesion"]
  
  pArmDLes_HL =my.wilcox.p.value(gArmADLesion, gArmHL, paired=FALSE)
  pArmDNL_HL =my.wilcox.p.value(gArmADNonlesion, gArmHL, paired=FALSE)

  pCheekDLes_HL =my.wilcox.p.value(gCheekADLesion, gCheekHL, paired=FALSE)
  pCheekDNL_HL =my.wilcox.p.value(gCheekADNonlesion, gCheekHL, paired=FALSE) 
  
  gArmAD= mydata2$gene_relative[mydata2$SiteDisease =="ArmAD"]
  pArmAD_HL =my.wilcox.p.value(gArmAD, gArmHL, paired=FALSE)

  gCheekAD= mydata2$gene_relative[mydata2$SiteDisease =="CheekAD"]
  pCheekAD_HL =my.wilcox.p.value(gCheekAD, gCheekHL, paired=FALSE)
  
  gArmNonlesion=mydata2$gene_relative[mydata2$SiteLesion =="ArmNonlesion"]
  pArmLes_NL=my.wilcox.p.value(gArmADLesion, gArmNonlesion, paired=FALSE)
  gCheekNonlesion=mydata2$gene_relative[mydata2$SiteLesion =="CheekNonlesion"]
  pCheekLes_NL=my.wilcox.p.value(gCheekADLesion, gCheekNonlesion, paired=FALSE)
  
  gCheek=mydata2$gene_relative[mydata2$Site =="Cheek"]
  gArm=mydata2$gene_relative[mydata2$Site =="Arm"]
  pCheek_Arm=my.wilcox.p.value(gCheek, gArm, paired=FALSE)
  pLesCheek_Arm=my.wilcox.p.value(gCheekADLesion, gArmADLesion, paired=FALSE)
  pADNLCheek_Arm=my.wilcox.p.value(gCheekADNonlesion, gArmADNonlesion, paired=FALSE)
  pHLCheek_Arm=my.wilcox.p.value(gCheekHL, gArmHL, paired=FALSE)
  
  mArm=mean(gArm)
  mArmAD=mean(gArmAD)
  mArmADLesion =mean(gArmADLesion)
  mArmADNonlesion = mean(gArmADNonlesion)
  mArmHL=mean(gArmHL)
  mArmNonlesion=mean(gArmNonlesion)
  
  mCheek=mean(gCheek)
  mCheekAD=mean(gCheekAD)
  mCheekADLesion =mean(gCheekADLesion)
  mCheekADNonlesion = mean(gCheekADNonlesion)
  mCheekHL=mean(gCheekHL)
  mCheekNonlesion=mean(gCheekNonlesion)
  
  gADLes=mydata2$gene_relative[mydata2$DiseaseLesion =="ADLesion"]
  gADNL=mydata2$gene_relative[mydata2$DiseaseLesion =="ADNonlesion"]
  gHL=mydata2$gene_relative[mydata2$DiseaseLesion =="HLNonlesion"]
  mADLes=mean(gADLes)
  mADNL=mean(gADNL)
  mHL=mean(gHL)
  pADLes_HL =my.wilcox.p.value(gADLes, gHL, paired=FALSE)
  pADLes_ADNL =my.wilcox.p.value(gADLes, gADNL, paired=FALSE)
  pADNL_HL =my.wilcox.p.value(gADNL, gHL, paired=FALSE)
  
  
  
  
  Aresult <- cor.test(mydata2$IGA[mydata2$Site=="Arm"], mydata2$gene_relative[mydata2$Site=="Arm"], method = "spearman", exact = FALSE)
  Aspearman_rho <- Aresult$estimate
  Aspearman_pval <- Aresult$p.value
  Aresult2 <- cor.test(mydata2$IGA[mydata2$Site=="Arm"], mydata2$gene_relative[mydata2$Site=="Arm"], method = "pearson", exact = FALSE)
  Apearson_rho <- Aresult2$estimate
  Apearson_pval <- Aresult2$p.value
  
  Cresult <- cor.test(mydata2$IGA[mydata2$Site=="Cheek"], mydata2$gene_relative[mydata2$Site=="Cheek"], method = "spearman", exact = FALSE)
  Cspearman_rho <- Cresult$estimate
  Cspearman_pval <- Cresult$p.value
  Cresult2 <- cor.test(mydata2$IGA[mydata2$Site=="Cheek"], mydata2$gene_relative[mydata2$Site=="Cheek"], method = "pearson", exact = FALSE)
  Cpearson_rho <- Cresult2$estimate
  Cpearson_pval <- Cresult2$p.value
  mydata2 <- mydata2 %>%
    mutate(
      present = if_else(gene_original > 0, 1, 0)  # 1 = present, 0 = absent
    )
  tab <- table(mydata2$Disease, mydata2$present)
  #tab
  #prop.table(tab, 1)
  if (ncol(tab) < 2 || nrow(tab) < 2) {pf_Disease=1}else{ pf_Disease=fisher.test(tab)$p}
  #for larger data use this : chisq.test(tab, correct = FALSE)
  tab1 <- table(mydata2$Lesion, mydata2$present)
  if (ncol(tab1) < 2 || nrow(tab1) < 2) {pf_Lesion=1}else{pf_Lesion=fisher.test(tab1)$p}
  tab2 <- table(mydata2$Site, mydata2$present)
  if (ncol(tab2) < 2 || nrow(tab2) < 2) {pf_Site=1}else{pf_Site=fisher.test(tab2)$p}
  myArm=mydata2[Site=="Arm",]
  rm(tab,tab1,tab2)
  tab <- table(myArm$Disease, myArm$present)
  if (ncol(tab) < 2 || nrow(tab) < 2) {pf_Arm_Disease=1}else{pf_Arm_Disease=fisher.test(tab)$p}
  tab1 <- table(myArm$Lesion, myArm$present)
  if (ncol(tab1) < 2 || nrow(tab1) < 2) {pf_Arm_Lesion=1}else{pf_Arm_Lesion=fisher.test(tab1)$p}
  myCheek=mydata2[Site=="Cheek",]
  rm(tab,tab1)
  tab <- table(myCheek$Disease, myCheek$present)
  if (ncol(tab) < 2 || nrow(tab) < 2) {pf_Cheek_Disease=1}else{pf_Cheek_Disease=fisher.test(tab)$p}
  tab1 <- table(myCheek$Lesion, myCheek$present)
  if (ncol(tab1) < 2 || nrow(tab1) < 2) {pf_Cheek_Lesion=1}else{pf_Cheek_Lesion=fisher.test(tab1)$p}
 
  present_ArmADLesion=sum(mydata2[SiteDiseaseLesion=="ArmADLesion",]$present)
  present_ArmADNonlesion=sum(mydata2[SiteDiseaseLesion=="ArmADNonlesion",]$present)
  present_CheekADLesion=sum(mydata2[SiteDiseaseLesion=="CheekADLesion",]$present)
  present_CheekADNonlesion=sum(mydata2[SiteDiseaseLesion=="CheekADNonlesion",]$present)
  present_ArmHLNonlesion=sum(mydata2[SiteDiseaseLesion=="ArmHLNonlesion",]$present)
  present_CheekHLNonlesion=sum(mydata2[SiteDiseaseLesion=="CheekHLNonlesion",]$present)
  word_KP=paste(genename, "|IGASpearmanCorrelation_rho =", round(spearman_rho, 2), 
                "|IGAPearsonCorrelation_rho =", round(pearson_rho, 2), 
                "|Arm-IGASpearmanCorrelation_rho =", round(Aspearman_rho, 2), 
                "|Arm-IGAPearsonCorrelation_rho =", round(Apearson_rho, 2),
                "|Cheek-IGASpearmanCorrelation_rho =", round(Cspearman_rho, 2), 
                "|Cheek-IGAPearsonCorrelation_rho =", round(Cpearson_rho, 2), 
                "|IGASpearmanCorrelation_pval =", sprintf("%.02f",spearman_pval),
                "|IGAPearsonCorrelation_pval =", sprintf("%.02f",pearson_pval),
                "|Arm-IGASpearmanCorrelation_pval =", sprintf("%.02f",Aspearman_pval),
                "|Arm-PearsonCorrelation_pval =", sprintf("%.02f",Apearson_pval),
                "|Cheek-IGASpearmanCorrelation_pval =", sprintf("%.02f",Cspearman_pval),
                "|Cheek-PearsonCorrelation_pval =", sprintf("%.02f",Cpearson_pval),
                "|ZeroNotDetectable%=",sprintf("%.02f", 100*percentZeroNA),
                "|KP.SiteDiseaseLesion=", sprintf("%.02f",KP_relative_SiteDiseaseLesion),
                "|KP.site=", sprintf("%.02f",KP_relative_Site),
                "|KP.disease=", sprintf("%.02f",KP_relative_Disease), 
                "|KP.lesion=", sprintf("%.02f",KP_relative_Lesion),
                "|KP.SID=", sprintf("%.02f",KP_relative_SID),
                "|KP.SiteDisease=", sprintf("%.02f",KP_relative_SiteDisease),
                "|KP.SiteLesion=", sprintf("%.02f",KP_relative_SiteLesion),
                "|KP.DiseaseLesion=", sprintf("%.02f",KP_relative_DiseaseLesion),
                "|Wilcox.ArmAD.Lesion_vs_Nonlesion=", sprintf("%.02f",pArm),
                "|Wilcox.Arm.ADLesion_vs_HL=", sprintf("%.02f",pArmDLes_HL),
                "|Wilcox.Arm.ADNonlesion_vs_HL=", sprintf("%.02f",pArmDNL_HL),
                "|Wilcox.Arm.AD_vs_HL=", sprintf("%.02f",pArmAD_HL),
                "|Wilcox.Arm.Lesion_vs_Nonlesion=", sprintf("%.02f",pArmLes_NL),
                "|Wilcox.CheekAD.Lesion_vs_Nonlesion=", sprintf("%.02f",pCheek),
                "|Wilcox.Cheek.ADLesion_vs_HL=", sprintf("%.02f",pCheekDLes_HL),
                "|Wilcox.Cheek.ADNonlesion_vs_HL=", sprintf("%.02f",pCheekDNL_HL),
                "|Wilcox.Cheek.AD_vs_HL=", sprintf("%.02f",pCheekAD_HL),
                "|Wilcox.Cheek.Lesion_vs_Nonlesion=", sprintf("%.02f",pCheekLes_NL),
                "|Wilcox.AllAD.Lesion_vs_Nonlesion=", sprintf("%.02f",pADLes_ADNL),
                "|Wilcox.All.ADLesion_vs_HL=", sprintf("%.02f",pADLes_HL),
                "|Wilcox.All.ADNonlesion_vs_HL=", sprintf("%.02f",pADNL_HL),
                
                "|Wilcox.Cheek_vs_Arm=", sprintf("%.02f",pCheek_Arm),
                "|Wilcox.Lesion.Cheek_vs_Arm=", sprintf("%.02f",pLesCheek_Arm),
                "|Wilcox.ADNonlesion.Cheek_vs_Arm=", sprintf("%.02f",pADNLCheek_Arm),
                "|Wilcox.Healthy.Cheek_vs_Arm=", sprintf("%.02f",pHLCheek_Arm),
                "|pFisher.AllDisease=", pf_Disease,
                "|pFisher.ArmDisease=", pf_Arm_Disease,
                "|pFisher.CheekDisease=", pf_Cheek_Disease,
                "|pFisher.AllLesion=", pf_Lesion,
                "|pFisher.ArmLesion=", pf_Arm_Lesion,
                "|pFisher.CheekLesion=", pf_Cheek_Lesion,
                "|pFisher.AllSite=", pf_Site,
                "|meanAlloriginal=",meanAllOriginal,
                "|meanAll=",meanAllRelative,
                "|meanArm=",mArm,  "|meanArmAD=",mArmAD, "|meanArmHealthy=",mArmHL,
                "|meanArmADNonlesion=",mArmADNonlesion,  
                "|meanArmADLesion=",mArmADLesion,
                "|meanArmNonlesion=",mArmNonlesion,
                "|meanCheek=",mCheek,  "|meanCheekAD=",mCheekAD, "|meanCheekHealthy=",mCheekHL,
                "|meanCheekADNonlesion=",mCheekADNonlesion,  
                "|meanCheekADLesion=",mCheekADLesion,
                "|meanCheekNonlesion=",mCheekNonlesion,
                "|meanADLesion=", mADLes,
                "|meanADNonlesion=", mADNL,
                "|meanHealthy=", mHL,
                "|Present_ArmADLesion=", present_ArmADLesion,
                "|Present_ArmADNonlesion=",present_ArmADNonlesion,
                "|Present_ArmHLNonlesion=",present_ArmHLNonlesion,
                "|Present_CheekADLesion=",present_CheekADLesion,
                "|Present_CheekADNonlesion=",present_CheekADNonlesion,
                "|Present_CheekHLNonlesion=",present_CheekHLNonlesion
  )
  print(word_KP)
  print(paste(genename, 
              spearman_rho, pearson_rho,
              Aspearman_rho,  Apearson_rho, Cspearman_rho, Cpearson_rho, 
              spearman_pval, pearson_pval, 
             Aspearman_pval,Apearson_pval,Cspearman_pval,Cpearson_pval,
              percentZeroNA,
              KP_relative_SiteDiseaseLesion,
              KP_relative_Site,
              KP_relative_Disease, 
              KP_relative_Lesion,
              KP_relative_SID,
              KP_relative_SiteDisease,
              KP_relative_SiteLesion,
              KP_relative_DiseaseLesion,
              pArm,
              pArmDLes_HL,
              pArmDNL_HL,
              pArmAD_HL,
              pArmLes_NL,
              pCheek,
              pCheekDLes_HL,
              pCheekDNL_HL,
              pCheekAD_HL,
              pCheekLes_NL,
             pADLes_ADNL,
            pADLes_HL,
             pADNL_HL,
              pCheek_Arm,
              pLesCheek_Arm,
              pADNLCheek_Arm,
              pHLCheek_Arm,
          pf_Disease,
         pf_Arm_Disease,
          pf_Cheek_Disease,
          pf_Lesion,
           pf_Arm_Lesion,
           pf_Cheek_Lesion,
           pf_Site,
              meanAllOriginal,
              meanAllRelative,
              mArm,  mArmAD, mArmHL,
              mArmADNonlesion,  
              mArmADLesion,
              mArmNonlesion,
              mCheek,  mCheekAD, mCheekHL,
              mCheekADNonlesion,  
              mCheekADLesion,
              mCheekNonlesion, 
              mADLes,
              mADNL,
              mHL,
         present_ArmADLesion,
         present_ArmADNonlesion,
         present_ArmHLNonlesion,
        present_CheekADLesion,
        present_CheekADNonlesion,
        present_CheekHLNonlesion,
              collapse = "||"
  )  )
  #print(word_KP)
  
  if((KP_relative_SiteDiseaseLesion <=0.1)||(spearman_pval <=0.1)||(pearson_pval<=0.1)) {
    plot_bar_with_p_wilcox_withDots(data = mydata2, var = gene_relative, group = SiteDiseaseLesion,
      titlename= genename, filename = paste0("SiteDiseaseLesion.", safe_genename),p_cutoff = 0.1)
  }
  # if(KP_relative_SiteLesion <=0.1){
  #   plot_bar_with_p_wilcox_withDots(data = mydata2, var = gene_relative, group = SiteLesion,
  #                                   titlename= genename, filename = paste0("SiteLesion.", safe_genename),p_cutoff = 0.1)
  # }
  # if(KP_relative_SiteDisease <=0.1){
  #   plot_bar_with_p_wilcox_withDots(data = mydata2, var = gene_relative, group = SiteDisease,
  #                                   titlename= genename, filename = paste0("SiteDisease.", safe_genename),p_cutoff = 0.1)
  # }
  # if(KP_relative_DiseaseLesion <=0.1){
  #   plot_bar_with_p_wilcox_withDots(data = mydata2, var = gene_relative, group = DiseaseLesion,
  #                                   titlename= genename, filename = paste0("DiseaseLesion.", safe_genename),p_cutoff = 0.1)
  # }
  # if(KP_relative_Site <=0.1){
  #   plot_bar_with_p_wilcox_withDots(data = mydata2, var = gene_relative, group = Site,
  #                                   titlename= genename, filename = paste0("Site.", safe_genename),p_cutoff = 0.1)
  # } 
  # if(KP_relative_Disease <=0.1){
  #   plot_bar_with_p_wilcox_withDots(data = mydata2, var = gene_relative, group = Disease,
  #                                   titlename= genename, filename = paste0("Disease.", safe_genename),p_cutoff = 0.1)
  # }
  # if(KP_relative_Lesion <=0.1){
  #   plot_bar_with_p_wilcox_withDots(data = mydata2, var = gene_relative, group = Lesion,
  #                                   titlename= genename, filename = paste0("Lesion.", safe_genename),p_cutoff = 0.1)
  # }
  # 

  
  img_dir <- "."
  img_rel <- file.path(img_dir, paste0("SiteDiseaseLesion.", safe_genename, ".wilcox.png"))
  if(file.exists(img_rel)){
  doc <- doc %>%
    add_slide(layout = "Title and Content", master = "Office Theme") %>%
    ph_with(
      fpar(ftext(genename, prop = fp_text(font.size = 20, bold = TRUE))),
      location = ph_location_type(type = "title")) %>%
    ph_with(external_img(img_rel, width =4.8 , height = 6),
            location = ph_location(left = 0.5, top = 1.5)) %>%
    ph_with(
      fpar(ftext(word_KP, prop = fp_text(font.size = 10))),
      location = ph_location(left = 5.5, top = 0.5, width = 4.6, height = 6)
    )
  }else{
    doc <- doc %>%
      add_slide(layout = "Title and Content", master = "Office Theme") %>%
      ph_with(value = genename, location = ph_location_type(type = "title")) %>%
      ph_with(
        fpar(ftext(word_KP, prop = fp_text(font.size = 12))),
        location = ph_location(left = 0.5, top = 2, width = 9, height = 1)
      )
  }

}
   
print(doc, target = paste0(filename, ".pptx"))
