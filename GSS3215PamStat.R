##########################################
#Ping Hu April 16, 2025
### Would like to add into the PPT summary
#### Need to filter out the summary graph with the %zero to be less than 75%
##########################################
rm(list=ls())
args <- commandArgs(trailingOnly = TRUE)
filename1="GSS3215Pam_meta.txt"
X<-read.table(filename1, sep="\t", header=TRUE)
#d <- dim(A);
#mydata <- inner_join(mydata0, X, by = c("Cname_clean" = "ID"))
#mydata$VisitTreat =paste0(mydata$Visit, ".",mydata$Trt)
#mydata$VisitTreat2 =mydata$VisitTreat
#mydata$VisitTreat2[mydata$Visit=="BL"]="BL"
filename=args[1]
#filename="count7.txt.no_outlier.relab10"
#filename="genefamily.community.uniref90_go.relab10.name.xls.no_outlier"
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
library(rstatix)
library(ggpubr)
library(dplyr)

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

Clen=length(Cname)

Treat=rep("NA", Clen)
Visit=rep("NA", Clen)
Side=rep("NA", Clen)
SID=rep("NA", Clen)
splitname<-strsplit(Cname, "[_]")

for(mm in  1:Clen ){
  Visit[mm]=splitname[[mm]][4]
  Side[mm]=splitname[[mm]][5]
  SID[mm]=splitname[[mm]][3]
  Treat[mm]=splitname[[mm]][7]
}

mydata=data.frame(Cname, Visit,Side, SID, Treat )
mydata$VisitTreat =paste0(mydata$Visit, ".",mydata$Treat)
mydata$VisitTreat2 =mydata$VisitTreat
mydata$VisitTreat2[mydata$Visit=="BL"]="BL"
order_levels <- c("BL", "W8.A", "W8.B",   "W8.C",  "W8.D",   "W8.E") 
mydata$VisitTreat2 <- factor(mydata$VisitTreat2, levels = order_levels, ordered = TRUE)
# Create PowerPoint
doc <- read_pptx()

for (i in 1:d[1]){
  genename=A0[i,1]
  species <- sub(".*\\.(s__[^.]+)$", "\\1", genename)
  #species
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
  mydata0=data.frame(gene_relative,gene_original,  Cname)
  mydata1 <- inner_join(mydata0, mydata, by = c("Cname" = "Cname"))
  mydata2 <-inner_join(mydata1, X, by = c("Cname" = "ID"))
  mydata2 <- mydata2[mydata2$gene_original != "NA", ]

  KP_VisitTreat2=kruskal.test(gene_relative ~ VisitTreat2, data = mydata2)$p.value
  KP_VisitTreat=kruskal.test(gene_relative ~ VisitTreat, data = mydata2)$p.value
  KP_SID=kruskal.test(gene_relative ~ SID, data = mydata2)$p.value
  KP_Visit=kruskal.test(gene_relative ~ Visit, data = mydata2)$p.value
  KP_Treat=kruskal.test(gene_relative ~ Treat, data = mydata2)$p.value
  KP_Side=kruskal.test(gene_relative ~ Side, data = mydata2)$p.value
  
  gBL=mydata2$gene_relative[mydata2$VisitTreat2 =="BL"]
  gW8A=mydata2$gene_relative[mydata2$VisitTreat2 =="W8.A"]
  gW8B=mydata2$gene_relative[mydata2$VisitTreat2 =="W8.B"]
  gW8C=mydata2$gene_relative[mydata2$VisitTreat2 =="W8.C"]
  gW8D=mydata2$gene_relative[mydata2$VisitTreat2 =="W8.D"]
  gW8E=mydata2$gene_relative[mydata2$VisitTreat2 =="W8.E"]
  
  pW8B_W8A=my.wilcox.p.value(gW8B, gW8A, paired=FALSE)
  pW8C_W8A=my.wilcox.p.value(gW8C, gW8A, paired=FALSE)
  pW8D_W8A=my.wilcox.p.value(gW8D, gW8A, paired=FALSE)
  pW8E_W8A=my.wilcox.p.value(gW8E, gW8A, paired=FALSE)
  pW8A_BL=my.wilcox.p.value(gW8A, gBL, paired=FALSE)
  pW8B_BL=my.wilcox.p.value(gW8B, gBL, paired=FALSE)
  pW8C_BL=my.wilcox.p.value(gW8C, gBL, paired=FALSE)
  pW8D_BL=my.wilcox.p.value(gW8D, gBL, paired=FALSE)
  pW8E_BL=my.wilcox.p.value(gW8E, gBL, paired=FALSE)
  
 
  mBL=mean(gBL)
  mW8A=mean(gW8A)
  mW8B=mean(gW8B)
  mW8C=mean(gW8C)
  mW8D=mean(gW8D)
  mW8E=mean(gW8E)
  
  mydataBL <-mydata2[mydata2$Visit.x == "BL",]
  mydataW8 <-mydata2[mydata2$Visit.x == "W8",]
  ############################################
  #https://www.datanovia.com/en/blog/ggpubr-how-to-add-p-values-generated-elsewhere-to-a-ggplot/
  ###Baseline Correlation
  result <- cor.test(mydataBL$BaselineValue.PostAcneSpotJSV, mydataBL$gene_relative, method = "spearman", exact = FALSE)
  spearman_rho <- result$estimate
  spearman_pval <- result$p.value
  result2 <- cor.test(mydataBL$BaselineValue.PostAcneSpotJSV, mydataBL$gene_relative, method = "pearson", exact = FALSE)
  pearson_rho <- result2$estimate
  pearson_pval <- result2$p.value
  minP=min(pW8A_BL,pW8B_BL,pW8C_BL,pW8D_BL,pW8E_BL,pW8B_W8A, pW8C_W8A, pW8D_W8A,pW8E_W8A,KP_VisitTreat2,pearson_pval,spearman_pval,na.rm=TRUE )
  word=paste(species," |", genename,
                     " |PostAcneSpotJSVSpearmanCorrelation_rho =", round(spearman_rho, 2), 
                     " |PostAcneSpotJSVPearsonCorrelation_rho =", round(pearson_rho, 2), 
                     " |PostAcneSpotJSVSpearmanCorrelation_pval =", sprintf("%.02f",spearman_pval),
                     " |PostAcneSpotJSVPearsonCorrelation_pval =", sprintf("%.02f",pearson_pval),
                     " |KP.Treat=", sprintf("%.02f",KP_Treat),
                     " |KP.Side=", sprintf("%.02f",KP_Side),
                     " |KP.SID=", sprintf("%.02f",KP_SID), 
                     " |KP.Visit=", sprintf("%.02f",KP_Visit),
                     " |KP.VisitTreat=", sprintf("%.02f",KP_VisitTreat),
                     " |KP.VisitTreat2=", sprintf("%.02f",KP_VisitTreat2),
                     " |Wilcox.BL_vs_W8A=", sprintf("%.02f",pW8A_BL),
                     " |Wilcox.BL_vs_W8B=", sprintf("%.02f",pW8B_BL),
                     " |Wilcox.BL_vs_W8C=", sprintf("%.02f",pW8C_BL),
                     " |Wilcox.BL_vs_W8D=", sprintf("%.02f",pW8D_BL),
                     " |Wilcox.BL_vs_W8E=", sprintf("%.02f",pW8E_BL),
                     " |Wilcox.W8A_vs_W8B=", sprintf("%.02f",pW8B_W8A),
                     " |Wilcox.W8A_vs_W8C=", sprintf("%.02f",pW8C_W8A),
                     " |Wilcox.W8A_vs_W8D=", sprintf("%.02f",pW8D_W8A),
                     " |Wilcox.W8A_vs_W8E=", sprintf("%.02f",pW8E_W8A),
                     " |meanBL=",mBL,  " |meanW8A=",mW8A, " |meanW8B=",mW8B,
                     " |meanW8C=",mW8C,  " |meanW8D=",mW8D, " |meanW8E=",mW8E
  )
  print(word)
  if(minP<=0.1){
  # 1) Clean data (this prevents most wilcox_test crashes)
  df <- mydata2 %>%
    mutate(
      gene_relative = as.numeric(gene_relative),
      VisitTreat2   = as.factor(VisitTreat2)
    ) %>%
    filter(is.finite(gene_relative)) %>%     # removes NA, NaN, Inf, -Inf
    drop_na(VisitTreat2)
  
 
  p <- ggplot(df, aes(x = VisitTreat2, y = gene_relative)) +
    stat_summary(fun = mean, geom = "col", width = 0.7) +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
    #stat_pvalue_manual(stat.gene, label = "p.signif", tip.length = 0.01) +
    #facet_wrap(~Gene, scales = "free_y") +   # remove if single gene
    labs(title = species) +
    theme_classic()
  png(filename=paste0(species, ".bar.png"), width=1000, height=1000, res=300)
  print(p)
  dev.off()
  # run t-test and keep only p <= 0.1
  # stat.gene <- mydata2 %>%
  #   wilcox_test(gene_relative ~ VisitTreat2) %>%
  #   filter(p <= 0.1) %>%   # <— filter here
  #   mutate(y.position = 1.05 * max(mydata2$gene_relative))
  
  # plot
  png(filename=paste0(species, ".png"), width=1600, height=1600, res=300)
  pobserved <- ggboxplot(mydata, x = "VisitTreat2", y = "gene_relative", color = "VisitTreat2", add = "jitter", add.params = list(alpha = 0.6, width = 0.2))+labs(title = species)+
  
  # pobserved <- pobserved +
  #   #stat_pvalue_manual(stat.observed, label = "p", tip.length = 0.01, step.increase = 0.1) +
  #   labs(x = NULL) +
     theme(
       legend.position = "none",
  #     axis.text.x = element_text(angle = 45, hjust = 1)
     ) 
  print(pobserved)
  dev.off()
  #addWorksheet(wb, "observed");             writeData(wb, "observed", stat.observed)
  doc <- read_pptx()
  
  # Add slide with both plots and caption
  doc <- doc %>%
    add_slide(layout = "Title and Content", master = "Office Theme") %>%
    ph_with(
      fpar(ftext(species, prop = fp_text(font.size = 20, bold = TRUE))),
      location = ph_location_type(type = "title")
    ) %>%
    ph_with(dml(ggobj = pobserved), location = ph_location(left = 0.5, top = 1.5, width = 2.4, height = 6)) %>%
    ph_with(dml(ggobj = p), location = ph_location(left = 3, top = 1.5, width = 2.4, height = 6)) %>%
    ph_with(
      fpar(ftext(word, prop = fp_text(font.size = 10))),
      location = ph_location(left = 5.5, top = 0.5, width = 4.6, height = 6)
    )
  }
}
   
print(doc, target = paste0(filename, ".pptx"))
