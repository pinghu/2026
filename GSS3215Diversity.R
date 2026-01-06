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
filename="count7.txt.no_outlier.10filter"
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
d=dim(test)
d
X=apply(test, 1, sum)
test_filter2=test[(X/d[2] >=10),]

test_d=data.frame(t(test_filter2))

min(apply(test_d, 1, sum)) ###0

test_a <- decostand(test_d, method = "total")

Cname=rownames(test_d)
Cname_clean <-sub(".metaphlan4", "", Cname)
Clen=length(Cname) ##there are 3 annotation columns

Clen=length(Cname_clean) ##there are 3 annotation columns

shannon<-diversity(test_d)
simp<-diversity(test_d, "simpson")
invsimp<-diversity(test_d, "inv")

observed<-apply(test_d>0,1,sum)
N <- apply(test_d,1,sum)
###Richness Index#####
Menhinick_index<-observed/sqrt(N)
Margalef_index <-(observed-1)/log(N)

mydata0=data.frame(shannon,simp,invsimp, observed, Menhinick_index, Margalef_index , Cname_clean)

mydata <- inner_join(mydata0, A, by = c("Cname_clean" = "ID"))

mydata$VisitTreat =paste0(mydata$Visit, ".",mydata$Trt)
mydata$VisitTreat2 =mydata$VisitTreat
mydata$VisitTreat2[mydata$Visit=="BL"]="BL"
# order_levels <- c("ArmADLesion", "ArmADNonlesion", "ArmHLNonlesion",   "CheekADLesion",  "CheekADNonlesion",   "CheekHLNonlesion") 
# mydata$SiteDiseaseLesion <- factor(mydata$SiteDiseaseLesion, levels = order_levels, ordered = TRUE)

write.table(t(mydata), file =paste0(filename, ".alphadiversity.RShort"), col.names=TRUE, row.names=TRUE, sep = "\t")
addWorksheet(wb, "aDiversity");             writeData(wb, "aDiversity", mydata)
mydataBL <-mydata[mydata$Visit == "BL",]
mydataW8 <-mydata[mydata$Visit == "W8",]
############################################
#https://www.datanovia.com/en/blog/ggpubr-how-to-add-p-values-generated-elsewhere-to-a-ggplot/
###Baseline Correlation
result <- cor.test(mydataBL$BaselineValue.PostAcneSpotJSV, mydataBL$shannon, method = "spearman", exact = FALSE)
spearman_rho <- result$estimate
spearman_pval <- result$p.value
result2 <- cor.test(mydataBL$BaselineValue.PostAcneSpotJSV, mydataBL$shannon, method = "pearson", exact = FALSE)
pearson_rho <- result2$estimate
pearson_pval <- result2$p.value

KP_Visit=kruskal.test(shannon ~ Visit, data = mydata)$p.value
KP_Treat=kruskal.test(shannon ~ Trt, data = mydata)$p.value
KP_Side=kruskal.test(shannon ~ Side, data = mydata)$p.value
KP_SID=kruskal.test(shannon ~ SID, data = mydata)$p.value
KP_VisitTreat=kruskal.test(shannon ~ VisitTreat2, data = mydata)$p.value

gBL=mydata$shannon[mydata$Visit =="BL"]
gW8A=mydata$shannon[mydata$VisitTreat =="W8.A"]
gW8B=mydata$shannon[mydata$VisitTreat =="W8.B"]
gW8C=mydata$shannon[mydata$VisitTreat =="W8.C"]
gW8D=mydata$shannon[mydata$VisitTreat =="W8.D"]
gW8E=mydata$shannon[mydata$VisitTreat =="W8.E"]
pBL_W8A=my.t.test.p.value(gBL, gW8A, paired=FALSE)
pBL_W8B=my.t.test.p.value(gBL, gW8B, paired=FALSE)
pBL_W8C=my.t.test.p.value(gBL, gW8C, paired=FALSE)
pBL_W8D=my.t.test.p.value(gBL, gW8D, paired=FALSE)
pBL_W8E=my.t.test.p.value(gBL, gW8E, paired=FALSE)
pW8A_W8B=my.t.test.p.value(gW8A, gW8B, paired=FALSE)
pW8A_W8C=my.t.test.p.value(gW8A, gW8C, paired=FALSE)
pW8A_W8D=my.t.test.p.value(gW8A, gW8D, paired=FALSE)
pW8A_W8E=my.t.test.p.value(gW8A, gW8E, paired=FALSE)

mBL=mean(gBL)
mW8A=mean(gW8A)
mW8B=mean(gW8B)
mW8C=mean(gW8C)
mW8D=mean(gW8D)
mW8E=mean(gW8E)

word_shannon=paste("shannon - ttest", 
              "|PostAcneSpotJSVSpearmanCorrelation_rho =", round(spearman_rho, 2), 
              "|PostAcneSpotJSVPearsonCorrelation_rho =", round(pearson_rho, 2), 
              "|PostAcneSpotJSVSpearmanCorrelation_pval =", sprintf("%.02f",spearman_pval),
              "|PostAcneSpotJSVPearsonCorrelation_pval =", sprintf("%.02f",pearson_pval),
              "|KP.Treat=", sprintf("%.02f",KP_Treat),
              "|KP.Side=", sprintf("%.02f",KP_Side),
              "|KP.SID=", sprintf("%.02f",KP_SID), 
              "|KP.Visit=", sprintf("%.02f",KP_Visit),
              "|KP.VisitTreat=", sprintf("%.02f",KP_VisitTreat),
              "|Wilcox.BL_vs_W8A=", sprintf("%.02f",pBL_W8A),
              "|Wilcox.BL_vs_W8B=", sprintf("%.02f",pBL_W8B),
              "|Wilcox.BL_vs_W8C=", sprintf("%.02f",pBL_W8C),
              "|Wilcox.BL_vs_W8D=", sprintf("%.02f",pBL_W8D),
              "|Wilcox.BL_vs_W8E=", sprintf("%.02f",pBL_W8E),
              "|Wilcox.W8A_vs_W8B=", sprintf("%.02f",pW8A_W8B),
              "|Wilcox.W8A_vs_W8C=", sprintf("%.02f",pW8A_W8C),
              "|Wilcox.W8A_vs_W8D=", sprintf("%.02f",pW8A_W8D),
              "|Wilcox.W8A_vs_W8E=", sprintf("%.02f",pW8A_W8E),
              "|meanBL=",mBL,  "|meanW8A=",mW8A, "|meanW8B=",mW8B,
              "|meanW8C=",mW8C,  "|meanW8D=",mW8D, "|meanW8E=",mW8E
)

#print(word_shannon)
print(paste("shannon -ttest",
            "BL.spearman_rho", "BL.pearson_rho","BL.spearman_pval", "BL.pearson_pval", 
            "KP_SID","KP_Side","KP_Treat", "KP_Visit","KP_VisitTreat",
            "pBL_W8A","pBL_W8B","pBL_W8C","pBL_W8D","pBL_W8E","pW8A_W8B","pW8A_W8C","pW8A_W8D","pW8A_W8E",
            "mBL", "mW8A", "mW8B","mW8C", "mW8D", "mW8E",collapse = "\t"))
print(paste("shannon -ttest",
            spearman_rho, pearson_rho,spearman_pval, pearson_pval, 
            KP_SID,
            KP_Side,
            KP_Treat, 
            KP_Visit,
            KP_VisitTreat,
            pBL_W8A,
            pBL_W8B,
            pBL_W8C,
            pBL_W8D,
            pBL_W8E,
            pW8A_W8B,
            pW8A_W8C,
            pW8A_W8D,
            pW8A_W8E,
            mBL,  mW8A, mW8B,mW8C, mW8D, mW8E,
            collapse = "\t"
))

######################Now for observed, do the same ###############################
result <- cor.test(mydataBL$BaselineValue.PostAcneSpotJSV, mydataBL$observed, method = "spearman", exact = FALSE)
spearman_rho <- result$estimate
spearman_pval <- result$p.value
result2 <- cor.test(mydataBL$BaselineValue.PostAcneSpotJSV, mydataBL$observed, method = "pearson", exact = FALSE)
pearson_rho <- result2$estimate
pearson_pval <- result2$p.value

KP_Visit=kruskal.test(observed ~ Visit, data = mydata)$p.value
KP_Treat=kruskal.test(observed ~ Trt, data = mydata)$p.value
KP_Side=kruskal.test(observed ~ Side, data = mydata)$p.value
KP_SID=kruskal.test(observed ~ SID, data = mydata)$p.value
KP_VisitTreat=kruskal.test(observed ~ VisitTreat2, data = mydata)$p.value

gBL=mydata$observed[mydata$Visit =="BL"]
gW8A=mydata$observed[mydata$VisitTreat =="W8.A"]
gW8B=mydata$observed[mydata$VisitTreat =="W8.B"]
gW8C=mydata$observed[mydata$VisitTreat =="W8.C"]
gW8D=mydata$observed[mydata$VisitTreat =="W8.D"]
gW8E=mydata$observed[mydata$VisitTreat =="W8.E"]
pBL_W8A=my.t.test.p.value(gBL, gW8A, paired=FALSE)
pBL_W8B=my.t.test.p.value(gBL, gW8B, paired=FALSE)
pBL_W8C=my.t.test.p.value(gBL, gW8C, paired=FALSE)
pBL_W8D=my.t.test.p.value(gBL, gW8D, paired=FALSE)
pBL_W8E=my.t.test.p.value(gBL, gW8E, paired=FALSE)
pW8A_W8B=my.t.test.p.value(gW8A, gW8B, paired=FALSE)
pW8A_W8C=my.t.test.p.value(gW8A, gW8C, paired=FALSE)
pW8A_W8D=my.t.test.p.value(gW8A, gW8D, paired=FALSE)
pW8A_W8E=my.t.test.p.value(gW8A, gW8E, paired=FALSE)

mBL=mean(gBL)
mW8A=mean(gW8A)
mW8B=mean(gW8B)
mW8C=mean(gW8C)
mW8D=mean(gW8D)
mW8E=mean(gW8E)

word_observed=paste("observed - ttest", 
                   "|PostAcneSpotJSVSpearmanCorrelation_rho =", round(spearman_rho, 2), 
                   "|PostAcneSpotJSVPearsonCorrelation_rho =", round(pearson_rho, 2), 
                   "|PostAcneSpotJSVSpearmanCorrelation_pval =", sprintf("%.02f",spearman_pval),
                   "|PostAcneSpotJSVPearsonCorrelation_pval =", sprintf("%.02f",pearson_pval),
                   "|KP.Treat=", sprintf("%.02f",KP_Treat),
                   "|KP.Side=", sprintf("%.02f",KP_Side),
                   "|KP.SID=", sprintf("%.02f",KP_SID), 
                   "|KP.Visit=", sprintf("%.02f",KP_Visit),
                   "|KP.VisitTreat=", sprintf("%.02f",KP_VisitTreat),
                   "|Wilcox.BL_vs_W8A=", sprintf("%.02f",pBL_W8A),
                   "|Wilcox.BL_vs_W8B=", sprintf("%.02f",pBL_W8B),
                   "|Wilcox.BL_vs_W8C=", sprintf("%.02f",pBL_W8C),
                   "|Wilcox.BL_vs_W8D=", sprintf("%.02f",pBL_W8D),
                   "|Wilcox.BL_vs_W8E=", sprintf("%.02f",pBL_W8E),
                   "|Wilcox.W8A_vs_W8B=", sprintf("%.02f",pW8A_W8B),
                   "|Wilcox.W8A_vs_W8C=", sprintf("%.02f",pW8A_W8C),
                   "|Wilcox.W8A_vs_W8D=", sprintf("%.02f",pW8A_W8D),
                   "|Wilcox.W8A_vs_W8E=", sprintf("%.02f",pW8A_W8E),
                   "|meanBL=",mBL,  "|meanW8A=",mW8A, "|meanW8B=",mW8B,
                   "|meanW8C=",mW8C,  "|meanW8D=",mW8D, "|meanW8E=",mW8E
)

#print(word_shannon)
print(paste("observed -ttest",
            "BL.spearman_rho", "BL.pearson_rho","BL.spearman_pval", "BL.pearson_pval", 
            "KP_SID","KP_Side","KP_Treat", "KP_Visit","KP_VisitTreat",
            "pBL_W8A","pBL_W8B","pBL_W8C","pBL_W8D","pBL_W8E","pW8A_W8B","pW8A_W8C","pW8A_W8D","pW8A_W8E",
            "mBL", "mW8A", "mW8B","mW8C", "mW8D", "mW8E",collapse = "\t"))
print(paste("observed -ttest",
            spearman_rho, pearson_rho,spearman_pval, pearson_pval, 
            KP_SID,
            KP_Side,
            KP_Treat, 
            KP_Visit,
            KP_VisitTreat,
            pBL_W8A,
            pBL_W8B,
            pBL_W8C,
            pBL_W8D,
            pBL_W8E,
            pW8A_W8B,
            pW8A_W8C,
            pW8A_W8D,
            pW8A_W8E,
            mBL,  mW8A, mW8B,mW8C, mW8D, mW8E,
            collapse = "\t"
))



###########################################
library(rstatix)
library(ggpubr)
library(dplyr)

# run t-test and keep only p <= 0.1
stat.observed <- mydata %>%
  t_test(observed ~ VisitTreat2) %>%
  filter(p <= 0.1) %>%   # <— filter here
  mutate(y.position = 1.05 * max(mydata$observed))

# plot
png(filename=paste0(filename, ".observed.png"), width=1600, height=1600, res=300)
pobserved <- ggboxplot(mydata, x = "VisitTreat2", y = "observed", color = "VisitTreat2", add = "jitter", add.params = list(alpha = 0.6, width = 0.2))

pobserved <- pobserved +
  stat_pvalue_manual(stat.observed, label = "p", tip.length = 0.01, step.increase = 0.1) +
  labs(x = NULL) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) 
print(pobserved)
dev.off()
addWorksheet(wb, "observed");             writeData(wb, "observed", stat.observed)

#######################
doc <- read_pptx()

# Add slide with both plots and caption
doc <- doc %>%
  add_slide(layout = "Title and Content", master = "Office Theme") %>%
  ph_with(
    fpar(ftext("Alpha Diversity: observed", prop = fp_text(font.size = 20, bold = TRUE))),
    location = ph_location_type(type = "title")
  ) %>%
  ph_with(dml(ggobj = pobserved), location = ph_location(left = 0.5, top = 1.5, width = 4.8, height = 6)) %>%
  ph_with(
    fpar(ftext(word_observed, prop = fp_text(font.size = 10))),
    location = ph_location(left = 5.5, top = 0.5, width = 4.6, height = 6)
  )
print(doc, target = paste0(filename, ".observed.pptx"))

#####################################################################
# run t-test and keep only p <= 0.1
stat.shannon <- mydata %>%
  t_test(shannon ~ VisitTreat2) %>%
  filter(p <= 0.1) %>%   # <— filter here
  mutate(y.position = 1.05 * max(mydata$shannon))

# plot
png(filename=paste0(filename, ".shannn.png"), width=1600, height=1600, res=300)
pshannon <- ggboxplot(mydata, x = "VisitTreat2", y = "shannon", color = "VisitTreat2", add = "jitter", add.params = list(alpha = 0.6, width = 0.2))

pshannon <- pshannon +
  stat_pvalue_manual(stat.shannon, label = "p", tip.length = 0.01, step.increase = 0.1) +
  labs(x = NULL) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) 
print(pshannon)
dev.off()
addWorksheet(wb, "shannon");             writeData(wb, "shannon", stat.shannon)

#######################
# Add slide with both plots and caption
doc <- doc %>%
  add_slide(layout = "Title and Content", master = "Office Theme") %>%
  ph_with(
    fpar(ftext("Alpha Diversity: shannon", prop = fp_text(font.size = 20, bold = TRUE))),
    location = ph_location_type(type = "title")
  ) %>%
  ph_with(dml(ggobj = pshannon), location = ph_location(left = 0.5, top = 1.5, width = 4.8, height = 6)) %>%
  ph_with(
    fpar(ftext(word_shannon, prop = fp_text(font.size = 10))),
    location = ph_location(left = 5.5, top = 0.5, width = 4.6, height = 6)
  )
print(doc, target = paste0(filename, ".shannon.pptx"))
########################################beta Diversity ###################
x<-test_d
x_clean <- x[rowSums(x) > 0, ]
beta_dist <- vegdist(x_clean, index = "bray")

#beta_dist <- vegdist(x, index = "bray")
mds <- metaMDS(beta_dist)
ttt<-as.data.frame(mds$points)

Rclean <-sub(".metaphlan4", "", rownames(ttt))
rownames(ttt) =Rclean
ttt$Cname_clean <- rownames(ttt)
mds_data <- merge(ttt, mydata, by = "Cname_clean", all.x = TRUE, sort = FALSE)
rownames(mds_data) <- mds_data$Cname_clean


#install.packages('devtools')
library(devtools)
#install_github('fawda123/ggord')
library(ggord)
library(ggplot2)

 
png(filename=paste0(filename,".beta.png"),  width=1800, height=800, res=300)
p2<-ggplot(mds_data, aes(x = MDS1, y = MDS2, color = VisitTreat2)) +geom_point(size=3, alpha = 0.7)+
  theme_bw()+stat_chull(aes(color=VisitTreat2, fill=VisitTreat2), alpha=0.01, geom="polygon")
print(p2)
dev.off() 
doc <- doc %>%
  add_slide(layout = "Title and Content", master = "Office Theme") %>%
  ph_with(value = "beta diversity by Group", location = ph_location_type(type = "title")) %>%
  # Add left plot (p1)
  ph_with(dml(ggobj = p2), location = ph_location(left = 0.5, top = 1.2, width = 8, height = 6))

png(filename=paste0(filename,".beta.TreatXVisit.png"), width=1800, height=800, res=300)
p2Treat2Visit<-ggplot(mds_data, aes(x = MDS1, y = MDS2, color = Trt)) +geom_point(size=3, alpha = 0.7)+
  theme_bw()+stat_chull(aes(color=Trt, fill=Trt), alpha=0.01, geom="polygon")+
  facet_wrap(~ Visit) 
print(p2Treat2Visit)
dev.off() 

doc <- doc %>%
  add_slide(layout = "Title and Content", master = "Office Theme") %>%
  ph_with(value = "beta diversity by Visit", location = ph_location_type(type = "title")) %>%
  # Add left plot (p1)
  ph_with(dml(ggobj = p2Treat2Visit), location = ph_location(left = 0.5, top = 1.2, width = 8, height = 6))



############################################################
print("all samples")
Rclean <-sub(".metaphlan4", "", rownames(x_clean))
idx <- match(Rclean, mydata$Cname_clean) 
new_names <- mydata$NewID[idx]
rownames(x_clean) <- new_names
mydata<-mydata[idx,]
# calculate Bray-Curtis distance among samples
bc.dist <- vegdist(x_clean, method = "bray")
# cluster communities using average-linkage algorithm
bc.clust <- hclust(bc.dist, method = "average")

# 1) Save the cluster as a crisp PNG (300 dpi, matches PPT size 9x4.5 in)
img_path <- paste0(filename, ".bray-Curtis-cluster.png")
png(filename = img_path, width = 9, height = 6, units = "in", res = 300)
plot(bc.clust, ylab = "Bray–Curtis dissimilarity")
dev.off()

# 2) Add the PNG to the slide
doc <- doc %>%
  add_slide(layout = "Title and Content", master = "Office Theme") %>%
  ph_with("bray curtis cluster", location = ph_location_type(type = "title")) %>%
  ph_with(
    external_img(img_path, width = 9, height = 5),
    location = ph_location(left = 0.5, top = 1.2)
  )

extract_adonis2 <- function(fit) {
  # find p-value column
  pcol <- grep("^Pr\\(>F\\)$", colnames(fit), value = TRUE)
  if (length(pcol) == 0) stop("Couldn't find p-value column in adonis2 result.")
  
  # raw numeric p
  p_raw <- as.numeric(fit[[pcol]][1])
  
  # star like adonis2: ***, **, *, ., ' '
  star <- if (p_raw <= 0.001) "***"
  else if (p_raw <= 0.01) "**"
  else if (p_raw <= 0.05) "*"
  else if (p_raw <= 0.1) "."
  else ""
  
  # string formatted exactly like the table (3 decimals)
  p_str <- formatC(p_raw, format = "f", digits = 3)
  
  list(p = p_raw, p_str = p_str, star = star)
}

# Example
#res <- adonis2(bc.dist ~ mydata$SiteDiseaseLesion)
#out <- extract_adonis2(res)
#out$p       # raw numeric, e.g. 0.037037...
#out$p_str   # "0.037"
#out$star    # "*"

# Taxonomic (Bray-Curtis) dissimilarity explained
t1t="adonis result bc.dist ~ VisitTreat2"
t1=adonis2(bc.dist ~ mydata$VisitTreat2)
t2t="adonis result bc.dist ~ Visit"
t2=adonis2(bc.dist ~ mydata$Visit) #0.001
#t3t="adonis result bc.dist ~ Visit treat"
#t3=adonis2(bc.dist ~ mydata$VisitTreat) #0.085
t4t="adonis result bc.dist ~ SID"
t4=adonis2(bc.dist ~ mydata$SID) #0.001
t5t="adonis result bc.dist ~ VisitTreat"
t5=adonis2(bc.dist ~ mydata$VisitTreat)
t6t="adonis result bc.dist ~ CheekSide"
t6=adonis2(bc.dist ~ mydata$Side) #0.122

myAdonis <- paste(t1t, extract_adonis2(t1)$p,"|",
                  t2t, extract_adonis2(t2)$p,"|",
                  #t3t, extract_adonis2(t3)$p,"|",
                  t4t, extract_adonis2(t4)$p,"|",
                  t5t, extract_adonis2(t5)$p,"|",
                  t6t, extract_adonis2(t6)$p
                  )

print(myAdonis)

doc <- doc %>%
  add_slide(layout = "Title and Content", master = "Office Theme") %>%
  ph_with("Adonis 2 Test P value", location = ph_location_type(type = "title")) %>%
  ph_with(
    value = fpar(ftext(myAdonis, fp_text(font.size = 12))),
    location = ph_location(left = 0.5, top = 6.2, width = 9, height = 1)
  )




# ######################################
# jaccard.dist <- vegdist(x, method = "jaccard")
# mds2 <- metaMDS(jaccard.dist)
# mds_data2 <- as.data.frame(cbind(mds2$points, mydata))
# mds_data2$SampleID <- rownames(mds_data2)
# 
# png(filename=paste0(filename,".beta.PH_Buttocks_Category.jaccard.png"),  width=1200, height=800, res=300)
# p1<-ggplot(mds_data2, aes(x = MDS1, y = MDS2, color = PH_Buttocks_Category)) +
#   geom_point(size=3, alpha=0.7)+theme_bw()+
#   stat_chull(aes(color=PH_Buttocks_Category, fill=PH_Buttocks_Category), alpha=0.1, geom="polygon") 
# print(p1)
# dev.off()
# 
# # test.adonis <- adonis2(bc.dist ~ RashPH)
# # test.adonis <- as.data.frame(test.adonis$aov.tab)
# # test.adonis
###################################################################
#PAIRWISE PERMANOVA


library(vegan)
cbn <- combn(x=unique(mds_data$VisitTreat2), m = 2)
p <- c()
for (i in 1:ncol(cbn)) {
  ps.subs <- x_clean[mds_data$VisitTreat2 %in% cbn[, i], ]
  metadata_sub <- mds_data[mds_data$VisitTreat2 %in% cbn[, i], ]
  permanova_pairwise <- adonis2(
    vegdist(ps.subs, method = "bray") ~ metadata_sub$VisitTreat2
  )
  p <- c(p, permanova_pairwise$`Pr(>F)`[1])
}
p
p.adj <- p.adjust(p, method = "BH")
p.table <- cbind.data.frame(cbn[1,],cbn[2,] ,p=p, p.adj=p.adj)
p.table
addWorksheet(wb, "betaDiversity_VisitTreat2");writeData(wb, "betaDiversity_VisitTreat2", p.table)

pairwise_ft <- flextable(p.table) %>%
  autofit() %>%
  set_table_properties(width = 1, layout = "autofit") %>%
  theme_booktabs() %>%
  align(align = "center", part = "all") %>%
  fontsize(size = 9, part = "all")

doc <- doc %>%
  add_slide(layout = "Title and Content", master = "Office Theme") %>%
  
  # Custom title with font size 36
  ph_with(
    value = fpar(ftext("Beta Diversity with Bray-Curtis Distance & Pairwise Permanova Test by VisitTreat2",
                       fp_text(font.size = 36))),
    location = ph_location(left = 0.5, top = 0.3, width = 9, height = 1)
  ) %>%
  # Add the pairwise permanova table next to the plot
  ph_with(pairwise_ft, 
          location = ph_location(left = 6.2, top = 2.0, width = 3.5, height = 5))
###############################################################
cbn <- combn(x=unique(mds_data$Visit), m = 2)
p <- c()
for (i in 1:ncol(cbn)) {
  ps.subs <- x_clean[mds_data$Visit %in% cbn[, i], ]
  metadata_sub <- mds_data[mds_data$Visit %in% cbn[, i], ]
  permanova_pairwise <- adonis2(
    vegdist(ps.subs, method = "bray") ~ metadata_sub$Visit
  )
  p <- c(p, permanova_pairwise$`Pr(>F)`[1])
}
p
p.adj <- p.adjust(p, method = "BH")
p.table <- cbind.data.frame(cbn[1,],cbn[2,] ,p=p, p.adj=p.adj)
p.table
addWorksheet(wb, "betaDiversity_Visit");writeData(wb, "betaDiversity_Visit", p.table)

pairwise_ft <- flextable(p.table) %>%
  autofit() %>%
  set_table_properties(width = 1, layout = "autofit") %>%
  theme_booktabs() %>%
  align(align = "center", part = "all") %>%
  fontsize(size = 9, part = "all")

doc <- doc %>%
  add_slide(layout = "Title and Content", master = "Office Theme") %>%
  
  # Custom title with font size 36
  ph_with(
    value = fpar(ftext("Beta Diversity with Bray-Curtis Distance & Pairwise Permanova Test by Visit",
                       fp_text(font.size = 36))),
    location = ph_location(left = 0.5, top = 0.3, width = 9, height = 1)
  ) %>%
  # Add the pairwise permanova table next to the plot
  ph_with(pairwise_ft, 
          location = ph_location(left = 6.2, top = 2.0, width = 3.5, height = 5))



############################################################
cbn <- combn(x=unique(mds_data$Side), m = 2)
p <- c()
for (i in 1:ncol(cbn)) {
  ps.subs <- x_clean[mds_data$Side %in% cbn[, i], ]
  metadata_sub <- mds_data[mds_data$Side %in% cbn[, i], ]
  permanova_pairwise <- adonis2(
    vegdist(ps.subs, method = "bray") ~ metadata_sub$Side
  )
  p <- c(p, permanova_pairwise$`Pr(>F)`[1])
}
p
p.adj <- p.adjust(p, method = "BH")
p.table <- cbind.data.frame(cbn[1,],cbn[2,] ,p=p, p.adj=p.adj)
p.table
addWorksheet(wb, "betaDiversity_Side");writeData(wb, "betaDiversity_Side", p.table)

pairwise_ft <- flextable(p.table) %>%
  autofit() %>%
  set_table_properties(width = 1, layout = "autofit") %>%
  theme_booktabs() %>%
  align(align = "center", part = "all") %>%
  fontsize(size = 9, part = "all")

doc <- doc %>%
  add_slide(layout = "Title and Content", master = "Office Theme") %>%
  
  # Custom title with font size 36
  ph_with(
    value = fpar(ftext("Beta Diversity with Bray-Curtis Side & Pairwise Permanova Test by Side",
                       fp_text(font.size = 36))),
    location = ph_location(left = 0.5, top = 0.3, width = 9, height = 1)
  ) %>%
  # Add the pairwise permanova table next to the plot
  ph_with(pairwise_ft, 
          location = ph_location(left = 6.2, top = 2.0, width = 3.5, height = 5))



######################################################################
# cbn <- combn(x=unique(mds_data$Trt), m = 2)
# p <- c()
# for (i in 1:ncol(cbn)) {
#   ps.subs <- x_clean[mds_data$Trt %in% cbn[, i], ]
#   metadata_sub <- mds_data[mds_data$Trt %in% cbn[, i], ]
#   permanova_pairwise <- adonis2(
#     vegdist(ps.subs, method = "bray") ~ metadata_sub$Trt
#   )
#   p <- c(p, permanova_pairwise$`Pr(>F)`[1])
# }
# p
# p.adj <- p.adjust(p, method = "BH")
# p.table <- cbind.data.frame(cbn[1,],cbn[2,] ,p=p, p.adj=p.adj)
# p.table
# addWorksheet(wb, "betaDiversity_Trt");writeData(wb, "betaDiversity_Trt", p.table)
# 
# pairwise_ft <- flextable(p.table) %>%
#   autofit() %>%
#   set_table_properties(width = 1, layout = "autofit") %>%
#   theme_booktabs() %>%
#   align(align = "center", part = "all") %>%
#   fontsize(size = 9, part = "all")
# 
# doc <- doc %>%
#   add_slide(layout = "Title and Content", master = "Office Theme") %>%
#   
#   # Custom title with font size 36
#   ph_with(
#     value = fpar(ftext("Beta Diversity with Bray-Curtis Distance & Pairwise Permanova Test by Treat",
#                        fp_text(font.size = 36))),
#     location = ph_location(left = 0.5, top = 0.3, width = 9, height = 1)
#   ) %>%
#   # Add the pairwise permanova table next to the plot
#   ph_with(pairwise_ft, 
#           location = ph_location(left = 6.2, top = 2.0, width = 3.5, height = 5))
# 
# 

########################################################
cbn <- combn(x=unique(mds_data$VisitTreat), m = 2)
p <- c()
for (i in 1:ncol(cbn)) {
  ps.subs <- x_clean[mds_data$VisitTreat %in% cbn[, i], ]
  metadata_sub <- mds_data[mds_data$VisitTreat %in% cbn[, i], ]
  permanova_pairwise <- adonis2(
    vegdist(ps.subs, method = "bray") ~ metadata_sub$VisitTreat
  )
  p <- c(p, permanova_pairwise$`Pr(>F)`[1])
}
p
p.adj <- p.adjust(p, method = "BH")
p.table <- cbind.data.frame(cbn[1,],cbn[2,] ,p=p, p.adj=p.adj)
p.table
addWorksheet(wb, "betaDiversity_VisitTreat");writeData(wb, "betaDiversity_VisitTreat", p.table)

pairwise_ft <- flextable(p.table) %>%
  autofit() %>%
  set_table_properties(width = 1, layout = "autofit") %>%
  theme_booktabs() %>%
  align(align = "center", part = "all") %>%
  fontsize(size = 9, part = "all")

doc <- doc %>%
  add_slide(layout = "Title and Content", master = "Office Theme") %>%
  
  # Custom title with font size 36
  ph_with(
    value = fpar(ftext("Beta Diversity with Bray-Curtis Distance & Pairwise Permanova Test by VisitTreat",
                       fp_text(font.size = 36))),
    location = ph_location(left = 0.5, top = 0.3, width = 9, height = 1)
  ) %>%
  # Add the pairwise permanova table next to the plot
  ph_with(pairwise_ft, 
          location = ph_location(left = 6.2, top = 2.0, width = 3.5, height = 5))



############################################################################
##################################################################
saveWorkbook(wb, "Diversity_analysis.xlsx", overwrite = TRUE)
print(doc, target = "diversity_summary.pptx")
