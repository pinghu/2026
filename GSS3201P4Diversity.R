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
filename1="cleanMeta.txt"

args <- commandArgs(trailingOnly = TRUE)
filename<- args[1]
filename="GSS3201_salmonHumanSkin.ReadsCount"
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
idx <- match(A$cleanMeta, colnames(test_filter2))

# Optional safety check
if (any(is.na(idx))) {
  stop("Some A$cleanMeta values are not found in colnames(test_filter2).")
}

test_filter2_mat <- as.matrix(test_filter2[, idx])
new_names <- as.character(A$NewID)  # use the real column, capital D
stopifnot(length(new_names) == ncol(test_filter2_mat))

colnames(test_filter2_mat) <- new_names

#colnames(test_filter2_mat)  # should now show your No.Arm..., Yes.Cheek..., etc.

relative_fraction_mat <- sweep(
  test_filter2_mat,
  MARGIN = 2,                    # operate over columns
  STATS  = A$P4.original.total.P4,
  FUN    = "/"
)
## 4. If you want *percentages* instead of fractions:
relative_percent_mat <- relative_fraction_mat * 100
# ---- Save as tab-delimited file ----
write.table(
  relative_percent_mat,
  file = "relative_percent_mat.txt",
  sep = "\t",
  quote = FALSE,
  row.names = TRUE,
  col.names = TRUE
)
write.table(
  test_filter2_mat,
  file = "test_filter2_mat.txt",
  sep = "\t",
  quote = FALSE,
  row.names = TRUE,
  col.names = TRUE
)
write.table(
  test_filter2,
  file = "test_filter2.txt",
  sep = "\t",
  quote = FALSE,
  row.names = TRUE,
  col.names = TRUE
)

# ---- Save as Excel file ----
library(openxlsx)

write.xlsx(
  relative_percent_mat,
  file = "relative_percent_mat.xlsx",
  rowNames = TRUE
)
write.xlsx(
  test_filter2_mat,
  file = "test_filter2_mat.xlsx",
  rowNames = TRUE
)
write.xlsx(
  test_filter2,
  file = "test_filter2.xlsx",
  rowNames = TRUE
)


test_d=data.frame(t(test))

min(apply(test_d, 1, sum)) ###162338

test_a <- decostand(test_d, method = "total")

Cname=rownames(test_d)
Cname_clean <-sub(".metaphlan4", "", Cname)
Clen=length(Cname) ##there are 3 annotation columns

Clen=length(Cname_clean) ##there are 3 annotation columns
Lesion=rep("NA", Clen)
Site=rep("NA", Clen)
SID=rep("NA", Clen)
splitname<-strsplit(Cname_clean, "[_]")


for(mm in  1:Clen ){
  Site[mm]=splitname[[mm]][4]
  SID[mm]=splitname[[mm]][2]
  Lesion[mm]=splitname[[mm]][3]
}

Site[Site=="LC"]="C"
Site[Site=="RC"]="C"
SiteLesion=paste0(Site, Lesion)

matching_dataset <- A %>%
  filter(Cname_clean %in% cleanMeta) %>%
  arrange(match(Cname_clean, cleanMeta))


shannon<-diversity(test_d)
simp<-diversity(test_d, "simpson")
invsimp<-diversity(test_d, "inv")

observed<-apply(test_d>0,1,sum)
N <- apply(test_d,1,sum)
###Richness Index#####
Menhinick_index<-observed/sqrt(N)
Margalef_index <-(observed-1)/log(N)

mydata0=data.frame(shannon,simp,invsimp, observed, Menhinick_index, Margalef_index ,SID, Site, Lesion, SiteLesion, Cname_clean)

mydata <- inner_join(mydata0, matching_dataset, by = c("Cname_clean" = "cleanMeta"))
mydata$Disease=mydata$AtopicDermatitis
mydata$Disease[mydata$AtopicDermatitis=="Yes"]="AD"
mydata$Disease[mydata$AtopicDermatitis=="No"]="HL"
mydata$Site[mydata$Site=="A"]="Arm"
mydata$Site[mydata$Site=="C"]="Cheek"
mydata$Lesion[mydata$Lesion=="Les"]="Lesion"
mydata$Lesion[mydata$Lesion=="NL"]="Nonlesion"
mydata$SID=mydata$SubjectID
mydata$SiteLesion=paste0(mydata$Site, mydata$Lesion)
mydata$SiteDisease=paste0(mydata$Site, mydata$Disease)
mydata$SiteDiseaseLesion=paste0(mydata$SiteDisease, mydata$Lesion)
mydata$DiseaseLesion=paste0(mydata$Disease, mydata$Lesion)
order_levels <- c("ArmADLesion", "ArmADNonlesion", "ArmHLNonlesion",   "CheekADLesion",  "CheekADNonlesion",   "CheekHLNonlesion") 
mydata$SiteDiseaseLesion <- factor(mydata$SiteDiseaseLesion, levels = order_levels, ordered = TRUE)

write.table(t(mydata), file =paste0(filename, ".alphadiversity.RShort"), col.names=TRUE, row.names=TRUE, sep = "\t")
addWorksheet(wb, "aDiversity");             writeData(wb, "aDiversity", mydata)

############################################
#https://www.datanovia.com/en/blog/ggpubr-how-to-add-p-values-generated-elsewhere-to-a-ggplot/
result <- cor.test(mydata$IGA, mydata$shannon, method = "spearman", exact = FALSE)
spearman_rho <- result$estimate
spearman_pval <- result$p.value
result2 <- cor.test(mydata$IGA, mydata$shannon, method = "pearson", exact = FALSE)
pearson_rho <- result2$estimate
pearson_pval <- result2$p.value

Aresult <- cor.test(mydata$IGA[mydata$Site=="Arm"], mydata$shannon[mydata$Site=="Arm"], method = "spearman", exact = FALSE)
Aspearman_rho <- Aresult$estimate
Aspearman_pval <- Aresult$p.value
Aresult2 <- cor.test(mydata$IGA[mydata$Site=="Arm"], mydata$shannon[mydata$Site=="Arm"], method = "pearson", exact = FALSE)
Apearson_rho <- Aresult2$estimate
Apearson_pval <- Aresult2$p.value

Cresult <- cor.test(mydata$IGA[mydata$Site=="Cheek"], mydata$shannon[mydata$Site=="Cheek"], method = "spearman", exact = FALSE)
Cspearman_rho <- Cresult$estimate
Cspearman_pval <- Cresult$p.value
Cresult2 <- cor.test(mydata$IGA[mydata$Site=="Cheek"], mydata$shannon[mydata$Site=="Cheek"], method = "pearson", exact = FALSE)
Cpearson_rho <- Cresult2$estimate
Cpearson_pval <- Cresult2$p.value


KP_SiteDiseaseLesion=kruskal.test(shannon ~ SiteDiseaseLesion, data = mydata)$p.value
KP_Site=kruskal.test(shannon ~ Site, data = mydata)$p.value
KP_SID=kruskal.test(shannon ~ SID, data = mydata)$p.value
KP_Lesion=kruskal.test(shannon ~ Lesion, data = mydata)$p.value
KP_Disease=kruskal.test(shannon ~ Disease, data = mydata)$p.value
KP_SiteLesion=kruskal.test(shannon ~ SiteLesion, data = mydata)$p.value
KP_SiteDisease=kruskal.test(shannon ~ SiteDisease, data = mydata)$p.value
KP_DiseaseLesion=kruskal.test(shannon ~ DiseaseLesion, data = mydata)$p.value
gArmADLesion=mydata$shannon[mydata$SiteDiseaseLesion =="ArmADLesion"]
gArmADNonlesion=mydata$shannon[mydata$SiteDiseaseLesion =="ArmADNonlesion"]
pArm=my.t.test.p.value(gArmADLesion, gArmADNonlesion, paired=FALSE)
gCheekADLesion=mydata$shannon[mydata$SiteDiseaseLesion =="CheekADLesion"]
gCheekADNonlesion=mydata$shannon[mydata$SiteDiseaseLesion =="CheekADNonlesion"]
pCheek=my.t.test.p.value(gCheekADLesion, gCheekADNonlesion, paired=FALSE)
gCheekHL =mydata$shannon[mydata$SiteDiseaseLesion =="CheekHLNonlesion"]
gArmHL =mydata$shannon[mydata$SiteDiseaseLesion =="ArmHLNonlesion"]
pArmDLes_HL =my.t.test.p.value(gArmADLesion, gArmHL, paired=FALSE)
pArmDNL_HL =my.t.test.p.value(gArmADNonlesion, gArmHL, paired=FALSE)
pCheekDLes_HL =my.t.test.p.value(gCheekADLesion, gCheekHL, paired=FALSE)
pCheekDNL_HL =my.t.test.p.value(gCheekADNonlesion, gCheekHL, paired=FALSE) 
gArmAD= mydata$shannon[mydata$SiteDisease =="ArmAD"]
pArmAD_HL =my.t.test.p.value(gArmAD, gArmHL, paired=FALSE)
gCheekAD= mydata$shannon[mydata$SiteDisease =="CheekAD"]
pCheekAD_HL =my.t.test.p.value(gCheekAD, gCheekHL, paired=FALSE)
gArmNonlesion=mydata$shannon[mydata$SiteLesion =="ArmNonlesion"]
pArmLes_NL=my.t.test.p.value(gArmADLesion, gArmNonlesion, paired=FALSE)
gCheekNonlesion=mydata$shannon[mydata$SiteLesion =="CheekNonlesion"]
pCheekLes_NL=my.t.test.p.value(gCheekADLesion, gCheekNonlesion, paired=FALSE)
gCheek=mydata$shannon[mydata$Site =="Cheek"]
gArm=mydata$shannon[mydata$Site =="Arm"]
pCheek_Arm=my.t.test.p.value(gCheek, gArm, paired=FALSE)
pLesCheek_Arm=my.t.test.p.value(gCheekADLesion, gArmADLesion, paired=FALSE)
pADNLCheek_Arm=my.t.test.p.value(gCheekADNonlesion, gArmADNonlesion, paired=FALSE)
pHLCheek_Arm=my.t.test.p.value(gCheekHL, gArmHL, paired=FALSE)
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
gADLes=mydata$shannon[mydata$DiseaseLesion =="ADLesion"]
gADNL=mydata$shannon[mydata$DiseaseLesion =="ADNonlesion"]
gHL=mydata$shannon[mydata$DiseaseLesion =="HLNonlesion"]
mADLes=mean(gADLes)
mADNL=mean(gADNL)
mHL=mean(gHL)
pADLes_HL =my.t.test.p.value(gADLes, gHL, paired=FALSE)
pADLes_ADNL =my.t.test.p.value(gADLes, gADNL, paired=FALSE)
pADNL_HL =my.t.test.p.value(gADNL, gHL, paired=FALSE)

word_shannon=paste("shannon - ttest", 
              "|IGASpearmanCorrelation_rho =", round(spearman_rho, 2), 
              "|IGAPearsonCorrelation_rho =", round(pearson_rho, 2), 
              "|Arm-IGASpearmanCorrelation_rho =", round(Aspearman_rho, 2), 
              "|Arm-IGAPearsonCorrelation_rho =", round(Apearson_rho, 2),
              "|Cheek-IGASpearmanCorrelation_rho =", round(Cspearman_rho, 2), 
              "|Cheek-IGAPearsonCorrelation_rho =", round(Cpearson_rho, 2), 
              "|IGASpearmanCorrelation_pval =", sprintf("%.02f",spearman_pval),
              "|PearsonCorrelation_pval =", sprintf("%.02f",pearson_pval),
              "|Arm-IGASpearmanCorrelation_pval =", sprintf("%.02f",Aspearman_pval),
              "|Arm-PearsonCorrelation_pval =", sprintf("%.02f",Apearson_pval),
              "|Cheek-IGASpearmanCorrelation_pval =", sprintf("%.02f",Cspearman_pval),
              "|Cheek-PearsonCorrelation_pval =", sprintf("%.02f",Cpearson_pval),
              "|KP.SiteDiseaseLesion=", sprintf("%.02f",KP_SiteDiseaseLesion),
              "|KP.site=", sprintf("%.02f",KP_Site),
              "|KP.disease=", sprintf("%.02f",KP_Disease), 
              "|KP.lesion=", sprintf("%.02f",KP_Lesion),
              "|KP.SID=", sprintf("%.02f",KP_SID),
              "|KP.SiteDisease=", sprintf("%.02f",KP_SiteDisease),
              "|KP.SiteLesion=", sprintf("%.02f",KP_SiteLesion),
              "|KP.DiseaseLesion=", sprintf("%.02f",KP_DiseaseLesion),
              "|Wilcox.ArmAD.Lesion_vs_Nonlesion=", sprintf("%.02f",pArm),
              "|ttest.Arm.ADLesion_vs_HL=", sprintf("%.02f",pArmDLes_HL),
              "|ttest.Arm.ADNonlesion_vs_HL=", sprintf("%.02f",pArmDNL_HL),
              "|ttest.Arm.AD_vs_HL=", sprintf("%.02f",pArmAD_HL),
              "|ttest.Arm.Lesion_vs_Nonlesion=", sprintf("%.02f",pArmLes_NL),
              "|ttest.CheekAD.Lesion_vs_Nonlesion=", sprintf("%.02f",pCheek),
              "|ttest.Cheek.ADLesion_vs_HL=", sprintf("%.02f",pCheekDLes_HL),
              "|ttest.Cheek.ADNonlesion_vs_HL=", sprintf("%.02f",pCheekDNL_HL),
              "|ttest.Cheek.AD_vs_HL=", sprintf("%.02f",pCheekAD_HL),
              "|ttest.Cheek.Lesion_vs_Nonlesion=", sprintf("%.02f",pCheekLes_NL),
              "|ttest.AllAD.Lesion_vs_Nonlesion=", sprintf("%.02f",pADLes_ADNL),
              "|ttest.All.ADLesion_vs_HL=", sprintf("%.02f",pADLes_HL),
              "|ttest.All.ADNonlesion_vs_HL=", sprintf("%.02f",pADNL_HL),
              "|ttest.Cheek_vs_Arm=", sprintf("%.02f",pCheek_Arm),
              "|ttest.Lesion.Cheek_vs_Arm=", sprintf("%.02f",pLesCheek_Arm),
              "|ttest.ADNonlesion.Cheek_vs_Arm=", sprintf("%.02f",pADNLCheek_Arm),
              "|ttest.Healthy.Cheek_vs_Arm=", sprintf("%.02f",pHLCheek_Arm),
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
              "|meanHealthy=", mHL
)

#print(word_shannon)
print(paste("testname",
            "IGASpearmanCorrelation_rho","IGAPearsonCorrelation_rho", "IGASpearmanCorrelation_pval","IGAPearsonCorrelation_pval",
            "Arm-IGASpearmanCorrelation_rho","Arm-IGAPearsonCorrelation_rho", "Arm-IGASpearmanCorrelation_pval","Arm-IGAPearsonCorrelation_pval",
            "Cheek-IGASpearmanCorrelation_rho","Cheek-IGAPearsonCorrelation_rho", "Cheek-IGASpearmanCorrelation_pval","Cheek-IGAPearsonCorrelation_pval",
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
            collapse = "\t"
))
print(paste("shannon -ttest",
            spearman_rho, pearson_rho,spearman_pval, pearson_pval, 
            Aspearman_rho, Apearson_rho, Aspearman_pval, Apearson_pval,
            Cspearman_rho, Cpearson_rho, Cspearman_pval, Cpearson_pval,
            KP_SiteDiseaseLesion,
            KP_Site,
            KP_Disease, 
            KP_Lesion,
            KP_SID,
            KP_SiteDisease,
            KP_SiteLesion,
            KP_DiseaseLesion,
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
            collapse = "\t"
)  )
######################Now for observed, do the same ###############################
result <- cor.test(mydata$IGA, mydata$observed, method = "spearman", exact = FALSE)
spearman_rho <- result$estimate
spearman_pval <- result$p.value
result2 <- cor.test(mydata$IGA, mydata$observed, method = "pearson", exact = FALSE)
pearson_rho <- result2$estimate
pearson_pval <- result2$p.value

Aresult <- cor.test(mydata$IGA[mydata$Site=="Arm"], mydata$observed[mydata$Site=="Arm"], method = "spearman", exact = FALSE)
Aspearman_rho <- Aresult$estimate
Aspearman_pval <- Aresult$p.value
Aresult2 <- cor.test(mydata$IGA[mydata$Site=="Arm"], mydata$observed[mydata$Site=="Arm"], method = "pearson", exact = FALSE)
Apearson_rho <- Aresult2$estimate
Apearson_pval <- Aresult2$p.value

Cresult <- cor.test(mydata$IGA[mydata$Site=="Cheek"], mydata$observed[mydata$Site=="Cheek"], method = "spearman", exact = FALSE)
Cspearman_rho <- Cresult$estimate
Cspearman_pval <- Cresult$p.value
Cresult2 <- cor.test(mydata$IGA[mydata$Site=="Cheek"], mydata$observed[mydata$Site=="Cheek"], method = "pearson", exact = FALSE)
Cpearson_rho <- Cresult2$estimate
Cpearson_pval <- Cresult2$p.value

KP_SiteDiseaseLesion=kruskal.test(observed ~ SiteDiseaseLesion, data = mydata)$p.value
KP_Site=kruskal.test(observed ~ Site, data = mydata)$p.value
KP_SID=kruskal.test(observed ~ SID, data = mydata)$p.value
KP_Lesion=kruskal.test(observed ~ Lesion, data = mydata)$p.value
KP_Disease=kruskal.test(observed ~ Disease, data = mydata)$p.value
KP_SiteLesion=kruskal.test(observed ~ SiteLesion, data = mydata)$p.value
KP_SiteDisease=kruskal.test(observed ~ SiteDisease, data = mydata)$p.value
KP_DiseaseLesion=kruskal.test(observed ~ DiseaseLesion, data = mydata)$p.value
gArmADLesion=mydata$observed[mydata$SiteDiseaseLesion =="ArmADLesion"]
gArmADNonlesion=mydata$observed[mydata$SiteDiseaseLesion =="ArmADNonlesion"]
pArm=my.t.test.p.value(gArmADLesion, gArmADNonlesion, paired=FALSE)
gCheekADLesion=mydata$observed[mydata$SiteDiseaseLesion =="CheekADLesion"]
gCheekADNonlesion=mydata$observed[mydata$SiteDiseaseLesion =="CheekADNonlesion"]
pCheek=my.t.test.p.value(gCheekADLesion, gCheekADNonlesion, paired=FALSE)
gCheekHL =mydata$observed[mydata$SiteDiseaseLesion =="CheekHLNonlesion"]
gArmHL =mydata$observed[mydata$SiteDiseaseLesion =="ArmHLNonlesion"]
pArmDLes_HL =my.t.test.p.value(gArmADLesion, gArmHL, paired=FALSE)
pArmDNL_HL =my.t.test.p.value(gArmADNonlesion, gArmHL, paired=FALSE)
pCheekDLes_HL =my.t.test.p.value(gCheekADLesion, gCheekHL, paired=FALSE)
pCheekDNL_HL =my.t.test.p.value(gCheekADNonlesion, gCheekHL, paired=FALSE) 
gArmAD= mydata$observed[mydata$SiteDisease =="ArmAD"]
pArmAD_HL =my.t.test.p.value(gArmAD, gArmHL, paired=FALSE)
gCheekAD= mydata$observed[mydata$SiteDisease =="CheekAD"]
pCheekAD_HL =my.t.test.p.value(gCheekAD, gCheekHL, paired=FALSE)
gArmNonlesion=mydata$observed[mydata$SiteLesion =="ArmNonlesion"]
pArmLes_NL=my.t.test.p.value(gArmADLesion, gArmNonlesion, paired=FALSE)
gCheekNonlesion=mydata$observed[mydata$SiteLesion =="CheekNonlesion"]
pCheekLes_NL=my.t.test.p.value(gCheekADLesion, gCheekNonlesion, paired=FALSE)
gCheek=mydata$observed[mydata$Site =="Cheek"]
gArm=mydata$observed[mydata$Site =="Arm"]
pCheek_Arm=my.t.test.p.value(gCheek, gArm, paired=FALSE)
pLesCheek_Arm=my.t.test.p.value(gCheekADLesion, gArmADLesion, paired=FALSE)
pADNLCheek_Arm=my.t.test.p.value(gCheekADNonlesion, gArmADNonlesion, paired=FALSE)
pHLCheek_Arm=my.t.test.p.value(gCheekHL, gArmHL, paired=FALSE)
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
gADLes=mydata$observed[mydata$DiseaseLesion =="ADLesion"]
gADNL=mydata$observed[mydata$DiseaseLesion =="ADNonlesion"]
gHL=mydata$observed[mydata$DiseaseLesion =="HLNonlesion"]
mADLes=mean(gADLes)
mADNL=mean(gADNL)
mHL=mean(gHL)
pADLes_HL =my.t.test.p.value(gADLes, gHL, paired=FALSE)
pADLes_ADNL =my.t.test.p.value(gADLes, gADNL, paired=FALSE)
pADNL_HL =my.t.test.p.value(gADNL, gHL, paired=FALSE)

word_observed=paste("observed - ttest", 
                   "|IGASpearmanCorrelation_rho =", round(spearman_rho, 2), 
                   "|IGAPearsonCorrelation_rho =", round(pearson_rho, 2), 
                   "|Arm-IGASpearmanCorrelation_rho =", round(Aspearman_rho, 2), 
                   "|Arm-IGAPearsonCorrelation_rho =", round(Apearson_rho, 2),
                   "|Cheek-IGASpearmanCorrelation_rho =", round(Cspearman_rho, 2), 
                   "|Cheek-IGAPearsonCorrelation_rho =", round(Cpearson_rho, 2), 
                   "|IGASpearmanCorrelation_pval =", sprintf("%.02f",spearman_pval),
                   "|PearsonCorrelation_pval =", sprintf("%.02f",pearson_pval),
                   "|Arm-IGASpearmanCorrelation_pval =", sprintf("%.02f",Aspearman_pval),
                   "|Arm-PearsonCorrelation_pval =", sprintf("%.02f",Apearson_pval),
                   "|Cheek-IGASpearmanCorrelation_pval =", sprintf("%.02f",Cspearman_pval),
                   "|Cheek-PearsonCorrelation_pval =", sprintf("%.02f",Cpearson_pval),
                   
                   "|KP.SiteDiseaseLesion=", sprintf("%.02f",KP_SiteDiseaseLesion),
                   "|KP.site=", sprintf("%.02f",KP_Site),
                   "|KP.disease=", sprintf("%.02f",KP_Disease), 
                   "|KP.lesion=", sprintf("%.02f",KP_Lesion),
                   "|KP.SID=", sprintf("%.02f",KP_SID),
                   "|KP.SiteDisease=", sprintf("%.02f",KP_SiteDisease),
                   "|KP.SiteLesion=", sprintf("%.02f",KP_SiteLesion),
                   "|KP.DiseaseLesion=", sprintf("%.02f",KP_DiseaseLesion),
                   "|Wilcox.ArmAD.Lesion_vs_Nonlesion=", sprintf("%.02f",pArm),
                   "|ttest.Arm.ADLesion_vs_HL=", sprintf("%.02f",pArmDLes_HL),
                   "|ttest.Arm.ADNonlesion_vs_HL=", sprintf("%.02f",pArmDNL_HL),
                   "|ttest.Arm.AD_vs_HL=", sprintf("%.02f",pArmAD_HL),
                   "|ttest.Arm.Lesion_vs_Nonlesion=", sprintf("%.02f",pArmLes_NL),
                   "|ttest.CheekAD.Lesion_vs_Nonlesion=", sprintf("%.02f",pCheek),
                   "|ttest.Cheek.ADLesion_vs_HL=", sprintf("%.02f",pCheekDLes_HL),
                   "|ttest.Cheek.ADNonlesion_vs_HL=", sprintf("%.02f",pCheekDNL_HL),
                   "|ttest.Cheek.AD_vs_HL=", sprintf("%.02f",pCheekAD_HL),
                   "|ttest.Cheek.Lesion_vs_Nonlesion=", sprintf("%.02f",pCheekLes_NL),
                   "|ttest.AllAD.Lesion_vs_Nonlesion=", sprintf("%.02f",pADLes_ADNL),
                   "|ttest.All.ADLesion_vs_HL=", sprintf("%.02f",pADLes_HL),
                   "|ttest.All.ADNonlesion_vs_HL=", sprintf("%.02f",pADNL_HL),
                   "|ttest.Cheek_vs_Arm=", sprintf("%.02f",pCheek_Arm),
                   "|ttest.Lesion.Cheek_vs_Arm=", sprintf("%.02f",pLesCheek_Arm),
                   "|ttest.ADNonlesion.Cheek_vs_Arm=", sprintf("%.02f",pADNLCheek_Arm),
                   "|ttest.Healthy.Cheek_vs_Arm=", sprintf("%.02f",pHLCheek_Arm),
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
                   "|meanHealthy=", mHL
)

#print(word_observed)
print(paste("testname",
            "IGASpearmanCorrelation_rho",
            "IGAPearsonCorrelation_rho",
            "IGASpearmanCorrelation_pval",
            "IGAPearsonCorrelation_pval",
            "Arm-IGASpearmanCorrelation_rho","Arm-IGAPearsonCorrelation_rho", "Arm-IGASpearmanCorrelation_pval","Arm-IGAPearsonCorrelation_pval",
            "Cheek-IGASpearmanCorrelation_rho","Cheek-IGAPearsonCorrelation_rho", "Cheek-IGASpearmanCorrelation_pval","Cheek-IGAPearsonCorrelation_pval",
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
            collapse = "\t"
))
print(paste("observed -ttest",
            spearman_rho, pearson_rho,spearman_pval, pearson_pval, 
            Aspearman_rho, Apearson_rho, Aspearman_pval, Apearson_pval,
            Cspearman_rho, Cpearson_rho, Cspearman_pval, Cpearson_pval,
            KP_SiteDiseaseLesion,
            KP_Site,
            KP_Disease, 
            KP_Lesion,
            KP_SID,
            KP_SiteDisease,
            KP_SiteLesion,
            KP_DiseaseLesion,
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
            collapse = "\t"
)  )





###########################################
library(rstatix)
library(ggpubr)
library(dplyr)

# run t-test and keep only p <= 0.1
stat.observed <- mydata %>%
  t_test(observed ~ SiteDiseaseLesion) %>%
  filter(p <= 0.1) %>%   # <— filter here
  mutate(y.position = 1.05 * max(mydata$observed))

# plot
png(filename=paste0(filename, ".observed.png"), width=1600, height=1600, res=300)
pobserved <- ggboxplot(mydata, x = "SiteDiseaseLesion", y = "observed", color = "Site", add = "jitter", add.params = list(alpha = 0.6, width = 0.2))

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
  t_test(shannon ~ SiteDiseaseLesion) %>%
  filter(p <= 0.1) %>%   # <— filter here
  mutate(y.position = 1.05 * max(mydata$shannon))

# plot
png(filename=paste0(filename, ".shannn.png"), width=1600, height=1600, res=300)
pshannon <- ggboxplot(mydata, x = "SiteDiseaseLesion", y = "shannon", color = "Site", add = "jitter", add.params = list(alpha = 0.6, width = 0.2))

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
p2<-ggplot(mds_data, aes(x = MDS1, y = MDS2, color = SiteDiseaseLesion)) +geom_point(size=3, alpha = 0.7)+
  theme_bw()+stat_chull(aes(color=SiteDiseaseLesion, fill=SiteDiseaseLesion), alpha=0.01, geom="polygon")
print(p2)
dev.off() 
doc <- doc %>%
  add_slide(layout = "Title and Content", master = "Office Theme") %>%
  ph_with(value = "beta diversity by Group", location = ph_location_type(type = "title")) %>%
  # Add left plot (p1)
  ph_with(dml(ggobj = p2), location = ph_location(left = 0.5, top = 1.2, width = 8, height = 6))

png(filename=paste0(filename,".beta.GroupxSite.png"),  width=1800, height=800, res=300)
p2Group2Site<-ggplot(mds_data, aes(x = MDS1, y = MDS2, color = SiteDiseaseLesion)) +geom_point(size=3, alpha = 0.7)+
  theme_bw()+stat_chull(aes(color=SiteDiseaseLesion, fill=SiteDiseaseLesion), alpha=0.01, geom="polygon")+
  facet_wrap(~ Site) 
print(p2Group2Site)
dev.off() 

doc <- doc %>%
  add_slide(layout = "Title and Content", master = "Office Theme") %>%
  ph_with(value = "beta diversity by Site", location = ph_location_type(type = "title")) %>%
  # Add left plot (p1)
  ph_with(dml(ggobj = p2Group2Site), location = ph_location(left = 0.5, top = 1.2, width = 8, height = 6))



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
t1t="adonis result bc.dist ~ SiteDiseaseLesion"
t1=adonis2(bc.dist ~ mydata$SiteDiseaseLesion)
t2t="adonis result bc.dist ~ Site"
t2=adonis2(bc.dist ~ mydata$Site) #0.001
t3t="adonis result bc.dist ~ Disease"
t3=adonis2(bc.dist ~ mydata$Disease) #0.085
t4t="adonis result bc.dist ~ SID"
t4=adonis2(bc.dist ~ mydata$SID) #0.001
t5t="adonis result bc.dist ~ Lesion"
t5=adonis2(bc.dist ~ mydata$Lesion) #0.122
t6t="adonis result bc.dist ~ SiteLesion"
t6=adonis2(bc.dist ~ mydata$SiteLesion) #0.004
t7t="adonis result bc.dist ~ Sex"
t7=adonis2(bc.dist ~ mydata$Sex) #0.508
t8t="adonis result bc.dist ~ SiteDisease"
t8=adonis2(bc.dist ~ mydata$SiteDisease) #0.004
t9t="adonis result bc.dist ~ DiseaseLesion"
t9=adonis2(bc.dist ~ mydata$DiseaseLesion) #0.00

myAdonis <- paste(t1t, extract_adonis2(t1)$p,"|",
                  t2t, extract_adonis2(t2)$p,"|",
                  t3t, extract_adonis2(t3)$p,"|",
                  t4t, extract_adonis2(t4)$p,"|",
                  t5t, extract_adonis2(t5)$p,"|",
                  t6t, extract_adonis2(t6)$p,"|",
                  t7t, extract_adonis2(t7)$p,"|",
                  t8t, extract_adonis2(t8)$p,"|",
                  t9t, extract_adonis2(t9)$p
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
cbn <- combn(x=unique(mds_data$SiteDiseaseLesion), m = 2)
p <- c()
for (i in 1:ncol(cbn)) {
  ps.subs <- x_clean[mds_data$SiteDiseaseLesion %in% cbn[, i], ]
  metadata_sub <- mds_data[mds_data$SiteDiseaseLesion %in% cbn[, i], ]
  permanova_pairwise <- adonis2(
    vegdist(ps.subs, method = "bray") ~ metadata_sub$SiteDiseaseLesion
  )
  p <- c(p, permanova_pairwise$`Pr(>F)`[1])
}
p
p.adj <- p.adjust(p, method = "BH")
p.table <- cbind.data.frame(cbn[1,],cbn[2,] ,p=p, p.adj=p.adj)
p.table
addWorksheet(wb, "betaDiversity_Group");writeData(wb, "betaDiversity_Group", p.table)

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
    value = fpar(ftext("Beta Diversity with Bray-Curtis Distance & Pairwise Permanova Test by group",
                       fp_text(font.size = 36))),
    location = ph_location(left = 0.5, top = 0.3, width = 9, height = 1)
  ) %>%
  # Add the pairwise permanova table next to the plot
  ph_with(pairwise_ft, 
          location = ph_location(left = 6.2, top = 2.0, width = 3.5, height = 5))
###############################################################
cbn <- combn(x=unique(mds_data$Site), m = 2)
p <- c()
for (i in 1:ncol(cbn)) {
  ps.subs <- x_clean[mds_data$Site %in% cbn[, i], ]
  metadata_sub <- mds_data[mds_data$Site %in% cbn[, i], ]
  permanova_pairwise <- adonis2(
    vegdist(ps.subs, method = "bray") ~ metadata_sub$Site
  )
  p <- c(p, permanova_pairwise$`Pr(>F)`[1])
}
p
p.adj <- p.adjust(p, method = "BH")
p.table <- cbind.data.frame(cbn[1,],cbn[2,] ,p=p, p.adj=p.adj)
p.table
addWorksheet(wb, "betaDiversity_Site");writeData(wb, "betaDiversity_Site", p.table)

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
    value = fpar(ftext("Beta Diversity with Bray-Curtis Distance & Pairwise Permanova Test by Site",
                       fp_text(font.size = 36))),
    location = ph_location(left = 0.5, top = 0.3, width = 9, height = 1)
  ) %>%
  # Add the pairwise permanova table next to the plot
  ph_with(pairwise_ft, 
          location = ph_location(left = 6.2, top = 2.0, width = 3.5, height = 5))



############################################################
cbn <- combn(x=unique(mds_data$Disease), m = 2)
p <- c()
for (i in 1:ncol(cbn)) {
  ps.subs <- x_clean[mds_data$Disease %in% cbn[, i], ]
  metadata_sub <- mds_data[mds_data$Disease %in% cbn[, i], ]
  permanova_pairwise <- adonis2(
    vegdist(ps.subs, method = "bray") ~ metadata_sub$Disease
  )
  p <- c(p, permanova_pairwise$`Pr(>F)`[1])
}
p
p.adj <- p.adjust(p, method = "BH")
p.table <- cbind.data.frame(cbn[1,],cbn[2,] ,p=p, p.adj=p.adj)
p.table
addWorksheet(wb, "betaDiversity_Disease");writeData(wb, "betaDiversity_Disease", p.table)

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
    value = fpar(ftext("Beta Diversity with Bray-Curtis Distance & Pairwise Permanova Test by disease",
                       fp_text(font.size = 36))),
    location = ph_location(left = 0.5, top = 0.3, width = 9, height = 1)
  ) %>%
  # Add the pairwise permanova table next to the plot
  ph_with(pairwise_ft, 
          location = ph_location(left = 6.2, top = 2.0, width = 3.5, height = 5))



######################################################################
cbn <- combn(x=unique(mds_data$Lesion), m = 2)
p <- c()
for (i in 1:ncol(cbn)) {
  ps.subs <- x_clean[mds_data$Lesion %in% cbn[, i], ]
  metadata_sub <- mds_data[mds_data$Lesion %in% cbn[, i], ]
  permanova_pairwise <- adonis2(
    vegdist(ps.subs, method = "bray") ~ metadata_sub$Lesion
  )
  p <- c(p, permanova_pairwise$`Pr(>F)`[1])
}
p
p.adj <- p.adjust(p, method = "BH")
p.table <- cbind.data.frame(cbn[1,],cbn[2,] ,p=p, p.adj=p.adj)
p.table
addWorksheet(wb, "betaDiversity_Lesion");writeData(wb, "betaDiversity_Lesion", p.table)

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
    value = fpar(ftext("Beta Diversity with Bray-Curtis Distance & Pairwise Permanova Test by Lesion",
                       fp_text(font.size = 36))),
    location = ph_location(left = 0.5, top = 0.3, width = 9, height = 1)
  ) %>%
  # Add the pairwise permanova table next to the plot
  ph_with(pairwise_ft, 
          location = ph_location(left = 6.2, top = 2.0, width = 3.5, height = 5))



########################################################
cbn <- combn(x=unique(mds_data$SiteLesion), m = 2)
p <- c()
for (i in 1:ncol(cbn)) {
  ps.subs <- x_clean[mds_data$SiteLesion %in% cbn[, i], ]
  metadata_sub <- mds_data[mds_data$SiteLesion %in% cbn[, i], ]
  permanova_pairwise <- adonis2(
    vegdist(ps.subs, method = "bray") ~ metadata_sub$SiteLesion
  )
  p <- c(p, permanova_pairwise$`Pr(>F)`[1])
}
p
p.adj <- p.adjust(p, method = "BH")
p.table <- cbind.data.frame(cbn[1,],cbn[2,] ,p=p, p.adj=p.adj)
p.table
addWorksheet(wb, "betaDiversity_SiteLesion");writeData(wb, "betaDiversity_SiteLesion", p.table)

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
    value = fpar(ftext("Beta Diversity with Bray-Curtis Distance & Pairwise Permanova Test by SiteLesion",
                       fp_text(font.size = 36))),
    location = ph_location(left = 0.5, top = 0.3, width = 9, height = 1)
  ) %>%
  # Add the pairwise permanova table next to the plot
  ph_with(pairwise_ft, 
          location = ph_location(left = 6.2, top = 2.0, width = 3.5, height = 5))

############################################################################
cbn <- combn(x=unique(mds_data$Sex), m = 2)
p <- c()
for (i in 1:ncol(cbn)) {
  ps.subs <- x_clean[mds_data$Sex %in% cbn[, i], ]
  metadata_sub <- mds_data[mds_data$Sex %in% cbn[, i], ]
  permanova_pairwise <- adonis2(
    vegdist(ps.subs, method = "bray") ~ metadata_sub$Sex
  )
  p <- c(p, permanova_pairwise$`Pr(>F)`[1])
}
p
p.adj <- p.adjust(p, method = "BH")
p.table <- cbind.data.frame(cbn[1,],cbn[2,] ,p=p, p.adj=p.adj)
p.table
addWorksheet(wb, "betaDiversity_Sex");writeData(wb, "betaDiversity_Sex", p.table)

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
    value = fpar(ftext("Beta Diversity with Bray-Curtis Distance & Pairwise Permanova Test by Male1Female2",
                       fp_text(font.size = 36))),
    location = ph_location(left = 0.5, top = 0.3, width = 9, height = 1)
  ) %>%
  # Add the pairwise permanova table next to the plot
  ph_with(pairwise_ft, 
          location = ph_location(left = 6.2, top = 2.0, width = 3.5, height = 5))
########################################################################
cbn <- combn(x=unique(mds_data$SiteDisease), m = 2)
p <- c()
for (i in 1:ncol(cbn)) {
  ps.subs <- x_clean[mds_data$SiteDisease %in% cbn[, i], ]
  metadata_sub <- mds_data[mds_data$SiteDisease %in% cbn[, i], ]
  permanova_pairwise <- adonis2(
    vegdist(ps.subs, method = "bray") ~ metadata_sub$SiteDisease
  )
  p <- c(p, permanova_pairwise$`Pr(>F)`[1])
}
p
p.adj <- p.adjust(p, method = "BH")
p.table <- cbind.data.frame(cbn[1,],cbn[2,] ,p=p, p.adj=p.adj)
p.table
addWorksheet(wb, "betaDiversity_SiteDisease");writeData(wb, "betaDiversity_SiteDisease", p.table)

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
    value = fpar(ftext("Beta Diversity with Bray-Curtis Distance & Pairwise Permanova Test by SiteDisease",
                       fp_text(font.size = 36))),
    location = ph_location(left = 0.5, top = 0.3, width = 9, height = 1)
  ) %>%
  # Add the pairwise permanova table next to the plot
  ph_with(pairwise_ft, 
          location = ph_location(left = 6.2, top = 2.0, width = 3.5, height = 5))

############################################################################
cbn <- combn(x=unique(mds_data$DiseaseLesion), m = 2)
p <- c()
for (i in 1:ncol(cbn)) {
  ps.subs <- x_clean[mds_data$DiseaseLesion %in% cbn[, i], ]
  metadata_sub <- mds_data[mds_data$DiseaseLesion %in% cbn[, i], ]
  permanova_pairwise <- adonis2(
    vegdist(ps.subs, method = "bray") ~ metadata_sub$DiseaseLesion
  )
  p <- c(p, permanova_pairwise$`Pr(>F)`[1])
}
p
p.adj <- p.adjust(p, method = "BH")
p.table <- cbind.data.frame(cbn[1,],cbn[2,] ,p=p, p.adj=p.adj)
p.table
addWorksheet(wb, "betaDiversity_DiseaseLesion");writeData(wb, "betaDiversity_DiseaseLesion", p.table)

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
    value = fpar(ftext("Beta Diversity with Bray-Curtis Distance & Pairwise Permanova Test by DiseaseLesion",
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
