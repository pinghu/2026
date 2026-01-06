rm(list=ls())
library(ggplot2)
library(ggpubr)
library(rstatix)
library(dplyr)
library(tidyverse)
library("plotrix")
library(vegan)

filename1="~/Desktop/Disk2/project/GSS3149OralBiofilm/GSS3149OralBiofilm.meta.xls"
args <- commandArgs(trailingOnly = TRUE)
#filename<- args[1]
filename="species_counts.txt.10filter"
A<-read.table(filename1, sep="\t", header=TRUE)
d <- dim(A);
order_levels <- c("NoTreat","NCColgCP","Formula1","Formula2","Formula4", "Formula5")
A$Treat <- factor(A$Treat, levels = order_levels, ordered = TRUE)


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

test_d=data.frame(t(test))

min(apply(test_d, 1, sum)) ###162338

test_a <- decostand(test_d, method = "total")

Cname=rownames(test_d)
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

#https://stackoverflow.com/questions/30057765/histogram-ggplot-show-count-label-for-each-bin-for-each-category

shannon<-diversity(test_d)
simp<-diversity(test_d, "simpson")
invsimp<-diversity(test_d, "inv")

observed<-apply(test_d>0,1,sum)
N <- apply(test_d,1,sum)
###Richness Index#####
Menhinick_index<-observed/sqrt(N)
Margalef_index <-(observed-1)/log(N)

mydata0=data.frame(shannon,simp,invsimp, observed, Menhinick_index, Margalef_index ,SID)

mydata <- inner_join(mydata0, matching_dataset, by = c("SID" = "SID"))

write.table(t(mydata), file =paste0(filename, ".alphadiversity.RShort"), col.names=TRUE, row.names=TRUE, sep = "\t")

#saveRDS(mydata,paste0(filename, ".diversitydata"))
#saveRDS(test_d,paste0(filename,".test_d"))
#https://www.datanovia.com/en/blog/ggpubr-how-to-add-p-values-generated-elsewhere-to-a-ggplot/
png(filename=paste0(filename, ".shannon.Treat.png"), width=1600, height=1600, res=300)
stat.test <- mydata %>%
  t_test(shannon ~ Treat) %>%
  mutate(y.position = 1.05*max(mydata$shannon))
bxp <- ggboxplot(mydata, x = "Treat", y = "shannon", color = "Treat")
bxp<- bxp + stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, step.increase = 0.1)+
  ggtitle(paste0(filename, " Treat"))+
  geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.7) + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) #+ scale_x_discrete(limits = order_levels)
print(bxp)
dev.off()



png(filename=paste0(filename, ".shannon.TreatSig.png"), width=1600, height=2000, res=300)
# Perform t-test and add y.position for annotation
stat.test <- mydata %>%
  t_test(shannon ~ Treat) %>%
  mutate(y.position = 1.05 * max(mydata$shannon)) %>%
  filter(p <= 0.05)  # Keep only significant p-values
# Create boxplot
bxp <- ggboxplot(mydata, x = "Treat", y = "shannon", color = "Treat")
# Add significant p-values only
if (nrow(stat.test) > 0) {  # Ensure there are significant values before adding
  bxp <- bxp + stat_pvalue_manual(stat.test, label = "p", tip.length = 0.01, step.increase = 0.1)
}
bxp <- bxp +
  ggtitle(paste0(filename, " Treat")) +
  geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.7) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 60, hjust = 1)) 
# Print and save
print(bxp)
dev.off()


########################observed#################################
png(filename=paste0(filename, ".observed.Treat.png"), width=1600, height=2000, res=300)
# stat.test <- mydata %>%
#   t_test(observed ~ Treat) %>%
#   mutate(y.position = 1.05*max(mydata$observed))

library(dplyr)
library(rstatix)

# Remove NAs
mydata <- mydata %>% filter(!is.na(observed), !is.na(Treat))

# Check if any group has zero variance
sds <- mydata %>%
  group_by(Treat) %>%
  summarise(sd_val = sd(observed, na.rm = TRUE), .groups = "drop") %>%
  pull(sd_val)

if (any(sds == 0) || n_distinct(mydata$Treat) != 2) {
  # Constant group -> assign p = 1
  stat.test <- tibble(
    .y. = "observed",
    group1 = unique(mydata$Treat)[1],
    group2 = unique(mydata$Treat)[2],
    p = 1,
    p.signif = "ns",
    y.position = 1.05 * max(mydata$observed, na.rm = TRUE)
  )
} else {
  # Normal t-test
  stat.test <- t_test(mydata, observed ~ Treat) %>%
    add_significance("p") %>%
    mutate(y.position = 1.05 * max(mydata$observed, na.rm = TRUE))
}




bxp <- ggboxplot(mydata, x = "Treat", y = "observed", color = "Treat")
bxp<- bxp + stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, step.increase = 0.1)+
  ggtitle(paste0(filename, " Treat"))+
  geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.7) + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) #+ scale_x_discrete(limits = order_levels)
print(bxp)
dev.off()









png(filename=paste0(filename, ".observed.TreatSig.png"), width=1600, height=2000, res=300)
# Perform t-test and add y.position for annotation
# stat.test <- mydata %>%
#   t_test(observed ~ Treat) %>%
#   mutate(y.position = 1.05 * max(mydata$observed)) %>%
#   filter(p <= 0.05)  # Keep only significant p-values

library(dplyr)
library(rstatix)

safe_t_test <- function(df) {
  # Try Welch t-test; if it errors (e.g., constant data), return a dummy row with p=1
  tryCatch(
    t_test(df, observed ~ Treat),
    error = function(e) {
      grps <- levels(factor(df$Treat))
      tibble(
        .y.     = "observed",
        group1  = if (length(grps) >= 1) grps[1] else NA_character_,
        group2  = if (length(grps) >= 2) grps[2] else NA_character_,
        p       = 1,
        p.signif = "ns"
      )
    }
  )
}

stat.test <- mydata %>%
  safe_t_test() %>%
  mutate(y.position = 1.05 * max(mydata$observed, na.rm = TRUE)) %>%
  filter(p <= 0.05)  # Keep only significant p-values




# Create boxplot
bxp <- ggboxplot(mydata, x = "Treat", y = "observed", color = "Treat")
# Add significant p-values only
if (nrow(stat.test) > 0) {  # Ensure there are significant values before adding
  bxp <- bxp + stat_pvalue_manual(stat.test, label = "p", tip.length = 0.01, step.increase = 0.1)
}
bxp <- bxp +
  ggtitle(paste0(filename, " Treat")) +
  geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.7) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 60, hjust = 1)) 
# Print and save
print(bxp)
dev.off()



######################################################33333
KP_Treat_shannon=kruskal.test(shannon ~ Treat, data = mydata)$p.value 
KP_Treat_observed=kruskal.test(observed ~ Treat, data = mydata)$p.value 
print(paste0("Shannon KP_Treat=", KP_Treat_shannon,  
             "; KP_Treat_observed=", KP_Treat_observed))
#"Shannon KP_Treat=0.0136995515171783; KP_Treat_observed=0.270642020292283"
########################################
x<-test_d
beta_dist <- vegdist(x, index = "bray")
mds <- metaMDS(beta_dist)
mds_data <- cbind(as.data.frame(mds$points), mydata)
mds_data$SampleID <- rownames(mds_data)

#install.packages('devtools')
library(devtools)
#install_github('fawda123/ggord')
library(ggord)
library(ggplot2)

 png(filename=paste0(filename,".beta.Treat.png"),  width=1800, height=800, res=300)
 p1<-ggplot(mds_data, aes(x = MDS1, y = MDS2, color = Treat)) +geom_point(size=3, alpha = 0.7)+
   theme_bw()+stat_chull(aes(color=Treat, fill=Treat), alpha=0.1, geom="polygon")
 print(p1)
 dev.off() 
 
############################################################
print("all samples")
# calculate Bray-Curtis distance among samples
bc.dist <- vegdist(test_d, method = "bray")
# cluster communities using average-linkage algorithm
bc.clust <- hclust(bc.dist, method = "average")
# plot cluster diagram
png(filename=paste0(filename,".bray-Curtis-cluster.png"), width=3000, height=1600)
plot(bc.clust, ylab = "Bray-Curtis dissimilarity")
dev.off()
# Taxonomic (Bray-Curtis) dissimilarity explained
print("adonis result bc.dist ~ Treat")
adonis2(bc.dist ~ mydata$Treat) # p<=0.001

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
cbn <- combn(x=unique(mds_data$Treat), m = 2)

library(vegan)

p <- c()
for (i in 1:ncol(cbn)) {
  ps.subs <- x[mds_data$Treat %in% cbn[, i], ]
  metadata_sub <- mds_data[mds_data$Treat %in% cbn[, i], ]
  
  # Use adonis2 instead of adonis
  permanova_pairwise <- adonis2(
    vegdist(ps.subs, method = "bray") ~ metadata_sub$Treat
  )
  
  # Extract p-value from adonis2 result
  p <- c(p, permanova_pairwise$`Pr(>F)`[1])
}

p






# p <- c()
# for(i in 1:ncol(cbn)){
#   ps.subs <- x [c(mds_data$Treat %in% cbn[,i]), ]
#   metadata_sub <- mds_data[c(mds_data$Treat %in% cbn[,i]), ]
#   permanova_pairwise <- adonis(vegdist(ps.subs, method = "bray") ~ metadata_sub$Treat)
#   p <- c(p, permanova_pairwise$aov.tab$`Pr(>F)`[1])
# }

p.adj <- p.adjust(p, method = "BH")
p.table <- cbind.data.frame(cbn[1,],cbn[2,] ,p=p, p.adj=p.adj)
p.table
# Output the results to a file
write.table(p.table, paste0(filename, ".beta_Treat.txt"), sep = "\t", quote = FALSE)

