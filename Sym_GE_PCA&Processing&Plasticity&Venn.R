#####Symbiont DTV RNAseq analysis based on Colleen Boves code (https://github.com/seabove7/BelizeRT_Castillo_Bove/blob/main/RT_Host_Sym_GE.R)

###sym data
setwd("~/path/Sym_GE")

library("DESeq2")
library("ggplot2")
library("ggh4x")
library("dplyr")
library("vegan")
library("tidyverse")
library("pheatmap")
library("VennDiagram")
library("adegenet") 
library("WGCNA")
library("effectsize")
library("ggvenn")

 # if(!requireNamespace("BiocManager", quietly = TRUE))
 #   install.packages("BiocManager")
 # BiocManager::install("preprocessCore")

##### First going to get the count of genes in all sym samples ####
#read in counts
countData <- read.table("~/Desktop/BU NSF Postdoc/Nicola_Ofav_chpt3/Realign to genomes/FinalCounts_13Nov2024/Dtrenchii_counts_fixedSNP.txt")
head(countData)
length(countData[,1]) #55804 genes

names(countData)=sub(".fastq.Sym.counts.txt","",names(countData))
names(countData)
head(countData)

totalCounts=colSums(countData)
barplot(totalCounts, col="coral")
totalCounts

#fixedSNP with clones
# KA1.1   KA3.1   KA4.2   KA5.2  KB13.2  KB28.1  KB39.1  KB39.2  KE13.1  KE14.1 
# 437430  546207  573632  193257  801905 1044171  233697  432257  358819  146307 
# KE31.2   KE9.2  KF19.1  KF19.2  KF23.2   KF9.1  KG12.1  KG12.2  KG31.1  KG31.2 
# 105942  570093  126680  535874   84708  157061  689974  506562  191830  407233 
# VB11.1  VB11.2   VB1.2   VB5.1  VD10.1  VD15.2  VD17.1  VD17.2  VE14.2  VE18.1 
# 825846  358819  644044  929338  433837  565770  575633  312162 1065778  462595 
# VE18.2  VE30.1  VF23.2  VF31.1  VF40.1  VF40.2   VG1.1   VG1.2  VG17.2   VG6.1 
# 338973 1169007  912359 1039138  266919  200483  261628  235589 1134249  456301 

min(totalCounts) #84,708
max(totalCounts)  #1,169,007
mean(totalCounts) #508,302.7

#Remove genes with low counts from dataset - count less than 2 in more than 90% of samples
countData$low = apply(countData[,1:40],1,function(x){sum(x<=2)})#making new column counting number of samples with counts <=2 within each gene (host) 

countData<-countData[-which(countData$low>36),] #36 is ~90% of 40 samples - get rid of count less than 2 in more than 90% of samples
countData <- countData[,1:40]

#check after removing genes
totalCounts=colSums(countData)
barplot(totalCounts, col="coral")
totalCounts
min(totalCounts) #84,599
max(totalCounts)  #1,166,369
mean(totalCounts) #507,196.1

#Count the number of genes in each sample
nrow(countData)#2,351 for clones included dataset
SampleGenes = colSums(countData != 0)
NumSampleGenes <- data.frame(Column = names(countData), NonZeroCount = SampleGenes)


##### Now do real analysis with clones removed ####
#read in counts
countData <- read.table("~/path/Dtrenchii_counts_fixedSNP_noClones.txt")
head(countData)
length(countData[,1]) #55799 genes

names(countData)=sub(".fastq.Sym.counts.txt","",names(countData))
names(countData)
head(countData)

totalCounts=colSums(countData)
barplot(totalCounts, col="coral")
totalCounts

#fixedSNP_noClones
# KB13.2 KB28.1 KB39.1 KB39.2 KE13.1 KE14.1 KE31.2  KE9.2 KF19.1 KF19.2 KF23.2 
# 51790  69189  16303  29483  25420   9304   5831  39061   6758  38256   4504 
# KF9.1 KG12.1 KG12.2 KG31.1 KG31.2 VD10.1 VD15.2 VD17.1 VD17.2 VE14.2 VE18.1 
# 10955  45973  36742  11184  28114  28196  37773  40353  20128  77997  27801 
# VE18.2 VE30.1 VF23.2 VF31.1 VF40.1 VF40.2  VG1.1  VG1.2 VG17.2  VG6.1 
# 21562  82138  63146  73822  19040  13994  16876  16558  73763  32659 
  
min(totalCounts) #4504
max(totalCounts)  #82,138
mean(totalCounts) #33,583.53

#Be careful here, sample names on rows are stripped from the metadata
#in these steps, so order of samples in count table MUST match the
#order of samples in metadata, otherwise the metadata will get jumbled
#rt <- read.csv("~/Desktop/BU NSF Postdoc/Nicola_Ofav_chpt3/Thermvar_metadata_Apr2025_ClonesFixed.csv") #this file makes the right figures
rt <- read.csv("~/path/THERMVAR_MAIN_METADATA_Apr25_clonesfixed.csv") #this one doesn't


rt$Timepoint <- gsub("Preheat", "Pre-heat", rt$Time)
rt$Timepoint <- gsub("Postheat", "Post-heat", rt$Time)

rt$TreatTime <- gsub("Control_Preheat", "Control_Pre-heat", rt$TreatTime)
rt$TreatTime <- gsub("Control_Postheat", "Control_Post-heat", rt$TreatTime)
rt$TreatTime <- gsub("DTV_Preheat", "DTV_Pre-heat", rt$TreatTime)
rt$TreatTime <- gsub("DTV_Postheat", "DTV_Post-heat", rt$TreatTime)

#Remove rows from metadata that do not have GE
rt <- rt %>%
  filter(GE_sample != "")

#Remove rows from metadata that are clones
rt <- subset(rt, !(Genotype_Original %in% c("KA", "VB")))

head(rt)
conditions=data.frame(rt$Species, rt$Treatment, rt$Time, rt$Genotype, rt$TreatTime)
nrow(conditions) #32
names(conditions)=c("Species", "Treatment", "Time", "Genotype", "TreatTime")

#Remove genes with low counts from dataset - count less than 2 in more than 90% of samples
countData$low = apply(countData[,1:32],1,function(x){sum(x<=2)})#making new column counting number of samples with counts <=2 within each gene (host) 

countData<-countData[-which(countData$low>28),] #28 is ~90% of 32 samples - get rid of count less than 2 in more than 90% of samples
countData <- countData[,1:32]

#check after removing genes
totalCounts=colSums(countData)
barplot(totalCounts, col="coral")
totalCounts
min(totalCounts) #4385
max(totalCounts)  #79,222
mean(totalCounts) #32,395.12

#check that same number samples across counts and metadata
#also count the nuber of genes in each sample
ncol(countData) #32
nrow(countData)#2188 - maybe not good enough for WGCNA
nrow(conditions) #32

#first ran with design of experiment (next line) for initial check, then ran without design
#dds<-DESeqDataSetFromMatrix(countData=countData, colData=conditions, design=~species+treatment*time) #can only test for the main effects species, treatment, and time

#to make venn diagrams used this one to compare among groups
#dds<-DESeqDataSetFromMatrix(countData=countData, colData=conditions, design=~species+treatTime) 
##ran above to look at effect of time and treatment as one variable, wasn't significant so used first design for analysis
# I did use the treatTime design to generate PCA for this comparison, no significance was found with adonis2.

#deseq without giving design, seems to give slightly stronger trends so we stuck with this design for PCA and WGCNA analyses.
dds<-DESeqDataSetFromMatrix(countData=countData, colData=conditions, design=~1) #can only test for the main effects species, treatment, and time


#one step DESeq#one step DESeqtime()
dds<-DESeq(dds)#together
# estimating size factors
# estimating dispersions
# gene-wise dispersion estimates
# mean-dispersion relationship
# final dispersion estimates
# fitting model and testing

res<- results(dds)
############################# use rlog to transform 
rlog=rlogTransformation(dds, blind=TRUE) 
#save(rlog, file="rlogCountsGreedy.RData") #save for WGCNA - didn't do WGCNA
rld=assay(rlog)
rld_t=t(rld)
colnames(rld_t)
head(rld_t)
pca <- prcomp(rld_t,center = TRUE)
head(pca)
bstick(pca) #evaluates pca to see what to keep
screeplot(pca, bstick = TRUE, type = "lines") #makes figure showing results
li <- pca$sdev^2 / sum(pca$sdev^2)
pc1v <- round(li[1] * 100, 1)
pc2v <- round(li[2] * 100, 1)
pca_s <- as.data.frame(pca$x)
head(pca_s)
pca_s <- pca_s[,c(1,2)]
pca_s$Samples = row.names(pca_s)
pca_s$Treatment=conditions$Treatment
pca_s$Species=conditions$Species
pca_s$Time=conditions$Time
pca_s$Genotype=conditions$Genotype
pca_s$TreatTime=conditions$TreatTime
head(pca_s)

#reorder time variable for figures
f=c('Pre-heat','Post-heat')
pca_s <- within(pca_s, Time <- factor(Time, levels=f))

#reorder time variable for figures
f=c('Control_Pre-heat','Control_Post-heat', 'DTV_Pre-heat', 'DTV_Post-heat')
pca_s <- within(pca_s, TreatTime <- factor(TreatTime, levels=f))

#all three variables on same PCA - this shows separation by treatTime (interesting)
#cbPalette <- c( "#2C2673","#744899", "#EC4176", "#FFA45E")#treatTime colors #"#EC4176", "#EE93B0")
cbPalette <- c("#2C2673", "#9187FF","#E24A13","#FFA45E") 
#cbPalette <- c( "#D89362", "#8F3C00")#time
#pdf("PCA_Host_allgenes_rlog.pdf",height=5,width=6)
ggplot(pca_s, aes(PC1, PC2, color = TreatTime, shape=Species, group=TreatTime)) +
  geom_point(size = 4) +
  theme_bw() +
  guides(color = guide_legend(title = "Treatment x Heat"), shape = guide_legend(title = "Host Species"))+
  scale_colour_manual(values=cbPalette, labels = c("Control - Preheat", "Control - Postheat", 
                                                   "DTV - Preheat", "DTV - Postheat"))+
  theme( legend.title=element_text(size=14, face = "bold", color = "black"), 
         legend.text=element_text(size=14, face = "bold", color="black"),
         axis.title.x = element_text(size = 16, color = "black", face = "bold"),
         axis.text = element_text(color = "black", size = 14, face = "bold"),
         axis.title.y = element_text(size = 16, color = "black", face = "bold"),
         axis.ticks = element_line(color = "black"))+
  annotate("text", x = min(-25), y = max(20), 
           label = "Symbiodiniaceae", hjust = 0, vjust = 1, 
           size = 6, fontface = "bold", color = "black")+
  stat_ellipse()+
  xlab(paste0("PC1: ",pc1v,"% variance")) +
  ylab(paste0("PC2: ",pc2v,"% variance"))
#dev.off()


#faceted by species -- still see trend of separation by treatTime
# New facet label names for Species variable
sp.labs <- c("O. faveolata", "O. franksi")
names(sp.labs) <- c("Ofav", "Ofra")
#all three variables on same PCA
#cbPalette <- c( "#B395DA","#4A2E61", "#D89362", "#8F3C00")#Treat v Time
cbPalette <- c("#2C2673", "#9187FF","#E24A13","#FFA45E") #treatTime
#pdf("PCA_Host_allgenes_rlog.pdf",height=5,width=6)
ggplot(pca_s, aes(PC1, PC2, color = TreatTime, shape=Genotype, group=TreatTime)) +
  facet_grid2(Species~., labeller = labeller(species = sp.labs), 
              strip = strip_themed(background_y = elem_list_rect(fill = "black"), 
                             text_y = elem_list_text(size = 14, face = "bold", color = "white")))+
  geom_point(size = 4) +
  scale_shape_manual(values=c(16,15,18,17, 21,22,23,24))+
  theme_bw() +
  guides(color = guide_legend(title = "Treatment x Heat"), shape = guide_legend(title = "Host Genotype"))+
  scale_colour_manual(values=cbPalette, labels = c("Control - Preheat", "Control - Postheat", 
                                                   "DTV - Preheat", "DTV - Postheat"))+
  theme( legend.title=element_text(size=14, face = "bold", color = "black"), 
         legend.text=element_text(size=14, face = "bold", color="black"),
         axis.title.x = element_text(size = 16, color = "black", face = "bold"),
         axis.text = element_text(color = "black", size = 14, face = "bold"),
         axis.title.y = element_text(size = 16, color = "black", face = "bold"),
         axis.ticks = element_line(color = "black"))+
  stat_ellipse()+
  xlab(paste0("PC1: ",pc1v,"% variance")) +
  ylab(paste0("PC2: ",pc2v,"% variance")) 
#dev.off()

#faceted by time, all species together -- no big trends
# New facet label names for time variable
time.labs <- c("Preheat", "Postheat")
names(time.labs) <- c("preheat", "postheat")
cbPalette <- c("#2C2673", "#E24A13") #treatment
#pdf("PCA_Host_allgenes_rlog.pdf",height=5,width=6)
ggplot(pca_s, aes(PC1, PC2, color = Treatment, shape=Species, group=Treatment)) +
  facet_grid2(Time~., labeller = labeller(time = time.labs), 
              strip = strip_themed(background_y = elem_list_rect(fill = "black"), 
                                  text_y = elem_list_text(size = 14, face = "bold", color = "white")))+
  geom_point(size = 4) +
  theme_bw() +
  guides(color = guide_legend(title = "Treatment"), shape = guide_legend(title = "Host Species"))+
  scale_colour_manual(values=cbPalette, labels = c("Control", "DTV"))+
  theme( legend.title=element_text(size=14, face = "bold", color = "black"), 
         legend.text=element_text(size=14, face = "bold", color="black"),
         axis.title.x = element_text(size = 16, color = "black", face = "bold"),
         axis.text = element_text(color = "black", size = 14, face = "bold"),
         axis.title.y = element_text(size = 16, color = "black", face = "bold"),
         axis.ticks = element_line(color = "black"))+
  stat_ellipse()+
  xlab(paste0("PC1: ",pc1v,"% variance")) +
  ylab(paste0("PC2: ",pc2v,"% variance")) 
#dev.off()

#faceted by species - divide within species by pre and post heating
# New facet label names for treatment variable
treat.labs <- c("Control", "DTV")
names(treat.labs) <- c("control", "therm_var")
cbPalette <- c( "#D89362", "#8F3C00")#time
#pdf("PCA_Host_allgenes_rlog.pdf",height=5,width=6)
ggplot(pca_s, aes(PC1, PC2, color = Time, shape=Treatment, group=Time)) +
  facet_grid2(Species~., labeller = labeller(species = sp.labs), 
              strip = strip_themed(background_y = elem_list_rect(fill = "black"), 
                                   text_y = elem_list_text(size = 14, face = "bold", color = "white")))+
  geom_point(size = 4) +
  theme_bw() +
  guides(color = guide_legend(title = "Heat"), shape = guide_legend(title = "Treatment"))+
  scale_colour_manual(values=cbPalette, labels = c("Preheat", "Postheat"))+
  #scale_shape(labels = c("Control", "DTV"))+
  scale_shape_manual(values=c(15,18), labels = c("Control", "DTV"))+
  theme( legend.title=element_text(size=14, face = "bold", color = "black"), 
         legend.text=element_text(size=14, face = "bold", color="black"),
         axis.title.x = element_text(size = 16, color = "black", face = "bold"),
         axis.text = element_text(color = "black", size = 14, face = "bold"),
         axis.title.y = element_text(size = 16, color = "black", face = "bold"),
         axis.ticks = element_line(color = "black"))+
  stat_ellipse()+
  xlab(paste0("PC1: ",pc1v,"% variance")) +
  ylab(paste0("PC2: ",pc2v,"% variance")) 
#dev.off()

#faceted by time and species
cbPalette <- c( "#2C2673", "#9187FF","#E24A13","#FFA45E")#treatment 
# New facet label names for time variable
time.labs <- c("Preheat", "Postheat")
names(time.labs) <- c("preheat", "postheat")
#pdf("PCA_Host_allgenes_rlog.pdf",height=5,width=6)
ggplot(pca_s, aes(PC1, PC2, color = TreatTime, shape= Genotype, group=TreatTime)) +
  facet_grid2(Time~Species, labeller = labeller(species = sp.labs, time =time.labs), 
              strip = strip_themed(background_x = elem_list_rect(fill = "black"),
                                   background_y = elem_list_rect(fill = c("preheat"="grey20", "postheat"="grey70")), 
                                   text_y = elem_list_text(size = 14, face = "bold", color = "white"),
                                   text_x = elem_list_text(size = 14, face = "bold", color = "white")))+
  geom_point(size = 4,stroke = 1) +
  scale_shape_manual(values=c(16,15,18,17, 21,22,23,24))+
  theme_bw() +
  guides(color = guide_legend(title = "Treatment x Heat"), shape = guide_legend(title = "Host Species"))+
  scale_colour_manual(values=cbPalette, labels = c("Control - Preheat", "Control - Postheat", 
                                                   "DTV - Preheat", "DTV - Postheat"))+
  theme(plot.title = element_text(hjust = 0.5),
        title = element_text(size=14, face = "bold", color = "black"),
        legend.title=element_text(size=14, face = "bold", color = "black"), 
        legend.text=element_text(size=14, face = "bold", color="black"),
        axis.title.x = element_text(size = 16, color = "black", face = "bold"),
        axis.text = element_text(color = "black", size = 14, face = "bold"),
        axis.title.y = element_text(size = 16, color = "black", face = "bold"),
        axis.ticks = element_line(color = "black"))+
  stat_ellipse()+
  ggtitle("Symbiodiniaceae GE")+
  xlab(paste0("PC1: ",pc1v,"% variance")) +
  ylab(paste0("PC2: ",pc2v,"% variance")) 
#dev.off()

#faceted by treatment and species -- shows more variability in controls pre/postheat and tighter exp in therm_var
cbPalette <- c( "#2C2673", "#9187FF","#E24A13","#FFA45E")#treatTime
#pdf("PCA_Host_allgenes_rlog.pdf",height=5,width=6)
ggplot(pca_s, aes(PC1, PC2, color = TreatTime, shape= Genotype, group=TreatTime)) +
  facet_grid2(Treatment~Species,labeller = labeller(Species = sp.labs, Treatment =treat.labs), 
              strip = strip_themed(background_x = elem_list_rect(fill = "black"),
                                                                    background_y = elem_list_rect(fill = c("Control"="#2C2673", "DTV"="#E24A13")), 
                                                                    text_y = elem_list_text(size = 14, face = "bold", color = "white"),
                                                                    text_x = elem_list_text(size = 14, face = "bold", color = "white")))+
  geom_point(size = 4, stroke = 1) +
  scale_shape_manual(values=c(16,15,18,17, 21,22,23,24))+
  theme_bw() +
  guides(color = guide_legend(title = "Treatment x Heat"), shape = guide_legend(title = "Host Species"))+
  scale_colour_manual(values=cbPalette, labels = c("Control - Preheat", "Control - Postheat", 
                                                   "DTV - Preheat", "DTV - Postheat"))+
  theme(plot.title = element_text(hjust = 0.5),
        title = element_text(size=14, face = "bold", color = "black"),
        legend.title=element_text(size=14, face = "bold", color = "black"), 
        legend.text=element_text(size=14, face = "bold", color="black"),
        axis.title.x = element_text(size = 16, color = "black", face = "bold"),
        axis.text = element_text(color = "black", size = 14, face = "bold"),
        axis.title.y = element_text(size = 16, color = "black", face = "bold"),
        axis.ticks = element_line(color = "black"))+
  stat_ellipse()+
  ggtitle("Symbiodiniaceae GE")+
  xlab(paste0("PC1: ",pc1v,"% variance")) +
  ylab(paste0("PC2: ",pc2v,"% variance")) 
#dev.off()

#### Stats PCA ####
adonis2(pca_s[,1:2] ~ Species*TreatTime, data = pca_s, by = "terms", method='eu', na.rm = TRUE)
# DESeq = ~1
# adonis2(formula = pca_s[, 1:2] ~ species * treatTime, data = pca_s, method = "eu", by = "terms", na.rm = TRUE)
# Df SumOfSqs      R2      F Pr(>F)  
# species            1     9.38 0.00359 0.1220  0.889  
# treatTime          3   724.08 0.27710 3.1384  0.012 *
# species:treatTime  3    33.85 0.01295 0.1467  0.990  
# Residual          24  1845.75 0.70636                
# Total             31  2613.05 1.00000 


#overall PERMANOVA - no difference by host species - so omitted from formulas
adonis2(pca_s[,1:2] ~ Species+Treatment*Time, data = pca_s, by = "terms", method='eu', na.rm = TRUE)
#fixedSNP_noClones
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = pca_s[, 1:2] ~ species + treatment * time, data = pca_s, method = "eu", by = "terms", na.rm = TRUE)
#                 Df SumOfSqs      R2      F Pr(>F)    
# species         1     9.38 0.00359 0.1348  0.856    
# treatment       1    19.63 0.00751 0.2820  0.759    
# time            1   649.86 0.24870 9.3352  0.001 ***
# treatment:time  1    54.58 0.02089 0.7840  0.465    
# Residual       27  1879.59 0.71931                  
# Total          31  2613.05 1.00000 


library(pairwiseAdonis) 
# install.packages('devtools') #To install package
# library(devtools)
# install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

#both coral species pairwise PERMANOVA by treatTime
pairwise.adonis2(pca_s[,1:2] ~ TreatTime, data = pca_s, by = "terms", method='eu', na.rm = TRUE)
# $parent_call
# [1] "pca_s[, 1:2] ~ treatTime , strata = Null , permutations 999"
# 
# $therm_var_preheat_vs_control_preheat
# Df SumOfSqs     R2      F Pr(>F)
# treatTime  1    46.88 0.0465 0.6828  0.465
# Residual  14   961.27 0.9535              
# Total     15  1008.15 1.0000              
# 
# $therm_var_preheat_vs_control_postheat
# Df SumOfSqs      R2      F Pr(>F)   
# treatTime  1   447.23 0.31573 6.4597  0.005 **
#   Residual  14   969.26 0.68427                 
# Total     15  1416.49 1.00000                 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# $therm_var_preheat_vs_therm_var_postheat
# Df SumOfSqs      R2      F Pr(>F)   
# treatTime  1   391.80 0.38406 8.7297  0.002 **
#   Residual  14   628.35 0.61594                 
# Total     15  1020.15 1.00000                 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# $control_preheat_vs_control_postheat
# Df SumOfSqs      R2     F Pr(>F)  
# treatTime  1   312.64 0.19872 3.472  0.043 *
#   Residual  14  1260.63 0.80128               
# Total     15  1573.27 1.00000               
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# $control_preheat_vs_therm_var_postheat
# Df SumOfSqs      R2      F Pr(>F)  
# treatTime  1   222.27 0.19463 3.3834  0.056 .
# Residual  14   919.72 0.80537                
# Total     15  1141.99 1.00000                
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# $control_postheat_vs_therm_var_postheat
# Df SumOfSqs      R2      F Pr(>F)
# treatTime  1    27.33 0.02861 0.4124  0.609
# Residual  14   927.71 0.97139              
# Total     15   955.04 1.00000              
# 
# attr(,"class")
# [1] "pwadstrata" "list"      
# 


#to increase power for within species comparisons tried subsetting PCA data
Ofav<-subset(pca_s, Species=="O. faveolata")
adonis2(Ofav[,1:2] ~ TreatTime, data = Ofav,by = "terms", method='eu', na.rm = TRUE)
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = Ofav[, 1:2] ~ treatTime, data = Ofav, method = "eu", by = "terms", na.rm = TRUE)
#           Df SumOfSqs      R2      F Pr(>F)
# treatTime  3   479.51 0.31989 1.8814  0.127
# Residual  12  1019.47 0.68011              
# Total     15  1498.98 1.00000  

adonis2(Ofav[,1:2] ~ Treatment*Time, data = Ofav,by = "terms", method='eu', na.rm = TRUE)
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = Ofav[, 1:2] ~ time * treatment, data = Ofav, method = "eu", by = "terms", na.rm = TRUE)
#                 Df SumOfSqs      R2      F Pr(>F)  
# treatment       1    10.43 0.00696 0.1227  0.896  
# time            1   438.87 0.29278 5.1659  0.015 *
# treatment:time  1    30.21 0.02015 0.3556  0.719  
# Residual       12  1019.47 0.68011                
# Total          15  1498.98 1.00000   

pairwise.adonis2(Ofav[,1:2] ~ TreatTime, data = Ofav, by = "terms", method='eu', na.rm = TRUE)
# $parent_call
# [1] "Ofav[, 1:2] ~ treatTime , strata = Null , permutations 999"
# 
# $control_preheat_vs_therm_var_preheat
# Df SumOfSqs      R2      F Pr(>F)
# treatTime  1    29.38 0.06784 0.4367  0.565
# Residual   6   403.66 0.93216              
# Total      7   433.04 1.00000              
# 
# $control_preheat_vs_control_postheat
# Df SumOfSqs      R2      F Pr(>F)
# treatTime  1   232.06 0.27462 2.2715  0.155
# Residual   6   612.97 0.72538              
# Total      7   845.03 1.00000              
# 
# $control_preheat_vs_therm_var_postheat
# Df SumOfSqs      R2     F Pr(>F)
# treatTime  1   165.76 0.29503 2.511  0.153
# Residual   6   396.07 0.70497             
# Total      7   561.82 1.00000             
# 
# $therm_var_preheat_vs_control_postheat
# Df SumOfSqs      R2     F Pr(>F)
# treatTime  1   283.54 0.31264 2.729  0.124
# Residual   6   623.40 0.68736             
# Total      7   906.95 1.00000             
# 
# $therm_var_preheat_vs_therm_var_postheat
# Df SumOfSqs      R2      F Pr(>F)  
# treatTime  1   237.02 0.36832 3.4985  0.053 .
# Residual   6   406.49 0.63168                
# Total      7   643.52 1.00000                
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# $control_postheat_vs_therm_var_postheat
# Df SumOfSqs      R2      F Pr(>F)
# treatTime  1    11.26 0.01795 0.1097  0.818
# Residual   6   615.81 0.98205              
# Total      7   627.06 1.00000              
# 
# attr(,"class")
# [1] "pwadstrata" "list"      


Ofra<-subset(pca_s, Species=="O. franksi")
adonis2(Ofra[,1:2] ~ TreatTime, data = Ofra,by = "terms", method='eu', na.rm = TRUE)
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = Ofra[, 1:2] ~ treatTime, data = Ofra, method = "eu", by = "terms", na.rm = TRUE)
# Df SumOfSqs      R2      F Pr(>F)
# treatTime  3   278.41 0.25203 1.3478  0.307
# Residual  12   826.28 0.74797              
# Total     15  1104.69 1.00000 

adonis2(Ofra[,1:2] ~ Time*Treatment, data = Ofra,by = "terms", method='eu', na.rm = TRUE)
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = Ofra[, 1:2] ~ time * treatment, data = Ofra, method = "eu", by = "terms", na.rm = TRUE)
# 
#                Df SumOfSqs      R2      F Pr(>F)  
# time            1   242.92 0.21990 3.5280  0.045 *
# treatment       1    10.83 0.00980 0.1572  0.859  
# time:treatment  1    24.66 0.02233 0.3582  0.699  
# Residual       12   826.28 0.74797                
# Total          15  1104.69 1.00000  

pairwise.adonis2(Ofra[,1:2] ~ TreatTime, data = Ofra, by = "terms", method='eu', na.rm = TRUE)
# $parent_call
# [1] "Ofra[, 1:2] ~ treatTime , strata = Null , permutations 999"
# 
# $therm_var_preheat_vs_control_preheat
# Df SumOfSqs     R2      F Pr(>F)
# treatTime  1    19.08 0.0351 0.2183  0.832
# Residual   6   524.41 0.9649              
# Total      7   543.49 1.0000              
# 
# $therm_var_preheat_vs_control_postheat
# Df SumOfSqs      R2      F Pr(>F)  
# treatTime  1   174.61 0.34511 3.1619  0.032 *
#   Residual   6   331.33 0.65489                
# Total      7   505.94 1.00000                
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# $therm_var_preheat_vs_therm_var_postheat
# Df SumOfSqs      R2      F Pr(>F)  
# treatTime  1   167.88 0.44781 4.8659  0.022 *
#   Residual   6   207.01 0.55219                
# Total      7   374.89 1.00000                
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# $control_preheat_vs_control_postheat
# Df SumOfSqs      R2     F Pr(>F)
# treatTime  1    99.70 0.13868 0.966  0.363
# Residual   6   619.27 0.86132             
# Total      7   718.97 1.00000             
# 
# $control_preheat_vs_therm_var_postheat
# Df SumOfSqs      R2      F Pr(>F)
# treatTime  1    79.14 0.13786 0.9594  0.409
# Residual   6   494.95 0.86214              
# Total      7   574.09 1.00000              
# 
# $control_postheat_vs_therm_var_postheat
# Df SumOfSqs      R2      F Pr(>F)
# treatTime  1    16.41 0.05157 0.3262   0.81
# Residual   6   301.87 0.94843              
# Total      7   318.28 1.00000              
# 
# attr(,"class")
# [1] "pwadstrata" "list"   

AllT1<-subset(pca_s, Time=="Pre-heat")
pairwise.adonis2(AllT1[,1:2] ~ Treatment, data = AllT1, by = "terms", method='eu', na.rm = TRUE)
#N.S.

AllT2<-subset(pca_s, Time=="Post-heat")
pairwise.adonis2(AllT2[,1:2] ~ Treatment, data = AllT2, by = "terms", method='eu', na.rm = TRUE)
#N.S.

AllTr1<-subset(pca_s, Treatment=="Control")
pairwise.adonis2(AllTr1[,1:2] ~ Time, data = AllTr1, by = "terms", method='eu', na.rm = TRUE)
# $parent_call
# [1] "AllTr1[, 1:2] ~ time , strata = Null , permutations 999"
# 
# $preheat_vs_postheat
#           Df SumOfSqs      R2     F Pr(>F)  
# time      1   312.64 0.19872 3.472  0.044 *
# Residual 14  1260.63 0.80128               
# Total    15  1573.27 1.00000 

AllTr2<-subset(pca_s, Treatment=="DTV")
pairwise.adonis2(AllTr2[,1:2] ~ Time, data = AllTr2, by = "terms", method='eu', na.rm = TRUE)
# $parent_call
# [1] "AllTr2[, 1:2] ~ time , strata = Null , permutations 999"
# 
# $preheat_vs_postheat
#           Df SumOfSqs      R2      F Pr(>F)   
# time      1   391.80 0.38406 8.7297  0.005 **
# Residual 14   628.35 0.61594                 
# Total    15  1020.15 1.00000 

######vst looks the same as rlog- stick with rlog 

#### Calculate plasticity of gene expression between groups ####

library(ggbiplot) # plotting the PCA
# if(!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
#   BiocManager::install("ggfortify")
library(ggfortify) # plotting the PCA
library(vegan) # running the PERMANOVA (adonis2())
library(ggpubr) # for arranging multiple plots into single figure
source("~/Desktop/BU NSF Postdoc/Nicola_Ofav_chpt3/RNA analysis (Genome Counts)/Sym_GE/PCAplast_function.R") # source the plasticity function
####this function written by Colleen Bove and is available from her github

#https://github.com/seabove7/RandomFun/blob/main/Plasticity_function/PCAplast_function.R 
plast_out = PCAplast(pca = pca_s[,1:2], # the PCA dataframe containing the PCA eigenvalues and only the eigenvalues, can't insert prcomp object
         data = rt, # the condition/treatment data corresponding to samples
         sample_ID = "Samples", # the name of column that provide unique ID per sample (if blank, will pull rownames for this)
         num_pca =  "2", # the number of PCAs to include in analysis (default is 'all', but you can specify another number with a minimum of 2 PCAs)
         control_col = "TreatTime", # what the 'treatment' column is called
         control_lvl = "Control_Pre-heat", # control level of the treatment. If blank, a control mean per control level is assumed (when I leave this blank it throws an error and doesn't run)
         #group = "treatment" # the grouping column (i.e., colony). If blank, will assume control level grouping only!
         ) 

## Plot the plasticity (PC distances): overlay mean and 1 standard deviation
#reorder time variable for figures
f=c('Control_Post-heat', 'DTV_Pre-heat', 'DTV_Post-heat')
plast_out2 <- within(plast_out, TreatTime <- factor(TreatTime, levels=f))
cbPalette <- c("#9187FF","#E24A13","#FFA45E")
plast_plot <- ggplot(data = plast_out2, aes(x = TreatTime, y = dist, color = TreatTime)) + 
  geom_point(alpha = 0.5, position = position_jitter(width = 0.1)) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "errorbar", width = 0, size=1, colour = cbPalette) +
  stat_summary(fun = "mean", size = 0.8, colour = cbPalette) +
  guides(color = guide_legend(title = "Treatment x Time"))+
  theme_bw() +
  scale_colour_manual(values=cbPalette, labels = c( "Control - Post-heat", "DTV - Pre-heat", "DTV - Post-heat"))+
  scale_x_discrete(name ="Treatment x Time", 
                   labels=c("Control - Postheat","DTV - Preheat","DTV - Postheat"))+
  theme(plot.title = element_text(hjust = 0.5),
        title = element_text(size=14, face = "bold", color = "black"),
        legend.title=element_text(size=14, face = "bold", color = "black"), 
        legend.text=element_text(size=14, face = "bold", color="black"),
        axis.title.x = element_text(size = 16, color = "black", face = "bold"),
        axis.text = element_text(color = "black", size = 14, face = "bold"),
        axis.title.y = element_text(size = 16, color = "black", face = "bold"),
        axis.ticks = element_line(color = "black"))+
  annotate("text", x = min(0), y = max(20), 
           label = "Symbiodiniaceae", hjust = 0, vjust = 1, 
           size = 6, fontface = "bold", color = "black")+
  ylab("Gene Expression Plasticity")
plast_plot

#### Check significance of plasticity among groups
#assumption checks for ANOVA 
library("performance")
m <- lm(plast_out2$dist ~ TreatTime, data = plast_out2)
check_model(m)

hist(plast_out2$dist)#looks like a bell curve
shapiro.test(plast_out2$dist)#p-value = 0.8061, normal

#### two-way ANOVA
ANOVAmodel <- lm(plast_out2$dist ~ TreatTime, data = plast_out2)
anova(ANOVAmodel)
# Analysis of Variance Table
# 
# Response: plast_out2$dist
#           Df Sum Sq Mean Sq F value  Pr(>F)  
# treatTime  2 105.51  52.755  3.5809 0.04591 *
# Residuals 21 309.38  14.732     

ANOVAmodel <- lm(plast_out2$dist ~ Treatment*Time, data = plast_out2)
anova(ANOVAmodel)

# Analysis of Variance Table
# 
# Response: plast_out2$dist
#           Df  Sum Sq Mean Sq F value  Pr(>F)  
# treat      1  69.161  69.161  4.6945 0.04191 *
# time       1  36.349  36.349  2.4673 0.13118  
# Residuals 21 309.382  14.732 

#Tukey HSD test for significant pairwise comparisons
ANOVAmodel <- lm(plast_out2$dist ~ TreatTime, data=plast_out2)
SBPaov <- aov(ANOVAmodel)
TukeyHSD(SBPaov)
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = ANOVAmodel)
# 
# $treatTime
#                                       diff       lwr        upr     p adj
# therm_var_preheat-control_postheat   -5.108324 -9.945662 -0.2709853 0.0372747
# therm_var_postheat-control_postheat  -2.093807 -6.931146  2.7435311 0.5299447
# therm_var_postheat-therm_var_preheat  3.014516 -1.822822  7.8518548 0.2799291

#Tukey HSD test for significant pairwise comparisons
ANOVAmodel <- lm(plast_out2$dist ~ Treatment*Time, data=plast_out2)
SBPaov <- aov(ANOVAmodel)
TukeyHSD(SBPaov)
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = ANOVAmodel)
# 
# $treat
#                     diff       lwr        upr     p adj
# therm_var-control -3.601066 -7.057441 -0.1446896 0.0419095
# 
# $time
#                   diff       lwr      upr     p adj
# preheat-postheat -2.260887 -5.717263 1.195489 0.1881535

#### Effect Sizes ####
#effect size plasticity by treatTime
model <- aov(plast_out2$dist ~ TreatTime, data=plast_out2)

eta_squared(model, partial = TRUE)
# For one-way between subjects designs, partial eta squared is equivalent to eta squared. Returning eta
# squared.
# Effect Size for ANOVA
# 
# Parameter | Eta2 |       95% CI
# ---
#   treatTime | 0.25 | [0.00, 1.00]
# 
# - One-sided CIs: upper bound fixed at [1.00].

#effect size with all data
model <- aov(plast_out2$dist ~ Treatment*Time, data=plast_out2)

eta_squared(model, partial = TRUE)
# # Effect Size for ANOVA (Type I)
# 
# Parameter | Eta2 (partial) |       95% CI
# ---
# treat     |           0.18 | [0.00, 1.00]
# time      |           0.11 | [0.00, 1.00]
# 
# - One-sided CIs: upper bound fixed at [1.00].


#### VennDiagram of Genes ####

#MAKE SURE TO Re-RUN DESeq with the ~Species+treatTime through to res before running this part of code!!!
### These are lines: 125, 134, and 142

#Original DESeq design to compar Heat, DTV, and Host Species
#this wasn't very informative of differences among subgroups so ran differently in next section
# VENN Diagram to include both up and down regulated genes in common for transplant
resTrans1 <- results(dds, contrast=c("Time","Pre-heat", "Post-heat"))
resTrans2 <- results(dds, contrast=c("Treatment","DTV","Control"))
resTrans3 <- results(dds, contrast=c("Species","O. faveolata","O. franksi"))


head(resTrans1)
head(resTrans2)
head(resTrans3)

pTrans1up=row.names(resTrans1[resTrans1$padj<0.1 & !is.na(resTrans1$padj) & resTrans1$log2FoldChange>0,])
length(pTrans1up) #7
pTrans1down=row.names(resTrans2[resTrans1$padj<0.1 & !is.na(resTrans1$padj) & resTrans1$log2FoldChange<0,])
length(pTrans1down) #16
pTrans2up=row.names(resTrans2[resTrans2$padj<0.1 & !is.na(resTrans2$padj) & resTrans2$log2FoldChange>0,])
length(pTrans2up) #0
pTrans2down=row.names(resTrans2[resTrans2$padj<0.1 & !is.na(resTrans2$padj) & resTrans2$log2FoldChange<0,])
length(pTrans2down) #0
pTrans3up=row.names(resTrans3[resTrans3$padj<0.1 & !is.na(resTrans3$padj) & resTrans3$log2FoldChange>0,])
length(pTrans3up) #0
pTrans3down=row.names(resTrans3[resTrans3$padj<0.1 & !is.na(resTrans3$padj) & resTrans3$log2FoldChange<0,])
length(pTrans3down) #0

pTrans1=row.names(resTrans1[resTrans1$padj<0.1 & !is.na(resTrans1$padj),])
pTrans2=row.names(resTrans2[resTrans2$padj<0.1 & !is.na(resTrans2$padj),])
pTrans3=row.names(resTrans3[resTrans3$padj<0.1 & !is.na(resTrans3$padj),])

#UP
pdegs1_up=union(pTrans1up,pTrans2up)
pdegs1_up=union(pdegs1_up,pTrans3up)
length(pdegs1_up)
#7

#DOWN
pdegs1_down=union(pTrans1down,pTrans2down)
pdegs1_down=union(pdegs1_down,pTrans3down)
length(pdegs1_down)
#16

#ALL
pdegs1=union(pTrans1,pTrans2)
pdegs1=union(pdegs1,pTrans3)
length(pdegs1)
#23

###do UP, DOWN, ALL
candidates=list("Time"=pTrans1down, "Treatment"=pTrans2down, "HostSpecies"=pTrans3down)
candidates=list("Time"=pTrans1up, "Treatment"=pTrans2up, "HostSpecies"=pTrans3up)


#### venndiagram
#devtools::install_github("nicolash2/ggvenn")

library("ggvenn")
quartz()
ggplot() + 
  geom_venn(candidates, textsize=8) +
  scale_fill_manual(values=c("grey44", "#2C2673","#E24A13")) +
  theme_void()+
  theme(legend.position = "none")


#MAKE SURE TO Re-RUN DESeq with the ~Species+treatTime through to res before running this part of code!!!
### These are lines: 128, 134, and 142

### Fig. S7
# Venn with treatTime groupings, using DESeq ~Hostspecies+treatTime
# VENN Diagram to include both up and down regulated genes in common for transplant
resTrans1 <- results(dds, contrast=c("TreatTime","Control_Pre-heat", "Control_Post-heat"))
resTrans2 <- results(dds, contrast=c("TreatTime","DTV_Pre-heat","control_preheat"))
resTrans3 <- results(dds, contrast=c("TreatTime","DTV_Post-heat","Control_Post-heat"))
resTrans4 <- results(dds, contrast=c("TreatTime","DTV_Pre-heat","DTV_Post-heat"))

head(resTrans1)
head(resTrans2)
head(resTrans3)
head(resTrans4)

pTrans1up=row.names(resTrans1[resTrans1$padj<0.1 & !is.na(resTrans1$padj) & resTrans1$log2FoldChange>0,])
length(pTrans1up) #7
pTrans1down=row.names(resTrans1[resTrans1$padj<0.1 & !is.na(resTrans1$padj) & resTrans1$log2FoldChange<0,])
length(pTrans1down) #16
pTrans2up=row.names(resTrans2[resTrans2$padj<0.1 & !is.na(resTrans2$padj) & resTrans2$log2FoldChange>0,])
length(pTrans2up) #0
pTrans2down=row.names(resTrans2[resTrans2$padj<0.1 & !is.na(resTrans2$padj) & resTrans2$log2FoldChange<0,])
length(pTrans2down) #0
pTrans3up=row.names(resTrans3[resTrans3$padj<0.1 & !is.na(resTrans3$padj) & resTrans3$log2FoldChange>0,])
length(pTrans3up) #0
pTrans3down=row.names(resTrans3[resTrans3$padj<0.1 & !is.na(resTrans3$padj) & resTrans3$log2FoldChange<0,])
length(pTrans3down) #0
pTrans4up=row.names(resTrans4[resTrans4$padj<0.1 & !is.na(resTrans4$padj) & resTrans4$log2FoldChange>0,])
length(pTrans4up) #11
pTrans4down=row.names(resTrans4[resTrans4$padj<0.1 & !is.na(resTrans4$padj) & resTrans4$log2FoldChange<0,])
length(pTrans4down) #29


pTrans1=row.names(resTrans1[resTrans1$padj<0.1 & !is.na(resTrans1$padj),])
pTrans2=row.names(resTrans2[resTrans2$padj<0.1 & !is.na(resTrans2$padj),])
pTrans3=row.names(resTrans3[resTrans3$padj<0.1 & !is.na(resTrans3$padj),])
pTrans4=row.names(resTrans4[resTrans4$padj<0.1 & !is.na(resTrans4$padj),])

#UP
pdegs1_up=union(pTrans1up,pTrans2up)
pdegs1_up=union(pdegs1_up,pTrans3up)
pdegs1_up=union(pdegs1_up,pTrans4up)
length(pdegs1_up)
#16

#DOWN
pdegs1_down=union(pTrans1down,pTrans2down)
pdegs1_down=union(pdegs1_down,pTrans3down)
pdegs1_down=union(pdegs1_down,pTrans4down)
length(pdegs1_down)
#41

#ALL
pdegs1=union(pTrans1,pTrans2)
pdegs1=union(pdegs1,pTrans3)
pdegs1=union(pdegs1,pTrans4)
length(pdegs1)
#57

###do UP, DOWN, ALL
###do UP, DOWN, ALL
candidatesUp=list("Pre-heat DTV-Control"=pTrans2up, "Control pre-post"=pTrans1up, "DTV pre-post"=pTrans4up, "Post-heat DTV-Control"=pTrans3up)#, "Host - OfavOfra"=pTrans5up)
candidatesDown=list("Pre-heat DTV-Control"=pTrans2down, "Control pre-post"=pTrans1down, "DTV pre-post"=pTrans4down, "Post-heat DTV-Control"=pTrans3down)#, "Host - OfavOfra"=pTrans5down)


ggvenn(
  candidatesUp,
  fill_color = c("grey25", "grey25", "grey25", "grey25"),
  stroke_size = 0,
  set_name_size = 5,
  text_size = 5,
  show_percentage = FALSE)

ggvenn(
  candidatesDown,
  fill_color = c("grey25", "grey25", "grey25", "grey25"),
  stroke_size = 0,
  set_name_size = 5,
  text_size = 5,
  show_percentage = FALSE)

###Decided not to do WGCNA for Symbionts due to limited genes availale for analysis
