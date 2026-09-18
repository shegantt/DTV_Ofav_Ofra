#####Ofav DTV RNAseq analysis based on Colleen Boves code (https://github.com/seabove7/BelizeRT_Castillo_Bove/blob/main/RT_Host_Sym_GE.R)

###host data
setwd("~/")

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
library("ggvenn")

 # if(!requireNamespace("BiocManager", quietly = TRUE))
 #   install.packages("BiocManager")
 # BiocManager::install("preprocessCore")

##### First going to get the count of genes in all host samples ####
#read in counts
countData <- read.table("~/Ofav_counts_fixedSNP.txt")
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
min(totalCounts) #451,076
max(totalCounts)  #1,501,208
mean(totalCounts) #1,018,040

#Count the number of genes in each sample
nrow(countData)#16,909 for clones included dataset
SampleGenes = colSums(countData != 0)
NumSampleGenes <- data.frame(Column = names(countData), NonZeroCount = SampleGenes)


##### Now do real analysis with clones removed ####
#read in counts
countData <- read.table("~/FinalCounts_13Nov2024/Ofav_countsfixed_fixedSNP_noClones.txt")
head(countData)
length(countData[,1]) #35826 genes

names(countData)=sub(".fastq.host.counts.txt","",names(countData))
names(countData)
head(countData)

totalCounts=colSums(countData)
barplot(totalCounts, col="coral")
totalCounts

#fixedSNP_noClones
# KB13.2  KB28.1  KB39.1  KB39.2  KE13.1  KE14.1  KE31.2   KE9.2  KF19.1  KF19.2  KF23.2   KF9.1 
# 1066048  983062  452141 1207688  932055 1035177 1158441 1038914 1123245 1176579 1051869  628007 
# KG12.1  KG12.2  KG31.1  KG31.2  VD10.1  VD15.2  VD17.1  VD17.2  VE14.2  VE18.1  VE18.2  VE30.1 
# 1051820 1099529 1222731  900488 1130491 1182197  984099 1504792  828850 1374223 1406824 1019343 
# VF23.2  VF31.1  VF40.1  VF40.2   VG1.1   VG1.2  VG17.2   VG6.1 
# 1210168  979962  811199  758409 1371576  743384 1089053  499116
  
min(totalCounts) #452,141
max(totalCounts)  #1,504,792
mean(totalCounts) #noClones- 1,031,921

#Be careful here, sample names on rows are stripped from the metadata
#in these steps, so order of samples in count table MUST match the
#order of samples in metadata, otherwise the metadata will get jumbled

#import metadata
rt <- read.csv("~/THERMVAR_MAIN_METADATA_Apr25_clonesfixed.csv")

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
min(totalCounts) #450921
max(totalCounts)  #1500543
mean(totalCounts) #1029131

ncol(countData) #32
nrow(conditions) #32
nrow(countData)#16,447 genes



#initially tried with design (and used to make Venn Diagram) but not very informative
#dds<-DESeqDataSetFromMatrix(countData=countData, colData=conditions, design=~Species+Treatment*Time) #can only test for the main effects species, treatment, and time

#for venn diagrams, to get species influence and treatment effects
#dds<-DESeqDataSetFromMatrix(countData=countData, colData=conditions, design=~Species+TreatTime) 
#ran above to look at effect of time and treatment as one variable, wasn't significant so used first design for analysis
# I did use the treatTime design to generate PCA for this comparison, no significance was found with adonis2.

#DESeq without design (~1), (used for all PCA and WGCNA analyses)
dds<-DESeqDataSetFromMatrix(countData=countData, colData=conditions, design=~1) 

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
#save(rlog, file="rlogCountsGreedy.RData") #save for WGCNA 
rld=assay(rlog)
rld_t=t(rld)
colnames(rld_t)
head(rld_t)
pca <- prcomp(rld_t,center = TRUE) # also ran with PCA() that uses top 500 DEGs, but got less significance
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

#Fig 6. Host PCA - all three variables on same PCA #6A418B"
cbPalette <- c(  "#2C2673", "#9187FF","#E24A13","#FFA45E")#treatment
#pdf("PCA_Host_allgenes_rlog.pdf",height=5,width=6)
ggplot(pca_s, aes(PC1, PC2, color = TreatTime, shape=Species, group=TreatTime)) +
  geom_point(size = 4) +
  theme_bw() +
  guides(color = guide_legend(title = "Treatment x Heat"), shape = guide_legend(title = "Host Species"))+
  scale_colour_manual(values=cbPalette, labels = c("Control - Pre-heat", "Control - Post-heat",
                                                   "DTV - Pre-heat", "DTV - Post-heat"))+
  theme( legend.title=element_text(size=14, face = "bold", color = "black"), 
        legend.text=element_text(size=14, face = "bold", color="black"),
        axis.title.x = element_text(size = 16, color = "black", face = "bold"),
        axis.text = element_text(color = "black", size = 14, face = "bold"),
        axis.title.y = element_text(size = 16, color = "black", face = "bold"),
        axis.ticks = element_line(color = "black"))+
  annotate("text", x = min(-65), y = max(60), 
             label = "Host", hjust = 0, vjust = 1, 
             size = 6, fontface = "bold", color = "black")+
  stat_ellipse()+
  xlab(paste0("PC1: ",pc1v,"% variance")) +
  ylab(paste0("PC2: ",pc2v,"% variance"))
#dev.off()


#faceted by species
# New facet label names for Species variable
sp.labs <- c("O. faveolata", "O. franksi")
names(sp.labs) <- c("Ofav", "Ofra")
#all three variables on same PCA
#cbPalette <- c( "#B395DA","#4A2E61", "#D89362", "#8F3C00")#Treat v Time
cbPalette <- c( "#241F5C","#744899", "#EC4176", "#FFA45E")#treatTime
#pdf("PCA_Host_allgenes_rlog.pdf",height=5,width=6)
ggplot(pca_s, aes(PC1, PC2, color = TreatTime, shape=Genotype, group=TreatTime)) +
  facet_grid2(Species~., labeller = labeller(species = sp.labs), 
              strip = strip_themed(background_y = elem_list_rect(fill = "black"), 
                             text_y = elem_list_text(size = 14, face = "bold", color = "white")))+
  geom_point(size = 4, stroke = 1) +
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


#Fig. S6 Host - faceted by treatment and species
cbPalette <- c( "#2C2673", "#9187FF","#E24A13","#FFA45E")#treatTime
# New facet label names for treatment variable
treat.labs <- c("Control", "DTV")
names(treat.labs) <- c("control", "therm_var")
#pdf("PCA_Host_allgenes_rlog.pdf",height=5,width=6)
ggplot(pca_s, aes(PC1, PC2, color = TreatTime, shape= Genotype, group=TreatTime)) +
  facet_grid2(Treatment~Species, labeller = labeller(species = sp.labs, treatment =treat.labs), 
              strip = strip_themed(background_x = elem_list_rect(fill = c("O. faveolata"= "#8F3C00", "O. franksi"="#D89362")),
                             background_y = elem_list_rect(fill = c("Control"= "#2C2673", "DTV"="#E24A13")), 
                             text_y = elem_list_text(size = 14, face = "bold", color = "white"),
                             text_x = elem_list_text(size = 14, face = "bold", color = "white")))+
  geom_point(size = 4, stroke = 1) +
  scale_shape_manual(values=c(16,15,18,17, 21,22,23,24))+
  theme_bw() +
  guides(color = guide_legend(title = "Treatment x Heat"), shape = guide_legend(title = "Host Genotype"))+
  scale_colour_manual(values=cbPalette, labels = c("Control - Pre-heat", "Control - Post-heat", 
                                                   "DTV - Pre-heat", "DTV - Post-heat"))+
  theme( plot.title = element_text(hjust = 0.5),
         title = element_text(size=14, face = "bold", color = "black"),
         legend.title=element_text(size=14, face = "bold", color = "black"), 
         legend.text=element_text(size=14, face = "bold", color="black"),
         axis.title.x = element_text(size = 16, color = "black", face = "bold"),
         axis.text = element_text(color = "black", size = 14, face = "bold"),
         axis.title.y = element_text(size = 16, color = "black", face = "bold"),
         axis.ticks = element_line(color = "black"))+
  stat_ellipse()+
  ggtitle("Host GE")+
  xlab(paste0("PC1: ",pc1v,"% variance")) +
  ylab(paste0("PC2: ",pc2v,"% variance")) 
#dev.off()

#Fig. S6 Host - faceted by treatment and species
cbPalette <- c( "#2C2673", "#9187FF","#E24A13","#FFA45E")#treatTime
# New facet label names for treatment variable
treat.labs <- c("Stable", "DTV")
names(treat.labs) <- c("control", "therm_var")
#pdf("PCA_Host_allgenes_rlog.pdf",height=5,width=6)
ggplot(pca_s, aes(PC1, PC2, color = TreatTime, shape= Species, group=TreatTime)) +
  facet_grid2(Treatment~Species, labeller = labeller(species = sp.labs, treatment =treat.labs), 
              strip = strip_themed(background_x = elem_list_rect(fill = c("O. faveolata"= "#8b7355", "O. franksi"="#cdab7d")),
                                   background_y = elem_list_rect(fill = c("Stable"= "#2C2673", "DTV"="#E24A13")), 
                                   text_y = elem_list_text(size = 14, face = "bold", color = "white"),
                                   text_x = elem_list_text(size = 14, face = "bold", color = "white")))+
  geom_point(size = 4, stroke = 1) +
  scale_shape_manual(values=c(19, 17))+
  #scale_shape_manual(values=c(16,15,18,17, 21,22,23,24))+ #genotype shapes
  theme_bw() +
  guides(color = guide_legend(title = "Treatment x Heat"), shape = guide_legend(title = "Host Genotype"))+
  scale_colour_manual(values=cbPalette, labels = c("Control - Pre-heat", "Control - Post-heat", 
                                                   "DTV - Pre-heat", "DTV - Post-heat"))+
  theme( plot.title = element_text(hjust = 0.5),
         title = element_text(size=14, face = "bold", color = "black"),
         legend.title=element_text(size=14, face = "bold", color = "black"), 
         legend.text=element_text(size=14, face = "bold", color="black"),
         axis.title.x = element_text(size = 16, color = "black", face = "bold"),
         axis.text = element_text(color = "black", size = 14, face = "bold"),
         axis.title.y = element_text(size = 16, color = "black", face = "bold"),
         axis.ticks = element_line(color = "black"))+
  stat_ellipse()+
  ggtitle("Host GE")+
  xlab(paste0("PC1: ",pc1v,"% variance")) +
  ylab(paste0("PC2: ",pc2v,"% variance")) 
#dev.off()


#faceted by time and species
cbPalette <- c( "#2C2673", "#9187FF","#E24A13","#FFA45E")#treatTime
# New facet label names for time variable
time.labs <- c("Pre-heat", "Post-heat")
names(time.labs) <- c("preheat", "postheat")
#pdf("PCA_Host_allgenes_rlog.pdf",height=5,width=6)
ggplot(pca_s, aes(PC1, PC2, color = TreatTime, shape= Genotype, group=TreatTime)) +
  facet_grid2(Time~Species, labeller = labeller(species = sp.labs, time =time.labs), 
              strip = strip_themed(background_x = elem_list_rect(fill = "black"),
                             background_y = elem_list_rect(fill = c("preheat"="grey20", "postheat"="grey70")), 
                             text_y = elem_list_text(size = 14, face = "bold", color = "white"),
                             text_x = elem_list_text(size = 14, face = "bold", color = "white")))+
  geom_point(size = 4, stroke = 1) +
  scale_shape_manual(values=c(16,15,18,17, 21,22,23,24))+
  theme_bw() +
  guides(color = guide_legend(title = "Treatment x Heat"), shape = guide_legend(title = "Host Species"))+
  scale_colour_manual(values=cbPalette, labels = c("Control - Pre-heat", "Control - Post-heat", 
                                                   "DTV - Pre-heat", "DTV - Post-heat"))+
  theme(plot.title = element_text(hjust = 0.5),
        title = element_text(size=14, face = "bold", color = "black"),
        legend.title=element_text(size=14, face = "bold", color = "black"), 
        legend.text=element_text(size=14, face = "bold", color="black"),
        axis.title.x = element_text(size = 16, color = "black", face = "bold"),
        axis.text = element_text(color = "black", size = 14, face = "bold"),
        axis.title.y = element_text(size = 16, color = "black", face = "bold"),
        axis.ticks = element_line(color = "black"))+
  stat_ellipse()+
  ggtitle("Host GE")+
  xlab(paste0("PC1: ",pc1v,"% variance")) +
  ylab(paste0("PC2: ",pc2v,"% variance")) 
#dev.off()


adonis2(pca_s[,1:2] ~ Species+Treatment*Time, data = pca_s, by = "terms", method='eu', na.rm = TRUE)
#DESeq = ~1
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = pca_s[, 1:2] ~ species + treatment * time, data = pca_s, method = "eu", by = "terms", na.rm = TRUE)
#                 Df SumOfSqs      R2       F Pr(>F)    
# species         1  12184.2 0.52931 34.7448  0.001 ***
# treatment       1      8.5 0.00037  0.0241  0.969    
# time            1   1354.4 0.05884  3.8622  0.027 *  
# treatment:time  1      3.9 0.00017  0.0111  0.988    
# Residual       27   9468.3 0.41132                   
# Total          31  23019.3 1.00000     

#PERMDISP:
dist_s <- vegdist(pca_s[,1:2], method = "euclidean")

# PERMDISP for species
bd_Species <- betadisper(dist_s, pca_s$Species)
anova(bd_Species)                      # or: permutest(bd_species, permutations = 999)
# Response: Distances
# Df Sum Sq Mean Sq F value Pr(>F)
# Groups     1  146.3  146.31  1.1187 0.2986
# Residuals 30 3923.4  130.78 

# PERMDISP for treatment
bd_Treatment <- betadisper(dist_s, pca_s$Treatment)
anova(bd_Treatment)  
# Response: Distances
# Df  Sum Sq Mean Sq F value Pr(>F)
# Groups     1   50.46  50.455  0.5077 0.4816
# Residuals 30 2981.21  99.374 

#PERMDIP for time
bd_Time <- betadisper(dist_s, pca_s$Time)
anova(bd_treatTime) 
# Response: Distances
# Df  Sum Sq Mean Sq F value Pr(>F)
# Groups     1    0.11   0.114  0.0012  0.973
# Residuals 30 2957.07  98.569  

#tried with combined variable for treatment and time (Heat)
adonis2(pca_s[,1:2] ~ Species+TreatTime, data = pca_s, by = "terms", method='eu', na.rm = TRUE)
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = pca_s[, 1:2] ~ species + treatTime, data = pca_s, method = "eu", by = "terms", na.rm = TRUE)
#           Df SumOfSqs      R2       F Pr(>F)    
# species    1  11807.3 0.52583 34.0831  0.001 ***
# treatTime  3   1294.0 0.05763  1.2451  0.305    
# Residual  27   9353.5 0.41655                   
# Total     31  22454.8 1.00000  

#subset by host species to reduce noise and see if any trends emerge
#Same trends and significance as analysis above
ofav<-subset(pca_s, Species=="O. faveolata")

adonis2(ofav[,1:2] ~ Genotype+Treatment*Time, data = ofav,by = "terms", method='eu', na.rm = TRUE)
#DESeq =~1
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = ofav[, 1:2] ~ genotype + treatment * time, data = ofav, method = "eu", by = "terms", na.rm = TRUE)
#                 Df SumOfSqs      R2       F Pr(>F)    
# genotype        3   2682.9 0.63306 12.4646  0.001 ***
# treatment       1      9.3 0.00220  0.1301  0.850    
# time            1    867.8 0.20478 12.0961  0.003 ** 
# treatment:time  1     32.2 0.00759  0.4485  0.605    
# Residual        9    645.7 0.15237                   
# Total          15   4237.9 1.00000  

ofra<-subset(pca_s, Species=="O. franksi")

adonis2(ofra[,1:2] ~ Genotype+Treatment*Time, data = ofra, by = "terms",method='eu', na.rm = TRUE)
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = ofra[, 1:2] ~ genotype + treatment * time, data = ofra, method = "eu", by = "terms", na.rm = TRUE)
#                 Df SumOfSqs      R2       F Pr(>F)    
# genotype        3   5410.9 0.82020 25.3180  0.001 ***
# treatment       1      1.3 0.00020  0.0189  0.957    
# time            1    510.4 0.07737  7.1644  0.028 *  
# treatment:time  1     33.3 0.00504  0.4668  0.508    
# Residual        9    641.2 0.09719                   
# Total          15   6597.1 1.00000   


#comparisons to pairwise.adonis2
#Ofav preheat by treatment
ofavT1<-subset(ofav, Time=="Pre-heat")
adonis2(ofavT1[,1:2] ~ Treatment, data = ofavT1,by = "terms", method='eu', na.rm = TRUE)
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# adonis2(formula = ofavT1[, 1:2] ~ treatment, data = ofavT1, method = "eu", by = "terms", na.rm = TRUE)
#           Df SumOfSqs      R2      F Pr(>F)
# treatment  1    27.53 0.02209 0.1355  0.921
# Residual   6  1218.79 0.97791              
# Total      7  1246.32 1.00000 

#Ofav postheat by treatment
ofavT2<-subset(ofav, Time=="Post-heat")
adonis2(ofavT2[,1:2] ~ Treatment, data = ofavT2,by = "terms", method='eu', na.rm = TRUE)
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# adonis2(formula = ofavT2[, 1:2] ~ treatment, data = ofavT2, method = "eu", by = "terms", na.rm = TRUE)
#           Df SumOfSqs      R2      F Pr(>F)
# treatment  1    13.98 0.00658 0.0398  0.881
# Residual   6  2109.79 0.99342              
# Total      7  2123.77 1.00000 

#Ofav controls by time
ofavTr1<-subset(ofav, Treatment=="Control")
adonis2(ofavTr1[,1:2] ~ Time, data = ofavTr1,by = "terms", method='eu', na.rm = TRUE)
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# adonis2(formula = ofavTr1[, 1:2] ~ time, data = ofavTr1, method = "eu", by = "terms", na.rm = TRUE)
#           Df SumOfSqs      R2      F Pr(>F)
# time      1   304.94 0.11883 0.8091  0.458
# Residual  6  2261.32 0.88117              
# Total     7  2566.25 1.00000 

#Ofav DTV by time
ofavTr2<-subset(ofav, Treatment=="DTV")
adonis2(ofavTr2[,1:2] ~ Time, data = ofavTr2,by = "terms", method='eu', na.rm = TRUE)
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# adonis2(ofavTr1[,1:2] ~ time, data = ofavTr1,by = "terms", method='eu', na.rm = TRUE)
# adonis2(formula = ofavTr2[, 1:2] ~ time, data = ofavTr2, method = "eu", by = "terms", na.rm = TRUE)
#           Df SumOfSqs      R2      F Pr(>F)
# time      1   595.09 0.35798 3.3455  0.128
# Residual  6  1067.26 0.64202              
# Total     7  1662.35 1.00000     


#ofra preheat by treatment
ofraT1<-subset(ofra, Time=="Pre-heat")
adonis2(ofraT1[,1:2] ~ Treatment, data = ofraT1,by = "terms", method='eu', na.rm = TRUE)
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = ofraT1[, 1:2] ~ treatment, data = ofraT1, method = "eu", by = "terms", na.rm = TRUE)
#           Df SumOfSqs      R2      F Pr(>F)
# treatment  1     14.6 0.00357 0.0215  0.879
# Residual   6   4081.6 0.99643              
# Total      7   4096.2 1.00000 

#ofra postheat by treatment
ofraT2<-subset(ofra, Time=="Post-heat")
adonis2(ofraT2[,1:2] ~ Treatment, data = ofraT2,by = "terms", method='eu', na.rm = TRUE)
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = ofraT2[, 1:2] ~ treatment, data = ofraT2, method = "eu", by = "terms", na.rm = TRUE)
#           Df SumOfSqs      R2      F Pr(>F)
# treatment  1    19.97 0.01003 0.0608  0.868
# Residual   6  1970.52 0.98997              
# Total      7  1990.49 1.00000 

#Ofra controls by time
ofraTr1<-subset(ofra, Treatment=="Control")
adonis2(ofraTr1[,1:2] ~ Time, data = ofraTr1,by = "terms", method='eu', na.rm = TRUE)
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = ofraTr1[, 1:2] ~ time, data = ofraTr1, method = "eu", by = "terms", na.rm = TRUE)
#           Df SumOfSqs     R2      F Pr(>F)
# time      1    401.9 0.1047 0.7017  0.446
# Residual  6   3436.5 0.8953              
# Total     7   3838.4 1.0000  

#ofra DTV by time
ofraTr2<-subset(ofra, Treatment=="DTV")
adonis2(ofraTr2[,1:2] ~ Time, data = ofraTr2,by = "terms", method='eu', na.rm = TRUE)
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = ofraTr2[, 1:2] ~ time, data = ofraTr2, method = "eu", by = "terms", na.rm = TRUE)
#           Df SumOfSqs      R2      F Pr(>F)
# time      1   141.75 0.05141 0.3252  0.597
# Residual  6  2615.60 0.94859              
# Total     7  2757.34 1.00000    



library(pairwiseAdonis) 
# install.packages('devtools') #To install package
# library(devtools)
# install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")


#both coral species pairwise PERMANOVA by treatTime
pairwise.adonis2(pca_s[,1:2] ~ TreatTime, data = pca_s, by = "terms", method='eu', na.rm = TRUE)

#DESeq =~1
# $parent_call
# [1] "pca_s[, 1:2] ~ species*treatTime , strata = Null , permutations 999"

# $therm_var_preheat_vs_control_preheat
#           Df SumOfSqs      R2      F Pr(>F)
# treatTime  1     11.6 0.00102 0.0143  0.987
# Residual  14  11396.4 0.99898              
# Total     15  11408.0 1.00000              
# 
# $therm_var_preheat_vs_control_postheat
#           Df SumOfSqs     R2      F Pr(>F)
# treatTime  1      682 0.0616 0.9189  0.421
# Residual  14    10391 0.9384              
# Total     15    11073 1.0000              
# 
# $therm_var_preheat_vs_therm_var_postheat
#           Df SumOfSqs      R2      F Pr(>F)
# treatTime  1    656.4 0.06337 0.9472  0.395
# Residual  14   9702.4 0.93663              
# Total     15  10358.8 1.00000              
# 
# $control_preheat_vs_control_postheat
#           Df SumOfSqs      R2      F Pr(>F)
# treatTime  1    701.9 0.05547 0.8222  0.464
# Residual  14  11950.1 0.94453              
# Total     15  12652.0 1.00000              
# 
# $control_preheat_vs_therm_var_postheat
#           Df SumOfSqs      R2      F Pr(>F)
# treatTime  1    680.8 0.05701 0.8463  0.454
# Residual  14  11261.5 0.94299              
# Total     15  11942.3 1.00000              
# 
# $control_postheat_vs_therm_var_postheat
#           Df SumOfSqs      R2     F Pr(>F)
# treatTime  1      0.7 0.00007 0.001  0.997
# Residual  14  10256.2 0.99993             
# Total     15  10256.9 1.00000  


#PERMDISP
### Check your TreatTime groups
unique(pca_s$TreatTime)
table(pca_s$TreatTime)


### Create the Euclidean distance matrix
dist_euAll <- dist(pca_s[, 1:2], method = "euclidean")


### Get the TreatTime groups
groups <- unique(pca_s$TreatTime)

permdisp_results <- data.frame(
  Group1 = character(),
  Group2 = character(),
  F = numeric(),
  p = numeric(),
  stringsAsFactors = FALSE
)

### Run pairwise PERMDISP
for (i in 1:(length(groups) - 1)) {
  
  for (j in (i + 1):length(groups)) {
    
    group1 <- groups[i]
    group2 <- groups[j]
    
    # Keep only samples from these two TreatTime groups
    keep <- pca_s$TreatTime %in% c(group1, group2)
    
    # Subset distance matrix
    d_pair <- as.dist(
      as.matrix(dist_euAll)[keep, keep]
    )
    
    # Subset grouping variable
    group_pair <- droplevels(
      factor(Ofav$TreatTime[keep])
    )
    
    # Run PERMDISP
    bd <- betadisper(
      d_pair,
      group_pair
    )
    
    # Test for differences in dispersion
    test <- anova(bd)
    
    # Add results
    permdisp_results <- rbind(
      permdisp_results,
      data.frame(
        Group1 = group1,
        Group2 = group2,
        F = test$`F value`[1],
        p = test$`Pr(>F)`[1]
      )
    )
  }
}

### Adjust p-values for multiple comparisons
permdisp_results$p_adjusted <- p.adjust(
  permdisp_results$p,
  method = "BH"
)


### View results
permdisp_results


pairwise.adonis2(ofav[,1:2] ~ TreatTime, data = ofav, by = "terms", method='eu', na.rm = TRUE)
#none were significant
# $parent_call
# [1] "ofav[, 1:2] ~ treatTime , strata = Null , permutations 999"
# 
# $control_preheat_vs_therm_var_preheat
# Df SumOfSqs      R2      F Pr(>F)
# treatTime  1    27.53 0.02209 0.1355  0.928
# Residual   6  1218.79 0.97791              
# Total      7  1246.32 1.00000              
# 
# $control_preheat_vs_control_postheat
# Df SumOfSqs      R2      F Pr(>F)
# treatTime  1   304.94 0.11883 0.8091  0.457
# Residual   6  2261.32 0.88117              
# Total      7  2566.25 1.00000              
# 
# $control_preheat_vs_therm_var_postheat
# Df SumOfSqs      R2      F Pr(>F)
# treatTime  1   449.16 0.21952 1.6876  0.283
# Residual   6  1596.90 0.78048              
# Total      7  2046.06 1.00000              
# 
# $therm_var_preheat_vs_control_postheat
# Df SumOfSqs      R2     F Pr(>F)
# treatTime  1   428.02 0.19819 1.483  0.298
# Residual   6  1731.68 0.80181             
# Total      7  2159.70 1.00000             
# 
# $therm_var_preheat_vs_therm_var_postheat
# Df SumOfSqs      R2      F Pr(>F)
# treatTime  1   595.09 0.35798 3.3455  0.125
# Residual   6  1067.26 0.64202              
# Total      7  1662.35 1.00000              
# 
# $control_postheat_vs_therm_var_postheat
# Df SumOfSqs      R2      F Pr(>F)
# treatTime  1    13.98 0.00658 0.0398  0.884
# Residual   6  2109.79 0.99342              
# Total      7  2123.77 1.00000   

#PERMDISP
### Check your TreatTime groups
unique(ofav$TreatTime)
table(ofav$TreatTime)


### Create the Euclidean distance matrix
dist_euOfav <- dist(0fav[, 1:2], method = "euclidean")


### Get the TreatTime groups
groups <- unique(ofav$TreatTime)

permdisp_results <- data.frame(
  Group1 = character(),
  Group2 = character(),
  F = numeric(),
  p = numeric(),
  stringsAsFactors = FALSE
)

### Run pairwise PERMDISP
for (i in 1:(length(groups) - 1)) {
  
  for (j in (i + 1):length(groups)) {
    
    group1 <- groups[i]
    group2 <- groups[j]
    
    # Keep only samples from these two TreatTime groups
    keep <- ofav$TreatTime %in% c(group1, group2)
    
    # Subset distance matrix
    d_pair <- as.dist(
      as.matrix(dist_euOfav)[keep, keep]
    )
    
    # Subset grouping variable
    group_pair <- droplevels(
      factor(ofav$TreatTime[keep])
    )
    
    # Run PERMDISP
    bd <- betadisper(
      d_pair,
      group_pair
    )
    
    # Test for differences in dispersion
    test <- anova(bd)
    
    # Add results
    permdisp_results <- rbind(
      permdisp_results,
      data.frame(
        Group1 = group1,
        Group2 = group2,
        F = test$`F value`[1],
        p = test$`Pr(>F)`[1]
      )
    )
  }
}

### Adjust p-values for multiple comparisons
permdisp_results$p_adjusted <- p.adjust(
  permdisp_results$p,
  method = "BH"
)


### View results
permdisp_results

pairwise.adonis2(ofra[,1:2] ~ TreatTime, data = ofra, by = "terms", method='eu', na.rm = TRUE)
#none were significant
# $parent_call
# [1] "ofra[, 1:2] ~ treatTime , strata = Null , permutations 999"
# 
# $therm_var_preheat_vs_control_preheat
# Df SumOfSqs      R2      F Pr(>F)
# treatTime  1     14.6 0.00357 0.0215  0.898
# Residual   6   4081.6 0.99643              
# Total      7   4096.2 1.00000              
# 
# $therm_var_preheat_vs_control_postheat
# Df SumOfSqs      R2      F Pr(>F)
# treatTime  1   264.97 0.08821 0.5805   0.35
# Residual   6  2738.86 0.91179              
# Total      7  3003.83 1.00000              
# 
# $therm_var_preheat_vs_therm_var_postheat
# Df SumOfSqs      R2      F Pr(>F)
# treatTime  1   141.75 0.05141 0.3252  0.602
# Residual   6  2615.60 0.94859              
# Total      7  2757.34 1.00000              
# 
# $control_preheat_vs_control_postheat
# Df SumOfSqs     R2      F Pr(>F)
# treatTime  1    401.9 0.1047 0.7017  0.409
# Residual   6   3436.5 0.8953              
# Total      7   3838.4 1.0000              
# 
# $control_preheat_vs_therm_var_postheat
# Df SumOfSqs      R2      F Pr(>F)
# treatTime  1    246.8 0.06932 0.4469  0.518
# Residual   6   3313.2 0.93068              
# Total      7   3560.0 1.00000              
# 
# $control_postheat_vs_therm_var_postheat
# Df SumOfSqs      R2      F Pr(>F)
# treatTime  1    19.97 0.01003 0.0608  0.856
# Residual   6  1970.52 0.98997              
# Total      7  1990.49 1.00000       

# ##try again with vst normalization
# vst_norm=vst(dds, blind=TRUE) 
# v=assay(vst_norm)
# v_t=t(v)
# colnames(v_t)
# head(v_t)
# pca <- prcomp(v_t,center = TRUE)
# head(pca)
# li <- pca$sdev^2 / sum(pca$sdev^2)
# pc1v <- round(li[1] * 100, 1)
# pc2v <- round(li[2] * 100, 1)
# pca_s <- as.data.frame(pca$x)
# head(pca_s)
# pca_s <- pca_s[,c(1,2)]
# pca_s$Samples = row.names(pca_s)
# pca_s$treat=conditions$treat
# pca_s$species=conditions$species
# pca_s$time=conditions$time
# head(pca_s)
# 
# 
# #all three variables on same PCA
# cbPalette <- c( "red","ivory4")
# #pdf("PCA_Host_allgenes_vst.pdf",height=5,width=6)
# ggplot(pca_s, aes(PC1, PC2, color = time, shape=species, fill=treat, group=treat)) +
#   geom_point(size = 4) +
#   #geom_label(aes(label = Samples))+
#   #geom_text_repel(aes(label=Samples)) +
#   scale_shape_manual(values = c(21, 24) )+
#   scale_colour_manual(values=cbPalette)+
#   scale_fill_manual(values=c("grey", "black"))+
#   guides(fill=guide_legend(override.aes=list(shape=21)))+
#   guides(color=guide_legend(override.aes=list(shape=21)))+
#   guides(fill=guide_legend(override.aes=list(shape=21,fill=c("grey", "black"))))+
#   theme_bw() +
#   stat_ellipse()+
#   # geom_density2d(alpha=.5)+
#   # geom_polygon(alpha=.2)+
#   xlab(paste0("PC1: ",pc1v,"% variance")) +
#   ylab(paste0("PC2: ",pc2v,"% variance")) 
# #dev.off()

#PERMDISP
### Check your TreatTime groups
unique(ofra$TreatTime)
table(ofra$TreatTime)


### Create the Euclidean distance matrix
dist_eu <- dist(ofra[, 1:2], method = "euclidean")


### Get the TreatTime groups
groups <- unique(ofra$TreatTime)

permdisp_results <- data.frame(
  Group1 = character(),
  Group2 = character(),
  F = numeric(),
  p = numeric(),
  stringsAsFactors = FALSE
)

### Run pairwise PERMDISP
for (i in 1:(length(groups) - 1)) {
  
  for (j in (i + 1):length(groups)) {
    
    group1 <- groups[i]
    group2 <- groups[j]
    
    # Keep only samples from these two TreatTime groups
    keep <- ofra$TreatTime %in% c(group1, group2)
    
    # Subset distance matrix
    d_pair <- as.dist(
      as.matrix(dist_eu)[keep, keep]
    )
    
    # Subset grouping variable
    group_pair <- droplevels(
      factor(Ofav$TreatTime[keep])
    )
    
    # Run PERMDISP
    bd <- betadisper(
      d_pair,
      group_pair
    )
    
    # Test for differences in dispersion
    test <- anova(bd)
    
    # Add results
    permdisp_results <- rbind(
      permdisp_results,
      data.frame(
        Group1 = group1,
        Group2 = group2,
        F = test$`F value`[1],
        p = test$`Pr(>F)`[1]
      )
    )
  }
}

### Adjust p-values for multiple comparisons
permdisp_results$p_adjusted <- p.adjust(
  permdisp_results$p,
  method = "BH"
)


### View results
permdisp_results

######looks the same as rlog- stick with rlog for analyes

#### VennDiagram of Genes ####
#MAKE SURE TO Re-RUN DESeq with the ~Species+treatment*time through to res before running this part of code!!!
### These are lines: 125, 134, and 142

#This works to make Venn but not super informative about subgroups
# VENN Diagram to include both up and down regulated genes in common for transplant
resTrans1 <- results(dds, contrast=c("Time","Pre-heat", "Post-heat"))
resTrans2 <- results(dds, contrast=c("Treatment","DTV","Control"))
resTrans3 <- results(dds, contrast=c("Species","O. faveolata","O. franksi"))


head(resTrans1)
head(resTrans2)
head(resTrans3)

pTrans1up=row.names(resTrans1[resTrans1$padj<0.1 & !is.na(resTrans1$padj) & resTrans1$log2FoldChange>0,])
length(pTrans1up) #322
pTrans1down=row.names(resTrans2[resTrans1$padj<0.1 & !is.na(resTrans1$padj) & resTrans1$log2FoldChange<0,])
length(pTrans1down) #170
pTrans2up=row.names(resTrans2[resTrans2$padj<0.1 & !is.na(resTrans2$padj) & resTrans2$log2FoldChange>0,])
length(pTrans2up) #0
pTrans2down=row.names(resTrans2[resTrans2$padj<0.1 & !is.na(resTrans2$padj) & resTrans2$log2FoldChange<0,])
length(pTrans2down) #0
pTrans3up=row.names(resTrans3[resTrans3$padj<0.1 & !is.na(resTrans3$padj) & resTrans3$log2FoldChange>0,])
length(pTrans3up) #1835
pTrans3down=row.names(resTrans3[resTrans3$padj<0.1 & !is.na(resTrans3$padj) & resTrans3$log2FoldChange<0,])
length(pTrans3down) #1720

pTrans1=row.names(resTrans1[resTrans1$padj<0.1 & !is.na(resTrans1$padj),])
pTrans2=row.names(resTrans2[resTrans2$padj<0.1 & !is.na(resTrans2$padj),])
pTrans3=row.names(resTrans3[resTrans3$padj<0.1 & !is.na(resTrans3$padj),])

#UP
pdegs1_up=union(pTrans1up,pTrans2up)
pdegs1_up=union(pdegs1_up,pTrans3up)
length(pdegs1_up)
#2111

#DOWN
pdegs1_down=union(pTrans1down,pTrans2down)
pdegs1_down=union(pdegs1_down,pTrans3down)
length(pdegs1_down)
#1849

#ALL
pdegs1=union(pTrans1,pTrans2)
pdegs1=union(pdegs1,pTrans3)
length(pdegs1)
#3916

###do UP, DOWN, ALL
candidates=list("Time"=pTrans1up, "Treatment"=pTrans2up, "HostSpecies"=pTrans3up)
candidates=list("Time"=pTrans1down, "Treatment"=pTrans2down, "HostSpecies"=pTrans3down)

#### venndiagram
#devtools::install_github("nicolash2/ggvenn")

library("ggvenn")
quartz()
ggplot() + 
  geom_venn(candidates, textsize=8) +
  scale_fill_manual(values=c("grey44", "#2C2673","#E24A13")) +
  theme_void()+
  theme(legend.position = "none")

#Colleen Bove's venndiagramcode
# quartz()
# prettyvenn=venn.diagram(
#  x = candidates,
#  filename=NULL,
#  col = "transparent",
#  fill = c("coral2", "darkgreen", "royalblue1"),
#  alpha = 0.5,
#  #label.col = c("darkred", "white", "white", "white", "blue4"), #"white", "white", "white",
#  cex = 2.5,
#  fontfamily = "sans",
#  fontface = "bold",
#  cat.default.pos = "text",
#  cat.col = c("darkred", "darkgreen","blue4"),
#  cat.cex = 2.5,
#  cat.fontfamily = "sans",
#  cat.dist = c(0.08, 0.08, 0.03),
#  cat.pos = 1
# );
# 
# 
# grid.draw(prettyvenn)

### Fig. S7
#MAKE SURE TO Re-RUN DESeq with the ~Species+treatTime through to res before running this part of code!!!
### These are lines: 128, 134, and 142

#Venn with treatTime groupings, using DESeq ~Hostspecies+treatTime
# VENN Diagram to include both up and down regulated genes in common for transplant
resTrans1 <- results(dds, contrast=c("TreatTime","Control_Pre-heat", "Control_Post-heat"))
resTrans2 <- results(dds, contrast=c("TreatTime","Control_Pre-heat", "DTV_Pre-heat"))
resTrans3 <- results(dds, contrast=c("TreatTime","Control_Post-heat", "DTV_Post-heat"))
resTrans4 <- results(dds, contrast=c("TreatTime","DTV_Pre-heat","DTV_Post-heat"))
#resTrans5 <- results(dds, contrast=c("species","Ofav", "Ofra"))

head(resTrans1)
head(resTrans2)
head(resTrans3)
head(resTrans4)
#head(resTrans5)


pTrans1up=row.names(resTrans1[resTrans1$padj<0.1 & !is.na(resTrans1$padj) & resTrans1$log2FoldChange>0,])
length(pTrans1up) #322
pTrans1down=row.names(resTrans1[resTrans1$padj<0.1 & !is.na(resTrans1$padj) & resTrans1$log2FoldChange<0,])
length(pTrans1down) #170
pTrans2up=row.names(resTrans2[resTrans2$padj<0.1 & !is.na(resTrans2$padj) & resTrans2$log2FoldChange>0,])
length(pTrans2up) #1
pTrans2down=row.names(resTrans2[resTrans2$padj<0.1 & !is.na(resTrans2$padj) & resTrans2$log2FoldChange<0,])
length(pTrans2down) #0
pTrans3up=row.names(resTrans3[resTrans3$padj<0.1 & !is.na(resTrans3$padj) & resTrans3$log2FoldChange>0,])
length(pTrans3up) #0
pTrans3down=row.names(resTrans3[resTrans3$padj<0.1 & !is.na(resTrans3$padj) & resTrans3$log2FoldChange<0,])
length(pTrans3down) #0
pTrans4up=row.names(resTrans4[resTrans4$padj<0.1 & !is.na(resTrans4$padj) & resTrans4$log2FoldChange>0,])
length(pTrans4up) #216
pTrans4down=row.names(resTrans4[resTrans4$padj<0.1 & !is.na(resTrans4$padj) & resTrans4$log2FoldChange<0,])
length(pTrans4down) #125
# pTrans5up=row.names(resTrans5[resTrans5$padj<0.1 & !is.na(resTrans5$padj) & resTrans5$log2FoldChange>0,])
# length(pTrans5up) #1835
# pTrans5down=row.names(resTrans5[resTrans5$padj<0.1 & !is.na(resTrans5$padj) & resTrans5$log2FoldChange<0,])
# length(pTrans5down) #1720


pTrans1=row.names(resTrans1[resTrans1$padj<0.1 & !is.na(resTrans1$padj),])
pTrans2=row.names(resTrans2[resTrans2$padj<0.1 & !is.na(resTrans2$padj),])
pTrans3=row.names(resTrans3[resTrans3$padj<0.1 & !is.na(resTrans3$padj),])
pTrans4=row.names(resTrans4[resTrans4$padj<0.1 & !is.na(resTrans4$padj),])
#pTrans5=row.names(resTrans5[resTrans5$padj<0.1 & !is.na(resTrans5$padj),])


#UP
pdegs1_up=union(pTrans1up,pTrans2up)
pdegs1_up=union(pdegs1_up,pTrans3up)
pdegs1_up=union(pdegs1_up,pTrans4up)
#pdegs1_up=union(pdegs1_up,pTrans5up)
length(pdegs1_up)
#2182

#DOWN
pdegs1_down=union(pTrans1down,pTrans2down)
pdegs1_down=union(pdegs1_down,pTrans3down)
pdegs1_down=union(pdegs1_down,pTrans4down)
#pdegs1_down=union(pdegs1_down,pTrans5down)
length(pdegs1_down)
#1909

#ALL
pdegs1=union(pTrans1,pTrans2)
pdegs1=union(pdegs1,pTrans3)
pdegs1=union(pdegs1,pTrans4)
#pdegs1=union(pdegs1,pTrans5)
length(pdegs1)
#4021

###do UP, DOWN, ALL
candidatesUp=list("Preheat DTV-Control"=pTrans2up, "Control pre-post"=pTrans1up, "DTV pre-post"=pTrans4up, "Postheat DTV-Control"=pTrans3up)#, "Host - OfavOfra"=pTrans5up)
candidatesDown=list("Preheat DTV-Control"=pTrans2down, "Control pre-post"=pTrans1down, "DTV pre-post"=pTrans4down, "Postheat DTV-Control"=pTrans3down)#, "Host - OfavOfra"=pTrans5down)


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

###moved to WGCNA analysis from here######
