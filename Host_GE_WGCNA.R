##This code is adapted from Yingqi et al. 2023, Global Change Biology

setwd("~/path/Host_GE")

library(WGCNA)
library(flashClust)
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(diagram)
library(ggplot2)
library(cowplot)
library(dplyr)

# For PCA analysis
library("FactoMineR")
library("factoextra")

#### check data counts ####
CombCts=read.delim("~/path/Ofav_countsfixed_fixedSNP_noClones.txt",row.names=1,check.names=FALSE)
head(CombCts)
names(CombCts)=sub(".fastq.host.counts.txt","",names(CombCts))
names(CombCts)


#names(CombCts)=sub(".fastq.host.counts.txt","",names(CombCts))
names(CombCts)
ncol(CombCts) #32 samples in total
nrow(CombCts) #35826 genes in total

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)
allowWGCNAThreads()

#### Normalize reads, check for outliers, cluster samples by expression #################
# Ran rlogTransform.R from Yingqi Zhang et al. 2024

load(file = "rlogCountsGreedy.RData") #from Host_GE_DeSeqPCAVenn.R file
head(assay(rlog))

dat=as.data.frame(assay(rlog))

nrow(dat) # 16447 genes after removing isogroups with count less than 2 in more than 90% of samples; started with 35826 genes

# Reorganize data so that row corresponds to sample and column corresponds to gene
datExpr0 = as.data.frame(t(dat[,1:32]))

# Check for genes outliers and samples with too many missing values
gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK # TRUE

#import trait data and format dataframe for WGCNA
datTraits = read.csv("~/path/THERMVAR_MAIN_METADATA_Apr25_clonesfixed.csv",row.names = 1)

datTraits$Timepoint <- gsub("Preheat", "Pre-heat", datTraits$Time)
datTraits$Timepoint <- gsub("Postheat", "Post-heat", datTraits$Time)

datTraits$TreatTime <- gsub("Control_Preheat", "Control_Pre-heat", datTraits$TreatTime)
datTraits$TreatTime <- gsub("Control_Postheat", "Control_Post-heat", datTraits$TreatTime)
datTraits$TreatTime <- gsub("DTV_Preheat", "DTV_Pre-heat", datTraits$TreatTime)
datTraits$TreatTime <- gsub("DTV_Postheat", "DTV_Post-heat", datTraits$TreatTime)


#Remove rows from metadata that do not have GE
datTraits <- datTraits %>%
  filter(GE_sample != "")

#Remove rows from metadata that are clones
datTraits <- subset(datTraits, !(Genotype_Original %in% c("KA", "VB")))

#Convert Treatment, Timepoint, and TreatTime into binary columns
# Treatment
treatment_bin <- model.matrix(~ Treatment - 1, data = datTraits)
colnames(treatment_bin) <- gsub("Treatment", "", colnames(treatment_bin))

# Timepoint
time_bin <- model.matrix(~ Timepoint - 1, data = datTraits)
colnames(time_bin) <- gsub("Timepoint", "", colnames(time_bin))

# TreatTime
treattime_bin <- model.matrix(~ TreatTime - 1, data = datTraits)
colnames(treattime_bin) <- gsub("TreatTime", "", colnames(treattime_bin))

# Combine into one dataframe with same rownames as datTraits_filt
TraitsBinary_df <- data.frame(treatment_bin, time_bin, treattime_bin)
rownames(TraitsBinary_df) <- rownames(datTraits)

# Check result
head(TraitsBinary_df)

dim(TraitsBinary_df)
names(TraitsBinary_df)

table(rownames(TraitsBinary_df[1:32,])==rownames(datExpr0[1:32,])) # samples have matching names


# Cluster samples by expression
A=adjacency(t(datExpr0),type="signed") # SELECT SIGNED OR UNSIGNED HERE # this calculates the whole network connectivity
k = as.numeric(apply(A,2,sum))-1 # standardized connectivity
Z.k = scale(k)
thresholdZ.k = -2.5 # often -2.5
outlierColor = ifelse(Z.k<thresholdZ.k,"red","black")
sampleTree = flashClust(as.dist(1-A), method = "average")

#quartz() #use to make new windo with clustering tree
#par(cex = 0.6);
#par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)

# Convert traits to a color representation: red means high
traitColors = data.frame(numbers2colors(TraitsBinary_df,signed=FALSE))
dimnames(traitColors)[[2]] = paste(names(TraitsBinary_df))
datColors = data.frame(outlier = outlierColor,traitColors) # Combine outlier color info with trait info
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree,groupLabels=names(datColors),colors=datColors,main="Sample Dendrogram and Trait Heatmap")

# Remove outlying samples #no outliers but ran anyways
remove.samples = Z.k<thresholdZ.k | is.na(Z.k)
datExprOut = datExpr0[!remove.samples,]
datTraitsOut = datTraits[!remove.samples,]

dim(datExprOut)
dim(TraitsBinary_df)

datExprOutOfav= datExprOut[17:32,]
datExprOutOfra= datExprOut[1:16,]
datTraitsOutOfav = TraitsBinary_df[17:32,]
datTraitsOutOfra = TraitsBinary_df[1:16,]

save(datExprOut, datTraitsOut, datExprOutOfav, datExprOutOfra, datTraitsOutOfav, datTraitsOutOfra, file="SamplesAndTraits_OutliersRemoved.RData")

#######additional PCA analysis - not used for MS########

datTraitsOfav = datTraits[,2:8]
datTraitsOfra = datTraits[,2:8]


larOfav.pca = PCA(datExprOutOfav, graph = FALSE)
larOfra.pca = PCA(datExprOutOfra, graph = FALSE)

#Ofav
ind <- get_pca_ind(larOfav.pca)
PCcoordOfav = data.frame(ind$coord[,1:3])
#PCcoordOfav = data.frame(ind$coord[,2:3])
names(PCcoordOfav) = c("PC1", "PC2", "PC3")
#names(PCcoordOfav) = c("PC2", "PC3")
table(row.names(PCcoordOfav) == row.names(datTraitsOfav)) # TURE
PCcoordOfav = cbind(PCcoordOfav, datTraitsOfav)
#PCcoordOfav$Time = factor(PCcoordOfav$Time, levels = c("preheat", "postheat"))
#PCcoordOfav$Trmt = as.factor(PCcoordOfav$Trmt)
names(PCcoordOfav) = c("PC1", "PC2", "PC3","Actual_Genotype","Frag","Species","Tank","Treatment","Time","TreatTime")
#names(PCcoordOfav) = c("PC2", "PC3", "Time")

#Ofav PCA by Time
cbPalette <- c( "#B395DA","#4A2E61")#treatment
cbPalette <- c( "#D89362", "#8F3C00")#time
#pdf("PCA_Host_allgenes_rlog.pdf",height=5,width=6)
ggplot(PCcoordOfav, aes(PC1, PC2, color = Time, group=Time)) +
  geom_point(size = 4) +
  #geom_label(aes(label = Samples))+
  #geom_text_repel(aes(label=Samples)) +
  #scale_shape_manual(values = c(21, 24) )+
  scale_colour_manual(values=cbPalette)+
  #scale_fill_manual(values=c("#D89362", "#8F3C00"))+
  #guides(fill=guide_legend(override.aes=list(shape=21)))+
  #guides(color=guide_legend(override.aes=list(shape=21)))+
  #guides(fill=guide_legend(override.aes=list(shape=21,fill=c("#D89362", "#8F3C00"))))+
  theme_bw() +
  stat_ellipse()+
  # geom_density2d(alpha=.5)+
  # geom_polygon(alpha=.2)+
  xlab(paste0("PC1: ",round(larOfav.pca$eig[1,2], 1),"% variance")) +
  ylab(paste0("PC2: ",round(larOfav.pca$eig[2,2], 1),"% variance")) 
#dev.off()

#Ofra
ind <- get_pca_ind(larOfra.pca)
PCcoordOfra = data.frame(ind$coord[,1:3])
names(PCcoordOfra) = c("PC1", "PC2", "PC3")
table(row.names(PCcoordOfra) == row.names(datTraitsOfra)) # TURE
PCcoordOfra = cbind(PCcoordOfra, datTraitsOfra) # can subset like so datTraitsOutOfra_PCA[, 2:4])
#PCcoordOfra$Type = factor(PCcoordOfra$Type, levels = c("Offshore", "Cross1", "Cross2"))
#PCcoordOfra$Trmt = as.factor(PCcoordOfra$Trmt)
names(PCcoordOfra) = c("PC1", "PC2", "PC3","Actual_Genotype","Frag","Species","Tank","Treatment","Time","TreatTime")

#Ofra PCA by Time
cbPalette <- c( "#B395DA","#4A2E61")#treatment
cbPalette <- c( "#D89362", "#8F3C00")#time
#pdf("PCA_Host_allgenes_rlog.pdf",height=5,width=6)
ggplot(PCcoordOfra, aes(PC1, PC2, color = Time, group=Time)) +
  geom_point(size = 4) +
  #geom_label(aes(label = Samples))+
  #geom_text_repel(aes(label=Samples)) +
  #scale_shape_manual(values = c(21, 24) )+
  scale_colour_manual(values=cbPalette)+
  scale_fill_manual(values=c("#D89362", "#8F3C00"))+
  #guides(fill=guide_legend(override.aes=list(shape=21)))+
  #guides(color=guide_legend(override.aes=list(shape=21)))+
  #guides(fill=guide_legend(override.aes=list(shape=21,fill=c("#D89362", "#8F3C00"))))+
  theme_bw() +
  stat_ellipse()+
  # geom_density2d(alpha=.5)+
  # geom_polygon(alpha=.2)+
  xlab(paste0("PC1: ",round(larOfra.pca$eig[1,2], 1),"% variance")) +
  ylab(paste0("PC2: ",round(larOfra.pca$eig[2,2], 1),"% variance")) 
#dev.off()

#quartz()
#plot_grid(g1, g2, labels = c("a", "b"))

#### Both species - Choose a soft threshold power ####
powers = c(seq(1,30,by=1)) #may need to adjust these power values to hone in on proper sft value
sft = pickSoftThreshold(datExprOut, powerVector=powers, verbose =5, networkType="signed") #call network topology analysis function

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers, col="red")
abline(h=0.80, col="red")
# Choose power of 5

#### Both species -  TOM and GeneTree ####
# below is TOMstep.R from Yingqi Zhang et al. 2023
library(WGCNA)
options(stringsAsFactors = FALSE)
allowWGCNAThreads()
library(flashClust)

lnames = load(file="SamplesAndTraits_OutliersRemoved.RData") 

softPower=5
adjacencyOfav = adjacency(datExprOutOfav, power=softPower,type="signed")
# translate the adjacency into topological overlap matrix (TOM) and calculate the corresponding dissimilarity
TOMOfav= TOMsimilarity(adjacencyOfav,TOMType = "signed")
dissTOMOfav= 1-TOMOfav
# Generate a clustered gene tree
geneTreeOfav= flashClust(as.dist(dissTOMOfav), method="average")

adjacencyOfra = adjacency(datExprOutOfra, power=softPower,type="signed")
TOMOfra= TOMsimilarity(adjacencyOfra,TOMType = "signed")
dissTOMOfra= 1-TOMOfra
geneTreeOfra= flashClust(as.dist(dissTOMOfra), method="average")

save(dissTOMOfav, dissTOMOfra, geneTreeOfav, geneTreeOfra, file="TOM_Output.RData")

#### Both species - geneTree and dendrogram, set Min module size ####
### This is when I started using Jenna Dilworth et al. 2024 WGCNA code
### open this file instead if adjacency and TOM run earlier
load(file = "TOM_Output.RData")
load(file = "SamplesAndTraits_OutliersRemoved.RData")

###Cutting up the tree into highly co-expressed gene modules (by branches)
#Ofav
minModuleSize=35 #we only want large modules, this is considered relatively high
dynamicMods= cutreeDynamic(dendro= geneTreeOfav, distM= dissTOMOfav, deepSplit=0, pamRespectsDendro= FALSE, minClusterSize= minModuleSize)
table(dynamicMods) #lists the modules and how many genes are in each one

#Ofav
# dynamicMods
# 1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18   19   20 
# 1318  860  819  779  735  716  623  618  603  600  570  548  522  494  418  378  375  359  323  297 
# 21   22   23   24   25   26   27   28   29   30   31   32   33   34   35   36   37   38   39   40 
# 288  282  256  256  254  250  240  230  218  216  206  199  197  177  175  168  163  153  147  142 
# 41   42 
# 139  136 

#Ofra
minModuleSize=35 #we only want large modules, this is considered relatively high
dynamicMods= cutreeDynamic(dendro= geneTreeOfra, distM= dissTOMOfra, deepSplit=0, pamRespectsDendro= FALSE, minClusterSize= minModuleSize)
table(dynamicMods) #lists the modules and how many genes are in each one

#Ofra
# dynamicMods
# 1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18   19   20 
# 1590 1263  748  728  728  677  640  576  516  481  449  422  412  410  404  373  363  353  309  306 
# 21   22   23   24   25   26   27   28   29   30   31   32   33   34   35   36   37   38   39   40 
# 298  284  278  275  256  253  248  231  195  191  190  189  173  171  151  150  143  140  140  135 
# 41   42   43   44   45 
# 135  124  120  120  109 

##Plot module assignments under the gene tree
dynamicColors= labels2colors(dynamicMods)

#plot dendrogram and colors underneath --- comment in and out for whichever species you are running
#Ofav
#plotDendroAndColors(geneTreeOfav, dynamicColors, "Dynamic Tree Cut", dendroLabels= FALSE, hang=0.03, addGuide= TRUE, guideHang= 0.05, main= "Gene dendrogram and module colors sft=4, min Mod size=35")
#Ofra
plotDendroAndColors(geneTreeOfra, dynamicColors, "Dynamic Tree Cut", dendroLabels= FALSE, hang=0.03, addGuide= TRUE, guideHang= 0.05, main= "Gene dendrogram and module colors sft=5, min Mod size=35")

#Merge modules whose expression profiles are very similar or choose not to merge
#calculate eigengenes, I did this by species since their trees looked very different

### Make sure matches softpower chosen above!!! ###
softPower=5
#Ofav --- comment species in and out as needed
# MEList= moduleEigengenes(datExprOutOfav, colors= dynamicColors,softPower = softPower)
# MEs= MEList$eigengenes

#Ofra
MEList= moduleEigengenes(datExprOutOfra, colors= dynamicColors,softPower = softPower)
MEs= MEList$eigengenes

#Calculate dissimilarity of module eigengenes
MEDiss= 1-cor(MEs)

#Cluster module eigengenes
METree= flashClust(as.dist(MEDiss), method= "average")

#plot
plot(METree, main= "Clustering of module eigengenes", xlab= "", sub= "")

#start with a 95% similarity merge to get an initial sense of module trait correlations
#can take as high as 0.45 if needed to merge modules with 55%  similarity
#run each here and then one of the color abline lines below to see how they compare on same fig.
# MEDissThres= 0 #original 
# MEDissThres= 0.10 
# MEDissThres=0.15 
# MEDissThres=0.2 
# MEDissThres=0.3
# MEDissThres=0.35#tried this 17Dec2024, more significance but the same genes for preheat as 0.45
# MEDissThres=0.4
MEDissThres=0.45 #using this right now 16Dec2024 for Ofav and 9Jan2025 Ofra

##the lower the height cut, the more the modules have to be similar (more conservative)
#plot the cut line into tree
abline(h=MEDissThres, col="red")
abline(h=MEDissThres, col="blue")
abline(h=MEDissThres, col="green")

#Ofav --- comment each species in or out as needed for run
#merge= mergeCloseModules(datExprOutOfav, dynamicColors, cutHeight= MEDissThres, verbose =3)

#Ofra
merge= mergeCloseModules(datExprOutOfra, dynamicColors, cutHeight= MEDissThres, verbose =3)

#merge module colors and find new eigengenes on the new merged modules
mergedColors= merge$colors
mergedMEs= merge$newMEs
length(unique(mergedColors)) #shows how many modules there are now
#Softpowern = 5, DESeq =~1
#Ofav -> MEDissThres = 0.45 -> 14
#Ofra -> MEDissThres = 0.45 -> 14
#don't shift much so going to go with 0.45 module merging, since the trends observed are pretty stable.

##plot new module colors on gene tree under previous colors to see how they change
#Ofav --- comment species in and out as needed
#plotDendroAndColors(geneTreeOfav, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels= FALSE, hang=0.03, addGuide= TRUE, guideHang=0.05)
#Ofra
plotDendroAndColors(geneTreeOfra, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels= FALSE, hang=0.03, addGuide= TRUE, guideHang=0.05)

moduleColors= mergedColors
#create numerical lables corresponding to the colors
colorOrder= c("grey", standardColors(50))
moduleLabels= match(moduleColors, colorOrder)-1
MEs=mergedMEs #a dataframe of module Eigengenes for each module for each sample

#save module colors and labels for use in subsequent parts
#Ofav --- comment species in and out as needed
#save(MEs, moduleLabels, moduleColors, geneTreeOfav, file= "OfavNetwork_rlog_signed_merged0.45_sft5_31Mar2025.RData")

#Ofra
save(MEs, moduleLabels, moduleColors, geneTreeOfra, file= "OfraNetwork_rlog_signed_merged0.45_sft5_31Mar2025.RData")

#################Relating modules to traits and finding important genes




#### Species considered separately from here down in code
#### Can start here if you restarted R session --- heatmap and module exploration


##Start Trait Analysis here!!!
#### Load Trait Data and Add Physiology rows - run before Analysis ####
#import trait data and format dataframe for WGCNA, both coral host species
load(file = "SamplesAndTraits_OutliersRemoved.RData")
datTraits = read.csv("~/path/THERMVAR_MAIN_METADATA_Apr25_clonesfixed.csv",row.names = 1)

# Preserve rownames as an explicit column
datTraits$SampleID <- rownames(datTraits)
#rownames(datTraits) <- NULL

datTraits$Timepoint <- gsub("Preheat", "Pre-heat", datTraits$Time)
datTraits$Timepoint <- gsub("Postheat", "Post-heat", datTraits$Time)

datTraits$TreatTime <- gsub("Control_Preheat", "Control_Pre-heat", datTraits$TreatTime)
datTraits$TreatTime <- gsub("Control_Postheat", "Control_Post-heat", datTraits$TreatTime)
datTraits$TreatTime <- gsub("DTV_Preheat", "DTV_Pre-heat", datTraits$TreatTime)
datTraits$TreatTime <- gsub("DTV_Postheat", "DTV_Post-heat", datTraits$TreatTime)

#Remove rows from metadata that are clones
datTraits <- subset(datTraits, !(Genotype_Original %in% c("KA", "VB")))

#Make BuoyantWeight Calculations
## % change in weight
datTraits$bw_t1t0 <- ((datTraits$bw_avg_t1-datTraits$bw_avg_t0)/datTraits$bw_avg_t0)*100
datTraits$bw_t2t1 <- ((datTraits$bw_avg_t2-datTraits$bw_avg_t1)/datTraits$bw_avg_t1)*100
datTraits$bw_t3t2 <- ((datTraits$bw_avg_t3-datTraits$bw_avg_t2)/datTraits$bw_avg_t2)*100
datTraits$bw_t4t3 <- ((datTraits$bw_avg_t4-datTraits$bw_avg_t3)/datTraits$bw_avg_t3)*100

head(datTraits)


# Subset to species dataframes
OfavPhys_all <- datTraits %>% 
  filter(Species == "O. faveolata")

OfraPhys_all <- datTraits %>% 
  filter(Species == "O. franksi")

# Keep only the physiology variables of interest
OfavPhys_vars <- OfavPhys_all %>% 
  select(SampleID,Genotype_Original, Species, Treatment, Timepoint, GE_sample,
         bw_t3t2, bw_t4t3, pam_avg_t3, pam_avg_t4, 
         red_chan_avg_t3, red_chan_avg_t4, 
         gross_photo_t3, gross_photo_t4, 
         resp_t3, resp_t4, 
         PR_t3, PR_t4)

OfraPhys_vars <- OfraPhys_all %>% 
  select(SampleID, Genotype_Original, Species, Treatment, Timepoint, GE_sample,
         bw_t3t2, bw_t4t3, pam_avg_t3, pam_avg_t4, 
         red_chan_avg_t3, red_chan_avg_t4, 
         gross_photo_t3, gross_photo_t4, 
         resp_t3, resp_t4, 
         PR_t3, PR_t4)

# Collapse to unique combinations of identifiers, keeping first non-NA per group
#this is so we can use Respirometry data that was completed on frags that were not sampled for GE
OfavPhys_matched <-  OfavPhys_all %>%
  group_by(Genotype_Original, Species, Treatment, Timepoint) %>%
  summarise_at(
    vars(gross_photo_t3, gross_photo_t4, 
         resp_t3, resp_t4, 
         PR_t3, PR_t4),
    ~ if (all(is.na(.))) NA_real_ else dplyr::first(na.omit(.))
  ) %>%
  ungroup()

OfraPhys_matched <-  OfraPhys_all %>%
  group_by(Genotype_Original, Species, Treatment, Timepoint) %>%
  summarise_at(
    vars(gross_photo_t3, gross_photo_t4, 
         resp_t3, resp_t4, 
         PR_t3, PR_t4),
    ~ if (all(is.na(.))) NA_real_ else dplyr::first(na.omit(.))
  ) %>%
  ungroup()


# Now merge this back to your GE-sample-only dataframe
OfavPhys_GE <- OfavPhys_vars %>% 
  filter(GE_sample != "")

OfraPhys_GE <- OfraPhys_vars %>% 
  filter(GE_sample != "")

# Keep all original columns in OfavPhys_GE
# but only add certain summarized variables from OfavPhys_matched
OfavPhys_final <- OfavPhys_GE %>%
  left_join(
    OfavPhys_matched %>%
      select(Genotype_Original, Species, Treatment, Timepoint,
             # only include the summarised traits you want to add
             gross_photo_t3, gross_photo_t4,
             resp_t3, resp_t4, PR_t3, PR_t4),
    by = c("Genotype_Original", "Species", "Treatment", "Timepoint")
  )

OfraPhys_final <- OfraPhys_GE %>%
  left_join(
    OfraPhys_matched %>%
      select(Genotype_Original, Species, Treatment, Timepoint,
             # only include the summarised traits you want to add
             gross_photo_t3, gross_photo_t4,
             resp_t3, resp_t4, PR_t3, PR_t4),
    by = c("Genotype_Original", "Species", "Treatment", "Timepoint")
  )

#Append physiology data to binary traits for GE
dim(datTraitsOutOfav) # 16 samples, 8 traits
dim(OfavPhys_final) # 16 samples, 24 traits
datTraitsOutOfav <- bind_cols(datTraitsOutOfav,
  OfavPhys_final %>% select(bw_t3t2, bw_t4t3, pam_avg_t3, pam_avg_t4, 
                      red_chan_avg_t3, red_chan_avg_t4, 
                      gross_photo_t3.y, gross_photo_t4.y, 
                      resp_t3.y, resp_t4.y, 
                      PR_t3.y, PR_t4.y))

dim(datTraitsOutOfra) # 16 samples, 8 traits
dim(OfraPhys_final) # 16 samples, 24 traits
datTraitsOutOfra <- bind_cols(datTraitsOutOfra,
  OfraPhys_final %>% select(bw_t3t2, bw_t4t3, pam_avg_t3, pam_avg_t4, 
                      red_chan_avg_t3, red_chan_avg_t4, 
                      gross_photo_t3.y, gross_photo_t4.y, 
                      resp_t3.y, resp_t4.y, 
                      PR_t3.y, PR_t4.y))

#reorganize columns for heatmaps
# Desired trait order
traitOrder <- c("Control", "DTV", "Pre.heat", "Post.heat", "Control_Pre.heat","DTV_Pre.heat","Control_Post.heat","DTV_Post.heat", 
                "bw_t3t2","bw_t4t3", "pam_avg_t3", "pam_avg_t4", 
                "red_chan_avg_t3", "red_chan_avg_t4", "gross_photo_t3.y", "gross_photo_t4.y", 
                "resp_t3.y", "resp_t4.y", "PR_t3.y", "PR_t4.y")

# Reorder columns in datTraits
datTraitsOutOfav<- datTraitsOutOfav[, traitOrder]
datTraitsOutOfra<- datTraitsOutOfra[, traitOrder]


### Ofav time comparisons -- WGCNA ####
library(WGCNA)
library(stringr)
library(flashClust)
library(tidyverse)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);

# Load network data saved in the second part. Change this based on the level of merge you want to work with
lnames = load(file = "OfavNetwork_rlog_signed_merged0.45_sft5_31Mar2025.RData")

# class(datTraitsOutOfav) #binary representation of categorical traits
# 
# #Load metadata with continuous data and categorical names for module comparisons
# #load metadata if not already loaded
# datTraits= read.csv("~/Desktop/BU NSF Postdoc/Nicola_Ofav_chpt3/nic_DTV_metadata_fixed_noClones_binary.csv",row.names = 1)
# 
# #subset for Ofav only
# datTraitsOfav = subset(datTraits, Ofav == "1")
# 
# # Confirm the rows match
# stopifnot(all(rownames(datTraitsOfav) == rownames(datTraitsOutOfav))) #doesn't stop with error, so TRUE
# 
# #Categorical Traits from same frag as gene expression 
# datTraitsOutOfav <- datTraitsOutOfav %>%
#   mutate(Time = ifelse(preheat == 1, "Preheat", "Postheat"))
# 
# datTraitsOutOfav <- datTraitsOutOfav %>%
#   mutate(Treatment = ifelse(control == 1, "Control", "DTV"))


####################### Replot module dendrogram
### Make sure matches softpower chosen above!!! ###
softPower=5
MEList= moduleEigengenes(datExprOutOfav, moduleColors, softPower=softPower)$eigengenes

MEs<- MEList

#Calculate dissimilarity of module eigenegenes
MEDiss= 1-cor(MEs)
#Cluster module eigengenes
METree= flashClust(as.dist(MEDiss), method= "average")
plot(METree, main= "Clustering of module eigengenes", xlab= "", sub= "")
abline(h=.25,col="blue")
abline(h=.21,col="green")
#all modules are already merged above these cutoffs

#now for module trait heatmap
#correlate eigengenes with external clinical traits to look for most significant associations
nGenes = ncol(datExprOutOfav)
nSamples = nrow(datExprOutOfav)

# datTraits1=datTraitsOutOfav[,c(3:10)] #ALL
# head(datTraits1)

datTraits1=datTraitsOutOfav[,c(5:20)] #Excludes species, treatment, and time (includes only treatxtime for WGCNA)
head(datTraits1)

# Recalculate MEs with color labels (MEs=module eigengenes)
MEs0 = moduleEigengenes(datExprOutOfav, moduleColors,softPower=softPower)$eigengenes 
MEs = orderMEs(MEs0)

#correlations of traits with eigengenes

# for overall heatmap of treatment
moduleTraitCor2 = cor(MEs, datTraits1, use = "p"); #p=pearsons, #changes to MEs2 for just top 15 in heat map
moduleTraitPvalue = corPvalueStudent(moduleTraitCor2, nSamples);
Colors=sub("ME","",names(MEs)) #just takes off "MEs" in front of the color names

# for pretty fig of dominant symbiont
# moduleTraitCor = cor(MEs, datTraits2, use = "p"); #p=pearsons
# moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
# Colors=sub("ME","",names(MEs)) #just takes off "MEs" in front of the color names

#correlations of genes with eigengenes 
##used in next step of correlation process, will get done later in tutorial (can be down now or later)
moduleGeneCor=cor(MEs,datExprOutOfav)#run for treatment
moduleGenePvalue = corPvalueStudent(moduleGeneCor, nSamples)

# moduleGeneCor2=cor(MEs2,datExpr) #run for dominant Symbiont
# moduleGenePvalue = corPvalueStudent(moduleGeneCor2, nSamples)

#represent module trait correlations as a heatmap
# module-trait correlations
library(RColorBrewer)
modLabels=sub("ME","",names(MEs))

ps=signif(moduleTraitPvalue,1)
cors=signif(moduleTraitCor2,2)
textMatrix = ps;
#dev.off()
#displays only significant p values
textMatrix[ps>0.05]="-"
dim(textMatrix) = dim(moduleTraitCor2)#change based on what you are running above

# Will display correlations and their p-values
#textMatrix = paste(signif(moduleTraitCor2, 2),  "\n(",
#signif(moduleTraitPvalue, 1), ")", sep = "")

# Display the correlation values within a heatmap plot - all modules and all traits
#par( 1, 1);
#quartz()
par(font = 2)
labeledHeatmap(Matrix = moduleTraitCor2,
               xLabels = names(datTraits1),
               ySymbols = modLabels,
               yLabels = modLabels,
               colorLabels = FALSE,
               colors = colorRampPalette(c("dodgerblue","lightblue","white","coral","red"))(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.9,
               cex.lab = 0.9,
               zlim = c(-0.7,0.7))
#dev.off()
# module size barplot
labelShift=300 # increase to move module size labels to the right
#quartz()
#par(mar = c(6, 8.5, 3, 3));
mct=table(moduleColors)
mct[modLabels]
x=barplot(mct[rev(modLabels)],horiz=T,las=1,xlim=c(0,2500),col=rev(modLabels))
text(mct[rev(modLabels)]+labelShift,y=x,mct[rev(modLabels)],cex=0.9) 

mct

#selecting the 15 largest modules (all modules, cause only 15)
#Ofav - 0.45, soft5
MEs2 = MEs%>%select("MEblack", "MEbrown", "MEdarkmagenta", "MEdarkolivegreen", "MEdarkorange", 
                      "MEdarkturquoise", "MEgreenyellow", "MElightgreen", "MElightyellow",
                      "MEmagenta", "MEpurple", "MEred", "MEsaddlebrown", "MEturquoise") 

# MEs2 = MEs%>%select("MEblack", "MEblue", "MEbrown","MEbrown4","MEcyan", "MEdarkgreen", "MEdarkmagenta", "MEdarkolivegreen",
#                     "MEdarkred","MEdarkturquoise", "MEgreenyellow", "MEgrey60", "MElightcyan1", "MEmediumpurple3", "MEmidnightblue", "MEorangered4",
#                     "MEpink", "MEplum1", "MEpurple", "MEskyblue", "MEskyblue3", "MEsteelblue")
                    

####### creating a clustered complex heatmap - code from Chille et al. 2021 BMC Genomics
### https://github.com/echille/Mcapitata_Developmental_Gene_Expression_Timeseries/blob/v1.0.0/2a-WGCNA/Developmental_WGCNA.Rmd

#BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap) # for seeing module relatedness
#Create list of pvalues for eigengene correlation with specific life stages
heatmappval <- signif(moduleTraitPvalue, 1)
#Make list of heatmap row colors
htmap.colors <- names(MEs2)
htmap.colors <- gsub("ME", "", htmap.colors)

mod.sizes<-as.data.frame(mct)%>% #why filter here if it doesn't change anything in the heat map???
  filter(Freq>163)%>% #163 is the size of modules we are filtering by (for my dataset 164 is the size of the 13th largest module)
  #arrange(moduleColors=c("black", "brown4", "cyan", "darkgreen", "darkred",
  #                       "darkturquoise","green", "greenyellow", "lightcyan1", "midnightblue",
   #                      "orangered4", "plum1", "salmon"))
  arrange(moduleColors=c("black", "brown", "darkmagenta", "darkolivegreen", "darkorange", 
                         "darkturquoise", "greenyellow", "lightgreen", "lightyellow",
                         "magenta", "purple", "red", "saddlebrown", "turquoise"))


rownames(moduleTraitCor2) = gsub("ME", "", rownames(moduleTraitCor2))#all colors heatmap
#rownames(moduleTraitCor) = paste(rownames(moduleTraitCor), sep= " ", mod.sizes$Freq)

#quartz()
ht=Heatmap(moduleTraitCor2, name = "Corr.",
           col = blueWhiteRed(50), 
           row_names_side = "left", row_dend_side = "left",
           #right_annotation = size.annot,
           #width = unit(4, "in"), height = unit(8.5, "in"), 
           cluster_columns = FALSE,
           #cluster_rows = METree, row_split = 6, row_gap = unit(2.5, "mm"), border = TRUE,
           cell_fun = function(j, i, x, y, w, h, col) {
             if(heatmappval[i, j] <= 0.05) {
               grid.text(sprintf("%s", heatmappval[i, j]), x, y, gp = gpar(fontsize = 8, fontface = "bold"))
             }
             else {
               grid.text(sprintf("-"), x, y, gp = gpar(fontsize = 8, fontface = "plain"))
             }},
           column_names_gp =  gpar(fontsize = 10),
           column_names_rot = 45,
           row_names_gp = gpar(fontsize = 10, border = FALSE))
draw(ht)
#dev.off()
#map.grob =  grid.grabExpr(draw(ht)) 


#### Gene relationship to trait and important modules ####
#### Ofav - Time comparison by Categorical treatment ####
# Define variable weight containing the weight column of datTrait - leave weight as variable, but change names in first 2 commands

#for preheat v. postheat comparisons (for phys comparisons using Chla_ug.cm2 see further below)
weight = as.data.frame(datTraitsOutOfav$Pre.heat); 
#change to your trait name of interest: Combined, Control, Heat, pH, etc
names(weight) = "Pre.heat"

# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExprOutOfav, MEs, use = "p")); 
#finds pearson correlations of GE and module eigengenes
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

#making dataframe of pvalues for gene module membership values in each module
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExprOutOfav, weight, use = "p")); 

# data frame of correlations between expression and trait of interest using pearson correlations
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples)); 
#p-values for correlations between expression and trait of interest
names(geneTraitSignificance) = paste("GS.", names(weight), sep="");
names(GSPvalue) = paste("p.GS.", names(weight), sep="");
#geneTraitSignificance is a dataframe of correlation/covariance of GE data and stage value
#GSPvalue is a list of correlation of gTS pvalue

#need to adjust list below for my treatments and analyses###
#modules of interest for preheat
moduleCols=c("brown", "saddlebrown","darkturquoise","darkolivegreen", "lightyellow",
  "lightgreen", "magenta", "red", "lightyellow", "purple", "turquoise", "darkmagenta", "black", "greenyellow") #pick modules of interest from heatmap above


#plot scatter plots of gene significance vs module membership for all of these modules of interest
#add correlation and p-value, use this to look at how strong the modules are.
#quartz()
# par(mfrow=c(1,1))
# par(mar = c(2, 2, 2, 2));
# par(bg = "white");
for (module in moduleCols) {
  column = match(module, modNames);
  moduleGenes = moduleColors==module;
  verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                     abs(geneTraitSignificance[moduleGenes, 1]),
                     xlab = paste("ModMem", module),
                     ylab = "Gene Sig for Heat",
                     main = paste("MM vs. GS\n"),
                     cex.main = 1, cex.lab = 1, cex.axis = 1.2, col = module)
} 

#All modules for treatments checked with code above have positive correlations for 0.45 threshold


#### Visualization of Gene Networks --- Ofav preheat to postheat comparisons ####

#heatmaps of module expression with bar plot of eigengene for all samples
#use this to look at the different samples and make sure the modules make sense
#we want blocks by treatment, no single sample driving the differences.

#start with modules that looked strongest

#this section is computationally heavy - will abort abruptly. So ran in SCC 
# Add time from external metadata
MEs$Time <- datTraitsOutOfav$Time[match(rownames(MEs), rownames(datTraitsOutOfav))]

MEs.sorted <- MEs %>%
  arrange(Time) %>%
  dplyr::select(-Time)

which.module="brown" #pick module of interest
ME=MEs.sorted[, paste("ME",which.module, sep="")]
genes=datExprOutOfav[,moduleColors==which.module ]

#sort genes by Time
# Add time from external metadata
genes$Time <- datTraitsOutOfav$Time[match(rownames(genes), rownames(datTraitsOutOfav))]

genes.sorted <- genes %>%
  arrange(Time) %>%
  dplyr::select(-Time)

#quartz()
#par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
#plot expression  in all to make sure one sample isn't driving significance
plotMat(t(scale(genes.sorted) ),nrgcols=30,rlabels=F, clabels=rownames(genes.sorted), rcols=which.module,)

#par(mar=c(5, 4.2, 0, 0.7))
barplot(ME, col=which.module, main="", cex.main=2,
        ylab="eigengene expression",xlab="sample")


##########To output ME by sample, to see expression by sample/treatment

which.module="saddlebrown" #pick module of interest
ME=MEs[, paste("ME",which.module, sep="")]
mesaddlebrown<-ME

which.module="darkturquoise" #pick module of interest
ME=MEs[, paste("ME",which.module, sep="")]
medarkturquoise<-ME

which.module="darkolivegreen"
ME=MEs[, paste("ME",which.module, sep="")]
medarkolivegreen<-ME

which.module="lightyellow"
ME=MEs[, paste("ME",which.module, sep="")]
melightyellow<-ME

#meout<-data.frame(cbind(rownames(datExprOutOfav), medarkturquoise, melightcyan1))

#write.csv(meout,"Ofav_MEbySample_modsofinterest_9Dec2024.csv",quote=F,row.names=F) 


# get MEs from mods of interest from previously created ME dataframe 
# sort by treatment for boxplots
MEs$sort=rownames(MEs)
MEs.sorted <- MEs%>%
  separate(sort, into=c("Time", "garbage"), sep= "-")%>%
  dplyr::select(-garbage)%>%
  arrange(Time)%>%
  select(MEsaddlebrown, ,MEdarkturquoise, MEbrown, MEdarkolivegreen, MElightyellow,
         Time)%>%
  mutate_at("Time", as.factor)

MEs.sorted$Time = c("preheat","preheat", "postheat", "postheat","preheat", "postheat",
                    "postheat", "preheat", "preheat", "preheat","postheat", "postheat", 
                    "postheat", "postheat", "preheat", "preheat")
MEs.sorted$Treatment = c("control", "therm_var", "control", "therm_var","therm_var", "control",
                               "therm_var", "control", "therm_var", "control","control","therm_var", 
                               "control", "therm_var", "therm_var", "control")


medians <- MEs.sorted%>%
  group_by(Time)%>%
  dplyr::summarise(brown = median(MEbrown, na.rm=TRUE), saddlebrown= median(MEsaddlebrown, na.rm=TRUE), 
                   darkturquoise = median(MEdarkturquoise, na.rm=TRUE), darkolivegreen = median(MEdarkolivegreen, na.rm=TRUE),
                   lightyellow = median(MElightyellow, na.rm=TRUE),
                   .groups="keep")

# medians <- MEs.sorted%>%
#   group_by(Symbiont)%>%
#   dplyr::summarise(plum2 = median(MEplum2, na.rm=TRUE), darkorange2 = median(MEdarkorange2, na.rm=TRUE),
#                    purple= median(MEpurple, na.rm=TRUE), brown= median(MEbrown, na.rm=TRUE), 
#                    lavenderblush3= median(MElavenderblush3, na.rm=TRUE),lightcyan= median(MElightcyan, na.rm=TRUE),
#                    .groups="keep")


brown<- ggplot(MEs.sorted, aes(x = Time, y = MEbrown))+
  geom_boxplot(fill= "brown")+
  theme_classic()#+
# theme(axis.title.y=element_blank())+
# theme(axis.title.x=element_blank())+
# scale_y_continuous(breaks=c(-0.1,0, 0.1))+
# geom_hline(yintercept =0.06039775, linetype = "dashed", color = "grey")+
# geom_segment(x= "Heat", xend = "Heat", y = 0.06039775, yend = -0.04583882,color = "#D22B2B",
#              arrow = arrow(length = unit(0.03, "npc"), ends = "last", type = "closed"))+
# geom_segment(x= "pH", xend = "pH", y = 0.06039775, yend =	0.07624769,color = "#6082B6",
#              arrow = arrow(length = unit(0.03, "npc"), ends = "last", type = "closed"))+  
# geom_segment(x= "Combined", xend = "Combined", y =0.06039775, yend = -0.04583882,color = "#D22B2B",
#              arrow = arrow(length = unit(0.03, "npc"), ends = "last", type="closed"))+
# geom_segment(x= "Combined", xend = "Combined", y = -0.04583882, yend = (-0.04583882+(0.07624769-0.06039775)),color = "#6082B6",
#              arrow = arrow(length = unit(0.03, "npc"), ends = "last", type="closed"))
brown

saddlebrown<- ggplot(MEs.sorted, aes(x = Time, y = MEsaddlebrown))+
  geom_boxplot(fill= "saddlebrown")+
  theme_classic()#+
# theme(axis.title.y=element_blank())+
# theme(axis.title.x=element_blank())+
# scale_y_continuous(breaks=c(-0.1,0, 0.1))+
# geom_hline(yintercept =0.06039775, linetype = "dashed", color = "grey")+
# geom_segment(x= "Heat", xend = "Heat", y = 0.06039775, yend = -0.04583882,color = "#D22B2B",
#              arrow = arrow(length = unit(0.03, "npc"), ends = "last", type = "closed"))+
# geom_segment(x= "pH", xend = "pH", y = 0.06039775, yend =	0.07624769,color = "#6082B6",
#              arrow = arrow(length = unit(0.03, "npc"), ends = "last", type = "closed"))+  
# geom_segment(x= "Combined", xend = "Combined", y =0.06039775, yend = -0.04583882,color = "#D22B2B",
#              arrow = arrow(length = unit(0.03, "npc"), ends = "last", type="closed"))+
# geom_segment(x= "Combined", xend = "Combined", y = -0.04583882, yend = (-0.04583882+(0.07624769-0.06039775)),color = "#6082B6",
#              arrow = arrow(length = unit(0.03, "npc"), ends = "last", type="closed"))
saddlebrown

darkturquoise <- ggplot(MEs.sorted, aes(x = Time, y = MEdarkturquoise))+
  geom_boxplot(fill= "darkturquoise")+
  theme_classic()#+
# theme(axis.title.y=element_blank())+
# theme(axis.title.x=element_blank())+
# scale_y_continuous(breaks=c(-0.1,0, 0.1))+
# geom_hline(yintercept =0.06039775, linetype = "dashed", color = "grey")+
# geom_segment(x= "Heat", xend = "Heat", y = 0.06039775, yend = -0.04583882,color = "#D22B2B",
#              arrow = arrow(length = unit(0.03, "npc"), ends = "last", type = "closed"))+
# geom_segment(x= "pH", xend = "pH", y = 0.06039775, yend =	0.07624769,color = "#6082B6",
#              arrow = arrow(length = unit(0.03, "npc"), ends = "last", type = "closed"))+  
# geom_segment(x= "Combined", xend = "Combined", y =0.06039775, yend = -0.04583882,color = "#D22B2B",
#              arrow = arrow(length = unit(0.03, "npc"), ends = "last", type="closed"))+
# geom_segment(x= "Combined", xend = "Combined", y = -0.04583882, yend = (-0.04583882+(0.07624769-0.06039775)),color = "#6082B6",
#              arrow = arrow(length = unit(0.03, "npc"), ends = "last", type="closed"))
darkturquoise

darkolivegreen <- ggplot(MEs.sorted, aes(x = Time, y = MEdarkolivegreen))+
  geom_boxplot(fill= "darkolivegreen")+
  theme_classic() #+
# theme(axis.title.y=element_blank())+
# theme(axis.title.x=element_blank()) +
# scale_y_continuous(breaks=c(-0.1,0, 0.1))+
# geom_hline(yintercept = -0.03785620, linetype = "dashed", color = "grey")+
# geom_segment(x= "Heat", xend = "Heat", y = -0.03785620, yend = 0.06115896,color = "grey",
#              arrow = arrow(length = unit(0.03, "npc"), ends = "last", type = "closed"))
#geom_segment(x= "pH", xend = "pH", y = -0.07594524, yend = -0.04762537,color = "#6082B6",
#             arrow = arrow(length = unit(0.03, "npc"), ends = "last", type = "closed"))
darkolivegreen

lightyellow <- ggplot(MEs.sorted, aes(x = Time, y = MElightyellow))+
  geom_boxplot(fill= "lightyellow")+
  theme_classic() #+
# theme(axis.title.y=element_blank())+
# theme(axis.title.x=element_blank()) +
# scale_y_continuous(breaks=c(-0.1,0, 0.1))+
# geom_hline(yintercept = -0.03785620, linetype = "dashed", color = "grey")+
# geom_segment(x= "Heat", xend = "Heat", y = -0.03785620, yend = 0.06115896,color = "grey",
#              arrow = arrow(length = unit(0.03, "npc"), ends = "last", type = "closed"))
#geom_segment(x= "pH", xend = "pH", y = -0.07594524, yend = -0.04762537,color = "#6082B6",
#             arrow = arrow(length = unit(0.03, "npc"), ends = "last", type = "closed"))
lightyellow

library(cowplot)
library(grid)
library(gridExtra)
#makes pretty labeled figure of all boxplots from above (run each for sym and treatment)
boxes <- plot_grid(brown, saddlebrown, darkturquoise, darkolivegreen, lightyellow, ncol = 3, align = "v", 
                   axis="b",rel_heights = c(1,1,1), labels = "AUTO")
boxes

#create common x and y labels
y.grob <- textGrob("Eigengene Expression", 
                   gp=gpar(fontsize=12), rot=90)

x.grob <- textGrob("Time", 
                   gp=gpar(fontsize=12))

#add to plot
allplots <-grid.arrange(arrangeGrob(boxes, left = y.grob, bottom = x.grob))

#combine into figure with the heatmap
fig1 <- plot_grid(map.grob, allplots, rel_widths = c(1,2), labels = c("A", ""))
fig1


#### plotting massive table of all information - module membership, genes, gene names, etc. ####
annot=read.table("~/path/Ofav_Made_iso2geneName_16Dec2024.tab",sep="\t",quote="")
iso2go=read.table("~/path/Ofav_Made_iso2go_9Dec2024.tab",sep="\t",quote="")

probes=colnames(datExprOutOfav)

probes2go=match(probes,iso2go$V1)
summary(probes2go) # GO that match with probes from our dataset
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  1    3335    6750    6673    9936   13318    7619 

probes2annot = match(probes,annot$V1)
summary(probes2annot) # 14058 genes annotated
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  1    3320    6716    6641    9890   13260    7639

datGS.Traits=data.frame(cor(datExprOutOfav,datTraits1,use="p"))
names(datGS.Traits)=paste("cor",names(datGS.Traits),sep=".")
datME=moduleEigengenes(datExprOutOfav,moduleColors)$eigengenes
datKME=signedKME(datExprOutOfav, datME, outputColumnName="MM.")
datOutput=data.frame(ProbeID=names(datExprOutOfav),annot[probes2annot,],moduleColors,datKME,datGS.Traits)
#write.table(datOutput,"Ofav_AnnotatedNetworkAnalysisResultsHost_time_rlog_signed_sft4_merge0.45_16Dec2024.csv",row.names=F,sep=",")


#### Prep of Ofav WGCNA comparisons for GO Analysis for Fisher's Exact test######
# Generate categorical data for GO analysis for preheat and Physiology modules
ModColors=c("brown", "saddlebrown","darkturquoise","darkolivegreen", "lightyellow",
            "lightgreen", "magenta", "red", "lightyellow", "purple", "turquoise", "darkmagenta", "black", "greenyellow")

# # #Categorical GeneOntology by Module (Fisher Exact Test) with 0 (absent) or 1 (present) for each module
# for (col in ModColors) {
#   tab=datOutput[,c(1,4)]
#   
#   tab$modulesColors=as.character(tab$moduleColors)
#   tab$moduleColors[tab$moduleColors!=col]<-0
#   tab$moduleColors[tab$moduleColors==col]<-1
#   tab$moduleColors=as.factor(tab$moduleColors) 
#   print(col)
#   print(summary(tab)) #do counts match table of module colors? YES!!!
#   print(head(tab))
#   
#   write.csv(tab,file=paste("GO_MM",col,"_categorical_phys.csv", sep=""),quote=F,row.names=F)
# }


#GeneOntology by Module (Fisher Exact Test) with 0 or KME value to retain information on strength of membership
for (col in ModColors) {
  tab <- datOutput[, c("ProbeID", "moduleColors", paste0("MM.", col))]
  
  # Replace KME values with 0 if gene is not in this module
  tab$KME_val <- ifelse(tab$moduleColors == col, tab[[paste0("MM.", col)]], 0)
  
  # Keep only the KME values of interest
  out <- tab[, c("ProbeID", "KME_val")]
  
  print(col)
  print(summary(out$KME_val))  # check distribution
  print(head(out))
  
  write.csv(out, 
            file = paste0("GO_MM_", col, "_KME.csv"), 
            quote = FALSE, row.names = FALSE)
}

###Now input these files into GO_MWU.R for analysis!!!#####


####(NOT USED) Ofav time+treatment comparisons - same results as just time, so only used time for other analyses ####
library(stringr)
library(flashClust)
library(tidyverse)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);

# Load the expression and trait data saved in the first part 
lnames = load(file="SamplesAndTraits_OutliersRemoved.RData")  
#The variable lnames contains the names of loaded variables

# Load network data saved in the second part. Change this based on the level of merge you want to work with
lnames = load(file = "OfavNetwork_rlog_signed_merged0.45_sft8_16Dec2024.RData")

class(datTraitsOutOfav)

#Traits from same frag as gene expression

datTraitsOutOfav$TimeTreat = c("control_preheat", "therm_var_preheat", "control_postheat", "therm_var_postheat","therm_var_preheat", "control_postheat",
                               "therm_var_postheat", "control_preheat", "therm_var_preheat", "control_preheat","control_postheat","therm_var_postheat", 
                               "control_postheat", "therm_var_postheat", "therm_var_preheat", "control_preheat")


#######################Replot module dendrogram
MEList= moduleEigengenes(datExprOutOfav, moduleColors, softPower=softPower)$eigengenes

MEs<- MEList

#Calculate dissimilarity of module eigenegenes
MEDiss= 1-cor(MEs)
#Cluster module eigengenes
METree= flashClust(as.dist(MEDiss), method= "average")
plot(METree, main= "Clustering of module eigengenes", xlab= "", sub= "")
abline(h=.25,col="blue")
abline(h=.21,col="green")


#now for module trait heatmap
#correlate eigengenes with external clinical traits to look for most significant associations
nGenes = ncol(datExprOutOfav)
nSamples = nrow(datExprOutOfav)

datTraits1=datTraitsOutOfav[,c(3:10)] #ALL
head(datTraits1)

# Recalculate MEs with color labels (MEs=module eigengenes)
MEs0 = moduleEigengenes(datExprOutOfav, moduleColors,softPower=softPower)$eigengenes 
MEs = orderMEs(MEs0)

#correlations of traits with eigengenes

# for overall heatmap of treatment
moduleTraitCor2 = cor(MEs, datTraits1, use = "p"); #p=pearsons, #changes to MEs2 for just top 15 in heat map
moduleTraitPvalue = corPvalueStudent(moduleTraitCor2, nSamples);
Colors=sub("ME","",names(MEs)) #just takes off "MEs" in front of the color names

# for pretty fig of dominant symbiont
# moduleTraitCor = cor(MEs, datTraits2, use = "p"); #p=pearsons
# moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
# Colors=sub("ME","",names(MEs)) #just takes off "MEs" in front of the color names

#correlations of genes with eigengenes 
##used in next step of correlation process, will get done later in tutorial (can be down now or later)
moduleGeneCor=cor(MEs,datExprOutOfav)#run for treatment
moduleGenePvalue = corPvalueStudent(moduleGeneCor, nSamples)

# moduleGeneCor2=cor(MEs2,datExpr) #run for dominant Symbiont
# moduleGenePvalue = corPvalueStudent(moduleGeneCor2, nSamples)

#represent module trait correlations as a heatmap
# module-trait correlations
library(RColorBrewer)
modLabels=sub("ME","",names(MEs))

ps=signif(moduleTraitPvalue,1)
cors=signif(moduleTraitCor2,2)
textMatrix = ps;
#dev.off()
#displays only significant p values
textMatrix[ps>0.05]="-"
dim(textMatrix) = dim(moduleTraitCor2)#change based on what you are running above

# Will display correlations and their p-values
#textMatrix = paste(signif(moduleTraitCor2, 2),  "\n(",
#signif(moduleTraitPvalue, 1), ")", sep = "")

# Display the correlation values within a heatmap plot - all modules and all traits
#par( 1, 1);
#quartz()
labeledHeatmap(Matrix = moduleTraitCor2,
               xLabels = names(datTraits1),
               ySymbols = modLabels,
               yLabels = modLabels,
               colorLabels = FALSE,
               colors = colorRampPalette(c("blue","lightblue","white","coral","red"))(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.7,
               cex.lab = 0.7,
               zlim = c(-0.7,0.7))
#dev.off()
# module size barplot
labelShift=300 # increase to move module size labels to the right
#quartz()
#par(mar = c(6, 8.5, 3, 3));
mct=table(moduleColors)
mct[modLabels]
x=barplot(mct[rev(modLabels)],horiz=T,las=1,xlim=c(0,2500),col=rev(modLabels))
text(mct[rev(modLabels)]+labelShift,y=x,mct[rev(modLabels)],cex=0.9) 

mct
#selecting the 15 largest modules
#0.45 MEDissthreshold
MEs2 = MEs%>%select("MEblack", "MEbrown4", "MEcyan", "MEdarkgreen", "MEdarkred",
                    "MEdarkturquoise", "MEgreen", "MEgreenyellow", "MElightcyan1",
                    "MEmidnightblue", "MEorangered4", "MEplum1", "MEsalmon")
#0.35 MEDissThreshold
#MEs2 = MEs%>%select("MEblack", "MEblue", "MEbrown","MEbrown4","MEcyan", "MEdarkgreen", "MEdarkmagenta", "MEdarkolivegreen",
#                    "MEdarkred","MEdarkturquoise", "MEgreenyellow", "MEgrey60", "MElightcyan1", "MEmediumpurple3", "MEmidnightblue", "MEorangered4",
#                    "MEpink", "MEplum1", "MEpurple", "MEskyblue", "MEskyblue3", "MEsteelblue")


### creating a clustered complex heatmap - code from Chille et al. 2021 BMC Genomics
### https://github.com/echille/Mcapitata_Developmental_Gene_Expression_Timeseries/blob/v1.0.0/2a-WGCNA/Developmental_WGCNA.Rmd

#BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
#Create list of pvalues for eigengene correlation with specific life stages
heatmappval <- signif(moduleTraitPvalue, 1)
#Make list of heatmap row colors
htmap.colors <- names(MEs2)
htmap.colors <- gsub("ME", "", htmap.colors)

mod.sizes<-as.data.frame(mct)%>% #why filter here if it doesn't change anything in the heat map???
  filter(Freq>163)%>% #163 is the size of modules we are filtering by (for my dataset 164 is the size of the 13th largest module)
  #arrange(moduleColors=c("black", "brown4", "cyan", "darkgreen", "darkred",
  #                       "darkturquoise","green", "greenyellow", "lightcyan1", "midnightblue",
  #                      "orangered4", "plum1", "salmon"))
  arrange(moduleColors=c("black", "brown4", "cyan", "darkgreen", "darkred",
                         "darkturquoise", "green", "greenyellow", "lightcyan1",
                         "midnightblue", "orangered4", "plum1", "salmon"))


#0.35 MEDissThreshold
  # arrange(moduleColors=c("black", "blue", "brown","brown4","cyan", "darkgreen", "darkmagenta", "darkolivegreen",
  #                        "darkred","darkturquoise", "greenyellow", "grey60", "lightcyan1", "mediumpurple3", "midnightblue"
  #                        , "orangered4","pink", "plum1", "purple", "skyblue", "skyblue3", "steelblue"))


rownames(moduleTraitCor2) = gsub("ME", "", rownames(moduleTraitCor2))#all colors heatmap
#rownames(moduleTraitCor) = paste(rownames(moduleTraitCor), sep= " ", mod.sizes$Freq)

#quartz()
ht=Heatmap(moduleTraitCor2, name = "Corr.",
           col = blueWhiteRed(50), 
           row_names_side = "left", row_dend_side = "left",
           #right_annotation = size.annot,
           #width = unit(4, "in"), height = unit(8.5, "in"), 
           cluster_columns = FALSE,
           #cluster_rows = METree, row_split = 6, row_gap = unit(2.5, "mm"), border = TRUE,
           cell_fun = function(j, i, x, y, w, h, col) {
             if(heatmappval[i, j] <= 0.05) {
               grid.text(sprintf("%s", heatmappval[i, j]), x, y, gp = gpar(fontsize = 8, fontface = "bold"))
             }
             else {
               grid.text(sprintf("-"), x, y, gp = gpar(fontsize = 8, fontface = "plain"))
             }},
           column_names_gp =  gpar(fontsize = 10),
           column_names_rot = 45,
           row_names_gp = gpar(fontsize = 10, border = FALSE))
draw(ht)

map.grob =  grid.grabExpr(draw(ht)) 

#Gene relationship to trait and important modules:
#############Treatment
# Define variable weight containing the weight column of datTrait - leave weight as variable, but change names in first 2 commands
weight = as.data.frame(datTraitsOutOfav$therm_var_preheat); 
#change to your trait name of interest: Combined, Control, Heat, pH, etc
names(weight) = "therm_var_preheat"

#run separately for other variables after running lines after
# weight = as.data.frame(datTraitsOutOfav$postheat); 
# names(weight) = "postheat"


# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExprOutOfav, MEs, use = "p")); 
#finds pearson correlations of GE and module eigengenes
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

#making dataframe of pvalues for gene module membership values in each module
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExprOutOfav, weight, use = "p")); 

# data frame of correlations between expression and trait of interest using pearson correlations
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples)); 
#p-values for correlations between expression and trait of interest
names(geneTraitSignificance) = paste("GS.", names(weight), sep="");
names(GSPvalue) = paste("p.GS.", names(weight), sep="");
#geneTraitSignificance is a dataframe of correlation/covariance of GE data and stage value
#GSPvalue is a list of correlation of gTS pvalue

#need to adjust list below for my treatments and analyses###
#modules of interest for preheat
#moduleCols=c("darkturquoise", "darkred", "lightcyan1")#0.45 threshold
#moduleCols=c("steelblue", "blue", "darkred", "lightcyan1", "darkturquoise","darkolivegreen")#0.35 threshold
moduleCols=c("lightcyan1", "darkred")

#plot scatter plots of gene significance vs module membership for all of these modules of interest
#add correlation and p-value, use this to look at how strong the modules are.
#quartz()
# par(mfrow=c(1,1))
# par(mar = c(2, 2, 2, 2));
# par(bg = "white");
for (module in moduleCols) {
  column = match(module, modNames);
  moduleGenes = moduleColors==module;
  verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                     abs(geneTraitSignificance[moduleGenes, 1]),
                     xlab = paste("ModMem", module),
                     ylab = "Gene Sig for Heat",
                     main = paste("MM vs. GS\n"),
                     cex.main = 1, cex.lab = 1, cex.axis = 1.2, col = module)
} 

#All modules for treatments checked with code above have positive correlations


###Visualization of Gene Networks


#heatmaps of module expression with bar plot of eigengene for all samples
#use this to look at the different samples and make sure the modules make sense
#we want blocks by treatment, no single sample driving the differences.

#start with modules that looked strongest: "green","purple","cyan", "tan"

#sort ME by Time  ---- there is an issue here, causes session to abort, when too much data, if keeps crashiing try on SCC--
MEs$sort=rownames(MEs)
MEs.sorted <- MEs%>%
  separate(sort, into=c("garbage", "TimeTreat", "garbage2"), sep= "-")%>%
  dplyr::select(-garbage, -garbage2)%>%
  arrange(TimeTreat)%>%
  dplyr::select(-TimeTreat)

which.module="darkred" #pick module of interest
ME=MEs.sorted[, paste("ME",which.module, sep="")]
genes=datExprOutOfav[,moduleColors==which.module ]

#sort genes by treatment
genes$sort=rownames(genes)
genes.sorted <- genes%>%
  separate(sort, into=c("garbage", "TimeTreat", "garbage2"), sep= "-")%>%
  dplyr::select(-garbage, -garbage2)%>%
  arrange(TimeTreat)%>%
  dplyr::select(-TimeTreat)

#quartz()
#par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
#plot expression  in all to make sure one sample isn't driving significance
plotMat(t(scale(genes.sorted) ),nrgcols=30,rlabels=F, clabels=rownames(genes.sorted), rcols=which.module,)

#par(mar=c(5, 4.2, 0, 0.7))
barplot(ME, col=which.module, main="", cex.main=2,
        ylab="eigengene expression",xlab="sample")


##########To output ME by sample, to see expression by sample/treatment

which.module="lightcyan1" #pick module of interest
ME=MEs[, paste("ME",which.module, sep="")]
melightcyan1<-ME

meout<-data.frame(cbind(rownames(datExprOutOfav), mesteelblue))

write.csv(meout,"Ofav_MEbySample_modsofinterest_0.45TimeTreat_20Dec2024.csv",quote=F,row.names=F) 



# get MEs from mods of interest from previously created ME dataframe 
# sort by treatment for boxplots
MEs$sort=rownames(MEs)
MEs.sorted <- MEs%>%
  separate(sort, into=c("TimeTreat", "garbage"), sep= "-")%>%
  dplyr::select(-garbage)%>%
  arrange(TimeTreat)%>%
  select(MEdarkred, MElightcyan1,
         TimeTreat)%>%
  mutate_at("TimeTreat", as.factor)

MEs.sorted$Time = c("preheat","preheat", "postheat", "postheat","preheat", "postheat",
                    "postheat", "preheat", "preheat", "preheat","postheat", "postheat", 
                    "postheat", "postheat", "preheat", "preheat")
MEs.sorted$Treatment = c("control", "therm_var", "control", "therm_var","therm_var", "control",
                         "therm_var", "control", "therm_var", "control","control","therm_var", 
                         "control", "therm_var", "therm_var", "control")

MEs.sorted$TimeTreat = c("control_preheat", "therm_var_preheat", "control_postheat", "therm_var_postheat","therm_var_preheat", "control_postheat",
                               "therm_var_postheat", "control_preheat", "therm_var_preheat", "control_preheat","control_postheat","therm_var_postheat", 
                               "control_postheat", "therm_var_postheat", "therm_var_preheat", "control_preheat")



medians <- MEs.sorted%>%
  group_by(TimeTreat)%>%
  dplyr::summarise(darkred = median(MEdarkred, na.rm=TRUE), lightcyan1= median(MElightcyan1, na.rm=TRUE), 
                   .groups="keep")


# medians <- MEs.sorted%>%
#   group_by(Symbiont)%>%
#   dplyr::summarise(plum2 = median(MEplum2, na.rm=TRUE), darkorange2 = median(MEdarkorange2, na.rm=TRUE),
#                    purple= median(MEpurple, na.rm=TRUE), brown= median(MEbrown, na.rm=TRUE), 
#                    lavenderblush3= median(MElavenderblush3, na.rm=TRUE),lightcyan= median(MElightcyan, na.rm=TRUE),
#                    .groups="keep")


lightcyan1 <- ggplot(MEs.sorted, aes(x = TimeTreat, y = MElightcyan1))+
  geom_boxplot(fill= "lightcyan1")+
  theme_classic()#+
# theme(axis.title.y=element_blank())+
# theme(axis.title.x=element_blank())+
# scale_y_continuous(breaks=c(-0.1,0, 0.1))+
# geom_hline(yintercept =0.06039775, linetype = "dashed", color = "grey")+
# geom_segment(x= "Heat", xend = "Heat", y = 0.06039775, yend = -0.04583882,color = "#D22B2B",
#              arrow = arrow(length = unit(0.03, "npc"), ends = "last", type = "closed"))+
# geom_segment(x= "pH", xend = "pH", y = 0.06039775, yend =	0.07624769,color = "#6082B6",
#              arrow = arrow(length = unit(0.03, "npc"), ends = "last", type = "closed"))+  
# geom_segment(x= "Combined", xend = "Combined", y =0.06039775, yend = -0.04583882,color = "#D22B2B",
#              arrow = arrow(length = unit(0.03, "npc"), ends = "last", type="closed"))+
# geom_segment(x= "Combined", xend = "Combined", y = -0.04583882, yend = (-0.04583882+(0.07624769-0.06039775)),color = "#6082B6",
#              arrow = arrow(length = unit(0.03, "npc"), ends = "last", type="closed"))
lightcyan1


darkred <- ggplot(MEs.sorted, aes(x = TimeTreat, y = MEdarkred))+
  geom_boxplot(fill= "darkred")+
  theme_classic() #+
# theme(axis.title.y=element_blank())+
# theme(axis.title.x=element_blank()) +
# scale_y_continuous(breaks=c(-0.1,0, 0.1))+
# geom_hline(yintercept = -0.03785620, linetype = "dashed", color = "grey")+
# geom_segment(x= "Heat", xend = "Heat", y = -0.03785620, yend = 0.06115896,color = "grey",
#              arrow = arrow(length = unit(0.03, "npc"), ends = "last", type = "closed"))
#geom_segment(x= "pH", xend = "pH", y = -0.07594524, yend = -0.04762537,color = "#6082B6",
#             arrow = arrow(length = unit(0.03, "npc"), ends = "last", type = "closed"))
darkred

library(cowplot)
library(grid)
library(gridExtra)
#makes pretty labeled figure of all boxplots from above (run each for sym and treatment)
boxes <- plot_grid(darkred, lightcyan1, ncol = 2, align = "v", 
                   axis="b",rel_heights = c(1,1), labels = "AUTO")
boxes

#create common x and y labels
y.grob <- textGrob("Eigengene Expression", 
                   gp=gpar(fontsize=12), rot=90)

x.grob <- textGrob("Time", 
                   gp=gpar(fontsize=12))

#add to plot
allplots <-grid.arrange(arrangeGrob(boxes, left = y.grob, bottom = x.grob))

#combine into figure with the heatmap
fig1 <- plot_grid(map.grob, allplots, rel_widths = c(1,2), labels = c("A", ""))
fig1


###########plotting massive table of all information - module membership, genes, gene names, etc. ####
annot=read.table("~/path/Ofav_Made_iso2geneName_16Dec2024.tab",sep="\t",quote="")
iso2go=read.table("~/path/Ofav_Made_iso2go_9Dec2024.tab",sep="\t",quote="")

probes=colnames(datExprOutOfav)

probes2go=match(probes,iso2go$V1)
summary(probes2go) # ? genes annotated <-how is this determined?
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 1    3335    6750    6673    9936   13318    7619 

probes2annot = match(probes,annot$V1)
summary(probes2annot) # 14058 genes annotated <-how is this determined?
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  1    3320    6716    6641    9890   13260    7639 

datGS.Traits=data.frame(cor(datExprOutOfav,datTraits1,use="p"))
names(datGS.Traits)=paste("cor",names(datGS.Traits),sep=".")
datME=moduleEigengenes(datExprOutOfav,moduleColors)$eigengenes
datKME=signedKME(datExprOutOfav, datME, outputColumnName="MM.")
datOutput=data.frame(ProbeID=names(datExprOutOfav),annot[probes2annot,],moduleColors,datKME,datGS.Traits)
#write.table(datOutput,"Ofav_AnnotatedNetworkAnalysisResultsHost_time_rlog_signed_sft4_merge0.45_16Dec2024.csv",row.names=F,sep=",")


####Prep time+treatment for GO Analysis Categorical for Fisher's Exact test######
# Generate categorical data for GO analysis 
ModColors=c( "darkred", "lightcyan1")

for (col in ModColors) {
  tab=datOutput[,c(1,4)]
  
  tab$modulesColors=as.character(tab$moduleColors)
  #Categorical GeneOntology by Module (Fisher Exact Test)
  tab$moduleColors[tab$moduleColors!=col]<-0
  tab$moduleColors[tab$moduleColors==col]<-1 
  tab$moduleColors=as.factor(tab$moduleColors) 
  print(col)
  print(summary(tab)) #do counts match table of module colors? YES!!!
  print(head(tab))
  
  write.csv(tab,file=paste("GO_MM",col,"_categorical.csv", sep=""),quote=F,row.names=F)
}

###Now input these files into GO_MWU.R for analysis!!!#####



#### Ofra Time comparisons - WGCNA ####
library(WGCNA)
library(stringr)
library(flashClust)
library(tidyverse)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);


# Load network data saved in the second part. Change this based on the level of merge you want to work with
lnames = load(file = "OfraNetwork_rlog_signed_merged0.45_sft5_31Mar2025.RData")

class(datTraitsOutOfra)

#Load metadata with continuous data and categorical names for module comparisons
# 
# #subset for Ofra only
# datTraitsOfra = subset(datTraits, Ofra == "1")
# 
# # Confirm the rows match
# stopifnot(all(rownames(datTraitsOfra) == rownames(datTraitsOutOfra))) 
# #doesn't stop with error, so TRUE
# 
# 
# #Categorical Traits from same frag as gene expression 
# datTraitsOutOfra <- datTraitsOutOfra %>%
#   mutate(Time = ifelse(preheat == 1, "Preheat", "Postheat"))

# datTraitsOutOfra <- datTraitsOutOfra %>%
#   mutate(Treatment = ifelse(control == 1, "Control", "DTV"))


#######################Replot module dendrogram
softPower = 5
MEList= moduleEigengenes(datExprOutOfra, moduleColors, softPower=softPower)$eigengenes

MEs<- MEList

#Calculate dissimilarity of module eigenegenes
MEDiss= 1-cor(MEs)
#Cluster module eigengenes
METree= flashClust(as.dist(MEDiss), method= "average")
plot(METree, main= "Clustering of module eigengenes", xlab= "", sub= "")
abline(h=.25,col="blue")
abline(h=.21,col="green")


#now for module trait heatmap
#correlate eigengenes with external clinical traits to look for most significant associations
nGenes = ncol(datExprOutOfra)
nSamples = nrow(datExprOutOfra)

datTraits1=datTraitsOutOfra[,c(5:20)] #ALL
head(datTraits1)

# Recalculate MEs with color labels (MEs=module eigengenes)
MEs0 = moduleEigengenes(datExprOutOfra, moduleColors,softPower=softPower)$eigengenes 
MEs = orderMEs(MEs0)

#correlations of traits with eigengenes

# for overall heatmap of treatment
moduleTraitCor2 = cor(MEs, datTraits1, use = "p"); #p=pearsons, #changes to MEs2 for just top 15 in heat map
moduleTraitPvalue = corPvalueStudent(moduleTraitCor2, nSamples);
Colors=sub("ME","",names(MEs)) #just takes off "MEs" in front of the color names

#correlations of genes with eigengenes 
##used in next step of correlation process, will get done later in tutorial (can be down now or later)
moduleGeneCor=cor(MEs,datExprOutOfra)#run for treatment
moduleGenePvalue = corPvalueStudent(moduleGeneCor, nSamples)

#represent module trait correlations as a heatmap
# module-trait correlations
library(RColorBrewer)
modLabels=sub("ME","",names(MEs))

ps=signif(moduleTraitPvalue,1)
cors=signif(moduleTraitCor2,2)
textMatrix = ps;
#dev.off()
#displays only significant p values
textMatrix[ps>0.05]="-"
dim(textMatrix) = dim(moduleTraitCor2)#change based on what you are running above

# Will display correlations and their p-values
#textMatrix = paste(signif(moduleTraitCor2, 2),  "\n(",
#signif(moduleTraitPvalue, 1), ")", sep = "")

# Display the correlation values within a heatmap plot - all modules and all traits
#par( 1, 1);
#quartz()
par(font = 2)
labeledHeatmap(Matrix = moduleTraitCor2,
               xLabels = names(datTraits1),
               ySymbols = modLabels,
               yLabels = modLabels,
               colorLabels = FALSE,
               colors = colorRampPalette(c("dodgerblue","lightblue","white","coral","red"))(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.9,
               cex.lab = 0.9,
               zlim = c(-0.7,0.7))
#dev.off()
# module size barplot
labelShift=300 # increase to move module size labels to the right
#quartz()
#par(mar = c(6, 8.5, 3, 3));
mct=table(moduleColors)
mct[modLabels]
x=barplot(mct[rev(modLabels)],horiz=T,las=1,xlim=c(0,4000),col=rev(modLabels))
text(mct[rev(modLabels)]+labelShift,y=x,mct[rev(modLabels)],cex=0.9) 

mct
#selecting the 15 largest modules
#MEs2 = MEs%>%select("MEblack", "MEbrown4", "MEcyan", "MEdarkgreen", "MEdarkred",
#                    "MEdarkturquoise","MEgreen", "MEgreenyellow", "MElightcyan1", "MEmidnightblue",
#                    "MEorangered4", "MEplum1", "MEsalmon")

MEs2 = MEs%>%select("MEblack", "MEblue", "MEbrown", "MEcyan", "MEdarkgrey","MEdarkmagenta", "MEdarkorange", "MEdarkred", "MEivory",
                    "MElightcyan", "MElightcyan1", "MEpaleturquoise", "MEpink", "MEsienna3")


### creating a clustered complex heatmap - code from Chille et al. 2021 BMC Genomics
### https://github.com/echille/Mcapitata_Developmental_Gene_Expression_Timeseries/blob/v1.0.0/2a-WGCNA/Developmental_WGCNA.Rmd

#BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
#Create list of pvalues for eigengene correlation with specific life stages
heatmappval <- signif(moduleTraitPvalue, 1)
#Make list of heatmap row colors
htmap.colors <- names(MEs2)
htmap.colors <- gsub("ME", "", htmap.colors)

mod.sizes<-as.data.frame(mct)%>% #why filter here if it doesn't change anything in the heat map???
  filter(Freq>123)%>% #123 is the size of modules we are filtering by (for my dataset 124 is the size of the 14th largest module)
  arrange(moduleColors=c("black", "blue", "brown", "cyan", "darkgrey","darkmagenta", "darkorange", "darkred", "ivory",
                         "lightcyan", "lightcyan1", "paleturquoise", "pink", "sienna3"))


rownames(moduleTraitCor2) = gsub("ME", "", rownames(moduleTraitCor2))#all colors heatmap
#rownames(moduleTraitCor) = paste(rownames(moduleTraitCor), sep= " ", mod.sizes$Freq)

#quartz()
ht=Heatmap(moduleTraitCor2, name = "Corr.",
           col = blueWhiteRed(50), 
           row_names_side = "left", row_dend_side = "left",
           #right_annotation = size.annot,
           #width = unit(4, "in"), height = unit(8.5, "in"), 
           cluster_columns = FALSE,
           #cluster_rows = METree, row_split = 6, row_gap = unit(2.5, "mm"), border = TRUE,
           cell_fun = function(j, i, x, y, w, h, col) {
             if(heatmappval[i, j] <= 0.05) {
               grid.text(sprintf("%s", heatmappval[i, j]), x, y, gp = gpar(fontsize = 8, fontface = "bold"))
             }
             else {
               grid.text(sprintf("-"), x, y, gp = gpar(fontsize = 8, fontface = "plain"))
             }},
           column_names_gp =  gpar(fontsize = 10),
           column_names_rot = 45,
           row_names_gp = gpar(fontsize = 10, border = FALSE))
draw(ht)

#map.grob =  grid.grabExpr(draw(ht)) 

#### Ofra preheat Gene relationship to trait and important modules ####
# Define variable weight containing the weight column of datTrait - leave weight as variable, but change names in first 2 commands
weight = as.data.frame(datTraitsOutOfra$Pre.heat); 
#change to your trait name of interest: Combined, Control, Heat, pH, etc
names(weight) = "Pre.heat"

# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExprOutOfra, MEs, use = "p")); 
#finds pearson correlations of GE and module eigengenes
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

#making dataframe of pvalues for gene module membership values in each module
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExprOutOfra, weight, use = "p")); 

# data frame of correlations between expression and trait of interest using pearson correlations
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples)); 
#p-values for correlations between expression and trait of interest
names(geneTraitSignificance) = paste("GS.", names(weight), sep="");
names(GSPvalue) = paste("p.GS.", names(weight), sep="");
#geneTraitSignificance is a dataframe of correlation/covariance of GE data and stage value
#GSPvalue is a list of correlation of gTS pvalue



###Ofra preheat Visualization of Gene Networks####

#need to adjust list below for my treatments and analyses###
#modules of interest for preheat
moduleCols=c("ivory", "brown", "darkgrey")#0.45 threshold with softpower 5


#plot scatter plots of gene significance vs module membership for all of these modules of interest
#add correlation and p-value, use this to look at how strong the modules are.
#quartz()
# par(mfrow=c(1,1))
# par(mar = c(2, 2, 2, 2));
#par(bg = "grey"); # so we can see ivory points
for (module in moduleCols) {
  column = match(module, modNames);
  moduleGenes = moduleColors==module;
  verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                     abs(geneTraitSignificance[moduleGenes, 1]),
                     xlab = paste("ModMem", module),
                     ylab = "Gene Sig for Heat",
                     main = paste("MM vs. GS\n"),
                     cex.main = 1, cex.lab = 1, cex.axis = 1.2, col = module)
} 

#All modules for time checked with code above have positive correlations

#heatmaps of module expression with bar plot of eigengene for all samples
#use this to look at the different samples and make sure the modules make sense
#we want blocks by treatment, no single sample driving the differences.

#start with modules that looked strongest

#sort ME by Time    
# Add time from external metadata
MEs$Time <- datTraitsOutOfra$Time[match(rownames(MEs), rownames(datTraitsOutOfra))]

MEs.sorted <- MEs %>%
  arrange(Time) %>%
  dplyr::select(-Time)

which.module="ivory" #pick module of interest
ME=MEs.sorted[, paste("ME",which.module, sep="")]
genes=datExprOutOfra[,moduleColors==which.module ] # check with each color

#sort genes by Time
# Add time from external metadata
genes$Time <- datTraitsOutOfra$Time[match(rownames(genes), rownames(datTraitsOutOfra))]

genes.sorted <- genes %>%
  arrange(Time) %>%
  dplyr::select(-Time)

#quartz()
#par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
#plot expression  in all to make sure one sample isn't driving significance
plotMat(t(scale(genes.sorted) ),nrgcols=30,rlabels=F, clabels=rownames(genes.sorted), rcols=which.module,)

#par(mar=c(5, 4.2, 0, 0.7))
barplot(ME, col=which.module, main="", cex.main=2,
        ylab="eigengene expression",xlab="sample")

##########To output ME by sample, to see expression by sample/treatment

which.module="darkgrey" #pick module of interest
ME=MEs[, paste("ME",which.module, sep="")]
medarkgrey<-ME

which.module="brown"
ME=MEs[, paste("ME",which.module, sep="")]
mebrown<-ME

#meout<-data.frame(cbind(rownames(datExprOutOfav), medarkturquoise, melightcyan1))

#write.csv(meout,"Ofav_MEbySample_modsofinterest_9Dec2024.csv",quote=F,row.names=F) 



# get MEs from mods of interest from previously created ME dataframe 
# sort by time for boxplots
# Add time from external metadata
MEs$Time <- datTraitsOutOfra$Time[match(rownames(MEs), rownames(datTraitsOutOfra))]

MEs.sorted <- MEs %>%
  arrange(Time) %>%
  select(MEbrown, ,MEdarkgrey, MEivory,
         Time)%>%
  dplyr::select(-Time)

MEs.sorted$Time = c("preheat","preheat", "postheat", "postheat","preheat", "postheat",
                    "postheat", "preheat", "preheat", "preheat","postheat", "postheat", 
                    "postheat", "postheat", "preheat", "preheat")
# MEs.sorted$Treatment = c("control", "therm_var", "control", "therm_var","therm_var", "control",
#                          "therm_var", "control", "therm_var", "control","control","therm_var", 
#                          "control", "therm_var", "therm_var", "control")


medians <- MEs.sorted%>%
  group_by(Time)%>%
  dplyr::summarise(brown = median(MEbrown, na.rm=TRUE), darkgrey= median(MEdarkgrey, na.rm=TRUE), 
                   ivory = median(MEivory, na.rm=TRUE),
                   .groups="keep")


# medians <- MEs.sorted%>%
#   group_by(Symbiont)%>%
#   dplyr::summarise(plum2 = median(MEplum2, na.rm=TRUE), darkorange2 = median(MEdarkorange2, na.rm=TRUE),
#                    purple= median(MEpurple, na.rm=TRUE), brown= median(MEbrown, na.rm=TRUE), 
#                    lavenderblush3= median(MElavenderblush3, na.rm=TRUE),lightcyan= median(MElightcyan, na.rm=TRUE),
#                    .groups="keep")

# Set the desired order of Time levels
MEs.sorted$Time <- factor(MEs.sorted$Time, levels = c("preheat", "postheat"))

brown <- ggplot(MEs.sorted, aes(x = Time, y = MEbrown))+
  geom_boxplot(fill= "brown")+
  theme_classic()#+
# theme(axis.title.y=element_blank())+
# theme(axis.title.x=element_blank())+
# scale_y_continuous(breaks=c(-0.1,0, 0.1))+
# geom_hline(yintercept =0.06039775, linetype = "dashed", color = "grey")+
# geom_segment(x= "Heat", xend = "Heat", y = 0.06039775, yend = -0.04583882,color = "#D22B2B",
#              arrow = arrow(length = unit(0.03, "npc"), ends = "last", type = "closed"))+
# geom_segment(x= "pH", xend = "pH", y = 0.06039775, yend =	0.07624769,color = "#6082B6",
#              arrow = arrow(length = unit(0.03, "npc"), ends = "last", type = "closed"))+  
# geom_segment(x= "Combined", xend = "Combined", y =0.06039775, yend = -0.04583882,color = "#D22B2B",
#              arrow = arrow(length = unit(0.03, "npc"), ends = "last", type="closed"))+
# geom_segment(x= "Combined", xend = "Combined", y = -0.04583882, yend = (-0.04583882+(0.07624769-0.06039775)),color = "#6082B6",
#              arrow = arrow(length = unit(0.03, "npc"), ends = "last", type="closed"))
brown

darkgrey <- ggplot(MEs.sorted, aes(x = Time, y = MEdarkgrey))+
  geom_boxplot(fill= "darkgrey")+
  theme_classic()#+
# theme(axis.title.y=element_blank())+
# theme(axis.title.x=element_blank())+
# scale_y_continuous(breaks=c(-0.1,0, 0.1))+
# geom_hline(yintercept =0.06039775, linetype = "dashed", color = "grey")+
# geom_segment(x= "Heat", xend = "Heat", y = 0.06039775, yend = -0.04583882,color = "#D22B2B",
#              arrow = arrow(length = unit(0.03, "npc"), ends = "last", type = "closed"))+
# geom_segment(x= "pH", xend = "pH", y = 0.06039775, yend =	0.07624769,color = "#6082B6",
#              arrow = arrow(length = unit(0.03, "npc"), ends = "last", type = "closed"))+  
# geom_segment(x= "Combined", xend = "Combined", y =0.06039775, yend = -0.04583882,color = "#D22B2B",
#              arrow = arrow(length = unit(0.03, "npc"), ends = "last", type="closed"))+
# geom_segment(x= "Combined", xend = "Combined", y = -0.04583882, yend = (-0.04583882+(0.07624769-0.06039775)),color = "#6082B6",
#              arrow = arrow(length = unit(0.03, "npc"), ends = "last", type="closed"))
darkgrey

ivory <- ggplot(MEs.sorted, aes(x = Time, y = MEivory))+
  geom_boxplot(fill= "ivory")+
  theme_classic()#+
# theme(axis.title.y=element_blank())+
# theme(axis.title.x=element_blank()) +
# scale_y_continuous(breaks=c(-0.1,0, 0.1))+
# geom_hline(yintercept = -0.03785620, linetype = "dashed", color = "grey")+
# geom_segment(x= "Heat", xend = "Heat", y = -0.03785620, yend = 0.06115896,color = "grey",
#              arrow = arrow(length = unit(0.03, "npc"), ends = "last", type = "closed"))
#geom_segment(x= "pH", xend = "pH", y = -0.07594524, yend = -0.04762537,color = "#6082B6",
#             arrow = arrow(length = unit(0.03, "npc"), ends = "last", type = "closed"))
ivory

library(cowplot)
library(grid)
library(gridExtra)
#makes pretty labeled figure of all boxplots from above (run each for sym and treatment)
boxes <- plot_grid(brown, darkgrey, ivory, ncol = 3, align = "v", 
                   axis="b",rel_heights = c(1,1,1), labels = "AUTO")
boxes

#create common x and y labels
y.grob <- textGrob("Eigengene Expression", 
                   gp=gpar(fontsize=12), rot=90)

x.grob <- textGrob("Time", 
                   gp=gpar(fontsize=12))

#add to plot
allplots <-grid.arrange(arrangeGrob(boxes, left = y.grob, bottom = x.grob))

#combine into figure with the heatmap
fig1 <- plot_grid(map.grob, allplots, rel_widths = c(1,2), labels = c("A", ""))
fig1

#####plotting massive table of all information - module membership, genes, gene names, etc. ####
annot=read.table("~/path/Ofav_Made_iso2geneName_16Dec2024.tab",sep="\t",quote="")
iso2go=read.table("~/path/Ofav_Made_iso2go_9Dec2024.tab",sep="\t",quote="")

probes=colnames(datExprOutOfav)

probes2go=match(probes,iso2go$V1)
summary(probes2go) # ? genes annotated <-how is this determined?
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  1    3335    6750    6673    9936   13318    7619 

probes2annot = match(probes,annot$V1)
summary(probes2annot) # 14058 genes annotated <-how is this determined?
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  1    3320    6716    6641    9890   13260    7639

datGS.Traits=data.frame(cor(datExprOutOfav,datTraits1,use="p"))
names(datGS.Traits)=paste("cor",names(datGS.Traits),sep=".")
datME=moduleEigengenes(datExprOutOfav,moduleColors)$eigengenes
datKME=signedKME(datExprOutOfav, datME, outputColumnName="MM.")
datOutput=data.frame(ProbeID=names(datExprOutOfav),annot[probes2annot,],moduleColors,datKME,datGS.Traits)
#write.table(datOutput,"Ofav_AnnotatedNetworkAnalysisResultsHost_time_rlog_signed_sft4_merge0.45_16Dec2024.csv",row.names=F,sep=",")


#### Prep for preheat GO Analysis Yingqi's Categorical for Fisher's Exact test######
# Generate categorical data for GO analysis 
ModColors=c("brown", "darkgrey", "ivory", "sienna3", "brown", "pink", "darkred", "black",
            "darkmagenta", "lightcyan1", "lightcyan", "blue", "darkorange", "cyan", "paleturquoise")

# # #Categorical GeneOntology by Module (Fisher Exact Test) with 0 (absent) or 1 (present) for each module
# for (col in ModColors) {
#   tab=datOutput[,c(1,4)]
#   
#   tab$modulesColors=as.character(tab$moduleColors)
#   tab$moduleColors[tab$moduleColors!=col]<-0
#   tab$moduleColors[tab$moduleColors==col]<-1
#   tab$moduleColors=as.factor(tab$moduleColors) 
#   print(col)
#   print(summary(tab)) #do counts match table of module colors? YES!!!
#   print(head(tab))
#   
#   write.csv(tab,file=paste("GO_MM",col,"_categorical_phys.csv", sep=""),quote=F,row.names=F)
# }


#GeneOntology by Module (Fisher Exact Test) with 0 or KME value to retain information on strength of membership
for (col in ModColors) {
  tab <- datOutput[, c("ProbeID", "moduleColors", paste0("MM.", col))]
  
  # Replace KME values with 0 if gene is not in this module
  tab$KME_val <- ifelse(tab$moduleColors == col, tab[[paste0("MM.", col)]], 0)
  
  # Keep only the KME values of interest
  out <- tab[, c("ProbeID", "KME_val")]
  
  print(col)
  print(summary(out$KME_val))  # check distribution
  print(head(out))
  
  write.csv(out, 
            file = paste0("GO_MM_", col, "_KME.csv"), 
            quote = FALSE, row.names = FALSE)
}

###Now input these files into GO_MWU.R for analysis!!!#####


