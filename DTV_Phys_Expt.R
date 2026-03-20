#DTV_tank_Experiment - figs and analyses

library(ggplot2)
library(ggh4x)
library(Rmisc)
library(dplyr)
library(cowplot)
library(car)
library(lme4)
library(emmeans)
library(bestNormalize)
library(nlme)
library(performance)

setwd("~/")
phys.all <- read.csv("THERMVAR_MAIN_METADATA_Apr25_clonesfixed.csv")

str(phys.all)
phys.all$Timepoint <- factor(phys.all$Timepoint,levels=c("Preheat","Postheat"))


##### Fig 2 PAM - Fv/Fm, preheat to postheat ####

#TreatTime
phys.less <- phys.all[complete.cases(phys.all$pam_avg_t4),]

## Group by genotype, species, treatment, and timepoint and average PAM values
phys.less.geno <- phys.less %>%
  group_by(Genotype, Species, Treatment, Timepoint) %>%
  summarize_at(c('pam_avg_t3', 'pam_avg_t4'), mean, na.rm = TRUE)
phys.less.geno.df <- data.frame(phys.less.geno)

## Reshape to long format
phys.less.geno.time <- reshape(phys.less.geno.df,
                               varying = c('pam_avg_t3', 'pam_avg_t4'),
                               v.names = 'pam',
                               timevar = 'sampling',
                               direction = "long")
phys.less.geno.time$sampling <- as.factor(phys.less.geno.time$sampling)

## Create combined factor for Treatment and Timepoint
phys.less.geno.time$TreatTime <- paste(phys.less.geno.time$Treatment,
                                       ifelse(phys.less.geno.time$sampling == 1, "Pre-heat", "Post-heat"),
                                       sep = "_")
#convert to factors
phys.less.geno.time$TreatTime <- as.factor(phys.less.geno.time$TreatTime)
phys.less.geno.time$Species <- as.factor(phys.less.geno.time$Species)

## Summarize data with standard error
phys.less.geno.time.se <- summarySE(data = phys.less.geno.time,
                                    measurevar = "pam",
                                    groupvars = c("TreatTime", "Treatment", "Species", "sampling"))
#convert to factor
phys.less.geno.time.se$Species <- factor(phys.less.geno.time.se$Species)

#reorder variables for figures
f=c('Control_Pre-heat','Control_Post-heat', 'DTV_Pre-heat', 'DTV_Post-heat')
phys.less.geno.time <- within(phys.less.geno.time, TreatTime<- factor(TreatTime, levels=f))


## Plot
gg.pam <- ggplot(phys.less.geno.time, aes(x = sampling, y = pam, shape = Treatment, fill = TreatTime)) +
  scale_y_continuous(expand = c(0,0), limits = c(0.45,0.65)) + #this is set for our data, may need to adjust if using with different data
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1), alpha = 0.3, size = 4) +
  geom_errorbar(stat = 'summary',
                position = position_dodge(width = 0.5), width = 0.2) +
  geom_point(data = phys.less.geno.time.se,
             aes(x = sampling, y = pam),
             position = position_dodge(width = 0.5), size = 5) +
  facet_grid2(~Species,
              strip = strip_themed(background_x = elem_list_rect(fill = c("#8B7355","#CDAB7D")),
                                   text_x = elem_list_text(size = 14, face = "bold", color = "white")))+
  scale_fill_manual(name = "Treatment x Time", 
                    values = c("Control_Pre-heat" = "#2C2673",
                               "Control_Post-heat" = "#9187FF",
                               "DTV_Pre-heat" = "#E24A13",
                               "DTV_Post-heat" = "#FFA45E")) +
  scale_shape_manual(name = "Treatment",
                     values = c("Control" = 22,
                                "DTV" = 23)) +
  theme_bw() +
  theme(legend.title=element_text(size=14, face = "bold", color = "black"), 
        legend.text=element_text(size=14, face = "bold", color="black"),
        axis.text = element_text(color = "black", size = 14, face = "bold"),
        axis.title.y = element_text(size = 16, color = "black", face = "bold"),
        axis.ticks = element_line(color = "black"))+
  xlab("") +
  ylab("Fv/Fm") +
  ggtitle("") +
  scale_x_discrete(labels = c("Preheat", "Postheat")) #+
  #theme(axis.text.x = element_text(angle = 45, hjust = 1))

gg.pam



#### Fig S2 PAM - All timepoints ####

## Group by genotype, species, treatment, and timepoint and average PAM values
phys.less3.geno <- phys.less %>%
  group_by(Genotype, Species, Treatment, Timepoint) %>%
  summarize_at(c('pam_avg_t0','pam_avg_t1', 'pam_avg_t2','pam_avg_t3', 'pam_avg_t4'), mean, na.rm = TRUE)
phys.less3.geno.df <- data.frame(phys.less3.geno)

## Reshape to long format
phys.less3.geno.time <- reshape(phys.less3.geno.df,
                               varying = c('pam_avg_t0','pam_avg_t1', 'pam_avg_t2','pam_avg_t3', 'pam_avg_t4'),
                               v.names = 'pam',
                               timevar = 'sampling',
                               direction = "long")
phys.less3.geno.time$sampling <- as.factor(phys.less3.geno.time$sampling)

## Create combined factor for Treatment and Timepoint
phys.less3.geno.time$TreatTime <- paste(phys.less3.geno.time$Treatment,
                                       ifelse(phys.less3.geno.time$sampling == 1, "Pre-heat", "Post-heat"),
                                       sep = "_")
#convert to factors
phys.less3.geno.time$TreatTime <- as.factor(phys.less3.geno.time$TreatTime)
phys.less3.geno.time$Species <- as.factor(phys.less3.geno.time$Species)

## Summarize data with standard error
phys.less3.geno.time.se <- summarySE(data = phys.less3.geno.time,
                                    measurevar = "pam",
                                    groupvars = c("TreatTime", "Treatment", "Species", "sampling"))
#convert to factor
phys.less3.geno.time.se$Species <- factor(phys.less3.geno.time.se$Species)


#Plot
gg.pamAll <- ggplot(phys.less3.geno.time, aes(x = sampling, y = pam, shape = Treatment, fill = Treatment)) +
  scale_y_continuous(expand = c(0,0), limits = c(0.45,0.70)) + #this is set for our data, may need to adjust if using with different data
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1), alpha = 0.3, size = 4) +
  geom_errorbar(stat = 'summary',
                position = position_dodge(width = 0.5), width = 0.2) +
  geom_point(data = phys.less3.geno.time.se,
             aes(x = sampling, y = pam),
             position = position_dodge(width = 0.5), size = 5) +
  facet_grid2(~Species,
              strip = strip_themed(background_x = elem_list_rect(fill = c("#8B7355","#CDAB7D")),
                                   text_x = elem_list_text(size = 14, face = "bold", color = "white")))+
  scale_fill_manual(name = "Treatment", 
                    values = c("Control" = "#2C2673",
                               "DTV" = "#E24A13")) +
  scale_shape_manual(name = "Treatment",
                     values = c("Control" = 22,
                                "DTV" = 23)) +
  theme_bw() +
  theme(legend.title=element_text(size=14, face = "bold", color = "black"), 
        legend.text=element_text(size=14, face = "bold", color="black"),
        axis.text = element_text(color = "black", size = 14, face = "bold"),
        axis.title.y = element_text(size = 16, color = "black", face = "bold"),
        axis.title.x = element_text(size = 16, color = "black", face = "bold"),
        axis.ticks = element_line(color = "black"))+
  xlab("Day") +
  ylab("Fv/Fm") +
  ggtitle("") +
  scale_x_discrete(labels = c("0","30", "60","90", "115")) #+
#theme(axis.text.x = element_text(angle = 45, hjust = 1))
gg.pamAll


#### PAM stats  ####
### preheat and postheat ####
##both species
leveneTest(pam~Treatment,data=phys.less.geno.time) #ns
shapiro.test(phys.less.geno.time$pam) #ns

# Genotype is a random effect
lmer.pam <- lmer(pam ~ Treatment * Species * sampling + (1 | Genotype), data = phys.less.geno.time)

check_model(lmer.pam) # check assumptions more thoroughly, looks decent!

# Summary of the model
summary(lmer.pam)

# Type III ANOVA table (if needed)
Anova(lmer.pam, type = 3)
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: pam
# Chisq Df Pr(>Chisq)    
# (Intercept)                1529.4958  1  < 2.2e-16 ***
# Treatment                     1.8086  1  0.1786713    
# Species                       0.8466  1  0.3575199    
# sampling                     14.4222  1  0.0001461 ***
# Treatment:Species             2.3937  1  0.1218261    
# Treatment:sampling            0.2951  1  0.5869725    
# Species:sampling              0.3239  1  0.5692906    
# Treatment:Species:sampling    1.1311  1  0.2875474    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# aov.pam <- aov(pam~Treatment*Species*sampling,data=phys.less.geno.time) #doesn't control for clones
# summary(aov.pam) 
#                            Df   Sum Sq  Mean Sq F value  Pr(>F)   
# Treatment                   1 0.003137 0.003137   4.320 0.04576 * 
# Species                     1 0.000349 0.000349   0.481 0.49287   
# sampling                    1 0.007022 0.007022   9.672 0.00391 **
# Treatment:Species           1 0.000062 0.000062   0.085 0.77210   
# Treatment:sampling          1 0.000007 0.000007   0.010 0.92065   
# Species:sampling            1 0.000010 0.000010   0.013 0.90834   
# Treatment:Species:sampling  1 0.000163 0.000163   0.224 0.63918   
# Residuals                  32 0.023234 0.000726                   
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##favs only
phys.less.geno.time.fav <- subset(phys.less.geno.time,Species=="O. faveolata")

#Check assumptions for LME
# Homogeneity of variances
leveneTest(pam~Treatment,data=phys.less.geno.time.fav) #ns
# Normality of residuals
shapiro.test(phys.less.geno.time.fav$pam) #ns

# control for Genotype is a random effect
lmer.pam <- lmer(pam ~ Treatment * sampling + (1 | Genotype), data = phys.less.geno.time.fav)

check_model(lmer.pam) # check assumptions more thoroughly, looks decent!

# Summary of the model
summary(lmer.pam)

# Type III ANOVA table (if needed)
Anova(lmer.pam, type = 3)
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: pam
# Chisq Df Pr(>Chisq)    
# (Intercept)        1360.9250  1  < 2.2e-16 ***
#   Treatment             5.0963  1    0.02398 *  
#   sampling             40.6381  1  1.832e-10 ***
#   Treatment:sampling    0.8315  1    0.36184    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Normal anova, not suitable because comparing repeated measures
# aov.pam.fav <- aov(pam~Treatment*sampling,data=phys.less.geno.time.fav)
# summary(aov.pam.fav)
# Df   Sum Sq   Mean Sq F value Pr(>F)
# Treatment           1 0.000760 0.0007601   0.812  0.385
# sampling            1 0.002963 0.0029627   3.165  0.101
# Treatment:sampling  1 0.000038 0.0000375   0.040  0.845
# Residuals          12 0.011234 0.0009362  

#Pairwise comparisons
emmeans(lmer.pam, pairwise  ~ Treatment * sampling, adjust = "tukey")
# contrast                              estimate      SE df t.ratio p.value
# Control sampling1 - DTV sampling1      -0.0107 0.00475  9  -2.257  0.1797
# Control sampling1 - Control sampling2   0.0303 0.00475  9   6.375  0.0006 ***
# Control sampling1 - DTV sampling2       0.0134 0.00475  9   2.828  0.0778
# DTV sampling1 - Control sampling2       0.0410 0.00475  9   8.632  0.0001 ***
# DTV sampling1 - DTV sampling2           0.0242 0.00475  9   5.085  0.0030 ***
# Control sampling2 - DTV sampling2      -0.0168 0.00475  9  -3.547  0.0264 **
# 
# Degrees-of-freedom method: kenward-roger 
# P value adjustment: tukey method for comparing a family of 4 estimates 


##franks
phys.less.geno.time.fra <- subset(phys.less.geno.time,Species=="O. franksi")

leveneTest(pam~Treatment,data=phys.less.geno.time.fra) #ns
shapiro.test(phys.less.geno.time.fra$pam) #ns

# Genotype is a random effect
lmer.pam <- lmer(pam ~ Treatment * sampling + (1 | Genotype), data = phys.less.geno.time.fra)

check_model(lmer.pam) # check assumptions more thoroughly, looks decent!

# Summary of the model
summary(lmer.pam)

# Type III ANOVA table (if needed)
Anova(lmer.pam, type = 3)
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: pam
# Chisq Df Pr(>Chisq)    
# (Intercept)        1631.4958  1    < 2e-16 ***
#   Treatment             7.5868  1    0.00588 ** 
#   sampling              5.4447  1    0.01963 *  
#   Treatment:sampling    0.5612  1    0.45380    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#Pairwise comparisons
emmeans(lmer.pam, pairwise ~ Treatment * sampling, adjust = "tukey")
# contrast                              estimate     SE df t.ratio p.value
# Control sampling1 - DTV sampling1     -0.02817 0.0102  9  -2.754  0.0868 
# Control sampling1 - Control sampling2  0.02386 0.0102  9   2.333  0.1612
# Control sampling1 - DTV sampling2      0.00653 0.0102  9   0.638  0.9170
# DTV sampling1 - Control sampling2      0.05203 0.0102  9   5.088  0.0030 **
# DTV sampling1 - DTV sampling2          0.03469 0.0102  9   3.393  0.0333 *
# Control sampling2 - DTV sampling2     -0.01733 0.0102  9  -1.695  0.3796
# 
# Degrees-of-freedom method: kenward-roger 
# P value adjustment: tukey method for comparing a family of 4 estimates 

#ggsave("pam.pdf",width=4,height=3)


### All Pam T0-T90 time points ####
##favs only
phys.less3.geno.time.fav <- subset(phys.less3.geno.time,Species=="O. faveolata")
phys.less3.geno.time.fav <- subset(phys.less3.geno.time.fav,sampling!="5") #remove heat time point (T115)

leveneTest(pam~Treatment,data=phys.less3.geno.time.fav) #ns
shapiro.test(phys.less3.geno.time.fav$pam) #ns

# Genotype is a random effect
lmer.pam <- lmer(pam ~ Treatment * sampling + (1 | Genotype), data = phys.less3.geno.time.fav)
# lme.red <- lme(
#   pam ~ Treatment * sampling,                # Fixed effects
#   random = ~1 | Genotype,                         # Random intercepts for repeated measures
#   #weights = varIdent(form = ~1 | Treatment),      # Allows unequal variances by Treatment group
#   data = phys.less3.geno.time.fav
# )

check_model(lmer.pam) # check assumptions more thoroughly, pretty okay


# Summary of the model w/ ANOVA like summaries
Anova(lmer.pam, type = "III") 
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: pam
# Chisq Df Pr(>Chisq)    
# (Intercept)        1326.5500  1    < 2e-16 ***
# Treatment             3.2068  1    0.07333 .  
# sampling              6.8836  3    0.07570 .  
# Treatment:sampling    3.8249  3    0.28100  
  

#Pairwise comparisons
emmeans(lmer.pam, pairwise ~ Treatment * sampling, adjust = "tukey")
# $contrasts
# contrast                               estimate      SE df t.ratio p.value
# Control sampling1 - DTV sampling1     -1.44e-02 0.00804 21  -1.791  0.6324
# Control sampling1 - Control sampling2 -1.43e-02 0.00804 21  -1.784  0.6367
# Control sampling1 - DTV sampling2     -4.58e-02 0.00804 21  -5.702  0.0003
# Control sampling1 - Control sampling3 -8.28e-03 0.00804 21  -1.030  0.9641
# Control sampling1 - DTV sampling3     -2.58e-02 0.00804 21  -3.213  0.0667
# Control sampling1 - Control sampling4  5.04e-03 0.00804 21   0.627  0.9980
# Control sampling1 - DTV sampling4     -5.68e-03 0.00804 21  -0.707  0.9958
# DTV sampling1 - Control sampling2      5.56e-05 0.00804 21   0.007  1.0000
# DTV sampling1 - DTV sampling2         -3.14e-02 0.00804 21  -3.912  0.0151 **
# DTV sampling1 - Control sampling3      6.11e-03 0.00804 21   0.761  0.9935
# DTV sampling1 - DTV sampling3         -1.14e-02 0.00804 21  -1.423  0.8370
# DTV sampling1 - Control sampling4      1.94e-02 0.00804 21   2.418  0.2829
# DTV sampling1 - DTV sampling4          8.71e-03 0.00804 21   1.084  0.9534
# Control sampling2 - DTV sampling2     -3.15e-02 0.00804 21  -3.919  0.0149 **
# Control sampling2 - Control sampling3  6.06e-03 0.00804 21   0.754  0.9939
# Control sampling2 - DTV sampling3     -1.15e-02 0.00804 21  -1.429  0.8338
# Control sampling2 - Control sampling4  1.94e-02 0.00804 21   2.411  0.2860
# Control sampling2 - DTV sampling4      8.65e-03 0.00804 21   1.077  0.9549
# DTV sampling2 - Control sampling3      3.75e-02 0.00804 21   4.672  0.0027
# DTV sampling2 - DTV sampling3          2.00e-02 0.00804 21   2.489  0.2528
# DTV sampling2 - Control sampling4      5.09e-02 0.00804 21   6.330  0.0001
# DTV sampling2 - DTV sampling4          4.01e-02 0.00804 21   4.995  0.0013 **
# Control sampling3 - DTV sampling3     -1.75e-02 0.00804 21  -2.183  0.3994
# Control sampling3 - Control sampling4  1.33e-02 0.00804 21   1.658  0.7123
# Control sampling3 - DTV sampling4      2.60e-03 0.00804 21   0.323  1.0000
# DTV sampling3 - Control sampling4      3.09e-02 0.00804 21   3.841  0.0177
# DTV sampling3 - DTV sampling4          2.01e-02 0.00804 21   2.506  0.2458
# Control sampling4 - DTV sampling4     -1.07e-02 0.00804 21  -1.334  0.8754

##franks
phys.less3.geno.time.fra <- subset(phys.less3.geno.time,Species=="O. franksi")
phys.less3.geno.time.fra <- subset(phys.less3.geno.time.fra,sampling!="5") #remove heat time point (T115)

leveneTest(pam~Treatment,data=phys.less3.geno.time.fra) #ns
shapiro.test(phys.less3.geno.time.fra$pam) #ns

# Genotype is a random effect
lmer.pam <- lmer(pam ~ Treatment * sampling + (1 | Genotype), data = phys.less3.geno.time.fra)
# lme.red <- lme(
#   pam ~ Treatment * sampling,                # Fixed effects
#   random = ~1 | Genotype,                         # Random intercepts for repeated measures
#   #weights = varIdent(form = ~1 | Treatment),      # Allows unequal variances by Treatment group
#   data = phys.less3.geno.time.fra
# )

check_model(lmer.pam) # check assumptions more thoroughly, looks pretty good


# Summary of the model w/ ANOVA like summaries
Anova(lmer.pam, type = "III") 
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: pam
# Chisq Df Pr(>Chisq)    
# (Intercept)        2313.4592  1  < 2.2e-16 ***
#   Treatment            12.1520  1  0.0004904 ***
#   sampling             11.5822  3  0.0089605 ** 
#   Treatment:sampling    0.8676  3  0.8332315 

#Pairwise comparisons
# emmeans(lmer.pam, pairwise ~ Treatment * sampling, adjust = "tukey")
# $contrasts
# contrast                              estimate     SE df t.ratio p.value
# Control sampling1 - DTV sampling1     -0.03994 0.0115 21  -3.486  0.0380 **
# Control sampling1 - Control sampling2 -0.03594 0.0115 21  -3.137  0.0777
# Control sampling1 - DTV sampling2     -0.07704 0.0115 21  -6.723  <.0001
# Control sampling1 - Control sampling3 -0.02893 0.0115 21  -2.525  0.2385
# Control sampling1 - DTV sampling3     -0.06142 0.0115 21  -5.360  0.0006
# Control sampling1 - Control sampling4 -0.01486 0.0115 21  -1.297  0.8901
# Control sampling1 - DTV sampling4     -0.04303 0.0115 21  -3.755  0.0213
# DTV sampling1 - Control sampling2      0.00400 0.0115 21   0.349  1.0000
# DTV sampling1 - DTV sampling2         -0.03710 0.0115 21  -3.237  0.0635
# DTV sampling1 - Control sampling3      0.01101 0.0115 21   0.961  0.9752
# DTV sampling1 - DTV sampling3         -0.02147 0.0115 21  -1.874  0.5814
# DTV sampling1 - Control sampling4      0.02508 0.0115 21   2.189  0.3962
# DTV sampling1 - DTV sampling4         -0.00308 0.0115 21  -0.269  1.0000
# Control sampling2 - DTV sampling2     -0.04110 0.0115 21  -3.587  0.0306 **
# Control sampling2 - Control sampling3  0.00701 0.0115 21   0.612  0.9983
# Control sampling2 - DTV sampling3     -0.02547 0.0115 21  -2.223  0.3780
# Control sampling2 - Control sampling4  0.02108 0.0115 21   1.840  0.6023
# Control sampling2 - DTV sampling4     -0.00708 0.0115 21  -0.618  0.9982
# DTV sampling2 - Control sampling3      0.04811 0.0115 21   4.199  0.0080
# DTV sampling2 - DTV sampling3          0.01562 0.0115 21   1.364  0.8633
# DTV sampling2 - Control sampling4      0.06218 0.0115 21   5.427  0.0005
# DTV sampling2 - DTV sampling4          0.03401 0.0115 21   2.968  0.1080
# Control sampling3 - DTV sampling3     -0.03249 0.0115 21  -2.835  0.1387
# Control sampling3 - Control sampling4  0.01407 0.0115 21   1.228  0.9143
# Control sampling3 - DTV sampling4     -0.01410 0.0115 21  -1.230  0.9135
# DTV sampling3 - Control sampling4      0.04656 0.0115 21   4.063  0.0108
# DTV sampling3 - DTV sampling4          0.01839 0.0115 21   1.605  0.7426
# Control sampling4 - DTV sampling4     -0.02817 0.0115 21  -2.458  0.2657


#### Fig 2 Red channel only, preheat to postheat ####

#pull subset of data w/o zeroes for red channel
phys.less <- phys.all[complete.cases(phys.all$red_chan_avg_t4),]


## Group by each combination of genotype, species, treatment, and timepoint and channel values and average
phys.less.geno <- phys.less %>% 
  group_by(Genotype,Species,Treatment,Timepoint) %>%
  summarize_at(c('red_chan_avg_t3','red_chan_avg_t4'),mean,na.rm=T)
phys.less.geno.df <- data.frame(phys.less.geno)

## Reshape to long format
phys.less.geno.time <- reshape(phys.less.geno.df,
                               varying=c('red_chan_avg_t3','red_chan_avg_t4'),
                               v.names='red_chan',
                               timevar='sampling',
                               direction="long")
phys.less.geno.time$sampling <- as.factor(phys.less.geno.time$sampling)


## Create combined factor for Treatment and Timepoint = TreatTime
phys.less.geno.time$TreatTime <- paste(phys.less.geno.time$Treatment,
                                       ifelse(phys.less.geno.time$sampling == 1, "Pre-heat", "Post-heat"),
                                       sep = "_")
#convert to factors
phys.less.geno.time$TreatTime <- as.factor(phys.less.geno.time$TreatTime)
phys.less.geno.time$Species <- as.factor(phys.less.geno.time$Species)

## Summarize data with standard error
phys.less.geno.time.se <- summarySE(data=phys.less.geno.time,
                                    measurevar="red_chan",
                                    groupvars=c("TreatTime", "Treatment","Species","sampling"))

#convert to factor
phys.less.geno.time.se$Species <- factor(phys.less.geno.time.se$Species)

#reorder variables for figures
f=c('Control_Pre-heat','Control_Post-heat', 'DTV_Pre-heat', 'DTV_Post-heat')
phys.less.geno.time <- within(phys.less.geno.time, TreatTime<- factor(TreatTime, levels=f))

#Plot
gg.red <- ggplot(phys.less.geno.time, aes(x = sampling, y = red_chan, shape = Treatment, fill = TreatTime)) +
  scale_y_reverse(expand = c(0, 0), limits = c(80, 25)) + #this is set for our data, may need to adjust if using with different data
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1), alpha = 0.3, size = 4) +
  geom_errorbar(stat = 'summary',
                position = position_dodge(width = 0.5), width = 0.2) +
  geom_point(data = phys.less.geno.time.se,
             aes(x = sampling, y = red_chan),
             position = position_dodge(width = 0.5), size = 5) +
  facet_grid2(~Species,
              strip = strip_themed(background_x = elem_list_rect(fill = c("#8B7355","#CDAB7D")),
                                   text_x = elem_list_text(size = 14, face = "bold", color = "black")))+
  scale_fill_manual(name = "Treatment x Time", 
                    values = c("Control_Pre-heat" = "#2C2673",
                               "Control_Post-heat" = "#9187FF",
                               "DTV_Pre-heat" = "#E24A13",
                               "DTV_Post-heat" = "#FFA45E")) +
  scale_shape_manual(name = "Treatment",
                     values = c("Control" = 22,
                                "DTV" = 23)) +
  theme_bw() +
  theme(legend.title=element_text(size=14, face = "bold", color = "black"), 
        legend.text=element_text(size=14, face = "bold", color="black"),
        axis.text = element_text(color = "black", size = 14, face = "bold"),
        axis.title.y = element_text(size = 16, color = "black", face = "bold"),
        axis.ticks = element_line(color = "black"))+
  xlab("") +
  ylab("Red Channel Intensity") +
  ggtitle("") +
  scale_x_discrete(labels = c("Pre-heat", "Post-heat")) #+
#theme(axis.text.x = element_text(angle = 45, hjust = 1))
gg.red



#### Fig S3 Red channel only for ALL timepoints ####

#pull subset of data w/o zeroes for red channel
phys.less2 <- phys.all[complete.cases(phys.all$red_chan_avg_t4),]


## Group by each combination of genotype, species, treatment, and timepoint and channel values and average
phys.less2.geno <- phys.less2 %>% 
  group_by(Genotype,Species,Treatment,Timepoint) %>%
  summarize_at(c('red_chan_avg_t0','red_chan_avg_t1','red_chan_avg_t3','red_chan_avg_t4'),mean,na.rm=T)
phys.less2.geno.df <- data.frame(phys.less2.geno)

## Reshape to long format
phys.less2.geno.time <- reshape(phys.less2.geno.df,
                                varying=c('red_chan_avg_t0','red_chan_avg_t1','red_chan_avg_t3','red_chan_avg_t4'),
                                v.names='red_chan',
                                timevar='sampling',
                                direction="long")
phys.less2.geno.time$sampling <- as.factor(phys.less2.geno.time$sampling)


## Create combined factor for Treatment and Timepoint = TreatTime
phys.less2.geno.time$TreatTime <- paste(phys.less2.geno.time$Treatment,
                                        ifelse(phys.less2.geno.time$sampling == 1, "Pre-heat", "Post-heat"),
                                        sep = "_")
#convert to factors
phys.less2.geno.time$TreatTime <- as.factor(phys.less2.geno.time$TreatTime)
phys.less2.geno.time$Species <- as.factor(phys.less2.geno.time$Species)

## Summarize data with standard error
phys.less2.geno.time.se <- summarySE(data=phys.less2.geno.time,
                                     measurevar="red_chan",
                                     groupvars=c("TreatTime", "Treatment","Species","sampling"))

#convert to factor
phys.less2.geno.time.se$Species <- factor(phys.less2.geno.time.se$Species)

#Plot
gg.red <- ggplot(phys.less2.geno.time, aes(x = sampling, y = red_chan, shape = Treatment, fill = Treatment)) +
  scale_y_reverse(expand = c(0, 0), limits = c(110, 20)) + #this is set for our data, may need to adjust if using with different data
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1), alpha = 0.3, size = 4) +
  geom_errorbar(stat = 'summary',
                position = position_dodge(width = 0.5), width = 0.2) +
  geom_point(data = phys.less2.geno.time.se,
             aes(x = sampling, y = red_chan),
             position = position_dodge(width = 0.5), size = 5) +
  facet_grid2(~Species,
              strip = strip_themed(background_x = elem_list_rect(fill = c("#8B7355","#CDAB7D")),
                                   text_x = elem_list_text(size = 14, face = "bold", color = "white")))+
  scale_fill_manual(name = "Treatment", 
                    values = c("Control" = "#2C2673",
                               "DTV" = "#E24A13")) +
  scale_shape_manual(name = "Treatment",
                     values = c("Control" = 22,
                                "DTV" = 23)) +
  theme_bw() +
  theme(legend.title=element_text(size=14, face = "bold", color = "black"), 
        legend.text=element_text(size=14, face = "bold", color="black"),
        axis.text = element_text(color = "black", size = 14, face = "bold"),
        axis.title.y = element_text(size = 16, color = "black", face = "bold"),
        axis.title.x = element_text(size = 16, color = "black", face = "bold"),
        axis.ticks = element_line(color = "black"))+
  xlab("Day") +
  ylab("Red Channel Intensity") +
  ggtitle("") +
  scale_x_discrete(labels = c("0","30", "90", "115")) #+
#theme(axis.text.x = element_text(angle = 45, hjust = 1))
gg.red


#### Red Chan Color stats ####

### Red Intensity Only ### this is what we used in the manuscript

### preheat and postheat ####
##both species
phys.less.geno.time$Treatment <- as.factor(phys.less.geno.time$Treatment)

lmer.red <- lmer(red_chan ~ Treatment * sampling + (1|Genotype),
               data = phys.less.geno.time)
plot(lmer.red)                # residuals vs fitted
check_model(lmer.red)         # looks good enough to continue

# Summary of the model w/ ANOVA like summaries
Anova(lmer.red, type = "III") 
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: red_chan
# Chisq Df Pr(>Chisq)    
# (Intercept)        433.0361  1  < 2.2e-16 ***
# Treatment            7.2456  1   0.007108 ** 
# sampling             2.7393  1   0.097909 .  
# Treatment:sampling   4.2404  1   0.039473 * 

#Pairwise comparisons
emmeans(lmer.red, pairwise ~ Treatment * sampling, adjust = "tukey")
# $contrasts
# contrast                              estimate   SE df t.ratio p.value
# Control sampling1 - DTV sampling1         7.45 2.77 21   2.692  0.0608
# Control sampling1 - Control sampling2    -4.58 2.77 21  -1.655  0.3712
# Control sampling1 - DTV sampling2        -5.19 2.77 21  -1.876  0.2683
# DTV sampling1 - Control sampling2       -12.04 2.77 21  -4.347  0.0015 **
# DTV sampling1 - DTV sampling2           -12.65 2.77 21  -4.567  0.0009 ***
# Control sampling2 - DTV sampling2        -0.61 2.77 21  -0.220  0.9961
# 
# Degrees-of-freedom method: kenward-roger 
# P value adjustment: tukey method for comparing a family of 4 estimates 



# leveneTest(red_chan~Treatment,data=phys.less.geno.time) #sig
# phys.less.geno.time$log_redchan <- log(phys.less.geno.time$red_chan) #transform - still not equal
# phys.less.geno.time$log_redchan <- sqrt(phys.less.geno.time$red_chan) #transform - still not equal, run with variance unequal in LME 
# leveneTest(log_redchan~Treatment,data=phys.less.geno.time) #sig
# shapiro.test(phys.less.geno.time$red_chan) #ns
#bestNormalize(phys.less.geno.time$red_chan,loo=T)
#red.ord <- orderNorm(phys.less.geno.timeo$red_chan)
#phys.less.geno.time$red.ord <- red.ord[["x.t"]]
#shapiro.test(phys.less.geno.time$red.ord) #ns

# # Genotype is a random effect
# lme.red <- lme(
#   red_chan ~ Treatment * sampling * Species,                # Fixed effects
#   random = ~1 | Genotype,                         # Random intercepts for repeated measures
#   weights = varIdent(form = ~1 | Treatment),      # Allows unequal variances by Treatment group
#   data = phys.less.geno.time
# )


# Summary of the model w/ ANOVA like summaries
# summary(lme.red)
# Fixed effects:  red_chan ~ Treatment * sampling * Species 
# Value Std.Error DF   t-value p-value
# (Intercept)                               54.61250  2.780839 18 19.638859  0.0000 ***
# TreatmentDTV                             -10.69583  3.559545 18 -3.004832  0.0076 **
# sampling2                                  3.57917  2.577531 18  1.388603  0.1819
# SpeciesO. franksi                          2.43333  3.932700  6  0.618744  0.5588
# TreatmentDTV:sampling2                     6.12917  5.033956 18  1.217565  0.2391
# TreatmentDTV:SpeciesO. franksi             6.48333  5.033956 18  1.287920  0.2141
# sampling2:SpeciesO. franksi                2.00833  3.645180 18  0.550956  0.5884
# TreatmentDTV:sampling2:SpeciesO. franksi   3.87083  7.119089 18  0.543726  0.5933


##favs
phys.less.geno.time.fav <- subset(phys.less.geno.time,Species=="O. faveolata")
str(phys.less.geno.time.fav)
phys.less.geno.time.fav$Treatment <- as.factor(phys.less.geno.time.fav$Treatment)
str(phys.less.geno.time.fav)

lmer.red <- lmer(red_chan ~ Treatment * sampling + (1|Genotype),
                 data = phys.less.geno.time.fav)
plot(lmer.red)                # residuals vs fitted
check_model(lmer.red)         # Homogeneity of Variance is not flat - maybe do lme instead

# Summary of the model w/ ANOVA like summaries
Anova(lmer.red, type = "III") 
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: red_chan
# Chisq Df Pr(>Chisq)    
# (Intercept)        240.9661  1  < 2.2e-16 ***
# Treatment           12.7609  1  0.0003539 ***
# sampling             1.4290  1  0.2319358    
# Treatment:sampling   2.0952  1  0.1477620 


#Pairwise comparisons
emmeans(lmer.red, pairwise ~ Treatment * sampling, adjust = "tukey")
# contrasts
# contrast                              estimate   SE df t.ratio p.value
# Control sampling1 - DTV sampling1       10.696 2.99  9   3.572  0.0255 *
# Control sampling1 - Control sampling2   -3.579 2.99  9  -1.195  0.6445
# Control sampling1 - DTV sampling2        0.988 2.99  9   0.330  0.9868
# DTV sampling1 - Control sampling2      -14.275 2.99  9  -4.768  0.0046 **
# DTV sampling1 - DTV sampling2           -9.708 2.99  9  -3.242  0.0417 *
# Control sampling2 - DTV sampling2        4.567 2.99  9   1.525  0.4628

# leveneTest(red_chan~Treatment,data=phys.less.geno.time.fav) #ns
# shapiro.test(phys.less.geno.time.fav$red_chan) #ns
# 
# # Genotype is a random effect
# lme.red <- lme(
#   red_chan ~ Treatment * sampling,                # Fixed effects
#   random = ~1 | Genotype,                         # Random intercepts for repeated measures
#   weights = varIdent(form = ~1 | Treatment),      # Allows unequal variances by Treatment group
#   data = phys.less.geno.time.fav
# )
# 
# check_model(lme.red) # check assumptions more thoroughly, appropriate to proceed
# 
# # Summary of the model
# summary(lme.red)
# Fixed effects:  red_chan ~ Treatment * sampling 
# Value Std.Error DF   t-value p-value
# (Intercept)             54.61250  3.258986  9 16.757513  0.0000
# TreatmentDTV           -10.69583  2.787250  9 -3.837415  0.0040 **
# sampling2                3.57917  1.275466  9  2.806163  0.0205 *
# TreatmentDTV:sampling2   6.12917  3.941766  9  1.554929  0.1544



#Pairwise comparisons
# emmeans(lme.red, pairwise ~ Treatment * sampling, adjust = "tukey")
# $contrasts
# contrast                              estimate   SE df t.ratio p.value
# Control sampling1 - DTV sampling1       10.696 2.99  9   3.572  0.0255 *
# Control sampling1 - Control sampling2   -3.579 2.99  9  -1.195  0.6445
# Control sampling1 - DTV sampling2        0.988 2.99  9   0.330  0.9868
# DTV sampling1 - Control sampling2      -14.275 2.99  9  -4.768  0.0046 **
# DTV sampling1 - DTV sampling2           -9.708 2.99  9  -3.242  0.0417 *
# Control sampling2 - DTV sampling2        4.567 2.99  9   1.525  0.4628

##franks
phys.less.geno.time.fra <- subset(phys.less.geno.time,Species=="O. franksi")
str(phys.less.geno.time.fra)
phys.less.geno.time.fra$Treatment <- as.factor(phys.less.geno.time.fra$Treatment)
str(phys.less.geno.time.fra)

leveneTest(red_chan~Treatment,data=phys.less.geno.time.fra) #ns
shapiro.test(phys.less.geno.time.fra$red_chan) #ns

lmer.red <- lmer(red_chan ~ Treatment * sampling + (1|Genotype),
                 data = phys.less.geno.time.fra)
plot(lmer.red)                # residuals vs fitted
check_model(lmer.red)         # Homogeneity of Variance is not flat - maybe do lme instead

# Summary of the model w/ ANOVA like summaries
Anova(lmer.red, type = "III") 
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: red_chan
# Chisq Df Pr(>Chisq)    
# (Intercept)        422.1226  1    < 2e-16 ***
# Treatment            1.1509  1    0.28336    
# sampling             2.0249  1    0.15474    
# Treatment:sampling   3.2429  1    0.07173 .  


#Pairwise comparisons
emmeans(lmer.red, pairwise ~ Treatment * sampling, adjust = "tukey")
# $contrasts
# contrast                              estimate   SE df t.ratio p.value
# Control sampling1 - DTV sampling1         4.21 3.93  9   1.073  0.7135
# Control sampling1 - Control sampling2    -5.59 3.93  9  -1.423  0.5170
# Control sampling1 - DTV sampling2       -11.38 3.93  9  -2.897  0.0701
# DTV sampling1 - Control sampling2        -9.80 3.93  9  -2.496  0.1273
# DTV sampling1 - DTV sampling2           -15.59 3.93  9  -3.970  0.0142 *
# Control sampling2 - DTV sampling2        -5.79 3.93  9  -1.474  0.4896


#lmer.red <- lmer(red_chan ~ Treatment * sampling + (1 | Genotype), data = phys.less.geno.time.fra) - could use but has high singularity, better to use lme
# lme.red <- lme(
#   red_chan ~ Treatment * sampling,                # Fixed effects
#   random = ~1 | Genotype,                         # Random intercepts for repeated measures
#   # weights = varIdent(form = ~1 | Treatment),      # Allows unequal variances by Treatment group
#   data = phys.less.geno.time.fra
# )
# 
# check_model(lme.red) # check assumptions more thoroughly, assumptions are reasonably well met.
# 
# # Summary of the model
# summary(lme.red)

# 
# #Pairwise comparisons
# emmeans(lme.red, pairwise ~ Treatment * sampling, adjust = "tukey")




### All Red chan T0-T90 time points ####

### both species 
phys.less2.geno.time$Treatment <- as.factor(phys.less2.geno.time$Treatment)

leveneTest(red_chan~Treatment,data=phys.less2.geno.time) #ns
shapiro.test(phys.less2.geno.time$red_chan) #sig   
phys.less2.geno.time$log_redchan <- log(phys.less2.geno.time$red_chan) #transform - still not equal
shapiro.test(phys.less2.geno.time$log_redchan) #ns

lmer.red <- lmer(red_chan ~ Treatment * sampling + (1|Genotype),
                 data = phys.less2.geno.time)
plot(lmer.red)                # residuals vs fitted
check_model(lmer.red)         # Linearity is sketch and homogeneity is rough - let's see what corals separated looks like

# Summary of the model w/ ANOVA like summaries
Anova(lmer.red, type = "III") 


#Pairwise comparisons
emmeans(lmer.red, pairwise ~ Treatment * sampling, adjust = "tukey")

# Genotype is a random effect
# lme.red <- lme(
#   log_redchan ~ Treatment * sampling * Species,                # Fixed effects
#   random = ~1 | Genotype,                         # Random intercepts for repeated measures
#   #weights = varIdent(form = ~1 | Treatment),      # Allows unequal variances by Treatment group
#   data = phys.less2.geno.time
# )
# 
# check_model(lme.red) # check assumptions more thoroughly, looks pretty good
# 
# 
# # Summary of the model w/ ANOVA like summaries
# Anova(lme.red, type = "III") 



##franks
phys.less2.geno.time.fra <- subset(phys.less2.geno.time,Species=="O. franksi")
phys.less2.geno.time.fra <- subset(phys.less2.geno.time.fra,sampling!="4") #No photos day 60, so only 4 timepoints remain
str(phys.less2.geno.time.fra)
#phys.int.geno.fra$Treatment <- as.factor(phys.int.geno.fra$Treatment)
#str(phys.int.geno.fra)

leveneTest(red_chan~Treatment,data=phys.less2.geno.time.fra) #ns
shapiro.test(phys.less2.geno.time.fra$red_chan) #not normal

lmer.red <- lmer(red_chan ~ Treatment * sampling + (1|Genotype), data = phys.less2.geno.time.fra)
#plot(lmer.red)                # residuals vs fitted
check_model(lmer.red)         # Looks decent

# Summary of the model w/ ANOVA like summaries
Anova(lmer.red, type = "III") 
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: red_chan
# Chisq Df Pr(>Chisq)    
# (Intercept)        77.3886  1  < 2.2e-16 ***
#   Treatment           0.4271  1     0.5134    
#   sampling           36.6623  2  1.094e-08 ***
#   Treatment:sampling  1.2681  2     0.5304   

#Pairwise comparisons
emmeans(lmer.red, pairwise ~ Treatment * sampling, adjust = "tukey")
# $contrasts
# contrast                              estimate   SE df t.ratio p.value
# Control sampling1 - DTV sampling1        -3.38 5.18 15  -0.654  0.9845
# Control sampling1 - Control sampling2   -31.10 5.18 15  -6.008  0.0003 **
# Control sampling1 - DTV sampling2       -27.91 5.18 15  -5.391  0.0009
# Control sampling1 - Control sampling3   -12.18 5.18 15  -2.352  0.2338
# Control sampling1 - DTV sampling3        -7.96 5.18 15  -1.538  0.6472
# DTV sampling1 - Control sampling2       -27.72 5.18 15  -5.354  0.0009
# DTV sampling1 - DTV sampling2           -24.53 5.18 15  -4.738  0.0029 **
# DTV sampling1 - Control sampling3        -8.79 5.18 15  -1.698  0.5529
# DTV sampling1 - DTV sampling3            -4.58 5.18 15  -0.884  0.9445
# Control sampling2 - DTV sampling2         3.19 5.18 15   0.616  0.9880
# Control sampling2 - Control sampling3    18.93 5.18 15   3.656  0.0233 **
# Control sampling2 - DTV sampling3        23.14 5.18 15   4.470  0.0049
# DTV sampling2 - Control sampling3        15.74 5.18 15   3.040  0.0733
# DTV sampling2 - DTV sampling3            19.95 5.18 15   3.853  0.0160 **
# Control sampling3 - DTV sampling3         4.21 5.18 15   0.814  0.9605

# Genotype is a random effect
# lme.red <- lme(
#   red_chan ~ Treatment * sampling,                # Fixed effects
#   random = ~1 | Genotype,                         # Random intercepts for repeated measures
#   #weights = varIdent(form = ~1 | Treatment),      # Allows unequal variances by Treatment group
#   data = phys.less2.geno.time.fra
# )
# 
# check_model(lme.red) # check assumptions more thoroughly, looks pretty good
# 
# # Summary of the model w/ ANOVA like summaries
# Anova(lme.red, type = "III") 

# #Pairwise comparisons
# emmeans(lme.red, pairwise ~ Treatment * sampling, adjust = "tukey")


##favs
phys.less2.geno.time.fav <- subset(phys.less2.geno.time,Species=="O. faveolata")
phys.less2.geno.time.fav  <- subset(phys.less2.geno.time.fav ,sampling!="4") #No photos day 60, so only 4 timepoints remain
#str(phys.less2.geno.time.fav)


leveneTest(red_chan~Treatment,data=phys.less2.geno.time.fav) #ns
shapiro.test(phys.less2.geno.time.fav$red_chan) #sig

lmer2.red <- lmer(red_chan ~ Treatment * sampling + (1|Genotype), data = phys.less2.geno.time.fav)
plot(lmer2.red)                # residuals vs fitted
check_model(lmer2.red)         # looks borderline, but going to lean into the robustness of lmer

# Summary of the model w/ ANOVA like summaries
Anova(lmer2.red, type = "III") 
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: red_chan
# Chisq Df Pr(>Chisq)    
# (Intercept)        69.6132  1  < 2.2e-16 ***
#   Treatment           0.0006  1     0.9805    
#   sampling           32.8204  2  7.467e-08 ***
#   Treatment:sampling  2.9974  2     0.2234  

#Pairwise comparisons
emmeans(lmer2.red, pairwise ~ Treatment * sampling, adjust = "tukey")
# $contrasts
# contrast                              estimate   SE df t.ratio p.value
# Control sampling1 - DTV sampling1        0.125 5.12 15   0.024  1.0000
# Control sampling1 - Control sampling2  -26.950 5.12 15  -5.266  0.0011 **
# Control sampling1 - DTV sampling2      -27.367 5.12 15  -5.347  0.0009
# Control sampling1 - Control sampling3   -3.475 5.12 15  -0.679  0.9817
# Control sampling1 - DTV sampling3        7.221 5.12 15   1.411  0.7205
# DTV sampling1 - Control sampling2      -27.075 5.12 15  -5.290  0.0010
# DTV sampling1 - DTV sampling2          -27.492 5.12 15  -5.372  0.0009 **
# DTV sampling1 - Control sampling3       -3.600 5.12 15  -0.703  0.9786
# DTV sampling1 - DTV sampling3            7.096 5.12 15   1.386  0.7342
# Control sampling2 - DTV sampling2       -0.417 5.12 15  -0.081  1.0000
# Control sampling2 - Control sampling3   23.475 5.12 15   4.587  0.0039 **
# Control sampling2 - DTV sampling3       34.171 5.12 15   6.677  0.0001 
# DTV sampling2 - Control sampling3       23.892 5.12 15   4.668  0.0034
# DTV sampling2 - DTV sampling3           34.587 5.12 15   6.758  0.0001 **
# Control sampling3 - DTV sampling3       10.696 5.12 15   2.090  0.3424

# Genotype is a random effect
# lme.red <- lme(
#   log_redchan ~ Treatment * sampling,                # Fixed effects
#   random = ~1 | Genotype,                         # Random intercepts for repeated measures
#   weights = varIdent(form = ~1 | Treatment),      # Allows unequal variances by Treatment group
#   data = phys.less2.geno.time.fav
# )
# 
# check_model(lme.red) # check assumptions more thoroughly, looks pretty good
# 
# 
# # Summary of the model w/ ANOVA like summaries
# Anova(lme.red, type = "III") 
 
# #Pairwise comparisons
# emmeans(lme.red, pairwise ~ Treatment * sampling, adjust = "tukey")


### Buoyant weight ####

#pull subset of data w/o zeroes for red channel
phys.less <- phys.all[complete.cases(phys.all$bw_avg_t4),]

## % change
phys.less$bw_t1t0 <- ((phys.less$bw_avg_t1-phys.less$bw_avg_t0)/phys.less$bw_avg_t0)*100
phys.less$bw_t2t1 <- ((phys.less$bw_avg_t2-phys.less$bw_avg_t1)/phys.less$bw_avg_t1)*100
phys.less$bw_t3t2 <- ((phys.less$bw_avg_t3-phys.less$bw_avg_t2)/phys.less$bw_avg_t2)*100
phys.less$bw_t4t3 <- ((phys.less$bw_avg_t4-phys.less$bw_avg_t3)/phys.less$bw_avg_t3)*100


#### % change in weight preheat to postheat ####
phys.less.weight <- phys.less %>% 
  group_by(Genotype,Species,Treatment,Timepoint) %>%
  summarize_at(c('bw_t3t2','bw_t4t3'),mean,na.rm=T)
phys.less.weight.df <- data.frame(phys.less.weight)

## Reshape to long format
phys.less.weight.time <- reshape(phys.less.weight.df,
                               varying=c('bw_t3t2','bw_t4t3'),
                               v.names='weightChan',
                               timevar='sampling',
                               direction="long")
phys.less.weight.time$sampling <- as.factor(phys.less.weight.time$sampling)

## Create combined factor for Treatment and Timepoint = TreatTime
phys.less.weight.time$TreatTime <- paste(phys.less.weight.time$Treatment,
                                       ifelse(phys.less.weight.time$sampling == 1, "Pre-heat", "Post-heat"),
                                       sep = "_")
#convert to factors
phys.less.weight.time$TreatTime <- as.factor(phys.less.weight.time$TreatTime)
phys.less.weight.time$Species <- as.factor(phys.less.weight.time$Species)

## Summarize data with standard error
phys.less.weight.time.se <- summarySE(data=phys.less.weight.time,
                                    measurevar="weightChan",
                                    groupvars=c("TreatTime", "Treatment","Species","sampling"))
#reorder variables for figures
f=c('Control_Pre-heat','Control_Post-heat', 'DTV_Pre-heat', 'DTV_Post-heat')
phys.less.weight.time <- within(phys.less.weight.time, TreatTime<- factor(TreatTime, levels=f))


##dot & error bars
bw.both <- ggplot(phys.less.weight.time,aes(x=sampling,y=weightChan,shape=Treatment,fill=TreatTime))+
  #scale_y_continuous(expand = c(0, 0), limits = c(80, 25)) + #this is set for our data, may need to adjust if using with different data
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1), alpha = 0.3, size = 4) +
  # geom_errorbar(stat = 'summary',
  #               position = position_dodge(width = 0.5), width = 0.2) +
  geom_point(data = phys.less.weight.time.se,
             aes(x = sampling, y = weightChan),
             position = position_dodge(width = 0.5), size = 5) +
  geom_errorbar(stat = 'summary',position = position_dodge(width = 0.5), width = 0.2) +  
  facet_grid2(~Species,
              strip = strip_themed(background_x = elem_list_rect(fill = c("#8B7355","#CDAB7D")),
                                   text_x = elem_list_text(size = 14, face = "bold", color = "white")))+
  scale_fill_manual(name = "Treatment x Time", 
                    values = c("Control_Pre-heat" = "#2C2673",
                               "Control_Post-heat" = "#9187FF",
                               "DTV_Pre-heat" = "#E24A13",
                               "DTV_Post-heat" = "#FFA45E")) +
  scale_shape_manual(name = "Treatment",
                     values = c("Control" = 22,
                                "DTV" = 23)) +
  theme_bw() +
  theme(legend.title=element_text(size=14, face = "bold", color = "black"), 
        legend.text=element_text(size=14, face = "bold", color="black"),
        axis.text = element_text(color = "black", size = 14, face = "bold"),
        axis.title.y = element_text(size = 16, color = "black", face = "bold"),
        axis.ticks = element_line(color = "black"))+
  xlab("") +
  ylab("∆ Weight (%)") +
  ggtitle("") +
  scale_x_discrete(labels = c("Preheat", "Postheat")) 
bw.both

### Fig S4 % Change in weight all time points ####
phys.less2.weight <- phys.less %>% 
  group_by(Genotype,Species,Treatment,Timepoint) %>%
  summarize_at(c("bw_t1t0", "bw_t2t1", 'bw_t3t2','bw_t4t3'),mean,na.rm=T)
phys.less2.weight.df <- data.frame(phys.less2.weight)

## Reshape to long format
phys.less2.weight.time <- reshape(phys.less2.weight.df,
                                 varying=c("bw_t1t0", "bw_t2t1", 'bw_t3t2','bw_t4t3'),
                                 v.names='weightChan',
                                 timevar='sampling',
                                 direction="long")
phys.less2.weight.time$sampling <- as.factor(phys.less2.weight.time$sampling)

## Create combined factor for Treatment and Timepoint = TreatTime
phys.less2.weight.time$TreatTime <- paste(phys.less2.weight.time$Treatment,
                                         ifelse(phys.less2.weight.time$sampling == 1, "Pre-heat", "Post-heat"),
                                         sep = "_")
#convert to factors
phys.less2.weight.time$TreatTime <- as.factor(phys.less2.weight.time$TreatTime)
phys.less2.weight.time$Species <- as.factor(phys.less2.weight.time$Species)

## Summarize data with standard error
phys.less2.weight.time.se <- summarySE(data=phys.less2.weight.time,
                                      measurevar="weightChan",
                                      groupvars=c("TreatTime", "Treatment","Species","sampling"))


##dot & error bars
bw.both <- ggplot(phys.less2.weight.time,aes(x=sampling,y=weightChan,shape=Treatment,fill=Treatment))+
  #scale_y_continuous(expand = c(0, 0), limits = c(80, 25)) + #this is set for our data, may need to adjust if using with different data
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1), alpha = 0.3, size = 4) +
  # geom_errorbar(stat = 'summary',
  #               position = position_dodge(width = 0.5), width = 0.2) +
  geom_point(data = phys.less2.weight.time.se,
             aes(x = sampling, y = weightChan),
             position = position_dodge(width = 0.5), size = 5) +
  geom_errorbar(stat = 'summary',position = position_dodge(width = 0.5), width = 0.2) +  
  facet_grid2(~Species,
              strip = strip_themed(background_x = elem_list_rect(fill = c("#8B7355","#CDAB7D")),
                                   text_x = elem_list_text(size = 14, face = "bold", color = "white")))+
  scale_fill_manual(name = "Treatment", 
                    values = c("Control" = "#2C2673",
                               "DTV" = "#E24A13")) +
  scale_shape_manual(name = "Treatment",
                     values = c("Control" = 22,
                                "DTV" = 23)) +
  theme_bw() +
  theme(legend.title=element_text(size=14, face = "bold", color = "black"), 
        legend.text=element_text(size=14, face = "bold", color="black"),
        axis.text = element_text(color = "black", size = 14, face = "bold"),
        axis.title.y = element_text(size = 16, color = "black", face = "bold"),
        axis.title.x = element_text(size = 16, color = "black", face = "bold"),
        axis.ticks = element_line(color = "black"))+
  xlab("Day") +
  ylab("∆ Weight (%)") +
  ggtitle("") +
  scale_x_discrete(labels = c("30", "60","90", "115")) 
bw.both


#### Buoyant Weight Stats ####
### preheat and postheat ####
##favs only
phys.less.weight.time.fav <- subset(phys.less.weight.time,Species=="O. faveolata")

#Check assumptions for LME
# Homogeneity of variances
leveneTest(weightChan~Treatment,data=phys.less.weight.time.fav) #ns
histogram(phys.less.weight.time.fav$weightChan) #normal but likely an outlier, transform?
# Normality of residuals
shapiro.test(phys.less.weight.time.fav$weightChan) #barely sig - checked check_model results and decide to trust lme robustness


# control for Genotype is a random effect
lmer.weightChan <- lmer(weightChan ~ Treatment * sampling + (1 | Genotype), data = phys.less.weight.time.fav)


check_model(lmer.weightChan) # check assumptions more thoroughly, looks good!

# Summary of the model
summary(lmer.weightChan)

# Type III ANOVA table (if needed)
Anova(lmer.weightChan, type = 3)
# Analysis of Deviance Table (Type III tests)
# 
# Response: weightChan    ###significant but led by one sample at the preheat time point
# Chisq Df Pr(>Chisq)  
# (Intercept)        2.0489  1    0.15232  
# Treatment          4.6297  1    0.03142 *
# sampling           0.0468  1    0.82873  
# Treatment:sampling 4.7428  1    0.02942 *


#Pairwise comparisons
emmeans(lmer.weightChan, pairwise ~ Treatment * sampling, adjust = "tukey")
#all n.s.


##franks
phys.less.weight.time.fra <- subset(phys.less.weight.time,Species=="O. franksi")

leveneTest(weightChan~Treatment,data=phys.less2.weight.time.fra) #ns
shapiro.test(phys.less.weight.time.fra$weightChan)  #barely sig - checked check_model results and decide to trust lme robustness

# Genotype is a random effect
lmer.weightChan <- lmer(weightChan ~ Treatment * sampling + (1 | Genotype), data = phys.less.weight.time.fra)


check_model(lmer.weightChan) # check assumptions more thoroughly, looks pretty okay

# Summary of the model
summary(lmer.weightChan)

# Type III ANOVA table (if needed)
Anova(lmer.weightChan, type = 3)
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: weightChan
# Chisq Df Pr(>Chisq)   
# (Intercept)        9.1146  1   0.002536 **
#   Treatment          1.5317  1   0.215862   
#   sampling           5.8524  1   0.015556 * 
#   Treatment:sampling 0.6303  1   0.427250 


#Pairwise comparisons
emmeans(lmer.weightChan, pairwise ~ Treatment * sampling, adjust = "tukey")
#all n.s.



### All Buoyant weight T0-T90 time points ####
##favs only
phys.less2.weight.time.fav <- subset(phys.less2.weight.time,Species=="O. faveolata")
phys.less2.weight.time.fav <- subset(phys.less2.weight.time.fav,sampling!="4") #remove heat time point for this analysis

leveneTest(weightChan~Treatment,data=phys.less2.weight.time.fav) #ns
shapiro.test(phys.less2.weight.time.fav$weightChan) #ns

# control for Genotype is a random effect
lmer.weightChan <- lmer(weightChan ~ Treatment * sampling + (1 | Genotype), data = phys.less2.weight.time.fav)

                      
check_model(lmer.weightChan) # check assumptions more thoroughly, looks good!

# Summary of the model
summary(lmer.weightChan)

# Type III ANOVA table (if needed)
Anova(lmer.weightChan, type = 3)
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: weightChan
# Chisq Df Pr(>Chisq)    
# (Intercept)        25.8070  1  3.773e-07 ***
#   Treatment           0.0016  1   0.967823    
#   sampling           17.7051  2   0.000143 ***
#   Treatment:sampling  3.0810  2   0.214274   


#Pairwise comparisons
emmeans(lmer.weightChan, pairwise ~ Treatment * sampling, adjust = "tukey")
# $contrasts
# contrast                              estimate   SE df t.ratio p.value
# Control sampling1 - DTV sampling1       0.0206 0.51 15   0.040  1.0000
# Control sampling1 - Control sampling2   0.3945 0.51 15   0.774  0.9678
# Control sampling1 - DTV sampling2       0.5832 0.51 15   1.145  0.8552
# Control sampling1 - Control sampling3   2.0223 0.51 15   3.969  0.0128 **
# Control sampling1 - DTV sampling3       1.0413 0.51 15   2.044  0.3646
# DTV sampling1 - Control sampling2       0.3740 0.51 15   0.734  0.9744
# DTV sampling1 - DTV sampling2           0.5627 0.51 15   1.104  0.8722
# DTV sampling1 - Control sampling3       2.0018 0.51 15   3.929  0.0138
# DTV sampling1 - DTV sampling3           1.0207 0.51 15   2.003  0.3847
# Control sampling2 - DTV sampling2       0.1887 0.51 15   0.370  0.9989
# Control sampling2 - Control sampling3   1.6278 0.51 15   3.195  0.0553
# Control sampling2 - DTV sampling3       0.6467 0.51 15   1.269  0.7965
# DTV sampling2 - Control sampling3       1.4391 0.51 15   2.824  0.1074
# DTV sampling2 - DTV sampling3           0.4581 0.51 15   0.899  0.9409
# Control sampling3 - DTV sampling3      -0.9811 0.51 15  -1.925  0.4250



##franks
phys.less2.weight.time.fra <- subset(phys.less2.weight.time,Species=="O. franksi")
phys.less2.weight.time.fra <- subset(phys.less2.weight.time.fra,sampling!="4")

leveneTest(weightChan~Treatment,data=phys.less2.weight.time.fra) #ns
shapiro.test(phys.less2.weight.time.fra$weightChan) #ns

# Genotype is a random effect
lmer.weightChan <- lmer(weightChan ~ Treatment * sampling + (1 | Genotype), data = phys.less2.weight.time.fra)


check_model(lmer.weightChan) # check assumptions more thoroughly, looks pretty okay

# Summary of the model
summary(lmer.weightChan)

# Type III ANOVA table (if needed)
Anova(lmer.weightChan, type = 3)
# Analysis of Deviance Table (Type III tests)
# 
# Response: weightChan
# Chisq Df Pr(>Chisq)    
# (Intercept)        32.0282  1  1.520e-08 ***
#   Treatment           0.0005  1     0.9813    
#   sampling           23.7161  2  7.081e-06 ***
#   Treatment:sampling  2.2490  2     0.3248    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#Pairwise comparisons
emmeans(lmer.weightChan, pairwise ~ Treatment * sampling, adjust = "tukey")
# $contrasts
# contrast                              estimate   SE df t.ratio p.value
# Control sampling1 - DTV sampling1     -0.00585 0.25 15  -0.023  1.0000
# Control sampling1 - Control sampling2  0.45472 0.25 15   1.820  0.4828
# Control sampling1 - DTV sampling2      0.33033 0.25 15   1.322  0.7691
# Control sampling1 - Control sampling3  1.20478 0.25 15   4.822  0.0025 **
# Control sampling1 - DTV sampling3      1.58696 0.25 15   6.351  0.0002
# DTV sampling1 - Control sampling2      0.46057 0.25 15   1.843  0.4697
# DTV sampling1 - DTV sampling2          0.33618 0.25 15   1.345  0.7566
# DTV sampling1 - Control sampling3      1.21063 0.25 15   4.845  0.0024
# DTV sampling1 - DTV sampling3          1.59280 0.25 15   6.375  0.0002 **
# Control sampling2 - DTV sampling2     -0.12439 0.25 15  -0.498  0.9955
# Control sampling2 - Control sampling3  0.75006 0.25 15   3.002  0.0784
# Control sampling2 - DTV sampling3      1.13224 0.25 15   4.532  0.0044
# DTV sampling2 - Control sampling3      0.87445 0.25 15   3.500  0.0313
# DTV sampling2 - DTV sampling3          1.25662 0.25 15   5.029  0.0017 **
# Control sampling3 - DTV sampling3      0.38217 0.25 15   1.530  0.6522

##### Figs 2c and S6 Respirometry preheat to postheat ###
#### Gross Photosynthesis ####
#TreatTime
phys.less <- phys.all

## Group by genotype, species, treatment, and timepoint and average PAM values
phys.less.geno <- phys.less %>%
  group_by(Genotype, Species, Treatment, Timepoint) %>%
  summarize_at(c('gross_photo_t3', 'gross_photo_t4'), mean, na.rm = TRUE)
phys.less.geno.df <- data.frame(phys.less.geno)

## Reshape to long format
phys.less.geno.time <- reshape(phys.less.geno.df,
                               varying = c('gross_photo_t3', 'gross_photo_t4'),
                               v.names = 'Pgross',
                               timevar = 'sampling',
                               direction = "long")
phys.less.geno.time <- phys.less.geno.time %>%
  filter(!is.na(Pgross) & !is.nan(Pgross))
phys.less.geno.time$sampling <- as.factor(phys.less.geno.time$sampling)

## Create combined factor for Treatment and Timepoint
phys.less.geno.time$TreatTime <- paste(phys.less.geno.time$Treatment,
                                       ifelse(phys.less.geno.time$sampling == 1, "Pre-heat", "Post-heat"),
                                       sep = "_")
#convert to factors
phys.less.geno.time$TreatTime <- as.factor(phys.less.geno.time$TreatTime)
phys.less.geno.time$Species <- as.factor(phys.less.geno.time$Species)

## Summarize data with standard error
phys.less.geno.time.se <- summarySE(data = phys.less.geno.time,
                                    measurevar = "Pgross",
                                    groupvars = c("TreatTime", "Treatment", "Species", "sampling"))
#convert to factor
phys.less.geno.time.se$Species <- factor(phys.less.geno.time.se$Species)

#reorder variables for figures
f=c('Control_Pre-heat','Control_Post-heat', 'DTV_Pre-heat', 'DTV_Post-heat')
phys.less.geno.time <- within(phys.less.geno.time, TreatTime<- factor(TreatTime, levels=f))


## Plot
gg.pgross <- ggplot(phys.less.geno.time, aes(x = sampling, y = Pgross, shape = Treatment, fill = TreatTime)) +
  #scale_y_continuous(expand = c(0,0), limits = c(0.45,0.65)) + #this is set for our data, may need to adjust if using with different data
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1), alpha = 0.3, size = 4) +
  geom_errorbar(stat = 'summary',
                position = position_dodge(width = 0.5), width = 0.2) +
  geom_point(data = phys.less.geno.time.se,
             aes(x = sampling, y = Pgross),
             position = position_dodge(width = 0.5), size = 5) +
  facet_grid2(~Species,
              strip = strip_themed(background_x = elem_list_rect(fill = c("#8B7355","#CDAB7D")),
                                   text_x = elem_list_text(size = 14, face = "bold", color = "white")))+
  scale_fill_manual(name = "Treatment x Time", 
                    values = c("Control_Pre-heat" = "#2C2673",
                               "Control_Post-heat" = "#9187FF",
                               "DTV_Pre-heat" = "#E24A13",
                               "DTV_Post-heat" = "#FFA45E")) +
  scale_shape_manual(name = "Treatment",
                     values = c("Control" = 22,
                                "DTV" = 23)) +
  theme_bw() +
  theme(legend.title=element_text(size=14, face = "bold", color = "black"), 
        legend.text=element_text(size=14, face = "bold", color="black"),
        axis.text = element_text(color = "black", size = 14, face = "bold"),
        axis.title.y = element_text(size = 16, color = "black", face = "bold"),
        axis.ticks = element_line(color = "black"))+
  xlab("") +
  ylab("Gross Photosynthesis") +
  ggtitle("") +
  scale_x_discrete(labels = c("Preheat", "Postheat")) #+
#theme(axis.text.x = element_text(angle = 45, hjust = 1))

gg.pgross


#### Gross Photosynthesis stats  ####
#preheat to postheat
##both species
leveneTest(Pgross~Treatment,data=phys.less.geno.time) #ns
shapiro.test(phys.less.geno.time$Pgross) #ns

# Genotype is a random effect
lmer.pgross <- lmer(Pgross ~ Treatment * Species * sampling + (1 | Genotype), data = phys.less.geno.time)

# Summary of the model
summary(lmer.pgross)

# Type III ANOVA table (if needed)
Anova(lmer.pgross, type = 3)
#ns

##favs only
phys.less.geno.time.fav <- subset(phys.less.geno.time,Species=="O. faveolata")

#Check assumptions for LME
# Homogeneity of variances
leveneTest(Pgross~Treatment,data=phys.less.geno.time.fav) #ns
# Normality of residuals
shapiro.test(phys.less.geno.time.fav$Pgross) #ns

# control for Genotype is a random effect
lmer.pgross <- lmer(Pgross ~ Treatment * sampling + (1 | Genotype), data = phys.less.geno.time.fav)

check_model(lmer.pgross) # check assumptions more thoroughly, looks decent!

# Summary of the model
summary(lmer.pgross)

# Type III ANOVA table (if needed)
Anova(lmer.pgross, type = 3)
#ns

#Pairwise comparisons
emmeans(lmer.pgross, pairwise  ~ Treatment * sampling, adjust = "tukey")
#ns

##franks
phys.less.geno.time.fra <- subset(phys.less.geno.time,Species=="O. franksi")

leveneTest(Pgross~Treatment,data=phys.less.geno.time.fra) #ns
shapiro.test(phys.less.geno.time.fra$Pgross) #ns

# Genotype is a random effect
lmer.pgross <- lmer(Pgross ~ Treatment * sampling + (1 | Genotype), data = phys.less.geno.time.fra)

check_model(lmer.pgross) # check assumptions more thoroughly, looks decent!

# Summary of the model
summary(lmer.pgross)

# Type III ANOVA table (if needed)
Anova(lmer.pgross, type = 3)
#ns, timepoint almost id


#Pairwise comparisons
# emmeans(lmer.pgross, pairwise ~ Treatment * sampling, adjust = "tukey")
# contrast                              estimate    SE df t.ratio p.value
# Control sampling1 - DTV sampling1        0.529 0.889  9   0.595  0.9311
# Control sampling1 - Control sampling2    1.630 0.889  9   1.832  0.3199
# Control sampling1 - DTV sampling2        3.101 0.889  9   3.487  0.0289 ***
# DTV sampling1 - Control sampling2        1.100 0.889  9   1.237  0.6208
# DTV sampling1 - DTV sampling2            2.572 0.889  9   2.892  0.0706 ~
# Control sampling2 - DTV sampling2        1.472 0.889  9   1.655  0.3983


### dark respiration ####
#TreatTime
phys.less <- phys.all

## Group by genotype, species, treatment, and timepoint and average PAM values
phys.less.geno <- phys.less %>%
  group_by(Genotype, Species, Treatment, Timepoint) %>%
  summarize_at(c('resp_t3', 'resp_t4'), mean, na.rm = TRUE)
phys.less.geno.df <- data.frame(phys.less.geno)

## Reshape to long format
phys.less.geno.time <- reshape(phys.less.geno.df,
                               varying = c('resp_t3', 'resp_t4'),
                               v.names = 'Resp',
                               timevar = 'sampling',
                               direction = "long")
phys.less.geno.time <- phys.less.geno.time %>%
  filter(!is.na(Resp) & !is.nan(Resp))
phys.less.geno.time$sampling <- as.factor(phys.less.geno.time$sampling)

## Create combined factor for Treatment and Timepoint
phys.less.geno.time$TreatTime <- paste(phys.less.geno.time$Treatment,
                                       ifelse(phys.less.geno.time$sampling == 1, "Pre-heat", "Post-heat"),
                                       sep = "_")
#convert to factors
phys.less.geno.time$TreatTime <- as.factor(phys.less.geno.time$TreatTime)
phys.less.geno.time$Species <- as.factor(phys.less.geno.time$Species)

## Summarize data with standard error
phys.less.geno.time.se <- summarySE(data = phys.less.geno.time,
                                    measurevar = "Resp",
                                    groupvars = c("TreatTime", "Treatment", "Species", "sampling"))
#convert to factor
phys.less.geno.time.se$Species <- factor(phys.less.geno.time.se$Species)

#reorder variables for figures
f=c('Control_Pre-heat','Control_Post-heat', 'DTV_Pre-heat', 'DTV_Post-heat')
phys.less.geno.time <- within(phys.less.geno.time, TreatTime<- factor(TreatTime, levels=f))


## Plot
gg.resp <- ggplot(phys.less.geno.time, aes(x = sampling, y = Resp, shape = Treatment, fill = TreatTime)) +
  #scale_y_continuous(expand = c(0,0), limits = c(0.45,0.65)) + #this is set for our data, may need to adjust if using with different data
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1), alpha = 0.3, size = 4) +
  geom_errorbar(stat = 'summary',
                position = position_dodge(width = 0.5), width = 0.2) +
  geom_point(data = phys.less.geno.time.se,
             aes(x = sampling, y = Resp),
             position = position_dodge(width = 0.5), size = 5) +
  facet_grid2(~Species,
              strip = strip_themed(background_x = elem_list_rect(fill = c("#8B7355","#CDAB7D")),
                                   text_x = elem_list_text(size = 14, face = "bold", color = "white")))+
  scale_fill_manual(name = "Treatment x Time", 
                    values = c("Control_Pre-heat" = "#2C2673",
                               "Control_Post-heat" = "#9187FF",
                               "DTV_Pre-heat" = "#E24A13",
                               "DTV_Post-heat" = "#FFA45E")) +
  scale_shape_manual(name = "Treatment",
                     values = c("Control" = 22,
                                "DTV" = 23)) +
  theme_bw() +
  theme(legend.title=element_text(size=14, face = "bold", color = "black"), 
        legend.text=element_text(size=14, face = "bold", color="black"),
        axis.text = element_text(color = "black", size = 14, face = "bold"),
        axis.title.y = element_text(size = 16, color = "black", face = "bold"),
        axis.ticks = element_line(color = "black"))+
  xlab("") +
  ylab("Dark Respiration (nmol/min⁻¹mm⁻²)") +
  ggtitle("") +
  scale_x_discrete(labels = c("Preheat", "Postheat")) #+
#theme(axis.text.x = element_text(angle = 45, hjust = 1))

gg.resp


#### Dark Respiration stats  ####
#preheat to postheat
##both species
leveneTest(Resp~Treatment,data=phys.less.geno.time) #ns
shapiro.test(phys.less.geno.time$Resp) #sig
histogram(phys.less.geno.time$Resp) #looks a lil right skewed but not bad, lmer is robust to this

# Genotype is a random effect
lmer.Resp <- lmer(Resp ~ Treatment * Species * sampling + (1 | Genotype), data = phys.less.geno.time)

# Summary of the model
summary(lmer.Resp)

# Type III ANOVA table (if needed)
Anova(lmer.Resp, type = 3)
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: Resp
# Chisq Df Pr(>Chisq)    
# (Intercept)                37.4439  1  9.408e-10 ***
#   Treatment                   0.5101  1    0.47508    
# Species                     0.7189  1    0.39649    
# sampling                    5.6806  1    0.01715 *  
#   Treatment:Species           0.1480  1    0.70050    
# Treatment:sampling          0.8144  1    0.36681    
# Species:sampling            2.2659  1    0.13225    
# Treatment:Species:sampling  0.1780  1    0.67306    

##favs only
phys.less.geno.time.fav <- subset(phys.less.geno.time,Species=="O. faveolata")

#Check assumptions for LME
# Homogeneity of variances
leveneTest(Resp~Treatment,data=phys.less.geno.time.fav) #ns
# Normality of residuals
shapiro.test(phys.less.geno.time.fav$Resp) #sig
histogram(phys.less.geno.time.fav$Resp) #looks right skewed, try transforming

transResp <- log(phys.less.geno.time.fav$Resp + abs(min(phys.less.geno.time.fav$Resp, na.rm=TRUE)) + 1)
shapiro.test(transResp) #still sig, lmer is robust to this
histogram(transResp) #looks better, use this for anakysis

# control for Genotype is a random effect
lmer.Resp <- lmer(transResp ~ Treatment * sampling + (1 | Genotype), data = phys.less.geno.time.fav)

check_model(lmer.Resp) # check assumptions more thoroughly, looks decent!

# Summary of the model
summary(lmer.Resp)

# Type III ANOVA table (if needed)
Anova(lmer.Resp, type = 3)
#Response: transResp
# Chisq Df Pr(>Chisq)    
# (Intercept)        902.5176  1    < 2e-16 ***
# Treatment            0.8404  1    0.35927    
# sampling             6.3540  1    0.01171 *  
# Treatment:sampling   0.6317  1    0.42673  

#Pairwise comparisons
emmeans(lmer.Resp, pairwise  ~ Treatment * sampling, adjust = "tukey")
# $contrasts
# contrast                              estimate    SE df t.ratio p.value
# Control sampling1 - DTV sampling1       0.0706 0.077  9   0.917  0.7968
# Control sampling1 - Control sampling2  -0.1940 0.077  9  -2.521  0.1228
# Control sampling1 - DTV sampling2      -0.0369 0.077  9  -0.480  0.9617
# DTV sampling1 - Control sampling2      -0.2646 0.077  9  -3.437  0.0311 **
# DTV sampling1 - DTV sampling2          -0.1075 0.077  9  -1.397  0.5314
# Control sampling2 - DTV sampling2       0.1571 0.077  9   2.041  0.2429

##franks
phys.less.geno.time.fra <- subset(phys.less.geno.time,Species=="O. franksi")

leveneTest(Resp~Treatment,data=phys.less.geno.time.fra) #ns
shapiro.test(phys.less.geno.time.fra$Resp) #ns

# Genotype is a random effect
lmer.Resp <- lmer(Resp ~ Treatment * sampling + (1 | Genotype), data = phys.less.geno.time.fra)

check_model(lmer.Resp) # check assumptions more thoroughly, looks good!

# Summary of the model
summary(lmer.Resp)

# Type III ANOVA table (if needed)
Anova(lmer.Resp, type = 3)
#ns

#Pairwise comparisons
emmeans(lmer.Resp, pairwise  ~ Treatment * sampling, adjust = "tukey")
#ns



#### P:R ####
#TreatTime
phys.less <- phys.all


## Group by genotype, species, treatment, and timepoint and average PAM values
# phys.less.geno <- phys.less %>%
#   group_by(Genotype, Species, Treatment, Timepoint) %>%
#   summarize_at(c('PR_t3', 'PR_t4'), mean, na.rm = TRUE)
phys.less.geno <- phys.less %>%
  group_by(Genotype, Species, Treatment, Timepoint) %>%
  dplyr::summarize(across(c(PR_t3, PR_t4), ~ mean(.x, na.rm = TRUE)), .groups = "drop")



phys.less.geno.df <- data.frame(phys.less.geno)


## Reshape to long format
phys.less.geno.time <- reshape(phys.less.geno.df,
                               varying = c('PR_t3', 'PR_t4'),
                               v.names = 'PR',
                               timevar = 'sampling',
                               direction = "long")

phys.less.geno.time <- phys.less.geno.time %>%
  filter(!is.na(PR) & !is.nan(PR))
phys.less.geno.time$sampling <- as.factor(phys.less.geno.time$sampling)

## Create combined factor for Treatment and Timepoint
phys.less.geno.time$TreatTime <- paste(phys.less.geno.time$Treatment,
                                       ifelse(phys.less.geno.time$sampling == 1, "Pre-heat", "Post-heat"),
                                       sep = "_")
#convert to factors
phys.less.geno.time$TreatTime <- as.factor(phys.less.geno.time$TreatTime)
phys.less.geno.time$Species <- as.factor(phys.less.geno.time$Species)
phys.less.geno.time$PR[is.nan(phys.less.geno.time$PR)] <- NA
phys.less.geno.time$PR <- as.numeric(phys.less.geno.time$PR)

## Summarize data with standard error
phys.less.geno.time.se <- summarySE(data = phys.less.geno.time,
                                    measurevar = "PR",
                                    groupvars = c("TreatTime", "Treatment", "Species", "sampling"))
#convert to factor
phys.less.geno.time.se$Species <- factor(phys.less.geno.time.se$Species)

#reorder variables for figures
f=c('Control_Pre-heat','Control_Post-heat', 'DTV_Pre-heat', 'DTV_Post-heat')
phys.less.geno.time <- within(phys.less.geno.time, TreatTime<- factor(TreatTime, levels=f))


## Plot
gg.resp <- ggplot(phys.less.geno.time, aes(x = sampling, y = PR, shape = Treatment, fill = TreatTime)) +
  #scale_y_continuous(expand = c(0,0), limits = c(0.45,0.65)) + #this is set for our data, may need to adjust if using with different data
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1), alpha = 0.3, size = 4) +
  geom_errorbar(stat = 'summary',
                position = position_dodge(width = 0.5), width = 0.2) +
  geom_point(data = phys.less.geno.time.se,
             aes(x = sampling, y = PR),
             position = position_dodge(width = 0.5), size = 5) +
  facet_grid2(~Species,
              strip = strip_themed(background_x = elem_list_rect(fill = c("#8B7355","#CDAB7D")),
                                   text_x = elem_list_text(size = 14, face = "bold", color = "white")))+
  scale_fill_manual(name = "Treatment x Time", 
                    values = c("Control_Pre-heat" = "#2C2673",
                               "Control_Post-heat" = "#9187FF",
                               "DTV_Pre-heat" = "#E24A13",
                               "DTV_Post-heat" = "#FFA45E")) +
  scale_shape_manual(name = "Treatment",
                     values = c("Control" = 22,
                                "DTV" = 23)) +
  theme_bw() +
  theme(legend.title=element_text(size=14, face = "bold", color = "black"), 
        legend.text=element_text(size=14, face = "bold", color="black"),
        axis.text = element_text(color = "black", size = 14, face = "bold"),
        axis.title.y = element_text(size = 16, color = "black", face = "bold"),
        axis.ticks = element_line(color = "black"))+
  xlab("") +
  ylab("P/R") +
  ggtitle("") +
  scale_x_discrete(labels = c("Preheat", "Postheat")) #+
#theme(axis.text.x = element_text(angle = 45, hjust = 1))

gg.resp


#### P/R stats  ####
#preheat to postheat
##both species
leveneTest(PR~Treatment,data=phys.less.geno.time) #ns
shapiro.test(phys.less.geno.time$PR) #sig, but close, maybe dependon robustness of lmer
histogram(phys.less.geno.time$PR) #looks left skewed, going to try transforming

transPR <- log(phys.less.geno.time$PR + abs(min(phys.less.geno.time$PR, na.rm=TRUE)) + 1)
shapiro.test(transPR) #worse than before - try different transformation
histogram(transPR) #still left skewed, but worse

transPR <- (phys.less.geno.time$PR)^2
shapiro.test(transPR) #still sig but barely
histogram(transPR) #looks better, going to use this for analysis


# Genotype is a random effect
lmer.PR <- lmer(transPR ~ Treatment * Species * sampling + (1 | Genotype), data = phys.less.geno.time)

# Summary of the model
summary(lmer.PR)

# Type III ANOVA table (if needed)
Anova(lmer.PR, type = 3)
#ns


##favs only
phys.less.geno.time.fav <- subset(phys.less.geno.time,Species=="O. faveolata")

#Check assumptions for LME
# Homogeneity of variances
leveneTest(PR~Treatment,data=phys.less.geno.time.fav) #ns
# Normality of residuals
shapiro.test(phys.less.geno.time.fav$PR) #ns

# control for Genotype is a random effect
lmer.PR <- lmer(PR ~ Treatment * sampling + (1 | Genotype), data = phys.less.geno.time.fav)

check_model(lmer.PR) # check assumptions more thoroughly, looks decent!

# Summary of the model
summary(lmer.PR)

# Type III ANOVA table (if needed)
Anova(lmer.PR, type = 3)
#Response: PR
# Chisq Df Pr(>Chisq)    
# (Intercept)        185.0564  1  < 2.2e-16 ***
# Treatment            0.0064  1     0.9363    
# sampling            31.1933  1  2.336e-08 ***
# Treatment:sampling   0.8392  1     0.3596 


#Pairwise comparisons
emmeans(lmer.PR, pairwise  ~ Treatment * sampling, adjust = "tukey")
# $contrasts
# contrast                              estimate    SE df t.ratio p.value
# Control sampling1 - DTV sampling1       0.0159 0.199  9   0.080  0.9998
# Control sampling1 - Control sampling2   1.1110 0.199  9   5.585  0.0016 **
# Control sampling1 - DTV sampling2       0.8692 0.199  9   4.369  0.0080 **
# DTV sampling1 - Control sampling2       1.0951 0.199  9   5.505  0.0017 **
# DTV sampling1 - DTV sampling2           0.8533 0.199  9   4.290  0.0090 **
# Control sampling2 - DTV sampling2      -0.2418 0.199  9  -1.216  0.6330

##franks
phys.less.geno.time.fra <- subset(phys.less.geno.time,Species=="O. franksi")

leveneTest(PR~Treatment,data=phys.less.geno.time.fra) #ns
shapiro.test(phys.less.geno.time.fra$PR) #sig
histogram(phys.less.geno.time.fra$PR) #looks left skewed, try transforming

transPR <- (phys.less.geno.time.fra$PR)^2
shapiro.test(transPR) #ns
histogram(transPR) #looks better, use this for analysis

# Genotype is a random effect
lmer.PR <- lmer(transPR ~ Treatment * sampling + (1 | Genotype), data = phys.less.geno.time.fra)

check_model(lmer.PR) # check assumptions more thoroughly, much better!

# Summary of the model
summary(lmer.PR)

# Type III ANOVA table (if needed)
Anova(lmer.PR, type = 3)
#ns, but timepoint almost is

#Pairwise comparisons
emmeans(lmer.PR, pairwise  ~ Treatment * sampling, adjust = "tukey")
#ns



