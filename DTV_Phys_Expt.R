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

setwd("~/path")
phys.all <- read.csv("THERMVAR_MAIN_METADATA_Apr25_clonesfixed.csv")

str(phys.all)
phys.all$Timepoint <- factor(phys.all$Timepoint,levels=c("Preheat","Postheat"))


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


#### Color stats ####

### Red Intensity Only ### this is what we used in the manuscript

##both species
phys.less.geno.time$Treatment <- as.factor(phys.less.geno.time$Treatment)

leveneTest(red_chan~Treatment,data=phys.less.geno.time) #sig
phys.less.geno.time$log_redchan <- log(phys.less.geno.time$red_chan) #transform - still not equal
phys.less.geno.time$log_redchan <- sqrt(phys.less.geno.time$red_chan) #transform - still not equal, run with variance unequal in LME 
leveneTest(log_redchan~Treatment,data=phys.less.geno.time) #sig
shapiro.test(phys.less.geno.time$red_chan) #ns
#bestNormalize(phys.less.geno.time$red_chan,loo=T)
#red.ord <- orderNorm(phys.less.geno.timeo$red_chan)
#phys.less.geno.time$red.ord <- red.ord[["x.t"]]
#shapiro.test(phys.less.geno.time$red.ord) #ns

# Genotype is a random effect
lme.red <- lme(
  red_chan ~ Treatment * sampling * Species,                # Fixed effects
  random = ~1 | Genotype,                         # Random intercepts for repeated measures
  weights = varIdent(form = ~1 | Treatment),      # Allows unequal variances by Treatment group
  data = phys.less.geno.time
)


# Summary of the model w/ ANOVA like summaries
summary(lme.red)
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


#due to repeated measures ANOVA is not quite right
# aov.red <- aov(red.ord~Treatment*Species,data=phys.t4.geno)
# summary(aov.red) 
#                   Df Sum Sq Mean Sq F value Pr(>F)  
# Treatment          1  5.635   5.635   7.057 0.0172 *
# Species            1  0.319   0.319   0.400 0.5363  
# Treatment:Species  1  0.041   0.041   0.051 0.8244  
# Residuals         16 12.776   0.799                 
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


##franks
phys.less.geno.time.fra <- subset(phys.less.geno.time,Species=="O. franksi")
str(phys.less.geno.time.fra)
phys.less.geno.time.fra$Treatment <- as.factor(phys.less.geno.time.fra$Treatment)
str(phys.less.geno.time.fra)

leveneTest(red_chan~Treatment,data=phys.less.geno.time.fra) #ns
shapiro.test(phys.less.geno.time.fra$red_chan) #ns

# Genotype is a random effect
#lmer.red <- lmer(red_chan ~ Treatment * sampling + (1 | Genotype), data = phys.less.geno.time.fra) - could use but has high singularity, better to use lme
lme.red <- lme(
  red_chan ~ Treatment * sampling,                # Fixed effects
  random = ~1 | Genotype,                         # Random intercepts for repeated measures
 # weights = varIdent(form = ~1 | Treatment),      # Allows unequal variances by Treatment group
  data = phys.less.geno.time.fra
)

check_model(lme.red) # check assumptions more thoroughly, assumptions are reasonably well met.

# Summary of the model
summary(lme.red)
# Fixed effects:  red_chan ~ Treatment * sampling 
# Value Std.Error DF   t-value p-value
# (Intercept)            57.04583  2.776544  9 20.545623  0.0000
# TreatmentDTV           -4.21250  3.926627  9 -1.072804  0.3113
# sampling2               5.58750  3.926627  9  1.422977  0.1885
# TreatmentDTV:sampling2 10.00000  5.553089  9  1.800800  0.1053


#Pairwise comparisons
emmeans(lme.red, pairwise ~ Treatment * sampling, adjust = "tukey")
# $contrasts
# contrast                              estimate   SE df t.ratio p.value
# Control sampling1 - DTV sampling1         4.21 3.93  9   1.073  0.7135
# Control sampling1 - Control sampling2    -5.59 3.93  9  -1.423  0.5170
# Control sampling1 - DTV sampling2       -11.38 3.93  9  -2.897  0.0701
# DTV sampling1 - Control sampling2        -9.80 3.93  9  -2.496  0.1273
# DTV sampling1 - DTV sampling2           -15.59 3.93  9  -3.970  0.0142 *
# Control sampling2 - DTV sampling2        -5.79 3.93  9  -1.474  0.4896

##favs
phys.less.geno.time.fav <- subset(phys.less.geno.time,Species=="O. faveolata")
str(phys.less.geno.time.fav)
phys.less.geno.time.fav$Treatment <- as.factor(phys.less.geno.time.fav$Treatment)
str(phys.less.geno.time.fav)

leveneTest(red_chan~Treatment,data=phys.less.geno.time.fav) #ns
shapiro.test(phys.less.geno.time.fav$red_chan) #ns

# Genotype is a random effect
lme.red <- lme(
  red_chan ~ Treatment * sampling,                # Fixed effects
  random = ~1 | Genotype,                         # Random intercepts for repeated measures
  #weights = varIdent(form = ~1 | Treatment),      # Allows unequal variances by Treatment group
  data = phys.less.geno.time.fav
)

check_model(lme.red) # check assumptions more thoroughly, appropriate to proceed

# Summary of the model
summary(lme.red)
# Fixed effects:  red_chan ~ Treatment * sampling 
# Value Std.Error DF   t-value p-value
# (Intercept)             54.61250  3.518148  9 15.523082  0.0000 ***
# TreatmentDTV           -10.69583  2.994148  9 -3.572246  0.0060 **
# sampling2                3.57917  2.994148  9  1.195387  0.2625
# TreatmentDTV:sampling2   6.12917  4.234365  9  1.447482  0.1817



#Pairwise comparisons
emmeans(lme.red, pairwise ~ Treatment * sampling, adjust = "tukey")
# $contrasts
# contrast                              estimate   SE df t.ratio p.value
# Control sampling1 - DTV sampling1       10.696 2.99  9   3.572  0.0255 *
# Control sampling1 - Control sampling2   -3.579 2.99  9  -1.195  0.6445
# Control sampling1 - DTV sampling2        0.988 2.99  9   0.330  0.9868
# DTV sampling1 - Control sampling2      -14.275 2.99  9  -4.768  0.0046 **
# DTV sampling1 - DTV sampling2           -9.708 2.99  9  -3.242  0.0417 *
# Control sampling2 - DTV sampling2        4.567 2.99  9   1.525  0.4628






### All timepoints red chan

##both species
phys.less2.geno.time$Treatment <- as.factor(phys.less2.geno.time$Treatment)

leveneTest(red_chan~Treatment,data=phys.less2.geno.time) #ns
shapiro.test(phys.less2.geno.time$red_chan) #sig   
phys.less2.geno.time$log_redchan <- log(phys.less2.geno.time$red_chan) #transform - still not equal
shapiro.test(phys.less2.geno.time$log_redchan) #ns

# Genotype is a random effect
lme.red <- lme(
  log_redchan ~ Treatment * sampling * Species,                # Fixed effects
  random = ~1 | Genotype,                         # Random intercepts for repeated measures
  #weights = varIdent(form = ~1 | Treatment),      # Allows unequal variances by Treatment group
  data = phys.less2.geno.time
)

check_model(lme.red) # check assumptions more thoroughly, looks pretty good


# Summary of the model w/ ANOVA like summaries
Anova(lme.red, type = "III") 
# Analysis of Deviance Table (Type III tests)
# 
# Response: log_redchan
# Chisq Df Pr(>Chisq)    
# (Intercept)                2563.8905  1  < 2.2e-16 ***
#   Treatment                     0.0320  1     0.8581    
# sampling                     31.2839  3  7.407e-07 ***
#   Species                       1.3819  1     0.2398    
# Treatment:sampling            5.3475  3     0.1481    
# Treatment:Species             0.5769  1     0.4475    
# sampling:Species              4.0940  3     0.2515    
# Treatment:sampling:Species    2.3964  3     0.4943  



##franks
phys.less2.geno.time.fra <- subset(phys.less2.geno.time,Species=="O. franksi")
str(phys.less2.geno.time.fra)
#phys.int.geno.fra$Treatment <- as.factor(phys.int.geno.fra$Treatment)
#str(phys.int.geno.fra)

leveneTest(red_chan~Treatment,data=phys.less2.geno.time.fra) #ns
shapiro.test(phys.less2.geno.time.fra$red_chan) #ns

# Genotype is a random effect
lme.red <- lme(
  red_chan ~ Treatment * sampling,                # Fixed effects
  random = ~1 | Genotype,                         # Random intercepts for repeated measures
  #weights = varIdent(form = ~1 | Treatment),      # Allows unequal variances by Treatment group
  data = phys.less2.geno.time.fra
)

check_model(lme.red) # check assumptions more thoroughly, looks pretty good

# Summary of the model w/ ANOVA like summaries
Anova(lme.red, type = "III") 
# Analysis of Deviance Table (Type III tests)
# 
# Response: red_chan
#                       Chisq Df Pr(>Chisq)    
# (Intercept)        95.2113  1  < 2.2e-16 ***
#   Treatment           0.3436  1     0.5578    
# sampling           29.9946  3  1.384e-06 ***
#   Treatment:sampling  2.1638  3     0.5391  

#Pairwise comparisons
emmeans(lme.red, pairwise ~ Treatment * sampling, adjust = "tukey")
# $contrasts
# contrast                              estimate   SE df t.ratio p.value
# Control sampling1 - DTV sampling1        -3.38 5.77 21  -0.586  0.9987
# Control sampling1 - Control sampling2   -31.10 5.77 21  -5.389  0.0005 **
# Control sampling1 - DTV sampling2       -27.91 5.77 21  -4.836  0.0019
# Control sampling1 - Control sampling3   -12.18 5.77 21  -2.109  0.4406
# Control sampling1 - DTV sampling3        -7.96 5.77 21  -1.379  0.8565
# Control sampling1 - Control sampling4   -17.76 5.77 21  -3.077  0.0875
# Control sampling1 - DTV sampling4       -23.55 5.77 21  -4.080  0.0104
# DTV sampling1 - Control sampling2       -27.72 5.77 21  -4.802  0.0020
# DTV sampling1 - DTV sampling2           -24.53 5.77 21  -4.250  0.0071 **
# DTV sampling1 - Control sampling3        -8.79 5.77 21  -1.523  0.7871
# DTV sampling1 - DTV sampling3            -4.58 5.77 21  -0.793  0.9917
# DTV sampling1 - Control sampling4       -14.38 5.77 21  -2.491  0.2519
# DTV sampling1 - DTV sampling4           -20.17 5.77 21  -3.494  0.0373 *
# Control sampling2 - DTV sampling2         3.19 5.77 21   0.553  0.9991
# Control sampling2 - Control sampling3    18.93 5.77 21   3.279  0.0583 ~
# Control sampling2 - DTV sampling3        23.14 5.77 21   4.009  0.0122 *
# Control sampling2 - Control sampling4    13.34 5.77 21   2.311  0.3329
# Control sampling2 - DTV sampling4         7.55 5.77 21   1.309  0.8856
# DTV sampling2 - Control sampling3        15.74 5.77 21   2.726  0.1688
# DTV sampling2 - DTV sampling3            19.95 5.77 21   3.456  0.0404 *
# DTV sampling2 - Control sampling4        10.15 5.77 21   1.758  0.6521 -
# DTV sampling2 - DTV sampling4             4.36 5.77 21   0.756  0.9938
# Control sampling3 - DTV sampling3         4.21 5.77 21   0.730  0.9949
# Control sampling3 - Control sampling4    -5.59 5.77 21  -0.968  0.9743 -
# Control sampling3 - DTV sampling4       -11.38 5.77 21  -1.971  0.5223
# DTV sampling3 - Control sampling4        -9.80 5.77 21  -1.698  0.6886
# DTV sampling3 - DTV sampling4           -15.59 5.77 21  -2.700  0.1768 -
# Control sampling4 - DTV sampling4        -5.79 5.77 21  -1.003  0.9689

##favs
phys.less2.geno.time.fav <- subset(phys.less2.geno.time,Species=="O. faveolata")
str(phys.less2.geno.time.fav)


leveneTest(red_chan~Treatment,data=phys.less2.geno.time.fav) #ns
shapiro.test(phys.less2.geno.time.fav$red_chan) #sig
phys.less2.geno.time.fav$log_redchan <- log(phys.less2.geno.time.fav$red_chan) #transform - still not equal
shapiro.test(phys.less2.geno.time$log_redchan) #ns

# Genotype is a random effect
lme.red <- lme(
  log_redchan ~ Treatment * sampling,                # Fixed effects
  random = ~1 | Genotype,                         # Random intercepts for repeated measures
  #weights = varIdent(form = ~1 | Treatment),      # Allows unequal variances by Treatment group
  data = phys.less2.geno.time.fav
)

check_model(lme.red) # check assumptions more thoroughly, looks pretty good


# Summary of the model w/ ANOVA like summaries
Anova(lme.red, type = "III") 
# Analysis of Deviance Table (Type III tests)
# 
# Response: log_redchan
#                        Chisq Df Pr(>Chisq)    
# (Intercept)        2108.9566  1  < 2.2e-16 ***
#   Treatment             0.0372  1      0.847    
# sampling             36.4454  3  6.028e-08 ***
#   Treatment:sampling    6.2298  3      0.101   

#Pairwise comparisons
emmeans(lme.red, pairwise ~ Treatment * sampling, adjust = "tukey")
# $contrasts
# contrast                              estimate     SE df t.ratio p.value
# Control sampling1 - DTV sampling1       0.0140 0.0724 21   0.193  1.0000
# Control sampling1 - Control sampling2  -0.4052 0.0724 21  -5.597  0.0003 ***
# Control sampling1 - DTV sampling2      -0.4172 0.0724 21  -5.762  0.0002 ***
# Control sampling1 - Control sampling3  -0.0636 0.0724 21  -0.878  0.9850
# Control sampling1 - DTV sampling3       0.1574 0.0724 21   2.173  0.4047
# Control sampling1 - Control sampling4  -0.1296 0.0724 21  -1.790  0.6330
# Control sampling1 - DTV sampling4      -0.0449 0.0724 21  -0.620  0.9982
# DTV sampling1 - Control sampling2      -0.4192 0.0724 21  -5.790  0.0002 ***
# DTV sampling1 - DTV sampling2          -0.4311 0.0724 21  -5.955  0.0002 ***
# DTV sampling1 - Control sampling3      -0.0775 0.0724 21  -1.071  0.9561
# DTV sampling1 - DTV sampling3           0.1434 0.0724 21   1.980  0.5164
# DTV sampling1 - Control sampling4      -0.1436 0.0724 21  -1.983  0.5150
# DTV sampling1 - DTV sampling4          -0.0588 0.0724 21  -0.813  0.9904
# Control sampling2 - DTV sampling2      -0.0119 0.0724 21  -0.165  1.0000
# Control sampling2 - Control sampling3   0.3417 0.0724 21   4.719  0.0025 **
# Control sampling2 - DTV sampling3       0.5626 0.0724 21   7.770  <.0001 ***
# Control sampling2 - Control sampling4   0.2756 0.0724 21   3.807  0.0190 *
# Control sampling2 - DTV sampling4       0.3604 0.0724 21   4.977  0.0014 **
# DTV sampling2 - Control sampling3       0.3536 0.0724 21   4.884  0.0017 **
# DTV sampling2 - DTV sampling3           0.5745 0.0724 21   7.935  <.0001 ***
# DTV sampling2 - Control sampling4       0.2876 0.0724 21   3.972  0.0132
# DTV sampling2 - DTV sampling4           0.3723 0.0724 21   5.142  0.0009 ***
# Control sampling3 - DTV sampling3       0.2209 0.0724 21   3.051  0.0920
# Control sampling3 - Control sampling4  -0.0660 0.0724 21  -0.912  0.9815
# Control sampling3 - DTV sampling4       0.0187 0.0724 21   0.258  1.0000
# DTV sampling3 - Control sampling4      -0.2870 0.0724 21  -3.963  0.0135 *
# DTV sampling3 - DTV sampling4          -0.2022 0.0724 21  -2.793  0.1498
# Control sampling4 - DTV sampling4       0.0847 0.0724 21   1.170  0.9318

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
#preheat to postheat
##both species
leveneTest(pam~Treatment,data=phys.less.geno.time) #ns
shapiro.test(phys.less.geno.time$pam) #ns

# Genotype is a random effect
lmer.pam <- lmer(pam ~ Treatment * Species * sampling + (1 | Genotype), data = phys.less.geno.time)

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


###All Pam timepoints
##favs only
phys.less3.geno.time.fav <- subset(phys.less3.geno.time,Species=="O. faveolata")

leveneTest(pam~Treatment,data=phys.less3.geno.time.fra) #ns
shapiro.test(phys.less3.geno.time.fra$pam) #ns

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
# (Intercept)        1405.3428  1  < 2.2e-16 ***
#   Treatment             3.0839  1    0.07907 .  
#   sampling             44.2042  4  5.819e-09 ***
#   Treatment:sampling    3.7123  4    0.44634  
  

#Pairwise comparisons
emmeans(lmer.pam, pairwise ~ Treatment * sampling, adjust = "tukey")
# $contrasts
# contrast                               estimate      SE df t.ratio p.value
# Control sampling1 - DTV sampling1     -1.44e-02 0.00819 27  -1.756  0.7550
# Control sampling1 - Control sampling2 -1.43e-02 0.00819 27  -1.749  0.7589
# Control sampling1 - DTV sampling2     -4.58e-02 0.00819 27  -5.592  0.0002 ***
# Control sampling1 - Control sampling3 -8.28e-03 0.00819 27  -1.010  0.9890
# Control sampling1 - DTV sampling3     -2.58e-02 0.00819 27  -3.151  0.0933
# Control sampling1 - Control sampling4  5.04e-03 0.00819 27   0.615  0.9997
# Control sampling1 - DTV sampling4     -5.68e-03 0.00819 27  -0.693  0.9993
# Control sampling1 - Control sampling5  3.53e-02 0.00819 27   4.311  0.0062 **
# Control sampling1 - DTV sampling5      1.85e-02 0.00819 27   2.254  0.4472
# DTV sampling1 - Control sampling2      5.56e-05 0.00819 27   0.007  1.0000
# DTV sampling1 - DTV sampling2         -3.14e-02 0.00819 27  -3.836  0.0200 *
# DTV sampling1 - Control sampling3      6.11e-03 0.00819 27   0.746  0.9988
# DTV sampling1 - DTV sampling3         -1.14e-02 0.00819 27  -1.395  0.9184
# DTV sampling1 - Control sampling4      1.94e-02 0.00819 27   2.371  0.3797
# DTV sampling1 - DTV sampling4          8.71e-03 0.00819 27   1.063  0.9845
# DTV sampling1 - Control sampling5      4.97e-02 0.00819 27   6.067  0.0001 ***
# DTV sampling1 - DTV sampling5          3.29e-02 0.00819 27   4.011  0.0131 *
# Control sampling2 - DTV sampling2     -3.15e-02 0.00819 27  -3.843  0.0197 *
# Control sampling2 - Control sampling3  6.06e-03 0.00819 27   0.739  0.9989
# Control sampling2 - DTV sampling3     -1.15e-02 0.00819 27  -1.402  0.9162
# Control sampling2 - Control sampling4  1.94e-02 0.00819 27   2.365  0.3835
# Control sampling2 - DTV sampling4      8.65e-03 0.00819 27   1.056  0.9851
# Control sampling2 - Control sampling5  4.97e-02 0.00819 27   6.060  0.0001 ***
# Control sampling2 - DTV sampling5      3.28e-02 0.00819 27   4.004  0.0133 *
# DTV sampling2 - Control sampling3      3.75e-02 0.00819 27   4.582  0.0031 **
# DTV sampling2 - DTV sampling3          2.00e-02 0.00819 27   2.441  0.3422
# DTV sampling2 - Control sampling4      5.09e-02 0.00819 27   6.207  <.0001 ***
# DTV sampling2 - DTV sampling4          4.01e-02 0.00819 27   4.899  0.0014 **
# DTV sampling2 - Control sampling5      8.11e-02 0.00819 27   9.903  <.0001 ***
# DTV sampling2 - DTV sampling5          6.43e-02 0.00819 27   7.847  <.0001 ***
# Control sampling3 - DTV sampling3     -1.75e-02 0.00819 27  -2.141  0.5169
# Control sampling3 - Control sampling4  1.33e-02 0.00819 27   1.626  0.8245
# Control sampling3 - DTV sampling4      2.60e-03 0.00819 27   0.317  1.0000
# Control sampling3 - Control sampling5  4.36e-02 0.00819 27   5.321  0.0005 ***
# Control sampling3 - DTV sampling5      2.67e-02 0.00819 27   3.265  0.0734
# DTV sampling3 - Control sampling4      3.09e-02 0.00819 27   3.766  0.0236 *
# DTV sampling3 - DTV sampling4          2.01e-02 0.00819 27   2.458  0.3334
# DTV sampling3 - Control sampling5      6.11e-02 0.00819 27   7.462  <.0001 ***
# DTV sampling3 - DTV sampling5          4.43e-02 0.00819 27   5.406  0.0004 ***
# Control sampling4 - DTV sampling4     -1.07e-02 0.00819 27  -1.309  0.9430
# Control sampling4 - Control sampling5  3.03e-02 0.00819 27   3.695  0.0279 *
# Control sampling4 - DTV sampling5      1.34e-02 0.00819 27   1.639  0.8177
# DTV sampling4 - Control sampling5      4.10e-02 0.00819 27   5.004  0.0011 **
# DTV sampling4 - DTV sampling5          2.42e-02 0.00819 27   2.948  0.1407
# Control sampling5 - DTV sampling5     -1.68e-02 0.00819 27  -2.056  0.5705

##franks
phys.less3.geno.time.fra <- subset(phys.less3.geno.time,Species=="O. franksi")

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
# Analysis of Deviance Table (Type III tests)
# 
# Response: pam
# Chisq Df Pr(>Chisq)    
# (Intercept)        2013.3019  1  < 2.2e-16 ***
#   Treatment            12.1091  1  0.0005018 ***
#   sampling             21.7071  4  0.0002292 ***
#   Treatment:sampling    2.8515  4  0.5829772 

#Pairwise comparisons
emmeans(lmer.pam, pairwise ~ Treatment * sampling, adjust = "tukey")
# $contrasts
# contrast                              estimate     SE df t.ratio p.value
# Control sampling1 - DTV sampling1     -0.03994 0.0115 27  -3.480  0.0457 *
# Control sampling1 - Control sampling2 -0.03594 0.0115 27  -3.131  0.0972
# Control sampling1 - DTV sampling2     -0.07704 0.0115 27  -6.712  <.0001 ***
# Control sampling1 - Control sampling3 -0.02893 0.0115 27  -2.520  0.3021
# Control sampling1 - DTV sampling3     -0.06142 0.0115 27  -5.350  0.0004 ***
# Control sampling1 - Control sampling4 -0.01486 0.0115 27  -1.295  0.9464
# Control sampling1 - DTV sampling4     -0.04303 0.0115 27  -3.748  0.0246 *
# Control sampling1 - Control sampling5  0.00900 0.0115 27   0.784  0.9983
# Control sampling1 - DTV sampling5     -0.00833 0.0115 27  -0.726  0.9990
# DTV sampling1 - Control sampling2      0.00400 0.0115 27   0.348  1.0000
# DTV sampling1 - DTV sampling2         -0.03710 0.0115 27  -3.232  0.0787
# DTV sampling1 - Control sampling3      0.01101 0.0115 27   0.959  0.9923
# DTV sampling1 - DTV sampling3         -0.02147 0.0115 27  -1.871  0.6872
# DTV sampling1 - Control sampling4      0.02508 0.0115 27   2.185  0.4894
# DTV sampling1 - DTV sampling4         -0.00308 0.0115 27  -0.269  1.0000
# DTV sampling1 - Control sampling5      0.04894 0.0115 27   4.264  0.0070 **
# DTV sampling1 - DTV sampling5          0.03161 0.0115 27   2.754  0.2027
# Control sampling2 - DTV sampling2     -0.04110 0.0115 27  -3.580  0.0364 *
# Control sampling2 - Control sampling3  0.00701 0.0115 27   0.611  0.9998
# Control sampling2 - DTV sampling3     -0.02547 0.0115 27  -2.219  0.4686
# Control sampling2 - Control sampling4  0.02108 0.0115 27   1.837  0.7078
# Control sampling2 - DTV sampling4     -0.00708 0.0115 27  -0.617  0.9997
# Control sampling2 - Control sampling5  0.04494 0.0115 27   3.915  0.0165 *
# Control sampling2 - DTV sampling5      0.02761 0.0115 27   2.405  0.3611
# DTV sampling2 - Control sampling3      0.04811 0.0115 27   4.191  0.0084 **
# DTV sampling2 - DTV sampling3          0.01562 0.0115 27   1.361  0.9287
# DTV sampling2 - Control sampling4      0.06218 0.0115 27   5.417  0.0004 ***
# DTV sampling2 - DTV sampling4          0.03401 0.0115 27   2.963  0.1365
# DTV sampling2 - Control sampling5      0.08604 0.0115 27   7.496  <.0001 ***
# DTV sampling2 - DTV sampling5          0.06871 0.0115 27   5.986  0.0001 ***
# Control sampling3 - DTV sampling3     -0.03249 0.0115 27  -2.830  0.1762
# Control sampling3 - Control sampling4  0.01407 0.0115 27   1.226  0.9613
# Control sampling3 - DTV sampling4     -0.01410 0.0115 27  -1.228  0.9609
# Control sampling3 - Control sampling5  0.03793 0.0115 27   3.304  0.0674
# Control sampling3 - DTV sampling5      0.02060 0.0115 27   1.794  0.7329
# DTV sampling3 - Control sampling4      0.04656 0.0115 27   4.056  0.0117 *
# DTV sampling3 - DTV sampling4          0.01839 0.0115 27   1.602  0.8359
# DTV sampling3 - Control sampling5      0.07042 0.0115 27   6.134  0.0001 ***
# DTV sampling3 - DTV sampling5          0.05308 0.0115 27   4.624  0.0028 **
# Control sampling4 - DTV sampling4     -0.02817 0.0115 27  -2.454  0.3355
# Control sampling4 - Control sampling5  0.02386 0.0115 27   2.079  0.5561
# Control sampling4 - DTV sampling5      0.00653 0.0115 27   0.569  0.9999
# DTV sampling4 - Control sampling5      0.05203 0.0115 27   4.532  0.0036 **
# DTV sampling4 - DTV sampling5          0.03469 0.0115 27   3.022  0.1214
# Control sampling5 - DTV sampling5     -0.01733 0.0115 27  -1.510  0.8766


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

### Fig S4 % Change in weight all timepoints ####
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
###Preheat to Postheat Time points buoyant weights
##favs only
phys.less.weight.time.fav <- subset(phys.less.weight.time,Species=="O. faveolata")

#Check assumptions for LME
# Homogeneity of variances
leveneTest(weightChan~Treatment,data=phys.less.weight.time.fav) #ns
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



###All Time points buoyant weights
##favs only
phys.less2.weight.time.fav <- subset(phys.less2.weight.time,Species=="O. faveolata")

#Check assumptions for LME
# Homogeneity of variances
leveneTest(weightChan~Treatment,data=phys.less2.weight.time.fav) #ns
# Normality of residuals
shapiro.test(phys.less2.weight.time.fav$weightChan) #ns

# control for Genotype is a random effect
lmer.weightChan <- lmer(weightChan ~ Treatment * sampling + (1 | Genotype), data = phys.less2.weight.time.fav)

                      
check_model(lmer.weightChan) # check assumptions more thoroughly, looks good!

# Summary of the model
summary(lmer.weightChan)

# Type III ANOVA table (if needed)
Anova(lmer.weightChan, type = 3)
# Analysis of Deviance Table (Type III tests)
# 
# Response: weightChan
# Chisq Df Pr(>Chisq)    
# (Intercept)        30.9750  1  2.614e-08 ***
#   Treatment           0.0019  1     0.9656    
#   sampling           28.5146  3  2.832e-06 ***
#   Treatment:sampling  5.0545  3     0.1678  


#Pairwise comparisons
emmeans(lmer.weightChan, pairwise ~ Treatment * sampling, adjust = "tukey")
# $contrasts
# contrast                              estimate    SE df t.ratio p.value
# Control sampling1 - DTV sampling1       0.0206 0.476 21   0.043  1.0000
# Control sampling1 - Control sampling2   0.3945 0.476 21   0.828  0.9893
# Control sampling1 - DTV sampling2       0.5832 0.476 21   1.224  0.9155
# Control sampling1 - Control sampling3   2.0223 0.476 21   4.245  0.0072
# Control sampling1 - DTV sampling3       1.0413 0.476 21   2.186  0.3980
# Control sampling1 - Control sampling4   1.9237 0.476 21   4.038  0.0114
# Control sampling1 - DTV sampling4       2.3469 0.476 21   4.926  0.0015
# DTV sampling1 - Control sampling2       0.3740 0.476 21   0.785  0.9922
# DTV sampling1 - DTV sampling2           0.5627 0.476 21   1.181  0.9287
# DTV sampling1 - Control sampling3       2.0018 0.476 21   4.202  0.0079
# DTV sampling1 - DTV sampling3           1.0207 0.476 21   2.143  0.4218
# DTV sampling1 - Control sampling4       1.9031 0.476 21   3.995  0.0126
# DTV sampling1 - DTV sampling4           2.3263 0.476 21   4.883  0.0017
# Control sampling2 - DTV sampling2       0.1887 0.476 21   0.396  0.9999
# Control sampling2 - Control sampling3   1.6278 0.476 21   3.417  0.0439
# Control sampling2 - DTV sampling3       0.6467 0.476 21   1.357  0.8659
# Control sampling2 - Control sampling4   1.5291 0.476 21   3.210  0.0672
# Control sampling2 - DTV sampling4       1.9523 0.476 21   4.098  0.0100
# DTV sampling2 - Control sampling3       1.4391 0.476 21   3.021  0.0976
# DTV sampling2 - DTV sampling3           0.4581 0.476 21   0.961  0.9752
# DTV sampling2 - Control sampling4       1.3405 0.476 21   2.814  0.1442
# DTV sampling2 - DTV sampling4           1.7637 0.476 21   3.702  0.0239
# Control sampling3 - DTV sampling3      -0.9811 0.476 21  -2.059  0.4695
# Control sampling3 - Control sampling4  -0.0986 0.476 21  -0.207  1.0000
# Control sampling3 - DTV sampling4       0.3246 0.476 21   0.681  0.9967
# DTV sampling3 - Control sampling4       0.8824 0.476 21   1.852  0.5947
# DTV sampling3 - DTV sampling4           1.3056 0.476 21   2.741  0.1646
# Control sampling4 - DTV sampling4       0.4232 0.476 21   0.888  0.9840


##franks
phys.less2.weight.time.fra <- subset(phys.less2.weight.time,Species=="O. franksi")

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
# (Intercept)        37.2256  1  1.052e-09 ***
#   Treatment           0.0004  1     0.9841    
#   sampling           51.4117  3  3.997e-11 ***
#   Treatment:sampling  1.6556  3     0.6469    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#Pairwise comparisons
emmeans(lmer.weightChan, pairwise ~ Treatment * sampling, adjust = "tukey")
# $contrasts
# contrast                              estimate    SE df t.ratio p.value
# Control sampling1 - DTV sampling1     -0.00585 0.293 21  -0.020  1.0000
# Control sampling1 - Control sampling2  0.45472 0.293 21   1.552  0.7719
# Control sampling1 - DTV sampling2      0.33033 0.293 21   1.127  0.9432
# Control sampling1 - Control sampling3  1.20478 0.293 21   4.111  0.0097
# Control sampling1 - DTV sampling3      1.58696 0.293 21   5.415  0.0005
# Control sampling1 - Control sampling4  1.95183 0.293 21   6.661  <.0001
# Control sampling1 - DTV sampling4      1.98729 0.293 21   6.782  <.0001
# DTV sampling1 - Control sampling2      0.46057 0.293 21   1.572  0.7610
# DTV sampling1 - DTV sampling2          0.33618 0.293 21   1.147  0.9381
# DTV sampling1 - Control sampling3      1.21063 0.293 21   4.131  0.0093
# DTV sampling1 - DTV sampling3          1.59280 0.293 21   5.435  0.0005
# DTV sampling1 - Control sampling4      1.95768 0.293 21   6.680  <.0001
# DTV sampling1 - DTV sampling4          1.99314 0.293 21   6.802  <.0001
# Control sampling2 - DTV sampling2     -0.12439 0.293 21  -0.424  0.9998
# Control sampling2 - Control sampling3  0.75006 0.293 21   2.560  0.2251
# Control sampling2 - DTV sampling3      1.13224 0.293 21   3.864  0.0168
# Control sampling2 - Control sampling4  1.49711 0.293 21   5.109  0.0010
# Control sampling2 - DTV sampling4      1.53258 0.293 21   5.230  0.0008
# DTV sampling2 - Control sampling3      0.87445 0.293 21   2.984  0.1048
# DTV sampling2 - DTV sampling3          1.25662 0.293 21   4.288  0.0065
# DTV sampling2 - Control sampling4      1.62150 0.293 21   5.533  0.0004
# DTV sampling2 - DTV sampling4          1.65696 0.293 21   5.654  0.0003
# Control sampling3 - DTV sampling3      0.38217 0.293 21   1.304  0.8873
# Control sampling3 - Control sampling4  0.74705 0.293 21   2.549  0.2290
# Control sampling3 - DTV sampling4      0.78251 0.293 21   2.670  0.1864
# DTV sampling3 - Control sampling4      0.36487 0.293 21   1.245  0.9086
# DTV sampling3 - DTV sampling4          0.40034 0.293 21   1.366  0.8622
# Control sampling4 - DTV sampling4      0.03547 0.293 21   0.121  1.0000

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

