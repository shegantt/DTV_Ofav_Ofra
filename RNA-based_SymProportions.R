setwd("~/Desktop/~directory path~")

all_data <- read.csv("zoox_counts_new_fixed.csv")
head(all_data)

library(ggplot2)
library(Rmisc)

#model whether fracD is different between treatments
model1 <- lm(fracD ~ treat*time+spp, data = all_data)
summary(model1)

#use dplyr to make one combine fracA, fracB, fracC, and fracD into one column
library(tidyr)
all_data <- all_data %>% 
  gather(key = "symtype", value = "fracD", fracA, fracB, fracC, fracD)

all_data

#change the species names
all_data$spp <- gsub("Ofra", "O. franksi", all_data$spp)
all_data$spp <- gsub("Ofav", "O. faveolata", all_data$spp)

#change the time names
all_data$time <- gsub("preheat", "Pre-heat", all_data$time)
all_data$time <- gsub("postheat", "Post-heat", all_data$time)

#change the treatment names
all_data$treat <- gsub("control", "Control", all_data$treat)
all_data$treat <- gsub("therm_var", "DTV", all_data$treat)

#change symtype frac A to A
all_data$symtype <- gsub("fracA", "Symbiodinium", all_data$symtype)
all_data$symtype <- gsub("fracB", "Breviolum", all_data$symtype)
all_data$symtype <- gsub("fracC", "Cladocopium", all_data$symtype)
all_data$symtype <- gsub("fracD", "Durusdinium", all_data$symtype)

#reorder a factor
all_data$time <- factor(all_data$time, levels = c("Pre-heat", "Post-heat"))
f=c("Symbiodinium", "Breviolum", "Cladocopium", "Durusdinium")
all_data <- within(all_data, symtype <- factor(symtype, levels=f))

head(all_data)
colnames(all_data)[colnames(all_data) == 'symtype'] <- 'Genus'

library(ggh4x)
#plot the data and export as a pdf size 4x4
#pdf("Nicola_sym_av.pdf", width=4, height=4)
ggplot(all_data, aes(x=treat, y=fracD, fill=Genus)) +
  geom_bar(stat="identity", position="stack") +
  theme_bw() +
  xlab("Thermal variability treatment") +
  ylab("Relative Algal Symbiont Proportion") +
  scale_fill_manual(values=c("gray1", "grey25", "gray50", "gray85")) +
  facet_grid2(time~spp, scales = 'free', strip = strip_themed(
               background_y = list(element_rect(fill = "black"),
                                   element_rect(fill = "black")),
               background_x = list(element_rect(fill = "#8B7355"),
                                   element_rect(fill = "#CDAB7D")),
               text_y = elem_list_text(face = "bold", color = "white", size = 17),
               text_x = elem_list_text(face = "bold", color = "white", size = 17)))+
  #theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+
  theme(axis.text = element_text(color = "black", size = 14, face = "bold"),
        legend.title=element_text(size=14, face = "bold", color = "black"), 
        legend.text=element_text(size=14, face = "bold", color="black"),
        axis.title.x = element_text(size = 16, color = "black", face = "bold"),
        axis.title.y = element_text(size = 16, color = "black", face = "bold"))+
  theme(strip.text.x = element_text(face = "italic"))
#dev.off()

#looking at the data by species within genets
frank=subset(all_data, spp=="O. franksi")
fav=subset(all_data, spp=="O. faveolata")

head(all_data)

#O. franksi symbiont proportions
#pdf("Nicola_sym_frank_indiv.pdf", width=6, height=4)
ggplot(frank, aes(x=actual_gen, y=fracD, fill=Genus)) + 
  geom_col(position = "fill")+
  theme_bw() +
  xlab("Coral Genotype") +
  ylab("Relative Algal Symbiont Proportion") +
  scale_fill_manual(values=c("gray1", "grey25", "gray50", "gray85")) +
  facet_grid2(time~treat, scales = 'free', strip = 
              strip_themed(background_y = list(element_rect(fill = "black"), 
                                               element_rect(fill = "black")),
               background_x = list(element_rect(fill = "#2C2673"), element_rect(fill = "#E24A13")), 
                                   text_y = elem_list_text(face = "bold", color = "white", size = 17),
                                   text_x = elem_list_text(face = "bold", color = "white", size = 17)))+
  theme(axis.text = element_text(color = "black", size = 14, face = "bold"),
        legend.title=element_text(size=14, face = "bold", color = "black"), 
        legend.text=element_text(size=14, face = "bold", color="black"),
        axis.title.x = element_text(size = 16, color = "black", face = "bold"),
        axis.title.y = element_text(size = 16, color = "black", face = "bold"))
#dev.off()

#O. faveolata symbiont proportions
#pdf("Nicola_sym_fav_indiv.pdf", width=6, height=4)
ggplot(fav, aes(x=actual_gen, y=fracD, fill=Genus)) + 
  geom_col(position = "fill")+
  theme_bw() +
  xlab("Coral Genotype") +
  ylab("Relative Algal Symbiont Proportion") +
  scale_fill_manual(values=c("gray1", "grey25", "gray50", "gray85")) +
  facet_grid2(time~treat, scales = 'free', 
              strip = strip_themed(background_y = list(element_rect(fill = "black"), 
                                                       element_rect(fill = "black")),
              background_x = list(element_rect(fill = "#2C2673"), element_rect(fill = "#E24A13")),
                         text_y = elem_list_text(face = "bold", color = "white", size = 17),
                         text_x = elem_list_text(face = "bold", color = "white", size = 17)))+
  theme(axis.text = element_text(color = "black", size = 14, face = "bold"),
        legend.title=element_text(size=14, face = "bold", color = "black"), 
        legend.text=element_text(size=14, face = "bold", color="black"),
        axis.title.x = element_text(size = 16, color = "black", face = "bold"),
        axis.title.y = element_text(size = 16, color = "black", face = "bold"))
#dev.off()

#old_pallete c("sienna4", "lightseagreen", "#17365e", "deeppink4")
