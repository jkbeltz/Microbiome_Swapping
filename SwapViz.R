library(dplyr)
library(tidyverse)
library(plyr)
library(cowplot)
library(nlme)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(scales)
library(ggbreak)
library(patchwork)
library(ggpubr)


setwd("/Users/jackbeltz/Desktop/swap/swap")

####IMPORT####

cagepheno=read_csv("allpheno22.cagesum.csv")
treatpheno=read_csv("allpheno22.treatsum.csv")

### make them long###
cagepheno_long <- melt(cagepheno ,  id.vars = c('cage.number','pheno.treatment', 'label', 'cage.treatment'), variable.name = 'Phenotypes')
cagepheno_long=cagepheno_long[!grepl("sd", cagepheno_long$Phenotypes),]
cagepheno_long=cagepheno_long[!grepl("se", cagepheno_long$Phenotypes),]
treatpheno_long <- melt(treatpheno ,  id.vars = c('pheno.treatment', 'label', 'cage.treatment'), variable.name = 'Phenotypes')
treatpheno_long=treatpheno_long[!grepl("sd", treatpheno_long$Phenotypes),]
treatpheno_long=treatpheno_long[!grepl("se", treatpheno_long$Phenotypes),]

view(cagepheno_long)

TEST <- ggplot(cagepheno_long, aes(x=as.factor(pheno.treatment), y=value, color=cage.treatment))+
  geom_boxplot()+ 
  geom_jitter()+
  scale_x_discrete(name = "",limits = c("F", "N", "I", "S"), labels = c("Founder (Lab)", "None", "Initiation (Field)", "Summer (Field)"))+
  #scale_fill_manual(name = "",breaks=c("F","L"), labels =c("Field","Lab"), values = c("darkgreen","grey"))+
  scale_color_manual(name = "",breaks=c("F","L"), labels =c("Field","Lab"), values = c("darkgreen","darkgrey"))+
  4ii
  theme_bw(base_size = 12 )+
  theme(strip.background = element_blank(),
        strip.placement = "outside")
TEST
 

TEST1<-lme(LW ~ cage.treatment*pheno.treatment, random=~1|cage.number/cage.treatment, data=cagepheno) ##
anova(TEST1) ##TP + TREATMENT SIG / INTERCEPT = 0.08


##LDF INTERACTION
##DESS ALMOST INTERACTION