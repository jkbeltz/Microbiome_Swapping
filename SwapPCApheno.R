library(reshape2)
library(plyr)
library(tidyverse)
library(RColorBrewer)
library(cowplot)
library(patchwork)
library(latex2exp)
library(ggrepel)
library(scales)
library(grid)
library(ggtext)
library(ggpubr)
setwd("/Users/jackbeltz/Documents/PENN/Dissertation/CH 4 (Swapping)/Microbiome_swapping")

swapPCA <- read.csv("allpheno22.cagesum.csv") ##read data, generated in swapANAL
swapPCA$group <- paste(swapPCA1$cage.treatment, swapPCA1$pheno.treatment) ##create grouping variable
swapPCA <- swapPCA[, -c(6,7,9,10,12,13,15,16,18,20,21,23,24)]## remove SE and SD
View(swapPCA)



#### LISTS AND COLOR ASSIGNMENTS####

swapColorListPheno <- c("N"="grey",
                      "F"="blue",
                      "I"= "yellow", 
                      "S"= "darkgreen")

swapColorListPopulation <- c("L"="red", "F" = "black")

swapColorListMDS <- c("L"="red", "F" = "black")

swapColorListTrait<-c(
  "LD_M"=brewer.pal(n=12, name="Set3")[2],
  "LD_F"=brewer.pal(n=12, name="Set3")[3],
  "SR"=brewer.pal(n=12, name="Set3")[4],
  "DESS"=brewer.pal(n=12, name="Set3")[5],
  "PIG"=brewer.pal(n=12,name="Set3")[6],
  "DW"=brewer.pal(n=12, name="Set3")[7],
  "LW"=brewer.pal(n=12, name="Set3")[8])

swapPHENOLIST=c("N","F","I","S")
swapPOPULATIONLIST=c("F","L")
swapTraitsList<-data.frame(ABBRV=c("LD_M","LD_F","SR","DESS","PIG","DW","LW"),LONG=c("Development Time (Males)","Development Time (Females)","Starvation Resistance","Dessication","Pigmentation","Body Size","Lipid Weight"))
swapSigTrait3List<-data.frame(ABBRV=c("LD_M","LD_F","SR","DESS","DW"),LONG=c("Development Time (Males)","Development Time (Females)","Starvation Resistance","Dessication","Body Size"))

swapSelectionLabels=c("None","Founder","Initiation", "Summer-Evolved")
swapPopulationLabels=c("Field","Lab")
swapPhenoShape=c("N"=0,"F"=1,"I"=2,"S"=8)

#SETTING A CONSISTENT IN ALL FIGURES THEME

{
  {
    BASESIZE=8
    TITLESIZE=7
    FONTFAMILY="sans"
    TEXTCOLOR="black"
    AXISTEXTSIZE=7
    AXISTEXTCOLOR="grey40"
    AXISTITLESIZE=7
    LEGENDTEXTSIZE=6
    LABELSIZE=11
  }
  
  theme_evo <- theme_set(theme_bw(base_size=BASESIZE,
                                  base_family=FONTFAMILY))
  theme_evo <- theme_update(panel.background = element_blank(),
                            title = element_text(size=TITLESIZE, face= "bold"),
                            axis.text = element_text(size= AXISTEXTSIZE, color = AXISTEXTCOLOR, family = FONTFAMILY, face= "bold"),
                            axis.title = element_text(size= TITLESIZE, color = TEXTCOLOR, family = FONTFAMILY, face= "bold"),
                            legend.background = element_blank(),
                            legend.box.background = element_blank(),
                            legend.title=element_text(face= "bold",size = LEGENDTEXTSIZE),
                            legend.text=element_text(size = LEGENDTEXTSIZE),
                            legend.key =element_blank()
  )
  
}

# EMPTY PLOT
pBL<-ggplot()+theme_nothing()

{
  
  #### PREP####
  {
    # FUNTION TO FIND HULLS
    find_hull<-function(df) df[chull(df$PC1,df$PC2),]
 
    # SORT THE DATA
    swapPCA$pheno <-factor(swapPCA$pheno.treatment,swapPHENOLIST)
    swapPCA$Population=factor(swapPCA$cage.treatment,swapPOPULATIONLIST)
    swapPCA<-swapPCA %>% arrange(cage.number,pheno,Population)
    
    #View(ABpca2)
    #### CALCULATING THE PRINCIPAL COMPONENTS####
    swap_Lab <- swapPCA %>% subset(Population=="L")
    swap_Field <- swapPCA %>% subset(Population=="F")
    
    swap_N<- swapPCA %>% subset(pheno=="N")
    swap_F<- swapPCA %>% subset(pheno=="F")
    swap_I<- swapPCA %>% subset(pheno=="I")
    swap_S<- swapPCA %>% subset(pheno=="S")
  
    ### all variables
    swap.pca<-prcomp(swapPCA[5:11],scale=TRUE)

    ### sig variables
    swap_s.pca<-prcomp(swapPCA[c(5:8,10)],scale=TRUE)
    
    # GETTING PCA COORDINATES TO A DATA FRAME
    ##All Vars
    swap.pca.df<-data.frame(Cage=swapPCA$cage.number,Population=swapPCA$Population,Pheno=swapPCA$pheno,Group=swapPCA$group,swap.pca$x)

    ## manova for pheno / pop effect om all variables
    swap.PCA.aov<-manova(cbind(PC1,PC2,PC3,PC4,PC5,PC6,PC7)~Pheno*Population,swap.pca.df)
    summary(swap.PCA.aov) ### 
    
    ##Sig Vars 
    swap.s.pca.df<-data.frame(Cage=swapPCA$cage.number,Population=swapPCA$Population,Pheno=swapPCA$pheno,Group=swapPCA$group,swap_s.pca$x)
    
    ## manova for pheno / pop effect om sig variables
    swap.s.PCA.aov<-manova(cbind(PC1,PC2,PC3,PC4,PC5)~Pheno*Population,swap.s.pca.df)
    summary(swap.s.PCA.aov) ### 
    
    #### GETTING VECTORS TO A DATA FRAME####
    ##All Vars 
    
    swap.pca.vectors<-as.data.frame(swap.pca$rotation)
    swap.pca.vectors$ABBRV<-rownames(swap.pca.vectors)
    #View(swap.pca.vectors)
    ##SIG VARS ##
    
    swap_s.pca.vectors<-as.data.frame(swap_s.pca$rotation)
    swap_s.pca.vectors$ABBRV<-rownames(swap_s.pca.vectors)
    #View(swap_s.pca.vectors)
  
    
    # CALCULATING VAR EXPLAINED
    ##All Vars 
    swap.pca.var.explained<-swap.pca$sdev^2/sum(swap.pca$sdev^2)

    ##Sig Vars 
    swap_s.pca.var.explained<-swap_s.pca$sdev^2/sum(swap_s.pca$sdev^2)
    View(swap.pca.df)
    
    ####calculating hulls####
    # SUBSETTING AND CALCULATING HULLS FOR EACH Population FOR PLOTTING
    
    swap.FN.pca.df <-subset(swap.pca.df,Group=="F N")
    swap.FN.pca.hull<- swap.FN.pca.df  %>% ddply("Pheno",find_hull)
    
    swap.FF.pca.df <-subset(swap.pca.df,Group=="F F")
    swap.FF.pca.hull<- swap.FF.pca.df  %>% ddply("Population",find_hull)
    View(swap.FF.pca.df)
    swap.FI.pca.df <-subset(swap.pca.df,Group=="F I")
    swap.FI.pca.hull<- swap.FI.pca.df  %>% ddply("Pheno",find_hull)
    
    swap.FS.pca.df <-subset(swap.pca.df,Group=="F S")
    swap.FS.pca.hull<- swap.FS.pca.df  %>% ddply("Pheno",find_hull)
    
    swap.LN.pca.df <-subset(swap.pca.df,Group=="L N")
    swap.LN.pca.hull<- swap.LN.pca.df  %>% ddply("Pheno",find_hull)
    
    swap.LF.pca.df <-subset(swap.pca.df,Group=="L F")
    swap.LF.pca.hull<- swap.LF.pca.df  %>% ddply("Pheno",find_hull)
    
    swap.LI.pca.df <-subset(swap.pca.df,Group=="L I")
    swap.LI.pca.hull<- swap.LI.pca.df  %>% ddply("Pheno",find_hull)
    
    swap.LS.pca.df <-subset(swap.pca.df,Group=="L S")
    swap.LS.pca.hull<- swap.LS.pca.df  %>% ddply("Pheno",find_hull)
    
    #### HUlls using just sig variables
    swap.s.FN.pca.df <-subset(swap.s.pca.df,Group=="F N")
    swap.s.FN.pca.hull<- swap.s.FN.pca.df  %>% ddply("Pheno",find_hull)
    
    swap.s.FF.pca.df <-subset(swap.s.pca.df,Group=="F F")
    swap.s.FF.pca.hull<- swap.s.FF.pca.df  %>% ddply("Pheno",find_hull)
    
    swap.s.FI.pca.df <-subset(swap.s.pca.df,Group=="F I")
    swap.s.FI.pca.hull<- swap.s.FI.pca.df  %>% ddply("Pheno",find_hull)
    
    swap.s.FS.pca.df <-subset(swap.s.pca.df,Group=="F S")
    swap.s.FS.pca.hull<- swap.s.FS.pca.df  %>% ddply("Pheno",find_hull)
    
    swap.s.LN.pca.df <-subset(swap.s.pca.df,Group=="L N")
    swap.s.LN.pca.hull<- swap.s.LN.pca.df  %>% ddply("Pheno",find_hull)
    
    swap.s.LF.pca.df <-subset(swap.s.pca.df,Group=="L F")
    swap.s.LF.pca.hull<- swap.s.LF.pca.df  %>% ddply("Pheno",find_hull)
    
    swap.s.LI.pca.df <-subset(swap.s.pca.df,Group=="L I")
    swap.s.LI.pca.hull<- swap.s.LI.pca.df  %>% ddply("Pheno",find_hull)
    
    swap.s.LS.pca.df <-subset(swap.s.pca.df,Group=="L S")
    swap.s.LS.pca.hull<- swap.s.LS.pca.df  %>% ddply("Pheno",find_hull)
    
  }
  
  #### PLOT Config####
  {
    # SIZES
    {
      F3_POINTSIZE=2
      F3_LINESIZE=1
      F3_FILLALPHA=0.3
      
      F3_XLLIM=-4.3
      F3_XULIM=3.15
      F3A_YLLIM=-3.0
      F3A_YULIM=4.45
      F3B_YLLIM=-3.45
      F3B_YULIM=4.0
      
      F3_LGDJUST=c(0,1)
      F3_LGDPOS="none" ###chnage if want legends on indivisual plots to c(0,0)
      F3_LGDSPCX=0.02
      F3_LGDSPCY=0.1
      F3_LGDKEYSIZE=0.6
      
      F3B_ARROWSIZE=.5
      F3_ARROWHEAD=0.2
      F3B_TEXTSIZE=2.5
      F3_VECTORSCALEX=1.2
      F3_VECTORSCALEY=1.2
      F3_VECTOR_ALPHA=0.4
      F3_CIRCLE_COLOR="grey90"
      
      F3_VECTORS_POS_L=0.0
      F3_VECTORS_POS_B=0.62
      F3_VECTORS_POS_T=1.02
      F3_VECTORS_POS_R=0.38
      
      F4_VECTORS_POS_L=1.0
      F4_VECTORS_POS_B=0.62
      F4_VECTORS_POS_T=1.02
      F4_VECTORS_POS_R=1.38
      
      
      F3_PDFHW=7.1
      F3_PDFHH=3.5
      F3_PDFVW=3.5
      F3_PDFVH=7.1
      
      # NEEDED FOR THE PCA VECTORS PLOT
      angle <- seq(-pi, pi, length = 50) 
      circle <- data.frame(x = sin(angle), y = cos(angle))
    }
    
    
    #### PLOT with all variables####
    #####On Bloomoington ####
    
    swapPHENOS<-ggplot()+
      scale_color_manual(values = swapColorListPopulation) +
      scale_shape_manual(values = swapPhenoShape) +
      scale_fill_manual(values = swapColorListPheno)+
      
      geom_polygon(data=swap.FN.pca.hull, alpha=0,aes(x=PC1,y=PC2,fill=as.factor(Pheno),color=as.factor(Population)),size=F3_LINESIZE) +
      geom_polygon(data=swap.FN.pca.hull, alpha=F3_FILLALPHA,aes(x=PC1,y=PC2,fill=as.factor(Pheno))) +
      
      geom_polygon(data=swap.FF.pca.hull, alpha=0,aes(x=PC1,y=PC2,fill=as.factor(Pheno),color=as.factor(Population)),size=F3_LINESIZE) +
      geom_polygon(data=swap.FF.pca.hull, alpha=F3_FILLALPHA,aes(x=PC1,y=PC2,fill=as.factor(Pheno))) +
      
      geom_polygon(data=swap.FI.pca.hull, alpha=0,aes(x=PC1,y=PC2,fill=as.factor(Pheno),color=as.factor(Population)),size=F3_LINESIZE) +
      geom_polygon(data=swap.FI.pca.hull, alpha=F3_FILLALPHA,aes(x=PC1,y=PC2,fill=as.factor(Pheno))) +
      
      geom_polygon(data=swap.FS.pca.hull, alpha=0,aes(x=PC1,y=PC2,fill=as.factor(Pheno),color=as.factor(Population)),size=F3_LINESIZE) +
      geom_polygon(data=swap.FS.pca.hull, alpha=F3_FILLALPHA,aes(x=PC1,y=PC2,fill=as.factor(Pheno))) +
      
      geom_polygon(data=swap.LN.pca.hull, alpha=0,aes(x=PC1,y=PC2,fill=as.factor(Pheno),color=as.factor(Population)),size=F3_LINESIZE) +
      geom_polygon(data=swap.LN.pca.hull, alpha=F3_FILLALPHA,aes(x=PC1,y=PC2,fill=as.factor(Pheno))) +
      
      geom_polygon(data=swap.LF.pca.hull, alpha=0,aes(x=PC1,y=PC2,fill=as.factor(Pheno),color=as.factor(Population)),size=F3_LINESIZE) +
      geom_polygon(data=swap.LF.pca.hull, alpha=F3_FILLALPHA,aes(x=PC1,y=PC2,fill=as.factor(Pheno))) +
      
      geom_polygon(data=swap.LI.pca.hull, alpha=0,aes(x=PC1,y=PC2,fill=as.factor(Pheno),color=as.factor(Population)),size=F3_LINESIZE) +
      geom_polygon(data=swap.LI.pca.hull, alpha=F3_FILLALPHA,aes(x=PC1,y=PC2,fill=as.factor(Pheno))) +
      
      geom_polygon(data=swap.LS.pca.hull, alpha=0,aes(x=PC1,y=PC2,fill=as.factor(Pheno),color=as.factor(Population)),size=F3_LINESIZE) +
      geom_polygon(data=swap.LS.pca.hull, alpha=F3_FILLALPHA,aes(x=PC1,y=PC2,fill=as.factor(Pheno))) +
      
      
      geom_point(swap.FN.pca.df,
                 mapping = aes(x=PC1,y=PC2,fill=as.factor(Pheno),shape=as.factor(Pheno),color=Population),size=F3_POINTSIZE)+
      geom_point(swap.FF.pca.df,
                 mapping = aes(x=PC1,y=PC2,fill=as.factor(Pheno),shape=as.factor(Pheno),color=Population),size=F3_POINTSIZE)+
      geom_point(swap.FI.pca.df,
                 mapping = aes(x=PC1,y=PC2,fill=as.factor(Pheno),shape=as.factor(Pheno),color=Population),size=F3_POINTSIZE)+
      geom_point(swap.FS.pca.df,
                 mapping = aes(x=PC1,y=PC2,fill=as.factor(Pheno),shape=as.factor(Pheno),color=Population),size=F3_POINTSIZE)+
      geom_point(swap.LN.pca.df,
                 mapping = aes(x=PC1,y=PC2,fill=as.factor(Pheno),shape=as.factor(Pheno),color=Population),size=F3_POINTSIZE)+
      geom_point(swap.LF.pca.df,
                 mapping = aes(x=PC1,y=PC2,fill=as.factor(Pheno),shape=as.factor(Pheno),color=Population),size=F3_POINTSIZE)+
      geom_point(swap.LI.pca.df,
                 mapping = aes(x=PC1,y=PC2,fill=as.factor(Pheno),shape=as.factor(Pheno),color=Population),size=F3_POINTSIZE)+
      geom_point(swap.LS.pca.df,
                 mapping = aes(x=PC1,y=PC2,fill=as.factor(Pheno),shape=as.factor(Pheno),color=Population),size=F3_POINTSIZE)+
      scale_x_continuous(limits = c(F3_XLLIM, F3_XULIM)) + 
      scale_y_continuous(limits = c(F3A_YLLIM, F3A_YULIM)) +
      coord_fixed() + 
      xlab(paste("PC1 -",percent(swap.pca.var.explained[1], accuracy = 0.1))) +
      ylab(paste("PC2 -",percent(swap.pca.var.explained[2], accuracy = 0.1))) +
      ggtitle("TITLE") +
      theme(legend.justification=F3_LGDJUST, 
            legend.position = F3_LGDPOS,
            legend.direction = "vertical", 
            legend.box = "horizontal", 
            legend.spacing.x = unit(F3_LGDSPCX,"cm"),
            legend.spacing.y = unit(F3_LGDSPCY,"cm"),
            legend.key.size = unit(F3_LGDKEYSIZE,"cm"))
    
    swapPHENOS

    swap_sig_PHENOS<-ggplot()+
      scale_color_manual(values = swapColorListPopulation) +
      scale_shape_manual(values = swapPhenoShape) +
      scale_fill_manual(values = swapColorListPheno)+
      
      geom_polygon(data=swap.s.FN.pca.hull, alpha=0,aes(x=PC1,y=PC2,fill=as.factor(Pheno),color=as.factor(Population)),size=F3_LINESIZE) +
      geom_polygon(data=swap.s.FN.pca.hull, alpha=F3_FILLALPHA,aes(x=PC1,y=PC2,fill=as.factor(Pheno))) +
      
      geom_polygon(data=swap.s.FF.pca.hull, alpha=0,aes(x=PC1,y=PC2,fill=as.factor(Pheno),color=as.factor(Population)),size=F3_LINESIZE) +
      geom_polygon(data=swap.s.FF.pca.hull, alpha=F3_FILLALPHA,aes(x=PC1,y=PC2,fill=as.factor(Pheno))) +
      
      geom_polygon(data=swap.s.FI.pca.hull, alpha=0,aes(x=PC1,y=PC2,fill=as.factor(Pheno),color=as.factor(Population)),size=F3_LINESIZE) +
      geom_polygon(data=swap.s.FI.pca.hull, alpha=F3_FILLALPHA,aes(x=PC1,y=PC2,fill=as.factor(Pheno))) +
      
      geom_polygon(data=swap.s.FS.pca.hull, alpha=0,aes(x=PC1,y=PC2,fill=as.factor(Pheno),color=as.factor(Population)),size=F3_LINESIZE) +
      geom_polygon(data=swap.s.FS.pca.hull, alpha=F3_FILLALPHA,aes(x=PC1,y=PC2,fill=as.factor(Pheno))) +
      
      geom_polygon(data=swap.s.LN.pca.hull, alpha=0,aes(x=PC1,y=PC2,fill=as.factor(Pheno),color=as.factor(Population)),size=F3_LINESIZE) +
      geom_polygon(data=swap.s.LN.pca.hull, alpha=F3_FILLALPHA,aes(x=PC1,y=PC2,fill=as.factor(Pheno))) +
      
      geom_polygon(data=swap.s.LF.pca.hull, alpha=0,aes(x=PC1,y=PC2,fill=as.factor(Pheno),color=as.factor(Population)),size=F3_LINESIZE) +
      geom_polygon(data=swap.s.LF.pca.hull, alpha=F3_FILLALPHA,aes(x=PC1,y=PC2,fill=as.factor(Pheno))) +
      
      geom_polygon(data=swap.s.LI.pca.hull, alpha=0,aes(x=PC1,y=PC2,fill=as.factor(Pheno),color=as.factor(Population)),size=F3_LINESIZE) +
      geom_polygon(data=swap.s.LI.pca.hull, alpha=F3_FILLALPHA,aes(x=PC1,y=PC2,fill=as.factor(Pheno))) +
      
      geom_polygon(data=swap.s.LS.pca.hull, alpha=0,aes(x=PC1,y=PC2,fill=as.factor(Pheno),color=as.factor(Population)),size=F3_LINESIZE) +
      geom_polygon(data=swap.s.LS.pca.hull, alpha=F3_FILLALPHA,aes(x=PC1,y=PC2,fill=as.factor(Pheno))) +
      
      
      geom_point(swap.s.FN.pca.df,
                 mapping = aes(x=PC1,y=PC2,fill=as.factor(Pheno),shape=as.factor(Pheno),color=Population),size=F3_POINTSIZE)+
      geom_point(swap.s.FF.pca.df,
                 mapping = aes(x=PC1,y=PC2,fill=as.factor(Pheno),shape=as.factor(Pheno),color=Population),size=F3_POINTSIZE)+
      geom_point(swap.s.FI.pca.df,
                 mapping = aes(x=PC1,y=PC2,fill=as.factor(Pheno),shape=as.factor(Pheno),color=Population),size=F3_POINTSIZE)+
      geom_point(swap.s.FS.pca.df,
                 mapping = aes(x=PC1,y=PC2,fill=as.factor(Pheno),shape=as.factor(Pheno),color=Population),size=F3_POINTSIZE)+
      geom_point(swap.s.LN.pca.df,
                 mapping = aes(x=PC1,y=PC2,fill=as.factor(Pheno),shape=as.factor(Pheno),color=Population),size=F3_POINTSIZE)+
      geom_point(swap.s.LF.pca.df,
                 mapping = aes(x=PC1,y=PC2,fill=as.factor(Pheno),shape=as.factor(Pheno),color=Population),size=F3_POINTSIZE)+
      geom_point(swap.s.LI.pca.df,
                 mapping = aes(x=PC1,y=PC2,fill=as.factor(Pheno),shape=as.factor(Pheno),color=Population),size=F3_POINTSIZE)+
      geom_point(swap.s.LS.pca.df,
                 mapping = aes(x=PC1,y=PC2,fill=as.factor(Pheno),shape=as.factor(Pheno),color=Population),size=F3_POINTSIZE)+
      scale_x_continuous(limits = c(F3_XLLIM, F3_XULIM)) + 
      scale_y_continuous(limits = c(F3A_YLLIM, F3A_YULIM)) +
      coord_fixed() + 
      xlab(paste("PC1 -",percent(swap.pca.var.explained[1], accuracy = 0.1))) +
      ylab(paste("PC2 -",percent(swap.pca.var.explained[2], accuracy = 0.1))) +
      ggtitle("TITLE") +
      theme(legend.justification=F3_LGDJUST, 
            legend.position = F3_LGDPOS,
            legend.direction = "vertical", 
            legend.box = "horizontal", 
            legend.spacing.x = unit(F3_LGDSPCX,"cm"),
            legend.spacing.y = unit(F3_LGDSPCY,"cm"),
            legend.key.size = unit(F3_LGDKEYSIZE,"cm"))
    
    swap_sig_PHENOS
    