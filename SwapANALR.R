library(dplyr)
library(tidyverse)
setwd("/Users/jackbeltz/Desktop/swap/swap")

####STARVATION DATA####
starv_swap=read.csv("starv.swap.csv")
s=starv_swap

#View(starv_swap)
s <- s %>%
  mutate(mean.starv = ((ifelse(`X16` - `X12` > 0, `X16` - `X12`, 0) * 16) +
                        (ifelse(`X20` - `X16` > 0, `X20` - `X16`, 0) * 20) + 
                         (ifelse(`X24` - `X20` > 0, `X24` - `X20`, 0) * 24) + 
                         (ifelse(`X40` - `X24` > 0, `X40` - `X24`, 0) * 40) + 
                         (ifelse(`X44` - `X40` > 0, `X44` - `X40`, 0) * 44) + 
                         (ifelse(`X48` - `X44` > 0, `X48` - `X44`, 0) * 48) + 
                         (ifelse(`X64` - `X48` > 0, `X64` - `X48`, 0) * 64) + 
                         (ifelse(`X68` - `X64` > 0, `X68` - `X64`, 0) * 68) + 
                         (ifelse(`X72` - `X68` > 0, `X72` - `X68`, 0) * 72) + 
                         (ifelse(`X88` - `X72` > 0, `X88` - `X72`, 0) * 88) + 
                         (ifelse(`X92` - `X88` > 0, `X92` - `X88`, 0) * 92) + 
                         (ifelse(`X96` - `X92` > 0, `X96` - `X92`, 0) * 96) + 
                         (ifelse(`X112` - `X96` > 0, `X112` - `X96`, 0) * 112) + 
                         (ifelse(`X116` - `X112` > 0, `X116` - `X112`, 0) * 116) + 
                         (ifelse(`X120` - `X116` > 0, `X120` - `X116`, 0) * 120) + 
                         (ifelse(`X136` - `X120` > 0, `X136` - `X120`, 0) * 136) + 
                         (ifelse(`X140` - `X136` > 0, `X140` - `X136`, 0) * 140) + 
                         (ifelse(`X144` - `X140` > 0, `X144` - `X140`, 0) * 144) + 
                         (ifelse(`X160` - `X144` > 0, `X160` - `X144`, 0) * 160) + 
                         (ifelse(`X164` - `X160` > 0, `X164` - `X160`, 0) * 164) + 
                         (ifelse(`X168` - `X164` > 0, `X168` - `X164`, 0) * 168)) / count)

starv.cagesum=s %>%
  group_by(cage.number,pheno.treatment) %>%
  dplyr::summarise(label= paste0(cage.number, pheno.treatment),
                   cage.treatment = as.factor(cage.treatment),
                   N=n(),
                   STARV=mean(mean.starv,na.rm = TRUE),
                   STARVsd=sd(mean.starv,na.rm = TRUE),
                   STARVse=STARVsd/sqrt(N)
                   
  )

starv22.cagesum=unique(starv.cagesum)
#View(starv.cagesum)

starv.treatsum=s %>%
  group_by(cage.treatment,pheno.treatment) %>%
  dplyr::summarise(label= paste0(cage.treatment, pheno.treatment),
                   cage.treatment = as.factor(cage.treatment),
                   N=n(),
                   STARV=mean(mean.starv,na.rm = TRUE),
                   STARVsd=sd(mean.starv,na.rm = TRUE),
                   STARVse=STARVsd/sqrt(N)
                   
  )

starv22.treatsum=unique(starv.treatsum)
#View(starv22.treatsum)


####MALE DEVELOPEMNT TIME####
LDm_swap=read.csv("LDm.swap.csv")
ldm=LDm_swap
#View(ldm)
ldm=ldm %>% 
  rowwise %>%
  mutate(flycount = sum(c_across(starts_with("X236"):ends_with("X428")))) 

Mean.ldm=ldm %>% 
  rowwise() %>%
  dplyr::summarise(label=as.factor(ldm$label),
            mean.ldm= ((ldm$'X236' * 236) +
           (ldm$'X240' * 240) +
           (ldm$'X244' * 244) +
           (ldm$'X260' * 260) +
           (ldm$'X264' * 264) +
           (ldm$'X268' * 268) +
           (ldm$'X284' * 284) +
           (ldm$'X288' * 288) +
           (ldm$'X292' * 292) +
           (ldm$'X308' * 308) +
           (ldm$'X312' * 312) +
           (ldm$'X316' * 316) +
           (ldm$'X332' * 332) +
           (ldm$'X336' * 336) +
           (ldm$'X340' * 340) +
           (ldm$'X356' * 356) +
           (ldm$'X364' * 364) +
           (ldm$'X380' * 380) +
           (ldm$'X388' * 388) +
           (ldm$'X404' * 404) +
           (ldm$'X412' * 412) +
           (ldm$'X428' * 428)) / ldm$flycount)
           
Mean.ldm=unique(Mean.ldm)           
ldm$mean.ldm = Mean.ldm$mean.ldm
ldm.cagesum=ldm %>%
  group_by(cage.number,pheno.treatment) %>%
  dplyr::summarise(label= paste0(cage.number, pheno.treatment),
                   cage.treatment = as.factor(cage.treatment),
                   N=n(),
                   LDM=mean(mean.ldm,na.rm = TRUE),
                   LDMsd=sd(mean.ldm,na.rm = TRUE),
                   LDMse=LDMsd/sqrt(N)
                   
  )

ldm22.cagesum=unique(ldm.cagesum)
#View(ldm22.cagesum)


ldm.treatsum=ldm %>%
  group_by(cage.treatment,pheno.treatment) %>%
  dplyr::summarise(label= paste0(cage.treatment, pheno.treatment),
                   cage.treatment = as.factor(cage.treatment),
                   N=n(),
                   LDM=mean(mean.ldm,na.rm = TRUE),
                   LDMsd=sd(mean.ldm,na.rm = TRUE),
                   LDMse=LDMsd/sqrt(N)
                   
  )

ldm22.treatsum=unique(ldm.treatsum)
#View(ldm22.treatsum)

####FEMALE DEVELOPEMNT TIME####
LDf_swap=read.csv("LDf.swap.csv")
ldf=LDf_swap
#View(ldf)
ldf=ldf %>% 
  rowwise %>%
  mutate(flycount = sum(c_across(starts_with("X236"):ends_with("X428")))) 

Mean.ldf=ldf %>% 
  rowwise() %>%
  dplyr::summarise(label=as.factor(ldf$label),
                   mean.ldf= ((ldf$'X236' * 236) +
                                (ldf$'X240' * 240) +
                                (ldf$'X244' * 244) +
                                (ldf$'X260' * 260) +
                                (ldf$'X264' * 264) +
                                (ldf$'X268' * 268) +
                                (ldf$'X284' * 284) +
                                (ldf$'X288' * 288) +
                                (ldf$'X292' * 292) +
                                (ldf$'X308' * 308) +
                                (ldf$'X312' * 312) +
                                (ldf$'X316' * 316) +
                                (ldf$'X332' * 332) +
                                (ldf$'X336' * 336) +
                                (ldf$'X340' * 340) +
                                (ldf$'X356' * 356) +
                                (ldf$'X364' * 364) +
                                (ldf$'X380' * 380) +
                                (ldf$'X388' * 388) +
                                (ldf$'X404' * 404) +
                                (ldf$'X412' * 412) +
                                (ldf$'X428' * 428)) / ldf$flycount)

Mean.ldf=unique(Mean.ldf)           
ldf$mean.ldf = Mean.ldf$mean.ldf
#View(ldf)
ldf.cagesum=ldf %>%
  group_by(cage.number,pheno.treatment) %>%
  dplyr::summarise(label= paste0(cage.number, pheno.treatment),
                   cage.treatment = as.factor(cage.treatment),
                   N=n(),
                   LDF=mean(mean.ldf,na.rm = TRUE),
                   LDFsd=sd(mean.ldf,na.rm = TRUE),
                   LDFse=LDFsd/sqrt(N)
                   
  )

ldf22.cagesum=unique(ldf.cagesum)
#View(ldf22.cagesum)


ldf.treatsum=ldf %>%
  group_by(cage.treatment,pheno.treatment) %>%
  dplyr::summarise(label= paste0(cage.treatment, pheno.treatment),
                   cage.treatment = as.factor(cage.treatment),
                   N=n(),
                   LDF=mean(mean.ldf,na.rm = TRUE),
                   LDFsd=sd(mean.ldf,na.rm = TRUE),
                   LDFse=LDFsd/sqrt(N)
                   
  )

ldf22.treatsum=unique(ldf.treatsum)
#View(ldf22.treatsum)


####DESSICATION DATA####
Dess_swap=read.csv("Dess.swap.csv")
d=Dess_swap
#View(d)
d <- d %>%
  mutate(mean.dess = ((ifelse(`X104` - `X0` > 0, `X104` - `X0`, 0) * 104) +
                         (ifelse(`X248` - `X104` > 0, `X248` - `X104`, 0) * 248) + 
                         (ifelse(`X348` - `X248` > 0, `X348` - `X248`, 0) * 348) + 
                         (ifelse(`X462` - `X348` > 0, `X462` - `X348`, 0) * 462) + 
                         (ifelse(`X582` - `X462` > 0, `X582` - `X462`, 0) * 582) + 
                         (ifelse(`X716` - `X582` > 0, `X716` - `X582`, 0) * 716) + 
                         (ifelse(`X825` - `X716` > 0, `X825` - `X716`, 0) * 825) + 
                         (ifelse(`X939` - `X825` > 0, `X939` - `X825`, 0) * 939) + 
                         (ifelse(`X1060` - `X939` > 0, `X1060` - `X939`, 0) * 1060) + 
                         (ifelse(`X1201` - `X1060` > 0, `X1201` - `X1060`, 0) * 1201) + 
                         (ifelse(`X1316` - `X1201` > 0, `X1316` - `X1201`, 0) * 1316) + 
                         (ifelse(`X1492` - `X1316` > 0, `X1492` - `X1316`, 0) * 1492) + 
                         (ifelse(`X1537` - `X1492` > 0, `X1537` - `X1492`, 0) * 1537)) / vial.count)

dess.cagesum=d %>%
  group_by(cage.number,pheno.treatment) %>%
  dplyr::summarise(label= paste0(cage.number, pheno.treatment),
                   cage.treatment = as.factor(cage.treatment),
                   N=n(),
                   DESS=mean(mean.dess,na.rm = TRUE),
                   DESSsd=sd(mean.dess,na.rm = TRUE),
                   DESSse=DESSsd/sqrt(N)
                   
  )

dess22.cagesum=unique(dess.cagesum)
#View(dess22.cagesum)


dess.treatsum=d %>%
  group_by(cage.treatment,pheno.treatment) %>%
  dplyr::summarise(label= paste0(cage.treatment, pheno.treatment),
                   cage.treatment = as.factor(cage.treatment),
                   N=n(),
                   DESS=mean(mean.dess,na.rm = TRUE),
                   DESSsd=sd(mean.dess,na.rm = TRUE),
                   DESSse=DESSsd/sqrt(N)
                   
  )

dess22.treatsum=unique(dess.treatsum)
#View(dess22.treatsum)

####PIGMENTATION DATA####
Pig_swap=read.csv("pig.swap.csv")
pig=Pig_swap
Mean.pig=pig %>% 
  rowwise() %>%
  dplyr::summarise(sd.pig= sd(c(`X1.sum`,`X2.sum`,`X3.sum`,`X4.sum`,`X5.sum`,`X6.sum`,`X7.sum`,`X8.sum`,`X9.sum`,`X10.sum`,
                                `X11.sum`,`X12.sum`,`X13.sum`,`X14.sum`,`X15.sum`,`X16.sum`,`X17.sum`,`X18.sum`,`X19.sum`,`X20.sum`)),
                   se.pig = sd.pig/(sqrt(15)),
                   mean.pig= mean( c(`X1.sum`,`X2.sum`,`X3.sum`,`X4.sum`,`X5.sum`,`X6.sum`,`X7.sum`,`X8.sum`,`X9.sum`,`X10.sum`,
                                     `X11.sum`,`X12.sum`,`X13.sum`,`X14.sum`,`X15.sum`,`X16.sum`,`X17.sum`,`X18.sum`,`X19.sum`,`X20.sum`)))

pig$mean.pig = Mean.pig$mean.pig   
pig$se.pig =Mean.pig$se.pig 
pig.cagesum=pig %>%
  group_by(cage.number,pheno.treatment) %>%
  dplyr::summarise(label= paste0(cage.number, pheno.treatment),
                   cage.treatment = as.factor(cage.treatment),
                   PIG=mean(mean.pig,na.rm = TRUE),
                   PIGse=mean(se.pig,na.rm = TRUE)
                   
  )

pig22.cagesum=unique(pig.cagesum)
#View(pig22.cagesum)

pig.treatsum=pig %>%
  group_by(cage.treatment,pheno.treatment) %>%
  dplyr::summarise(label= paste0(cage.treatment, pheno.treatment),
                   cage.treatment = as.factor(cage.treatment),
                   PIG=mean(mean.pig,na.rm = TRUE),
                   PIGse=mean(se.pig,na.rm = TRUE)
                   
  )

pig22.treatsum=unique(pig.treatsum)
#View(pig22.treatsum)


####DRY WEIGHT / LIPID DATA ####
DW_swap=read.csv("DW.swap.csv")
dw=DW_swap

dw.cagesum=dw %>%
  group_by(cage.number,pheno.treatment) %>%
  dplyr::summarise(label= paste0(cage.number, pheno.treatment),
                   cage.treatment = as.factor(cage.treatment),
                   N=n(),
                   DW=mean(dry.weight,na.rm = TRUE),
                   DWsd=sd(dry.weight,na.rm = TRUE),
                   DWse=DWsd/sqrt(N),
                   LW=mean(lipid.weight,na.rm = TRUE),
                   LWsd=sd(lipid.weight,na.rm = TRUE),
                   LWse=LWsd/sqrt(N),
                   
  )

dw22.cagesum=unique(dw.cagesum)
#View(dw22.cagesum)


dw.treatsum=dw %>%
  group_by(cage.treatment,pheno.treatment) %>%
  dplyr::summarise(label= paste0(cage.treatment, pheno.treatment),
                   cage.treatment = as.factor(cage.treatment),
                   N=n(),
                   DW=mean(dry.weight,na.rm = TRUE),
                   DWsd=sd(dry.weight,na.rm = TRUE),
                   DWse=DWsd/sqrt(N),
                   LW=mean(lipid.weight,na.rm = TRUE),
                   LWsd=sd(lipid.weight,na.rm = TRUE),
                   LWse=LWsd/sqrt(N),
                   
  )

dw22.treatsum=unique(dw.treatsum)
#View(dw22.treatsum)


####ALL PHENOTYPES####
cagesum_list <- list(starv22.cagesum, ldm22.cagesum, ldf22.cagesum, dess22.cagesum, pig22.cagesum, dw22.cagesum)
treatsum_list <- list(starv22.treatsum, ldm22.treatsum, ldf22.treatsum, dess22.treatsum, pig22.treatsum, dw22.treatsum)

allpheno.cagesum=cagesum_list %>% reduce(full_join, by=c("label","cage.number","cage.treatment","pheno.treatment"))
allpheno22.cagesum <- allpheno.cagesum[ -c(5,9,13,17,23) ]
View(allpheno22.cagesum) ### grouped by cage #
write.csv(allpheno22.cagesum, '/Users/jackbeltz/Desktop/swap/swap/allpheno22.cagesum.csv', row.names = FALSE)



allpheno.treatsum=treatsum_list %>% reduce(full_join, by=c("label","cage.treatment","pheno.treatment"))
allpheno22.treatsum <- allpheno.treatsum[ -c(4,8,12,16,22) ]
View(allpheno22.treatsum) ### grouped by treatment group
write.csv(allpheno22.treatsum, '/Users/jackbeltz/Desktop/swap/swap/allpheno22.treatsum.csv', row.names = FALSE)
