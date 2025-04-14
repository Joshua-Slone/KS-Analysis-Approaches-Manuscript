.libPaths("~/R/rlib-4.2.1")

library(haven)
library(lubridate)
library(dplyr)
library(dtplyr)
library(data.table)
library(splines)
library(tidyr)
library(mvtnorm)
library(zoo)
library(RcppRoll)
library(survival)
library(stringr)


######################################
args   = commandArgs(TRUE)
TASKID = as.numeric(args[2])
njob   = TASKID

set.seed(1992 + njob)

#rm(list=ls())


##################################################################################################################
##################################################################################################################
### read in base dataset for bootstrap and MI

d_c=fread("dat_base_mi.csv") 

###########################################################################################################################
###########################################################################################################################
### bootstrap by strata phase 1 and phase 2

samp_strata=c()

for(i in 1:length(unique(d_c[d_c$V==0,]$strata))){
  stratum=unique(d_c[d_c$V==0,]$strata)[i]
  samp <- sample(d_c[d_c$strata==stratum & d_c$V==0,]$patient, length(d_c[d_c$strata==stratum & d_c$V==0,]$patient), replace=TRUE)
  samp_strata=c(samp_strata,samp)
}

samp_strata_v=c()

for(i in 1:length(unique(d_c[d_c$V==1,]$strata))){
  stratum=unique(d_c[d_c$V==1,]$strata)[i]
  samp <- sample(d_c[d_c$strata==stratum & d_c$V==1,]$patient, length(d_c[d_c$strata==stratum & d_c$V==1,]$patient), replace=TRUE)
  samp_strata_v=c(samp_strata_v,samp)
}


samp1=c(samp_strata,samp_strata_v)
df_table_strata <- as.data.frame(table(samp1))
a <- 1
dboot <- data.frame()

while(a <= max(df_table_strata$Freq)){
  
  id_boot <- df_table_strata[df_table_strata$Freq >= a, 1]
  dboot <- rbind(dboot, d_c[d_c$patient %in% id_boot, ])
  a <- a + 1
  
}

#### Honduras KS is fully validated so V=0 has 33 strata and V=1 has 34

##################################################################################################################
##################################################################################################################

##################################################################################################################
##################################################################################################################

#source("boot_mi_og_fxn_ks_site.R")   
source("boot_mi_og_fxn.R")   

imp_analysis <- mi_fxn(dat=dboot,n=10)
coef_prev1 <- imp_analysis$coef_prev
coef_inc1 <- imp_analysis$coef_inc
coef_tte1 <- imp_analysis$coef_tte

save(coef_prev1,coef_inc1,coef_tte1,file=paste0("boot_mi_coefs_",njob,".Rda")) 

