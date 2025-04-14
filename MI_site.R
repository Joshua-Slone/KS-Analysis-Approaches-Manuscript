mi_fxn <- function(dat, n){
  
  betas_mi_prev=betas_mi_inc=betas_mi_tte= vector("list", n)
  
  for(j in 1:n){
    
    ## Impute exclusion of patients with KS diagnosis before 2009-12-01; 1 patient from EA should be excluded
    d = dat |>mutate(exclusion_flag=if_else(is.na(ks_y_final) & V==0 & ks_y==1 & program==1,rbinom(length(program),1,1/265),0))|> filter(exclusion_flag==0) |>mutate(exclusion_flag=NULL)
    
    
    
    #####################################################################################################################
    #####################################################################################################################
    ### birth_d_error binomial draws by region
    

    d = d |>mutate(birth_d_error_imp_gr=if_else(program==1,rbinom(length(program),1,65/825),rbinom(length(program),1,1/114))) |> mutate(birth_d_error_imp=if_else(V==0, birth_d_error_imp_gr, birth_d_error)) 


    #####################################################################################################################
    #####################################################################################################################
    
    
    
    
    # ### birth_d_error full model
    # 
    # mod1 = glm(birth_d_error ~ death_y + ks_y + male + program, data= d[d$V==1,], family = binomial) 
    # # mod1 = glm(birth_d_error ~ death_y +  enrol_yr + ks_y + age_at_start + male + program + art_y, data= d[d$V==1,], family = binomial) 
    # 
    # ## to easily get design matrix for the full cohort  
    # d_1 = d |> mutate(birth_d_error=1,birth_d_error_days=1,start_error=1,stop_error=1,start_error_days=1,stop_error_days=1,enrol_d_error=1,enrol_d_error_days=1) 
    # v1=model.matrix(terms(mod1),d_1)
    # beta1 = mod1$coefficients
    # vcov1 = vcov(mod1)
    # rmvs1=rmvnorm(n=1,beta1,vcov1)
    # rmvs1 = as.vector(rmvs1)
    # ## predicted values
    # Xbeta1 = v1 %*% rmvs1
    # ## to get binary death_y, transform log-odds to probability
    # prob1=exp(Xbeta1)/(1+exp(Xbeta1))
    # birth_d_error_new=rbinom(length(prob1),1,prob1)
    # d$birth_d_error_imp_gr=birth_d_error_new
    # d$birth_d_error_imp=if_else(!is.na(d$birth_d_error) & d$V==1,d$birth_d_error,d$birth_d_error_imp_gr)
    
    ### birth_d_error days
    
    #### CCASAnet only has 1 birth_d in error so program was removed as a covariate, ART status as well given it should not relate to the errors themselves
    
    mod3 = lm(birth_d_error_days ~ death_y + ks_y + male, data= d[d$V==1 & d$birth_d_error==1,]) 
    #mod3 = lm(birth_d_error_days ~ death_y +  ns(enrol_yr,df=3) + ks_y + ns(age_at_start,df=3) + ns(total_days,df=3) + male + program + art_y, data= d[d$V==1 & d$birth_d_error==1,]) 
    d_1 = d |> mutate(birth_d_error=1,birth_d_error_days=1,start_error=1,stop_error=1,start_error_days=1,stop_error_days=1,enrol_d_error=1,enrol_d_error_days=1) 
    
    
    ## to easily get design matrix for the full cohort 
    v2=model.matrix(terms(mod3),d_1)
    beta2 = mod3$coefficients
    vcov2 = vcov(mod3)
    rmvs2=rmvnorm(n=1,beta2,vcov2)
    rmvs2 = as.vector(rmvs2)
    birth_d_error_days_new = as.vector((v2 %*% rmvs2) + sample(resid(mod3),1))
    
    d$birth_d_error_days_imp_gr=if_else(d$birth_d_error_imp_gr==1,birth_d_error_days_new,0)
    
    d$birth_d_error_days_imp=if_else(!is.na(d$birth_d_error_days) & d$V==1,d$birth_d_error_days,if_else(d$birth_d_error_imp==0,0,birth_d_error_days_new))
    
    d = d |> mutate(birth_d_imp=if_else(V==1,as.Date(birth_d_final),as.Date(birth_d) + birth_d_error_days_imp))
    
    ### New for updated enrol_d1
    d = d |> mutate(id=row_number(),start_d1 = if_else(V==1,start_d_final,start_d),stop_d1 = if_else(V==1,l_alive_d_final,stop_d), enrol_d1=  if_else(V==1,enrol_d_final,enrol_d))
    
    ###############################################################################
    
    #### expanded monthly data
    
    #add rows for each month
    setDT(d)
    # floor start and end month
    d[, start_month := lubridate::floor_date(start_d1, unit = "month")]
    d[, end_month := lubridate::floor_date(stop_d1, unit = "month")]
    
    #create sequence by month by id
    d1<-d[, .(year_month = format(seq(start_month, end_month, by = "1 month"), "%Y-%m")), by = .(patient,id,birth_d,birth_d_final,birth_d_imp,enrol_d,enrol_d_final,male, male_final,dis_d,dis_year_month= format(as.Date(dis_d), "%Y-%m"), canc_d_final,canc_final_year_month= format(as.Date(canc_d_final), "%Y-%m"),ks_y,ks_y_final,art_sd,art_year_month= format(as.Date(art_sd), "%Y-%m"), recart_d_final,art_final_year_month= format(as.Date(recart_d_final), "%Y-%m"),death_y,death_y_final,death_d, death_year_month= format(as.Date(death_d), "%Y-%m"),death_d_final,death_final_year_month= format(as.Date(death_d_final), "%Y-%m"),l_alive_d_final, start_d,start_d_final, stop_d, start_d1,stop_d1,enrol_d1,enrol_d1_year_month= format(as.Date(enrol_d1), "%Y-%m") , program,V)]
    
    ### merge cd4 values in by year_month from phase 1 and phase 2 long datasets and create variables below, carry forward cd4 12 months and then cutoff at 2010 forward
    load("all_cd4_1.Rda") 
    load("all_val_long4.Rda")
    d2<-merge(x=d1,y=all_val_long4[,c("patient","year_month","month_avg_cd4_final")],by=c("patient","year_month"),all.x = TRUE)
    
    d3<-merge(x=d2,y=all_cd4_1[,c("patient","year_month","month_avg_cd4")],by=c("patient","year_month"),all.x = TRUE)
    
    rm(d1,d2); invisible(gc())
    
    ### CD4 carry forward 12 months and back w/ in 6 months of enrollment, cd4 is value with backfill and carry forward
    d3=d3 |> lazy_dt() %>% group_by(patient,id) %>% mutate(new_v = if_else(is.na(month_avg_cd4_final),NA,row_number())) |> fill(new_v) |> mutate(row_diff = row_number()-new_v) |> fill(month_avg_cd4_final) |> mutate(month_avg_cd4_final=if_else(row_diff <= 12, month_avg_cd4_final,NA))  |> mutate(row=row_number(),cd4_1=month_avg_cd4_final) |> fill(cd4_1,.direction = "up") |> fill(new_v,.direction = "up") |> mutate(cd4_final=if_else(row < 7 & new_v <=7, cd4_1, month_avg_cd4_final)) |> ungroup() |> as.data.frame()
    
    d3=d3 |> lazy_dt() %>% group_by(patient,id) %>% mutate(new_v = if_else(is.na(month_avg_cd4),NA,row_number())) |> fill(new_v) |> mutate(row_diff = row_number()-new_v) |> fill(month_avg_cd4) |> mutate(month_avg_cd4=if_else(row_diff <= 12, month_avg_cd4,NA)) |> filter(year_month>="2010-01",year_month<="2019-12") |> mutate(row=row_number(),cd4_1=month_avg_cd4) |> fill(cd4_1,.direction = "up") |> fill(new_v,.direction = "up") |> mutate(cd4=if_else(row < 7 & new_v <=7, cd4_1, month_avg_cd4)) |> ungroup() |> as.data.frame()
    
    
    ### variables needed for imputation: death_y male program age art_status ks_ind
    
    ## ART has started indicator
    d3$art_ind<-ifelse(is.na(d3$art_year_month) | d3$art_year_month > d3$year_month,0,1)
    
    ## Death indicator (w.r.t. end of follow-up)
    d3$death_ind<-ifelse(is.na(d3$death_year_month) | d3$death_year_month != d3$year_month,0,1)
    
    ## Age each month
    d3$age_imp<-as.numeric(round((as.Date(paste(d3$year_month, "-01", sep=""))-d3$birth_d_imp)/365.25,1))
    
    d3$enrol_ind_1<-if_else(d3$year_month>d3$enrol_d1_year_month | d3$year_month==d3$enrol_d1_year_month,1,0)
    
    d3 = d3  |> lazy_dt() |> mutate(enrol_ind_month_d1=if_else(cumsum(enrol_ind_1)==1,1,0)) |> ungroup() |> as.data.frame()

    ## ART has started indicator
    d3$art_ind_final<-ifelse(d3$V==1 & (d3$art_final_year_month > d3$year_month | is.na(d3$art_final_year_month)),0,ifelse(d3$V==1 & (d3$art_final_year_month <= d3$year_month),1,NA))
    
    ## Death indicator (w.r.t. end of follow-up)
    d3$death_ind_final<-ifelse(d3$V==1 & (d3$death_final_year_month != d3$year_month | is.na(d3$death_final_year_month)),0,ifelse(d3$V==1 & d3$death_final_year_month == d3$year_month,1,NA))
    
    ## Age each month
    d3$age_final<-as.numeric(round((as.Date(paste(d3$year_month, "-01", sep=""))-as.Date(d3$birth_d_final))/365.25,1))
    
    ## KS indicator
    #d3 = d3 |> mutate(ks_ind = if_else(!is.na(dis_d) & dis_year_month <= year_month & as.Date(dis_d) >= as.Date("2010-01-01"),1,0))
    d3 = d3 |> lazy_dt() |> mutate(ks_ind = if_else(!is.na(dis_d) & dis_year_month <= year_month,1,0)) |> ungroup() |> as.data.frame()
    
    ## KS indicator validated (one patient had canc_d_final in 7/2009, make not KS diagnosis)
    d3 = d3 |> lazy_dt() |> mutate(ks_ind_final = if_else(V==0,NA, if_else(V==1 & !is.na(canc_d_final) & canc_final_year_month <= year_month & as.Date(canc_d_final) >= as.Date("2009-12-01"),1,0))) |> ungroup() |> as.data.frame()
    
    
    ## phase 1 and phase 2 indicating rows that need imputation for cd4 
    
    d3=d3 |> lazy_dt() |> group_by(patient,id) |> mutate(flag=if_else((!is.na(lag(cd4)) | row==1) & is.na(cd4),1,0)) |> mutate(phase1_imp_ind=if_else(flag==1 | (is.na(cd4) & (lag(flag,13)==1 & roll_sumr(cd4,13,na.rm=TRUE) == 0)) | (is.na(cd4) & (lag(flag,26)==1 & roll_sumr(cd4,26,na.rm=TRUE) == 0)) | (is.na(cd4) & (lag(flag,39)==1 & roll_sumr(cd4,39,na.rm=TRUE) == 0)) | (is.na(cd4) & (lag(flag,52)==1 & roll_sumr(cd4,52,na.rm=TRUE) == 0)) | (is.na(cd4) & (lag(flag,65)==1 & roll_sumr(cd4,65,na.rm=TRUE) == 0)) | (is.na(cd4) & (lag(flag,78)==1 & roll_sumr(cd4,78,na.rm=TRUE) == 0)) | (is.na(cd4) & (lag(flag,91)==1 & roll_sumr(cd4,91,na.rm=TRUE) == 0)) | (is.na(cd4) & (lag(flag,104)==1 & roll_sumr(cd4,104,na.rm=TRUE) == 0))| (is.na(cd4) & (lag(flag,117)==1 & roll_sumr(cd4,117,na.rm=TRUE) == 0)) ,1,0)) |> mutate(phase1_imp_ind=if_else(!is.na(phase1_imp_ind),phase1_imp_ind,0))  |> ungroup() |> as.data.frame()
    
    d4=d3 |> lazy_dt() |> filter(V==1) |> group_by(patient,id) |> mutate(flag=if_else((!is.na(lag(cd4_final)) | row==1) & is.na(cd4_final),1,0)) |> mutate(phase2_imp_ind=if_else(flag==1 | (is.na(cd4_final) & (lag(flag,13)==1 & roll_sumr(cd4_final,13,na.rm=TRUE) == 0)) | (is.na(cd4_final) & (lag(flag,26)==1 & roll_sumr(cd4_final,26,na.rm=TRUE) == 0)) | (is.na(cd4_final) & (lag(flag,39)==1 & roll_sumr(cd4_final,39,na.rm=TRUE) == 0)) | (is.na(cd4_final) & (lag(flag,52)==1 & roll_sumr(cd4_final,52,na.rm=TRUE) == 0)) | (is.na(cd4_final) & (lag(flag,65)==1 & roll_sumr(cd4_final,65,na.rm=TRUE) == 0)) | (is.na(cd4_final) & (lag(flag,78)==1 & roll_sumr(cd4_final,78,na.rm=TRUE) == 0)) | (is.na(cd4_final) & (lag(flag,91)==1 & roll_sumr(cd4_final,91,na.rm=TRUE) == 0)) | (is.na(cd4_final) & (lag(flag,104)==1 & roll_sumr(cd4_final,104,na.rm=TRUE) == 0))| (is.na(cd4_final) & (lag(flag,117)==1 & roll_sumr(cd4_final,117,na.rm=TRUE) == 0)) ,1,0)) |> mutate(phase2_imp_ind=if_else(!is.na(phase2_imp_ind) ,phase2_imp_ind,0))  |> ungroup() |> as.data.frame()
    
    d5=d3 |> lazy_dt() |> filter(V!=1) |> mutate(phase2_imp_ind=0)  |> ungroup() |> as.data.frame()
    
    d6 = rbind(d4,d5) |> arrange(patient,id)
    
    rm(d4,d5); invisible(gc())
    
    d6 = d6 |> lazy_dt() |> mutate(cd4_rt = sqrt(cd4), cd4_final_rt = sqrt(cd4_final)) |> ungroup() |> as.data.frame()
    
    d6$age<-as.numeric(round((as.Date(paste(d6$year_month, "-01", sep=""))-as.Date(d6$birth_d))/365.25,1))
    d6$year<-as.numeric(substr(d6$year_month,1,4)) 
    
    #####################################################################################################################
    
    ### imputation for phase 1 cd4 values based on error prone variables
    mod5 = lm(cd4_rt ~  male +program + ns(age_imp,df=4) +art_ind +ks_ind+ death_ind, data= d6[!is.na(d6$cd4_rt),]) 
    #mod5 = lm(cd4_rt ~  male +program + ns(age_imp,df=4) +art_ind + death_ind +ks_ind, data= d6[!is.na(d6$cd4_rt),]) 
    #ch=d6[!is.na(d6$cd4_rt),]
    
    ## to easily get design matrix for the full cohort
    d6_1 = d6 |> mutate(cd4_rt=1)
    v3=model.matrix(terms(mod5),d6_1)
    beta2 = mod5$coefficients
    vcov2 = vcov(mod5)
    rmvs2=rmvnorm(n=1,beta2,vcov2)
    rmvs2 = as.vector(rmvs2)
    cd4_rt_new = as.vector((v3 %*% rmvs2) + sample(resid(mod5),1))
    d6$cd4_rt_init_imp=if_else(d6$phase1_imp_ind==1,cd4_rt_new,d6$cd4_rt)
    ### any imputed CD4<0 assign a 1
    d6$cd4_rt_init_imp=if_else(d6$cd4_rt_init_imp<0,1,d6$cd4_rt_init_imp)
    d7=d6 %>% group_by(patient,id) |> fill(cd4_rt_init_imp)
    
    
    ## imputation for sex
    
    ### imputation for sex (male)
    #### using PPV and NPV for male indicator from phase II
    #### Stratified by program
    
    #ch= d |> filter(V==1,program==1) 
    #ch1= d |> filter(V==1,program==0) 
    
    #table(ch$male,ch$male_final) 
    #2/825 for program==1
    #PPV=329/331
    #NPV=494/494
    
    #ch1= d |> filter(V==1,program==0) 
    #table(ch1$male,ch1$male_final) 
    #0/103 for program==0
    
    
    # mod7 = glm(male_final ~ program + art_ind + ks_ind + male+ ns(age_imp,df=3) + death_ind  + ns(year,df=3)  + ns(cd4_rt_init_imp,df=3), data= d7[d7$V==1,], family = binomial) 
    # 
    # mod7_1=summary(update(mod7, method = "brglm_fit"))
    
    ch = d7 |> lazy_dt() |> count(patient,id,male,program,V,male_final) |>mutate(male_imp_gr=if_else(male==1 & program==1,rbinom(length(program),1,329/331),male)) |> mutate(male_imp=if_else(V==0, male_imp_gr, male_final)) |> ungroup() |> as.data.frame()
    
    d7 = d7 |> left_join(ch[,c('patient','id','male_imp','male_imp_gr')],by=c('patient','id'))
    
    #####################################################################################################################
    
    ### imputation for phase 2 cd4 values
    
    mod6 = lm(cd4_final_rt ~  male_final +program+ ns(age_final,df=4) + ns(cd4_rt_init_imp,df=4) +art_ind_final + death_ind_final +ks_ind_final, data= d7[d7$V==1,])
    ## to easily get design matrix for phase 2
    d7_1 = d7[d7$V==1,] |> mutate(cd4_final_rt=1)
    v4=model.matrix(terms(mod6),d7_1)
    beta2 = mod6$coefficients
    vcov2 = vcov(mod6)
    rmvs2=rmvnorm(n=1,beta2,vcov2)
    rmvs2 = as.vector(rmvs2)
    cd4_final_rt_new = as.vector((v4 %*% rmvs2) + sample(resid(mod6),1))
    d7_phase2 = d7[d7$V==1,]
    d7_phase2$cd4_rt_imp_final=if_else(d7_phase2$phase2_imp_ind==1,cd4_final_rt_new,d7_phase2$cd4_final_rt)
    
    d7_phase1 = d7 |> lazy_dt() |> filter(V!=1) |> mutate(cd4_rt_imp_final=NA) |> ungroup() |> as.data.frame()
    
    d8 = rbind(d7_phase1,d7_phase2) |> lazy_dt() |> arrange(patient,id) |> group_by(patient,id) |> fill(cd4_rt_imp_final) |> ungroup() |> as.data.frame()
    
    rm(d6,d6_1,d7,d7_phase1); invisible(gc())
    
    ### imputation for phase 1 cd4 values based on error free phase 2 variables
    
    mod6 = lm(cd4_rt_imp_final ~ male_imp +program + ns(age,df=4) + ns(age_imp,df=4)  + ns(cd4_rt_init_imp,df=4) + art_ind + death_ind +ks_ind, data= d8[d8$V==1,])
    ## to easily get design matrix for phase 2
    d8_1 = d8 |> mutate(cd4_rt_imp_final=1)
    v4=model.matrix(terms(mod6),d8_1)
    beta2 = mod6$coefficients
    vcov2 = vcov(mod6)
    rmvs2=as.vector(rmvnorm(n=1,beta2,vcov2))
    d8$cd4_rt_imp_gr = as.vector((v4 %*% rmvs2) + sample(resid(mod6),1))
    d8$cd4_rt_imp = if_else(!is.na(d8$cd4_rt_imp_final) & d8$V==1,d8$cd4_rt_imp_final, d8$cd4_rt_imp_gr)
    
    rm(d8_1); invisible(gc())
    
    ###############################################################################################################################################
    ###############################################################################################################################################
    
    #### Imputing KS indicator by site
    
    ### imputation for KS indicator
    #### using PPV and NPV for KS indicator from phase II
    #### Stratified by program: using PPV and NPV for KS indicator from phase II
    ### Bring in the strata to get the sites, d3 is now the dataset with the strata
    load("wgts_by_hand.Rda")
    d3$site=word(d3$strata,1)
    d3$site=str_extract(d3$site, "[^-]+")
    d8 = d8 |> left_join(d3[,c('patient','site')])
    
    load("ks_probs.Rda") 
    
    d8$npv <- npv[match(d8$site, names(npv))]
    d8$ppv <- ppv[match(d8$site, names(ppv))]
    
    d8 = d8 |> lazy_dt() |>mutate(ks_ind_imp_gr=if_else(ks_ind==1,rbinom(length(ppv),1,ppv),rbinom(length(npv),1,1-npv))) |> mutate(ks_ind_imp=if_else(V==0, ks_ind_imp_gr, ks_ind_final)) |> ungroup() |> as.data.frame()
     
    
    #####################################################################################################################
    
    
    
    d8 = d8 |> lazy_dt() |> group_by(patient,id) |> mutate(last_month=if_else(row_number()==n(),1,0), ks_ind_imp=if_else(cumsum(ks_ind_imp)>=1,1,0)) |> ungroup() |> as.data.frame()
    
    
    #####################################################################################################################
    #####################################################################################################################
    #####################################################################################################################
    
    
    ## imputation for art indicator
    mod7 = glm(art_ind_final ~ art_ind + death_ind +ks_ind_imp + male_imp + program  + ns(age_imp,df=3) + ns(cd4_rt_imp,df=3), data= d8[d8$V==1,], family = binomial) 

    #mod7 = glm(art_ind_final ~ art_ind + death_ind +ks_ind_imp + male_imp + program  + ns(age_imp,df=4) + year + ns(cd4_rt_imp,df=4) + ns(cd4_rt_init_imp,df=4), data= d8[d8$V==1,], family = binomial) 
    d8_1 = d8 |> mutate(art_ind_final=1)
    v4=model.matrix(terms(mod7),d8_1)
    beta3 = mod7$coefficients
    vcov3 = vcov(mod7)
    rmvs3=rmvnorm(n=1,beta3,vcov3)
    rmvs3 = as.vector(rmvs3)
    ## predicted values
    Xbeta3 = v4 %*% rmvs3
    prob2=exp(Xbeta3)/(1+exp(Xbeta3))
    d8$art_ind_imp_gr=rbinom(length(prob2),1,prob2)
    d8$art_ind_imp=if_else(!is.na(d8$art_ind_final) & d8$V==1,d8$art_ind_final,d8$art_ind_imp_gr)
    
    d8 = d8 |> lazy_dt() |> group_by(patient,id) |>  mutate(art_ind_imp=if_else(cumsum(art_ind_imp)>=1,1,0)) |>  ungroup() |> as.data.frame()
    
    rm(d8_1); invisible(gc())
    

    ## imputation for death indicator by region NPV and PPV
    
     d8 = d8 |>mutate(death_ind_imp_gr=if_else(last_month==1 & program==1 & death_ind==1,rbinom(length(program),1,75/83),if_else(last_month==1 & program==1 & death_ind==0,rbinom(length(program),1,1-(740/742)),if_else(last_month==1 & program==0 & death_ind==1,rbinom(length(program),1,17/21),if_else(last_month==1 & program==0 & death_ind==0,0,0))))) |> mutate(death_ind_imp=if_else(V==0, death_ind_imp_gr, death_ind_final)) 
    
    
      ###############################################################################################################################################
    ###############################################################################################################################################
    ### MI KS prevalence
    
    d_prev = d8 |> lazy_dt() |> filter(as.Date(enrol_d1)>=as.Date("2010-01-01") | (as.Date(enrol_d1)>=as.Date("2009-12-01") & program==1)) |>
      group_by(patient,id) |> mutate(ks_d_imp=if_else(cumsum(ks_ind_imp)==1 & is.na(canc_d_final),as.Date(paste(year_month, "-01", sep="")),if_else(!is.na(canc_d_final),as.Date(canc_d_final),NA))) |> filter(as.Date(ks_d_imp) >= as.Date(enrol_d1) - 60) |>
      mutate(prevalent = if_else(is.na(ks_d_imp),0,if_else(as.Date(ks_d_imp) <= as.Date(enrol_d1) + 60 & as.Date(ks_d_imp) >= as.Date(enrol_d1) - 60,1,0))) |>
      filter(enrol_d1_year_month==year_month | (row_number() ==1 & enrol_d1_year_month=="2009-12")) |> mutate(year_imp = as.numeric(substr(year_month,1,4)) ) |>
      mutate(age1= age_imp ,year1=year_imp,art_ind1=if_else(is.na(recart_d_final),0,if_else(as.Date(recart_d_final)<as.Date(enrol_d_final),1,0)),cd4_rt_imp_final1=cd4_rt_imp_final) |> ungroup() |> as.data.frame()
    
    fit_mi_prev = glm(prevalent ~ male_imp + program + I(age_imp/10) + year_imp + art_ind_imp + cd4_rt_imp, data= d_prev, family = binomial) 
    
    ###############################################################################################################################################
    ###############################################################################################################################################
    ### MI KS incidence
    
    d11 = d8 |> select(patient,id,program,year_month,art_ind_imp,male_imp,enrol_d1_year_month,cd4_rt_imp,birth_d_imp,ks_ind_imp)
    
    d11$year<-as.numeric(substr(d11$year_month,1,4)) 
    
    d11$cd4_rt_imp1= round(d11$cd4_rt_imp)
    
    ### follow up starts at enrollment
    suppressWarnings({d_inc_a =   d11 |> lazy_dt() |> group_by(patient,id)  |> filter(year_month>=enrol_d1_year_month & year_month>="2010-01") |> mutate(ct=rleid(cd4_rt_imp1),date=as.Date(paste(year_month, "-01", sep="")),enrol_d1_year_month_date=as.Date(paste(enrol_d1_year_month, "-01", sep=""))) |> arrange(patient,id,year_month) |>  filter(row_number() <= min(which(ks_ind_imp == 1))) |> ungroup() |> as.data.frame()})
    
    d_inc_a1 = d_inc_a |> lazy_dt() |> mutate(prevalent = if_else(ks_ind_imp==1 & date <= enrol_d1_year_month_date + 60,1,0))|> group_by(patient,id) |> filter(all(prevalent == 0)) |> ungroup() |> as.data.frame()
    
    d_inc_a2=d_inc_a1 |> lazy_dt() |> group_by(patient,id,program,year,art_ind_imp,male_imp,cd4_rt_imp1,ct,birth_d_imp)|> summarise(ks_ind_imp=max(ks_ind_imp),min_date=min(year_month),max_date=max(year_month)) |>  arrange(patient,id,min_date)  |> ungroup() |> as.data.frame()
    
    
    #d_inc_a2$age_imp<-as.numeric(as.Date(paste(d_inc_a2$max_date, "-01", sep="")) - as.Date(d_inc_a2$birth_d_imp))/365.25 

    d_inc_a2= d_inc_a2 |> group_by(patient,id) |> mutate(age_imp = dplyr::first(as.numeric(as.Date(paste(min_date, "-01", sep="")) - as.Date(birth_d_imp))/365.25))
    
        
    d_inc_a2$n = as.numeric(substr(d_inc_a2$max_date,6,7))-as.numeric(substr(d_inc_a2$min_date,6,7)) + 1
    d_inc_a2$pys=d_inc_a2$n/12
    d_inc_a2$year<-as.numeric(substr(d_inc_a2$min_date,1,4))
    
    fit_mi_inc = glm(ks_ind_imp ~ male_imp + program + I(age_imp/10) + I(year-2010) + art_ind_imp + cd4_rt_imp1, offset = log(pys/1000), family = poisson, data = d_inc_a2)
    
    
    ###############################################################################################################################################
    ###############################################################################################################################################
    ### MI time after KS diagnosis to death
    
    d_tte = d8 |> group_by(patient,id) |> filter(any(ks_ind_imp==1)) |> mutate(ks_d_imp=if_else(cumsum(ks_ind_imp)==1 & is.na(canc_d_final),as.Date(paste(year_month, "-01", sep="")),if_else(!is.na(canc_d_final),as.Date(canc_d_final),NA))) |>  filter(as.Date(ks_d_imp) >= as.Date(enrol_d1) - 60) |> mutate(censor_d_imp=if_else(row_number()==n() & death_ind_imp !=1,as.Date(paste(year_month, "-01", sep="")) %m+% months(1),NA),death_d_imp=if_else(death_ind_imp==1,as.Date(paste(year_month, "-01", sep="")) %m+% months(1),NA)) |> fill(c(death_d_imp,censor_d_imp),.direction = "downup") |> filter(cumsum(ks_ind_imp)==1) |> mutate(ks_d_imp=if_else(is.na(canc_d_final),as.Date(paste(year_month, "-01", sep="")),as.Date(canc_d_final)), time_imp=if_else(is.na(death_d_imp),as.numeric(as.Date(censor_d_imp)-as.Date(ks_d_imp)+1),as.numeric(as.Date(death_d_imp)-as.Date(ks_d_imp)+1)),canc_yr_imp=as.numeric(format(ks_d_imp,'%Y')), death_ind_imp1=if_else(is.na(death_d_imp),0,1))
    
    fit_mi_tte = coxph(Surv(time_imp,death_ind_imp) ~ male_imp+program+I(age_imp/10)+canc_yr_imp+art_ind_imp + cd4_rt_imp, data = d_tte)
    
    
    ###############################################################################################################################################
    ###############################################################################################################################################
    ### MI stored results for each MI run for each bootstrap replicate
    
    betas_mi_prev[[j]] = coef(fit_mi_prev)
    betas_mi_inc[[j]] = coef(fit_mi_inc)
    betas_mi_tte[[j]] = coef(fit_mi_tte)
    
    
    
    
  }
  
  ddest <- list("coef_prev" = betas_mi_prev, "coef_inc" = betas_mi_inc,"coef_tte" = betas_mi_tte)
  
  return(ddest) 
  
}







###########################################################################################################################################################################################################################################
###########################################################################################################################################################################################################################################
###########################################################################################################################################################################################################################################
###########################################################################################################################################################################################################################################
###########################################################################################################################################################################################################################################
###########################################################################################################################################################################################################################################
###########################################################################################################################################################################################################################################
###########################################################################################################################################################################################################################################
###########################################################################################################################################################################################################################################
###########################################################################################################################################################################################################################################
###########################################################################################################################################################################################################################################
###########################################################################################################################################################################################################################################
###########################################################################################################################################################################################################################################
###########################################################################################################################################################################################################################################





