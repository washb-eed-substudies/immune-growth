rm(list=ls())

source(here::here("0-config.R"))
#source(here::here("src/0-gam-functions.R"))

d <- readRDS(paste0(dropboxDir,"Data/Cleaned/Audrie/bangladesh-immune-growth-analysis-dataset.rds"))

#Set list of adjustment variables
#Make vectors of adjustment variable names
Wvars<-c("sex","birthord", "momage","momheight","momedu", 
         "hfiacat", "Nlt18","Ncomp", "watmin", "walls", 
         "floor", "roof", "HHwealth", "tr", "cesd_sum_t2", 
         "ari7d_t2", "diar7d_t2", "nose7d_t2", "life_viol_any_t3")

Wvars[!(Wvars %in% colnames(d))]



#Add in time varying covariates:

#NOTES
#Does monsoon_ut2 need to be replaced with monsoon_ht2 for growth measures? (and agemth_ut2 with agedays_ht2?)
Wvars2<-c("ageday_bt2", "ageday_at2",  "month_bt2", "month_at2") 
Wvars3<-c("ageday_bt3", "ageday_at3", "month_bt3", "month_at3", 
          "laz_t2", "waz_t2", "cesd_sum_ee_t3", "pss_sum_mom_t3", 
          "ari7d_t3", "diar7d_t3", "nose7d_t3") 
Wvars23<-c("ageday_bt2", "ageday_at3", "month_bt2", "month_at3", 
           "laz_t2", "waz_t2", "cesd_sum_ee_t3", "pss_sum_mom_t3", 
           "ari7d_t3", "diar7d_t3", "nose7d_t3")
Wvars_anthro23<-c("ageday_bt2", "ageday_at2", "ageday_at3", "month_bt2", "month_at2", "month_at3", 
                  "cesd_sum_ee_t3", "pss_sum_mom_t3", "ari7d_t3", "diar7d_t3", "nose7d_t3")

W2_immmune.W2_anthro <- c(Wvars, Wvars2) %>% unique(.)
W3_immune.W3_anthro <- c(Wvars, Wvars3) %>% unique(.)
W2_immune.W3_anthro <- c(Wvars, Wvars23) %>% unique(.)
W2_immune.W23_anthro <- c(Wvars, Wvars_anthro23) %>% unique(.)

add_hcz <- function(j, W){
  if (j=="hcz_t3"){Wset=c(W, "hcz_t2")}
  else {Wset=W}
  return(Wset)
}

#Loop over exposure-outcome pairs

#### Hypothesis 1: immune status associated with concurrent child growth ####
# all immune ratios at Y1 v. growth outcomes at Y1
Xvars <- c("t2_ratio_pro_il10", "t2_ratio_il2_il10", "t2_ratio_gmc_il10", "t2_ratio_th1_il10", "t2_ratio_th2_il10",     
           "t2_ratio_th17_il10", "t2_ratio_th1_th2", "t2_ratio_th1_th17", "t2_ln_agp", "t2_ln_crp", "sumscore_t2_Z")
Yvars <- c("laz_t2", "waz_t2", "whz_t2" ,"hcz_t2") 

#Fit models
H1_adj_nofever_models <- NULL
for(i in Xvars){
  for(j in Yvars){
    print(i)
    print(j)
    Wset <- W2_immmune.W2_anthro
    res_adj <- fit_RE_gam(d=d, X=i, Y=j,  W=Wset, forcedW=c("ageday_bt2", "ageday_at2"))
    res <- data.frame(X=i, Y=j, fit=I(list(res_adj$fit)), dat=I(list(res_adj$dat)))
    H1_adj_nofever_models <- bind_rows(H1_adj_nofever_models, res)
  }
}

# all immune outcomes at y2 and growth outcomes at y2
Xvars <- c("t3_ratio_pro_il10", "t3_ratio_il2_il10", "t3_ratio_gmc_il10", "t3_ratio_th1_il10", "t3_ratio_th2_il10",     
           "t3_ratio_th17_il10", "t3_ratio_th1_th2", "t3_ratio_th1_th17", "sumscore_t3_Z")            
Yvars <- c("laz_t3", "waz_t3", "whz_t3" ,"hcz_t3") 

for(i in Xvars){
  for(j in Yvars){
    print(i)
    print(j)
    Wset <- add_hcz(j, W3_immune.W3_anthro)
    res_adj <- fit_RE_gam(d=d, X=i, Y=j,  W=Wset, forcedW=c("ageday_bt3", "ageday_at3"))
    res <- data.frame(X=i, Y=j, fit=I(list(res_adj$fit)), dat=I(list(res_adj$dat)))
    H1_adj_nofever_models <- bind_rows(H1_adj_nofever_models, res)
  }
}

#Get primary contrasts
H1_adj_nofever_res <- NULL
for(i in 1:nrow(H1_adj_nofever_models)){
  res <- data.frame(X=H1_adj_nofever_models$X[i], Y=H1_adj_nofever_models$Y[i])
  preds <- predict_gam_diff(fit=H1_adj_nofever_models$fit[i][[1]], d=H1_adj_nofever_models$dat[i][[1]], quantile_diff=c(0.25,0.75), Xvar=res$X, Yvar=res$Y)
  H1_adj_nofever_res <-  bind_rows(H1_adj_nofever_res , preds$res)
}

#Make list of plots
H1_adj_nofever_plot_list <- NULL
H1_adj_nofever_plot_data <- NULL
for(i in 1:nrow(H1_adj_nofever_models)){
  res <- data.frame(X=H1_adj_nofever_models$X[i], Y=H1_adj_nofever_models$Y[i])
  simul_plot <- gam_simul_CI(H1_adj_nofever_models$fit[i][[1]], H1_adj_nofever_models$dat[i][[1]], xlab=res$X, ylab=res$Y, title="")
  H1_adj_nofever_plot_list[[i]] <-  simul_plot$p
  H1_adj_nofever_plot_data <-  rbind(H1_adj_nofever_plot_data, data.frame(Xvar=res$X, Yvar=res$Y, adj=0, simul_plot$pred %>% subset(., select = c(Y,X,id,fit,se.fit,uprP, lwrP,uprS,lwrS))))
}


#Save models
#saveRDS(H1_adj_nofever_models, paste0(dropboxDir,"results/stress-growth-models/models/H1_adj_nofever_models.RDS"))

#Save results
saveRDS(H1_adj_nofever_res, here("results/adjusted/H1_adj_nofever_res.RDS"))

#Save plots
#saveRDS(H1_adj_nofever_plot_list, paste0(dropboxDir,"results/stress-growth-models/figure-objects/H1_adj_nofever_splines.RDS"))

#Save plot data
saveRDS(H1_adj_nofever_plot_data, here('figure-data/H1_adj_nofever_spline_data.RDS'))



#### Hypothesis 2: immune status and subsequent growth ####
# all immune outcomes at y1 v. growth at y2
Xvars <- c("t2_ratio_pro_il10", "t2_ratio_il2_il10", "t2_ratio_gmc_il10", "t2_ratio_th1_il10", "t2_ratio_th2_il10",     
           "t2_ratio_th17_il10", "t2_ratio_th1_th2", "t2_ratio_th1_th17", "t2_ln_agp", "t2_ln_crp", "sumscore_t2_Z")            
Yvars <- c("laz_t3", "waz_t3", "whz_t3" ,"hcz_t3") 


#Fit models
H2_adj_nofever_models <- NULL
for(i in Xvars){
  for(j in Yvars){
    print(i)
    print(j)
    Wset <- add_hcz(j, W2_immune.W3_anthro)
    res_adj <- fit_RE_gam(d=d, X=i, Y=j,  W=Wset, forcedW=c("ageday_bt2", "ageday_at3"))
    res <- data.frame(X=i, Y=j, fit=I(list(res_adj$fit)), dat=I(list(res_adj$dat)))
    H2_adj_nofever_models <- bind_rows(H2_adj_nofever_models, res)
  }
}

#Get primary contrasts
H2_adj_nofever_res <- NULL
for(i in 1:nrow(H2_adj_nofever_models)){
  res <- data.frame(X=H2_adj_nofever_models$X[i], Y=H2_adj_nofever_models$Y[i])
  preds <- predict_gam_diff(fit=H2_adj_nofever_models$fit[i][[1]], d=H2_adj_nofever_models$dat[i][[1]], quantile_diff=c(0.25,0.75), Xvar=res$X, Yvar=res$Y)
  H2_adj_nofever_res <-  bind_rows(H2_adj_nofever_res , preds$res)
}

#Make list of plots
H2_plot_list <- NULL
H2_plot_data <- NULL
for(i in 1:nrow(H2_adj_nofever_models)){
  res <- data.frame(X=H2_adj_nofever_models$X[i], Y=H2_adj_nofever_models$Y[i])
  simul_plot <- gam_simul_CI(H2_adj_nofever_models$fit[i][[1]], H2_adj_nofever_models$dat[i][[1]], xlab=res$X, ylab=res$Y, title="")
  H2_plot_list[[i]] <-  simul_plot$p
  H2_plot_data <-  rbind(H2_plot_data, data.frame(Xvar=res$X, Yvar=res$Y, adj=0, simul_plot$pred%>% subset(., select = c(Y,X,id,fit,se.fit,uprP, lwrP,uprS,lwrS))))
}


#Save models
#saveRDS(H2_adj_nofever_models, paste0(dropboxDir,"results/stress-growth-models/models/adj_nofever_H2_adj_nofever_models.RDS"))

#Save results
saveRDS(H2_adj_nofever_res, here("results/adjusted/H2_adj_nofever_res.RDS"))


#Save plots
#saveRDS(H2_plot_list, paste0(dropboxDir,"results/stress-growth-models/figure-objects/H2_adj_nofever_splines.RDS"))

#Save plot data
saveRDS(H2_plot_data, here('figure-data/H2_adj_nofever_spline_data.RDS'))



#### Hypothesis 3: immune status and child growth velocity ####
# immune ratios at y1 and growth velocity outcomes between y1 and y2
Xvars <- c("t2_ratio_pro_il10", "t2_ratio_il2_il10", "t2_ratio_gmc_il10", "t2_ratio_th1_il10", "t2_ratio_th2_il10",     
           "t2_ratio_th17_il10", "t2_ratio_th1_th2", "t2_ratio_th1_th17", "t2_ln_agp", "t2_ln_crp", "sumscore_t2_Z")            
Yvars <- c("len_velocity_t2_t3", "wei_velocity_t2_t3", "hc_velocity_t2_t3")

#Fit models
H3_adj_nofever_models <- NULL
for(i in Xvars){
  for(j in Yvars){
    print(i)
    print(j)
    res_adj <- fit_RE_gam(d=d, X=i, Y=j,  W=W2_immune.W23_anthro, forcedW=c("ageday_bt2", "ageday_at2", "ageday_at3"))
    res <- data.frame(X=i, Y=j, fit=I(list(res_adj$fit)), dat=I(list(res_adj$dat)))
    H3_adj_nofever_models <- bind_rows(H3_adj_nofever_models, res)
  }
}

#Get primary contrasts
H3_adj_nofever_res <- NULL
for(i in 1:nrow(H3_adj_nofever_models)){
  res <- data.frame(X=H3_adj_nofever_models$X[i], Y=H3_adj_nofever_models$Y[i])
  preds <- predict_gam_diff(fit=H3_adj_nofever_models$fit[i][[1]], d=H3_adj_nofever_models$dat[i][[1]], quantile_diff=c(0.25,0.75), Xvar=res$X, Yvar=res$Y)
  H3_adj_nofever_res <-  bind_rows(H3_adj_nofever_res , preds$res)
}

#Make list of plots
H3_plot_list <- NULL
H3_plot_data <- NULL
for(i in 1:nrow(H3_adj_nofever_models)){
  res <- data.frame(X=H3_adj_nofever_models$X[i], Y=H3_adj_nofever_models$Y[i])
  simul_plot <- gam_simul_CI(H3_adj_nofever_models$fit[i][[1]], H3_adj_nofever_models$dat[i][[1]], xlab=res$X, ylab=res$Y, title="")
  H3_plot_list[[i]] <-  simul_plot$p
  H3_plot_data <-  rbind(H3_plot_data, data.frame(Xvar=res$X, Yvar=res$Y, adj=0, simul_plot$pred%>% subset(., select = c(Y,X,id,fit,se.fit,uprP, lwrP,uprS,lwrS))))
}


#Save models
#saveRDS(H3_adj_nofever_models, paste0(dropboxDir,"results/stress-growth-models/models/adj_nofever_H3_adj_nofever_models.RDS"))

#Save results
saveRDS(H3_adj_nofever_res, here("results/adjusted/H3_adj_nofever_res.RDS"))


#Save plots
#saveRDS(H3_plot_list, paste0(dropboxDir,"results/stress-growth-models/figure-objects/H3_adj_nofever_splines.RDS"))

#Save plot data
saveRDS(H3_plot_data, here('figure-data/H3_adj_nofever_spline_data.RDS'))


#### Hypothesis ####
# immune ratios at y1 v. change in growth outcomes between y1 and y2
Xvars <- c("t2_ratio_pro_il10", "t2_ratio_il2_il10", "t2_ratio_gmc_il10", "t2_ratio_th1_il10", "t2_ratio_th2_il10",     
           "t2_ratio_th17_il10", "t2_ratio_th1_th2", "t2_ratio_th1_th17", "t2_ln_agp", "t2_ln_crp", "sumscore_t2_Z")            
Yvars <- c("delta_laz_t2_t3", "delta_waz_t2_t3", "delta_whz_t2_t3", "delta_hcz_t2_t3")

#Fit models
delta_growth_adj_nofever_models <- NULL
for(i in Xvars){
  for(j in Yvars){
    print(i)
    print(j)
    res_adj <- fit_RE_gam(d=d, X=i, Y=j,  W=W2_immune.W23_anthro, forcedW=c("ageday_bt2", "ageday_at2", "ageday_at3"))
    res <- data.frame(X=i, Y=j, fit=I(list(res_adj$fit)), dat=I(list(res_adj$dat)))
    delta_growth_adj_nofever_models <- bind_rows(delta_growth_adj_nofever_models, res)
  }
}

#Get primary contrasts
delta_growth_adj_nofever_res <- NULL
for(i in 1:nrow(delta_growth_adj_nofever_models)){
  res <- data.frame(X=delta_growth_adj_nofever_models$X[i], Y=delta_growth_adj_nofever_models$Y[i])
  preds <- predict_gam_diff(fit=delta_growth_adj_nofever_models$fit[i][[1]], d=delta_growth_adj_nofever_models$dat[i][[1]], quantile_diff=c(0.25,0.75), Xvar=res$X, Yvar=res$Y)
  delta_growth_adj_nofever_res <-  bind_rows(delta_growth_adj_nofever_res , preds$res)
}

#Make list of plots
delta_growth_adj_nofever_plot_list <- NULL
delta_growth_adj_nofever_plot_data <- NULL
for(i in 1:nrow(delta_growth_adj_nofever_models)){
  res <- data.frame(X=delta_growth_adj_nofever_models$X[i], Y=delta_growth_adj_nofever_models$Y[i])
  simul_plot <- gam_simul_CI(delta_growth_adj_nofever_models$fit[i][[1]], delta_growth_adj_nofever_models$dat[i][[1]], xlab=res$X, ylab=res$Y, title="")
  delta_growth_adj_nofever_plot_list[[i]] <-  simul_plot$p
  delta_growth_adj_nofever_plot_data <-  rbind(delta_growth_adj_nofever_plot_data, data.frame(Xvar=res$X, Yvar=res$Y, adj=0, simul_plot$pred%>% subset(., select = c(Y,X,id,fit,se.fit,uprP, lwrP,uprS,lwrS))))
}


#Save models
#saveRDS(delta_growth_adj_nofever_models, paste0(dropboxDir,"results/stress-growth-models/models/adj_nofever_delta_growth_adj_nofever_models.RDS"))

#Save results
saveRDS(delta_growth_adj_nofever_res, here("results/adjusted/delta_growth_adj_nofever_res.RDS"))


#Save plots
#saveRDS(delta_growth_adj_nofever_plot_list, paste0(dropboxDir,"results/stress-growth-models/figure-objects/delta_growth_adj_nofever_adj_nofever_splines.RDS"))

#Save plot data
saveRDS(delta_growth_adj_nofever_plot_data, here("figure-data/delta_growth_adj_nofever_spline_data.RDS"))

