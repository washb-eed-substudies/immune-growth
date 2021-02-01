rm(list=ls())

source(here::here("0-config.R"))
#source(here::here("src/0-gam-functions.R"))

d <- readRDS(paste0(dropboxDir,"Data/Cleaned/Audrie/bangladesh-immune-growth-analysis-dataset.rds"))

#Set list of adjustment variables
#Make vectors of adjustment variable names
Wvars<-c("sex","birthord", "momage","momheight","momedu", 
         "hfiacat", "Nlt18","Ncomp", "watmin", "walls", "floor", "roof", "elec", "asset_wardrobe",
         "asset_table", "asset_chair", "asset_clock","asset_khat", "asset_chouki", 
         "asset_radio", "asset_tv", "asset_refrig", "asset_bike", "asset_moto", "asset_sewmach", 
         "asset_mobile", "n_cattle", "n_goat", "n_chicken")

Wvars[!(Wvars %in% colnames(d))]

# # #test new package
# X="t3_ratio_pro_il10"
# Y="laz_t3"
# d$momheight2 <- d$momheight+0.01
# d$hfiacat2 <- d$hfiacat
# W=c("birthord", "momage","momheight","momedu", "ageday_at1", "ageday_at2", "momheight2","hfiacat2")
# forcedW=W[grepl("age_", W)|grepl("agedays_", W)|grepl("ageday_", W)]
# V="sex"
# id="clusterid"
# family = "gaussian"
# pval = 0.2
# print=TRUE
# res_adj <- fit_RE_gam(d=d, X="t3_ratio_pro_il10", Y="laz_t3",  W=c(Wvars, "ageday_at1", "ageday_at2", "momheight2","hfiacat2"), V="sex")


#Add in time varying covariates:

#NOTES
#Does monsoon_ut2 need to be replaced with monsoon_ht2 for growth measures? (and agemth_ut2 with agedays_ht2?)
Wvars2<-c("monsoon_at2", "ageday_at2", "tr", "cesd_sum_t2", "fever7d_t2", "ari7d_t2", "diar7d_t2", "nose7d_t2") 
Wvars3<-c("lenhei_med_t2", "weight_med_t2", "monsoon_at2", "monsoon_at3", "ageday_at2", "ageday_at3", "tr",
          "cesd_sum_t2", "cesd_sum_ee_t3", "pss_sum_mom_t3", "life_viol_any_t3", "fever7d_t2", "fever7d_t3", 
          "ari7d_t2", "ari7d_t3", "diar7d_t2", "diar7d_t3", "nose7d_t2", "nose7d_t3") 
Wvars23<-c("anthro_days_btwn_t2_t3")

W2_F2.W2_anthro <- c(Wvars, Wvars2) %>% unique(.)
W2_F2.W3_anthro <- c(Wvars, Wvars3) %>% unique(.)
W2_F2.W23_anthro <- c(Wvars, Wvars3, Wvars23) %>% unique(.)


#Loop over exposure-outcome pairs

#### Hypothesis 1: immune status associated with concurrent child growth ####
# all immune ratios at Y1 v. growth outcomes at Y1
Xvars <- c("t2_ratio_pro_il10", "t2_ratio_il2_il10", "t2_ratio_gmc_il10", "t2_ratio_th1_il10", "t2_ratio_th2_il10",     
           "t2_ratio_th17_il10", "t2_ratio_th1_th2", "t2_ratio_th1_th17", "t2_ln_agp", "t2_ln_crp")            
Yvars <- c("laz_t2", "waz_t2", "whz_t2" ,"hcz_t2") 

pick_covariates_H1 <- function(j){
  if(grepl("_t2", j)){Wset = W2_F2.W2_anthro}
  if(grepl("_t3", j)){Wset = W2_F2.W3_anthro}
  return(Wset)
}
#Fit models
H1_adj_models <- NULL
for(i in Xvars){
  for(j in Yvars){
    print(i)
    print(j)
    Wset <- pick_covariates_H1(j)
    res_adj <- fit_RE_gam(d=d, X=i, Y=j,  W=Wset)
    res <- data.frame(X=i, Y=j, fit=I(list(res_adj$fit)), dat=I(list(res_adj$dat)))
    H1_adj_models <- bind_rows(H1_adj_models, res)
  }
}

# all immune outcomes at y2 and growth outcomes at y2
Xvars <- c("t3_ratio_pro_il10", "t3_ratio_il2_il10", "t3_ratio_gmc_il10", "t3_ratio_th1_il10", "t3_ratio_th2_il10",     
           "t3_ratio_th17_il10", "t3_ratio_th1_th2", "t3_ratio_th1_th17")            
Yvars <- c("laz_t3", "waz_t3", "whz_t3" ,"hcz_t3") 

for(i in Xvars){
  for(j in Yvars){
    print(i)
    print(j)
    Wset <- pick_covariates_H1(j)
    res_adj <- fit_RE_gam(d=d, X=i, Y=j,  W=Wset)
    res <- data.frame(X=i, Y=j, fit=I(list(res_adj$fit)), dat=I(list(res_adj$dat)))
    H1_adj_models <- bind_rows(H1_adj_models, res)
  }
}

#Get primary contrasts
H1_adj_res <- NULL
for(i in 1:nrow(H1_adj_models)){
  res <- data.frame(X=H1_adj_models$X[i], Y=H1_adj_models$Y[i])
  preds <- predict_gam_diff(fit=H1_adj_models$fit[i][[1]], d=H1_adj_models$dat[i][[1]], quantile_diff=c(0.25,0.75), Xvar=res$X, Yvar=res$Y)
  H1_adj_res <-  bind_rows(H1_adj_res , preds$res)
}

#Make list of plots
H1_adj_plot_list <- NULL
H1_adj_plot_data <- NULL
for(i in 1:nrow(H1_adj_models)){
  res <- data.frame(X=H1_adj_models$X[i], Y=H1_adj_models$Y[i])
  simul_plot <- gam_simul_CI(H1_adj_models$fit[i][[1]], H1_adj_models$dat[i][[1]], xlab=res$X, ylab=res$Y, title="")
  H1_adj_plot_list[[i]] <-  simul_plot$p
  H1_adj_plot_data <-  rbind(H1_adj_plot_data, data.frame(Xvar=res$X, Yvar=res$Y, adj=0, simul_plot$pred %>% subset(., select = c(Y,X,id,fit,se.fit,uprP, lwrP,uprS,lwrS))))
}


#Save models
#saveRDS(H1_adj_models, paste0(dropboxDir,"results/stress-growth-models/models/H1_adj_models.RDS"))

#Save results
#saveRDS(H1_adj_res, here("results/adjusted/H1_adj_res.RDS"))


#Save plots
#saveRDS(H1_adj_plot_list, paste0(dropboxDir,"results/stress-growth-models/figure-objects/H1_adj_splines.RDS"))

#Save plot data
saveRDS(H1_adj_plot_data, here("figure-data/H1_adj_spline_data.RDS"))



#### Hypothesis 2: immune status and subsequent growth ####
# all immune outcomes at y1 v. growth at y2
Xvars <- c("t2_ratio_pro_il10", "t2_ratio_il2_il10", "t2_ratio_gmc_il10", "t2_ratio_th1_il10", "t2_ratio_th2_il10",     
           "t2_ratio_th17_il10", "t2_ratio_th1_th2", "t2_ratio_th1_th17", "t2_ln_agp", "t2_ln_crp")            
Yvars <- c("laz_t3", "waz_t3", "whz_t3" ,"hcz_t3") 


#Fit models
H2_adj_models <- NULL
for(i in Xvars){
  for(j in Yvars){
    print(i)
    print(j)
    res_adj <- fit_RE_gam(d=d, X=i, Y=j,  W=W2_F2.W3_anthro)
    res <- data.frame(X=i, Y=j, fit=I(list(res_adj$fit)), dat=I(list(res_adj$dat)))
    H2_adj_models <- bind_rows(H2_adj_models, res)
  }
}

#Get primary contrasts
H2_adj_res <- NULL
for(i in 1:nrow(H2_adj_models)){
  res <- data.frame(X=H2_adj_models$X[i], Y=H2_adj_models$Y[i])
  preds <- predict_gam_diff(fit=H2_adj_models$fit[i][[1]], d=H2_adj_models$dat[i][[1]], quantile_diff=c(0.25,0.75), Xvar=res$X, Yvar=res$Y)
  H2_adj_res <-  bind_rows(H2_adj_res , preds$res)
}

#Make list of plots
H2_plot_list <- NULL
H2_plot_data <- NULL
for(i in 1:nrow(H2_adj_models)){
  res <- data.frame(X=H2_adj_models$X[i], Y=H2_adj_models$Y[i])
  simul_plot <- gam_simul_CI(H2_adj_models$fit[i][[1]], H2_adj_models$dat[i][[1]], xlab=res$X, ylab=res$Y, title="")
  H2_plot_list[[i]] <-  simul_plot$p
  H2_plot_data <-  rbind(H2_plot_data, data.frame(Xvar=res$X, Yvar=res$Y, adj=0, simul_plot$pred%>% subset(., select = c(Y,X,id,fit,se.fit,uprP, lwrP,uprS,lwrS))))
}


#Save models
#saveRDS(H2_adj_models, paste0(dropboxDir,"results/stress-growth-models/models/adj_H2_adj_models.RDS"))

#Save results
#saveRDS(H2_adj_res, here("results/adjusted/H2_adj_res.RDS"))


#Save plots
#saveRDS(H2_plot_list, paste0(dropboxDir,"results/stress-growth-models/figure-objects/H2_adj_splines.RDS"))

#Save plot data
saveRDS(H2_plot_data, here("figure-data/H2_adj_spline_data.RDS"))



#### Hypothesis 3: immune status and child growth velocity ####
# immune ratios at y1 and growth velocity outcomes between y1 and y2
Xvars <- c("t2_ratio_pro_il10", "t2_ratio_il2_il10", "t2_ratio_gmc_il10", "t2_ratio_th1_il10", "t2_ratio_th2_il10",     
           "t2_ratio_th17_il10", "t2_ratio_th1_th2", "t2_ratio_th1_th17", "t2_ln_agp", "t2_ln_crp")            
Yvars <- c("len_velocity_t2_t3", "wei_velocity_t2_t3", "hc_velocity_t2_t3")

#Fit models
H3_adj_models <- NULL
for(i in Xvars){
  for(j in Yvars){
    res_adj <- fit_RE_gam(d=d, X=i, Y=j,  W=W2_F2.W23_anthro)
    res <- data.frame(X=i, Y=j, fit=I(list(res_adj$fit)), dat=I(list(res_adj$dat)))
    H3_adj_models <- bind_rows(H3_adj_models, res)
  }
}

#Get primary contrasts
H3_adj_res <- NULL
for(i in 1:nrow(H3_adj_models)){
  res <- data.frame(X=H3_adj_models$X[i], Y=H3_adj_models$Y[i])
  preds <- predict_gam_diff(fit=H3_adj_models$fit[i][[1]], d=H3_adj_models$dat[i][[1]], quantile_diff=c(0.25,0.75), Xvar=res$X, Yvar=res$Y)
  H3_adj_res <-  bind_rows(H3_adj_res , preds$res)
}

#Make list of plots
H3_plot_list <- NULL
H3_plot_data <- NULL
for(i in 1:nrow(H3_adj_models)){
  res <- data.frame(X=H3_adj_models$X[i], Y=H3_adj_models$Y[i])
  simul_plot <- gam_simul_CI(H3_adj_models$fit[i][[1]], H3_adj_models$dat[i][[1]], xlab=res$X, ylab=res$Y, title="")
  H3_plot_list[[i]] <-  simul_plot$p
  H3_plot_data <-  rbind(H3_plot_data, data.frame(Xvar=res$X, Yvar=res$Y, adj=0, simul_plot$pred%>% subset(., select = c(Y,X,id,fit,se.fit,uprP, lwrP,uprS,lwrS))))
}


#Save models
#saveRDS(H3_adj_models, paste0(dropboxDir,"results/stress-growth-models/models/adj_H3_adj_models.RDS"))

#Save results
#saveRDS(H3_adj_res, here("results/adjusted/H3_adj_res.RDS"))


#Save plots
#saveRDS(H3_plot_list, paste0(dropboxDir,"results/stress-growth-models/figure-objects/H3_adj_splines.RDS"))

#Save plot data
saveRDS(H3_plot_data, here("figure-data/H3_adj_spline_data.RDS"))


#### Hypothesis ####
# immune ratios at y1 v. change in growth outcomes between y1 and y2
Xvars <- c("t2_ratio_pro_il10", "t2_ratio_il2_il10", "t2_ratio_gmc_il10", "t2_ratio_th1_il10", "t2_ratio_th2_il10",     
           "t2_ratio_th17_il10", "t2_ratio_th1_th2", "t2_ratio_th1_th17", "t2_ln_agp", "t2_ln_crp")            
Yvars <- c("delta_laz_t2_t3", "delta_waz_t2_t3", "delta_whz_t2_t3", "delta_hcz_t2_t3")

#Fit models
delta_growth_adj_models <- NULL
for(i in Xvars){
  for(j in Yvars){
    print(i)
    print(j)
    res_adj <- fit_RE_gam(d=d, X=i, Y=j,  W=W2_F2.W23_anthro)
    res <- data.frame(X=i, Y=j, fit=I(list(res_adj$fit)), dat=I(list(res_adj$dat)))
    delta_growth_adj_models <- bind_rows(delta_growth_adj_models, res)
  }
}

#Get primary contrasts
delta_growth_adj_res <- NULL
for(i in 1:nrow(delta_growth_adj_models)){
  res <- data.frame(X=delta_growth_adj_models$X[i], Y=delta_growth_adj_models$Y[i])
  preds <- predict_gam_diff(fit=delta_growth_adj_models$fit[i][[1]], d=delta_growth_adj_models$dat[i][[1]], quantile_diff=c(0.25,0.75), Xvar=res$X, Yvar=res$Y)
  delta_growth_adj_res <-  bind_rows(delta_growth_adj_res , preds$res)
}

#Make list of plots
delta_growth_adj_plot_list <- NULL
delta_growth_adj_plot_data <- NULL
for(i in 1:nrow(delta_growth_adj_models)){
  res <- data.frame(X=delta_growth_adj_models$X[i], Y=delta_growth_adj_models$Y[i])
  simul_plot <- gam_simul_CI(delta_growth_adj_models$fit[i][[1]], delta_growth_adj_models$dat[i][[1]], xlab=res$X, ylab=res$Y, title="")
  delta_growth_adj_plot_list[[i]] <-  simul_plot$p
  delta_growth_adj_plot_data <-  rbind(delta_growth_adj_plot_data, data.frame(Xvar=res$X, Yvar=res$Y, adj=0, simul_plot$pred%>% subset(., select = c(Y,X,id,fit,se.fit,uprP, lwrP,uprS,lwrS))))
}


#Save models
#saveRDS(delta_growth_adj_models, paste0(dropboxDir,"results/stress-growth-models/models/adj_delta_growth_adj_models.RDS"))

#Save results
#saveRDS(delta_growth_adj_res, here("results/adjusted/delta_growth_adj_res.RDS"))


#Save plots
#saveRDS(delta_growth_adj_plot_list, paste0(dropboxDir,"results/stress-growth-models/figure-objects/delta_growth_adj_adj_splines.RDS"))

#Save plot data
saveRDS(delta_growth_adj_plot_data, here("figure-data/delta_growth_adj_spline_data.RDS"))


# Adjust Pvalues with Benjamini-Hochberg procedure
full_res <- rbind(H1_adj_res, H2_adj_res, H3_adj_res, delta_growth_adj_res)
full_res$corrected.Pval <- p.adjust(full_res[['Pval']], method="BH")

H1_corr_res<-full_res[1:nrow(H1_adj_res),]
H2_corr_res<-full_res[(nrow(H1_adj_res)+1):(nrow(H1_adj_res)+nrow(H2_adj_res)),]
H3_corr_res<-full_res[(nrow(H1_adj_res)+nrow(H2_adj_res)+1):(nrow(H1_adj_res)+nrow(H2_adj_res)+nrow(H3_adj_res)),]
delta_growth_corr_res<-full_res[(nrow(full_res)-nrow(delta_growth_adj_res)+1):nrow(full_res),]

#Save results
saveRDS(H1_corr_res, here("results/adjusted/H1_adj_res.RDS"))
saveRDS(H2_corr_res, here("results/adjusted/H2_adj_res.RDS"))
saveRDS(H3_corr_res, here("results/adjusted/H3_adj_res.RDS"))
saveRDS(delta_growth_corr_res, here("results/adjusted/delta_growth_adj_res.RDS"))
