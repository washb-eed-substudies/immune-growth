###### RUN ADJUSTED CONCURRENT ANALYSES WITHOUT ADJUSTING FOR PRIOR GROWTH ######
rm(list=ls())

source(here::here("0-config.R"))
#source(here::here("src/0-gam-functions.R"))

dfull <- readRDS(paste0(dropboxDir,"Data/Cleaned/Audrie/bangladesh-immune-growth-analysis-dataset.rds"))
d <- dfull %>% filter(tr %in% c("Control", "Nutrition + WSH"))


#Set list of adjustment variables
#Make vectors of adjustment variable names
#removed roof because of low variability
Wvars<-c("sex","birthord", "momage","momheight","momedu", 
         "hfiacat", "Nlt18","Ncomp", "watmin", "walls", 
         "floor", "HHwealth", "tr", "cesd_sum_t2", 
         "ari7d_t2", "diar7d_t2", "nose7d_t2", "life_viol_any_t3")

Wvars[!(Wvars %in% colnames(d))]

Wvars2<-c("ageday_bt2", "ageday_at2",  "month_bt2", "month_at2") 
Wvars3<-c("ageday_bt3", "ageday_at3", "month_bt3", "month_at3", 
          "cesd_sum_ee_t3", "pss_sum_mom_t3", 
          "ari7d_t3", "diar7d_t3", "nose7d_t3") 

W2_immmune.W2_anthro <- c(Wvars, Wvars2) %>% unique(.)
W3_immune.W3_anthro <- c(Wvars, Wvars3) %>% unique(.)

#Loop over exposure-outcome pairs

#### Hypothesis 1: immune status associated with concurrent child growth ####
# all immune ratios at Y1 v. growth outcomes at Y1
Xvars <- c("t2_ratio_pro_il10", "t2_ratio_il2_il10", "t2_ratio_gmc_il10", "t2_ratio_th1_il10", "t2_ratio_th2_il10",     
           "t2_ratio_th17_il10", "t2_ratio_th1_th2", "t2_ratio_th1_th17", "t2_ln_agp", "t2_ln_crp", "sumscore_t2_Z", "t2_ln_ifn")
Yvars <- c("laz_t2", "waz_t2", "whz_t2" ,"hcz_t2") 

#Fit models
H1_adj_nofever_models <- NULL
for(i in Xvars){
  for(j in Yvars){
    print(i)
    print(j)
    Wset <- W2_immmune.W2_anthro
    res_adj <- fit_RE_gam(d=d, X=i, Y=j,  W=Wset)
    res <- data.frame(X=i, Y=j, fit=I(list(res_adj$fit)), dat=I(list(res_adj$dat)))
    H1_adj_nofever_models <- bind_rows(H1_adj_nofever_models, res)
  }
}

# all immune outcomes at y2 and growth outcomes at y2
Xvars <- c("t3_ratio_pro_il10", "t3_ratio_il2_il10", "t3_ratio_gmc_il10", "t3_ratio_th1_il10", "t3_ratio_th2_il10",     
           "t3_ratio_th17_il10", "t3_ratio_th1_th2", "t3_ratio_th1_th17", "t3_ln_crp", "t3_ln_agp", "sumscore_t3_Z", "t3_ln_ifn")            
Yvars <- c("laz_t3", "waz_t3", "whz_t3" ,"hcz_t3") 

for(i in Xvars){
  for(j in Yvars){
    print(i)
    print(j)
    Wset <- W3_immune.W3_anthro
    if (i %in% c("t3_ln_crp", "t3_ln_agp")){dfunc <- dfull}
    else {dfunc <- d}
    res_adj <- fit_RE_gam(d=dfunc, X=i, Y=j,  W=Wset)
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
saveRDS(H1_adj_nofever_res, here("results/adjusted/concurrent-growth-nopriorgrowth-res.RDS"))

#Save plots
#saveRDS(H1_adj_nofever_plot_list, paste0(dropboxDir,"results/stress-growth-models/figure-objects/H1_adj_nofever_splines.RDS"))

#Save plot data
saveRDS(H1_adj_nofever_plot_data, here('figure-data/concurrent-growth-nopriorgrowth-spline-data.RDS'))

# load results with adjustment for prior growth
with_growth <- readRDS("C:/Users/Sophia/Documents/immune-growth/results/adjusted/H1_adj_nofever_res.RDS")
view(with_growth)
# view results without adjusting for prior growth
view(H1_adj_nofever_res)
combined <- merge(with_growth, H1_adj_nofever_res, by = c("Y", "X"))
combined$CI.x <- paste(combined$point.diff.x %>% round(4), "(", combined$lb.diff.x %>% round(4), ", ", combined$ub.diff.x %>% round(4), ")", sep = "")
combined$CI.y <- paste(combined$point.diff.y %>% round(4), "(", combined$lb.diff.y %>% round(4), ", ", combined$ub.diff.y %>% round(4), ")", sep = "")
selected <- combined %>% select(X, Y, q1.x, q1.y, q3.x, q3.y, CI.x, CI.y, Pval.x, Pval.y)
selected
names(selected) <- gsub("\\.y", " no growth", gsub("\\.x", " w/ growth", names(selected)))
write.csv(selected, "C:/Users/Sophia/Documents/WASH/WASH Immune and Growth/growth-adjustment-effects.csv", row.names = F)
