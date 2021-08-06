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



#Add in time varying covariates:
Wvars2<-c("ageday_bt2", "ageday_at2",  "month_bt2", "month_at2", "laz_t1_cat", "waz_t1_cat") 
Wvars3<-c("ageday_bt3", "ageday_at3", "month_bt3", "month_at3", 
          "laz_t2_cat", "waz_t2_cat", "cesd_sum_ee_t3", "pss_sum_mom_t3", 
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

add_hcz <- function(i, j, W){
  if (j=="hcz_t3"){
    if(grepl("t2", i)){Wset=c(W, "hcz_t2")}
    else {Wset=c(W, "hcz_t2_cat")}}
  else if (j=="hcz_t2"){Wset=c(W, "hcz_t1_cat")}
  else {Wset=W}
  return(Wset)
}

#Loop over exposure-outcome pairs

#### Hypothesis 1: immune status associated with concurrent child growth ####
# all immune ratios at Y1 v. growth outcomes at Y1
Xvars <- c("t3_ln_crp", "t3_ln_agp")
Yvars <- c("laz_t3", "waz_t3", "whz_t3" ,"hcz_t3") 

#Fit models
H1_adj_nofever_models <- NULL
for(i in Xvars){
  for(j in Yvars){
    print(i)
    print(j)
    Wset <- add_hcz(i, j, W3_immune.W3_anthro)
    res_adj <- fit_RE_gam(d=d, X=i, Y=j,  W=Wset)
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



#Save results
saveRDS(H1_adj_nofever_res, here("results/adjusted/crp_agp_t3subset_res.RDS"))

#Save plot data
saveRDS(H1_adj_nofever_plot_data, here('figure-data/crp_agp_t3subset_spline_data.RDS'))

H1_adj_nofever_res$CI <- paste(H1_adj_nofever_res$point.diff %>% round(4), "(", H1_adj_nofever_res$lb.diff %>% round(4), ", ", H1_adj_nofever_res$ub.diff %>% round(4), ")", sep = "")
H1_adj_nofever_res %>% select(X, Y, N, q1, q3, pred.q1, pred.q3, CI, Pval) %>% 
  write.csv("C:/Users/Sophia/Documents/WASH/WASH Immune and Growth/crp_agp_t3_subset_results.csv", row.names = F)


## EMM analyses CRP and AGP ##
#Analysis
gam.analysis <- function(save, Xvar = NULL, Yvar = NULL, data = d, Wvar = NULL, Vvar = NULL){
  for(i in Xvar){
    for(j in Yvar){
      for(k in Vvar){
        print(i)
        print(j)
        print(k)
        res_adj <- fit_RE_gam(d=data, X=i, Y=j,  W=Wvar, V=k)
        res <- data.frame(X=i, Y=j,V=k, int.p =res_adj$int.p, fit=I(list(res_adj$fit)), dat=I(list(res_adj$dat)))
        save <- bind_rows(save, res)
      }
    }
  }
  return(save)
}

crp_agp_EMM_models <- NULL

crp_agp_t2 <- c("t2_ln_crp", "t2_ln_agp")
crp_agp_t3 <- c("t3_ln_crp", "t3_ln_agp")
growth_t2 <- c("laz_t2", "waz_t2", "whz_t2", "hcz_t2")
growth_t3 <- c("laz_t3", "waz_t3", "whz_t3", "hcz_t3")
growth_t2_t3 <- c("delta_laz_t2_t3", "delta_waz_t2_t3", "delta_whz_t2_t3", "delta_hcz_t2_t3",
                  "len_velocity_t2_t3", "wei_velocity_t2_t3", "hc_velocity_t2_t3")

W2_immmune.W2_anthro <- c(Wvars, Wvars2) %>% unique(.)
W3_immune.W3_anthro <- c(Wvars, Wvars3) %>% unique(.)
W2_immune.W3_anthro <- c(Wvars, Wvars23) %>% unique(.)
W2_immune.W23_anthro <- c(Wvars, Wvars_anthro23) %>% unique(.)

W2_immmune.W2_anthro <- W2_immmune.W2_anthro[W2_immmune.W2_anthro != "tr"]
W3_immune.W3_anthro <- W3_immune.W3_anthro[W3_immune.W3_anthro != "tr"]
W2_immune.W3_anthro <- W2_immune.W3_anthro[W2_immune.W3_anthro != "tr"]
W2_immune.W23_anthro <- W2_immune.W23_anthro[W2_immune.W23_anthro != "tr"]

#analysis 1
crp_agp_EMM_models <- gam.analysis(crp_agp_EMM_models, crp_agp_t2, growth_t2, d, W2_immmune.W2_anthro, c("tr"))
crp_agp_EMM_models <- gam.analysis(crp_agp_EMM_models, crp_agp_t3, growth_t3, dfull, W3_immune.W3_anthro, c("tr"))
crp_agp_EMM_models <- gam.analysis(crp_agp_EMM_models, crp_agp_t2, growth_t3, d, W2_immune.W3_anthro, c("tr"))
crp_agp_EMM_models <- gam.analysis(crp_agp_EMM_models, crp_agp_t2, growth_t2_t3, d, W2_immune.W23_anthro, c("tr"))


gam.results <- function(models, save){
  for(i in 1:nrow(models)){
    preds <- predict_gam_emm(fit=models$fit[i][[1]], d=models$dat[i][[1]], quantile_diff=c(0.25,0.75), Xvar=models$X[i], Yvar=models$Y[i])
    gamm_diff_res <- data.frame(V=models$V[i] , preds$res) 
    gamm_diff_res <- gamm_diff_res %>% 
      mutate(int.Pval = c(rep(NA, nrow(gamm_diff_res)-1), models$int.p[[i]]))
    save <-  bind_rows(save, gamm_diff_res)
  }
  save
}

res <- NULL

res <- gam.results(crp_agp_EMM_models, res) 

saveRDS(res, here("results/adjusted/crp_agp_EMM_res.RDS"))

res$CI <- paste(res$point.diff %>% round(4), "(", res$lb.diff %>% round(4), ", ", res$ub.diff %>% round(4), ")", sep = "")
res %>% select(X, Y, Vlevel, N, q1, q3, pred.q1, pred.q3, point.diff, lb.diff, ub.diff, Pval, int.Pval) %>% 
  write.csv("C:/Users/Sophia/Documents/WASH/WASH Immune and Growth/crp_agp_EMM_results.csv", row.names = F)

