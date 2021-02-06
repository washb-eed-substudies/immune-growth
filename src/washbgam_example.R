
rm(list=ls())

source(here::here("0-config.R"))
library(washbgam)

d <- readRDS(paste0(dropboxDir,"Data/Cleaned/Audrie/bangladesh-immune-growth-analysis-dataset.rds"))

#Example:

#Fit GAM model with random effects for childid
res_unadj <- fit_RE_gam(d=d, X="t3_ratio_pro_il10", Y="laz_t3",  W=NULL)

#Get predictions of differences from the 25th percentile of exposure
preds_unadj <- predict_gam_diff(fit=res_unadj$fit, d=res_unadj$dat, quantile_diff=c(0.25,0.75), Xvar="t3_cort_z01", Yvar="laz_t3")


#Primary parameter we are estimating: difference between 25th and 75th percentile of the exposure
preds_unadj$res

#Fit spline with simultaneous confidence intervals
simul_plot <- gam_simul_CI(res_unadj$fit, res_unadj$dat, xlab="t3_ratio_pro_il10", ylab="laz_t3", title="example title")
simul_plot$p


#Adjusted
Wvars<-c("sex","birthord", "momage","momheight","momedu")

Wvars[!(Wvars %in% colnames(d))]
d$tr <- "control"
res_adj <- fit_RE_gam(d=d, X="t3_ratio_pro_il10", Y="laz_t3",  W=Wvars)




#### Loop over exposure-outcome pairs ####

#### Hypothesis 1: immune status associated with concurrent child growth ####
# all immune ratios at Y1 v. growth outcomes at Y1
Xvars <- c("t2_ratio_pro_il10", "t2_ratio_il2_il10", "t2_ratio_gmc_il10", "t2_ratio_th1_il10", "t2_ratio_th2_il10",     
           "t2_ratio_th17_il10", "t2_ratio_th1_th2", "t2_ratio_th1_th17", "t2_ln_agp", "t2_ln_crp", "sumscore_t2_Z")            
Yvars <- c("laz_t2", "waz_t2", "whz_t2" ,"hcz_t2") 

#Fit models
H1_models <- NULL
for(i in Xvars){
  for(j in Yvars){
    res_unadj <- fit_RE_gam(d=d, X=i, Y=j,  W=NULL)
    res <- data.frame(X=i, Y=j,fit=I(list(res_unadj$fit)), dat=I(list(res_unadj$dat)))
    H1_models <- bind_rows(H1_models, res)
  }
}

# all immune outcomes at y2 and growth outcomes at y2
Xvars <- c("t3_ratio_pro_il10", "t3_ratio_il2_il10", "t3_ratio_gmc_il10", "t3_ratio_th1_il10", "t3_ratio_th2_il10",     
           "t3_ratio_th17_il10", "t3_ratio_th1_th2", "t3_ratio_th1_th17", "sumscore_t3_Z")            
Yvars <- c("laz_t3", "waz_t3", "whz_t3" ,"hcz_t3") 

#Fit models
for(i in Xvars){
  for(j in Yvars){
    res_unadj <- fit_RE_gam(d=total_d, X=i, Y=j,  W=NULL)
    res <- data.frame(X=i, Y=j, N=res_unadj$n, fit=I(list(res_unadj$fit)), dat=I(list(res_unadj$dat)))
    H1_models <- bind_rows(H1_models, res)
  }
}

#Get primary contrasts
H1_res <- NULL
for(i in 1:nrow(H1_models)){
  res <- data.frame(X=H1_models$X[i], Y=H1_models$Y[i])
  preds <- predict_gam_diff(fit=H1_models$fit[i][[1]], d=H1_models$dat[i][[1]], H1_models$N[i], quantile_diff=c(0.25,0.75), Xvar=res$X, Yvar=res$Y)
  H1_res <-  bind_rows(H1_res , preds$res)
}

#Make list of plots
H1_plot_list <- NULL
H1_plot_data <- NULL
for(i in 1:nrow(H1_models)){
  res <- data.frame(X=H1_models$X[i], Y=H1_models$Y[i])
  simul_plot <- gam_simul_CI(H1_models$fit[i][[1]], H1_models$dat[i][[1]], xlab=res$X, ylab=res$Y, title="")
  H1_plot_list[[i]] <-  simul_plot$p
  H1_plot_data <-  rbind(H1_plot_data, data.frame(Xvar=res$X, Yvar=res$Y, adj=0, simul_plot$pred))
}
