

rm(list=ls())

source(here::here("0-config.R"))

d <- readRDS(paste0(dropboxDir,"Data/Cleaned/Audrie/bangladesh-immune-growth-analysis-dataset.rds"))


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

#Example:

#Fit GAM model with random effects for childid
res_unadj <- fit_RE_gam(d=d, X="t3_ln_agp", Y="hcz_t3",  W=NULL)
preds_unadj <- predict_gam_diff(fit=res_unadj$fit, d=res_unadj$dat, quantile_diff=c(0.25,0.75), Xvar="t3_ln_agp", Yvar="hcz_t3")
preds_unadj$res
simul_plot_unadj <- gam_simul_CI(res_unadj$fit, res_unadj$dat, xlab="delta_TS", ylab="laz_t3", title="unadj")
simul_plot_unadj$p



res_unadj_emm <- fit_RE_gam(d=d, X="t3_ln_agp", Y="hcz_t3",  W=NULL, V="tr")
res_unadj_emm$int.p
preds_unadj_emm <- predict_gam_emm(fit=res_unadj_emm$fit, d=res_unadj_emm$dat, quantile_diff=c(0.25,0.75), Xvar="t3_ln_agp", Yvar="hcz_t3")
preds_unadj_emm$res



res_adj <- fit_RE_gam(d=d, X="t3_ln_agp", Y="hcz_t3",  W=W3_immune.W3_anthro)
preds_adj <- predict_gam_diff(fit=res_adj$fit, d=res_adj$dat, quantile_diff=c(0.25,0.75), Xvar="t3_ln_agp", Yvar="hcz_t3")
preds_adj$res
simul_plot <- gam_simul_CI(res_adj$fit, res_adj$dat, xlab="delta_TS", ylab="laz_t3", title="adj")
simul_plot$p


res_adj_emm <- fit_RE_gam(d=d, X="t3_ln_agp", Y="hcz_t3",  W=W2_immune.W3_anthro[!(W2_immune.W3_anthro=="tr")], V="tr")
res_adj_emm$int.p
preds_adj_emm <- predict_gam_emm(fit=res_adj_emm$fit, d=res_adj_emm$dat, quantile_diff=c(0.25,0.75), Xvar="t3_ln_agp", Yvar="hcz_t3")
preds_adj_emm$res

preds <- predict_gam_emm(fit=models$fit[i][[1]], d=models$dat[i][[1]], quantile_diff=c(0.25,0.75), Xvar=models$X[i], Yvar=models$Y[i])


ggplot(d) +
  geom_point(aes(x=t3_ln_agp, y=hcz_t3), alpha=0.1) +
  geom_smooth(aes(x=t3_ln_agp, y=hcz_t3, group=tr, color=tr)) +
  geom_smooth(aes(x=t3_ln_agp, y=hcz_t3)) +
    theme(legend.position = "right")
  

  
simul_plot <- gam_simul_CI(res_adj_emm$fit, res_adj_emm$dat[res_adj_emm$dat$V=="Control",], xlab="delta_TS", ylab="laz_t3", title="Control")
simul_plot$p

simul_plot <- gam_simul_CI(res_adj_emm$fit, res_adj_emm$dat[res_adj_emm$dat$V=="WSH",], xlab="delta_TS", ylab="laz_t3", title="WSH")
simul_plot$p



res <- glm(hcz_t3~t3_ln_agp, data=d)
summary(res)
res <- glm(hcz_t3~t3_ln_agp +tr, data=d)
summary(res)

res <- glm(hcz_t3~t3_ln_agp*tr, data=d)
summary(res)


d %>% group_by(tr) %>% summarise(N=n(), MN=mean(t3_ln_agp, na.rm=T))
