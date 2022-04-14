


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



m = res_adj_emm$fit
newdata = res_adj_emm$dat
nreps = 10000
xlab = ""
ylab = ""

emm_simul_CI <- function (m, newdata, nreps = 10000, xlab = "", ylab = "", 
          title = ""#, gam_diff = NULL
          ){
  set.seed(12345)
  require(mgcv)
  require(dplyr)
  newdata <- newdata %>% mutate(dummy = 0)
  Wvars <- colnames(newdata)[!(colnames(newdata) %in% c("Y", "X", "V", "id", "dummy"))]
  for (i in Wvars) {
    if (class(newdata[, i]) == "character" | class(newdata[, 
                                                           i]) == "factor") {
      newdata[, i] <- Mode(newdata[, i])
    }
    else {
      newdata[, i] <- median(newdata[, i])
    }
  }
  newdata <- newdata[order(newdata$X), ]
  newdata$V <- as.character(newdata$V)
  newdata2 <- newdata1 <- newdata
  newdata1$V <- unique(newdata$V)[1]
  newdata2$V <- unique(newdata$V)[2]
  
  Vb <- vcov(m, unconditional = TRUE)
  pred1 <- predict(m, newdata1, se.fit = TRUE)
  pred2 <- predict(m, newdata2, se.fit = TRUE)
  fit1 <- pred1$fit
  fit2 <- pred2$fit
  se.fit1 <- pred1$se.fit
  se.fit2 <- pred2$se.fit
  BUdiff <- MASS::mvrnorm(n = nreps, mu = rep(0, nrow(Vb)), Sigma = Vb)
  Cg1 <- predict(m, newdata1, type = "lpmatrix")
  Cg2 <- predict(m, newdata2, type = "lpmatrix")
  simDev1 <- Cg1 %*% t(BUdiff)
  simDev2 <- Cg2 %*% t(BUdiff)
  absDev1 <- abs(sweep(simDev1, 1, se.fit1, FUN = "/"))
  absDev2 <- abs(sweep(simDev2, 1, se.fit2, FUN = "/"))
  masd1 <- apply(absDev1, 2L, max)
  masd2 <- apply(absDev2, 2L, max)
  crit1 <- quantile(masd1, prob = 0.95, type = 8)
  crit2 <- quantile(masd2, prob = 0.95, type = 8)
  pred1 <- data.frame(newdata1, fit = pred1$fit, se.fit = pred1$se.fit)
  pred2 <- data.frame(newdata2, fit = pred2$fit, se.fit = pred2$se.fit)
  pred1 <- mutate(pred1, uprP = fit + (2 * se.fit1), lwrP = fit - 
                   (2 * se.fit1), uprS = fit + (crit1 * se.fit1), lwrS = fit - 
                   (crit1 * se.fit1)) %>% arrange(X) %>% mutate(Vlevel=unique(newdata$V)[1])
  pred2 <- mutate(pred2, uprP = fit + (2 * se.fit2), lwrP = fit - 
                   (2 * se.fit2), uprS = fit + (crit2 * se.fit2), lwrS = fit - 
                   (crit2 * se.fit2)) %>% arrange(X) %>% mutate(Vlevel=unique(newdata$V)[2])
  pred <- bind_rows(pred1, pred2) 
  p <- ggplot(pred) + geom_ribbon(aes(x = X, ymin = lwrS, ymax = uprS, group=V, fill=V, color=V), 
                                  alpha = 0.5) + 
    # geom_path(aes(x = X, y = lwrS, group=V, color=V)) + 
    # geom_path(aes(x = X, y = uprS, group=V, color=V)) + 
    geom_path(aes(x = X, y = fit, group=V, color=V)) + 
    xlab(xlab) + ylab(ylab) + ggtitle(title)

  return(list(p = p, pred = pred))
}



d$agp_t2
d$delta_waz_t2_t3



# wei_velocity_t2_t3
# t2_ln_agp

res_adj_emm <- fit_RE_gam(d=d, X="t2_ln_agp", Y="wei_velocity_t2_t3",  W=W2_immune.W23_anthro, V="tr")
res_adj_emm$int.p

preds_adj_emm <- predict_gam_emm(fit=res_adj_emm$fit, d=res_adj_emm$dat, quantile_diff=c(0.25,0.75), Xvar="t2_ln_agp", Yvar="hcz_t3")
preds_adj_emm$res


res_adj_emm <- fit_RE_gam(d=d, X="t2_ln_agp", Y="wei_velocity_t2_t3",  W=W2_immune.W23_anthro, V="tr")
res_adj_emm$int.p

preds_adj_emm <- predict_gam_emm(fit=res_adj_emm$fit, d=res_adj_emm$dat, quantile_diff=c(0.25,0.75), Xvar="t3_ln_agp", Yvar="wei_velocity_t2_t3")
preds_adj_emm$res

emm_simul_plot <- emm_simul_CI(res_adj_emm$fit, res_adj_emm$dat, xlab="t2_ln_agp", ylab="wei_velocity_t2_t3", title="adj")
emm_simul_plot$p
