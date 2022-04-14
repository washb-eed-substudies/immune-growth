

rm(list=ls())

source(here::here("0-config.R"))

d <- readRDS(paste0(dropboxDir,"Data/Cleaned/Audrie/bangladesh-immune-growth-analysis-dataset.rds"))

d <- d %>% filter(!is.na(hcz_t2), !is.na(t2_ln_agp))

Wvars<-c("sex","birthord", "momage","momheight","momedu", 
         "hfiacat", "Nlt18","Ncomp", "watmin", "walls", 
         "floor", "HHwealth", "tr", "cesd_sum_t2", 
         "ari7d_t2", "diar7d_t2", "nose7d_t2", "life_viol_any_t3")
growth_t2 <- c("laz_t2", "waz_t2", "whz_t2", "hcz_t2")

Wvars[!(Wvars %in% colnames(d))]



#Add in time varying covariates:
Wvars2<-c("ageday_bt2", "ageday_at2",  "month_bt2", "month_at2", "laz_t1_cat", "waz_t1_cat", "hcz_t1_cat") 
Wvars2_cont<-c("ageday_bt2", "ageday_at2",  "month_bt2", "month_at2", "laz_t1", "waz_t1", "hcz_t1") 

W2_immmune.W2_anthro <- c(Wvars, Wvars2) %>% unique(.)
W2_immmune.W2_anthro_cont <- c(Wvars, Wvars2_cont) %>% unique(.)



#---------------------------
# unadjusted
#---------------------------

#Fit GAM model with random effects for childid
res_unadj <- fit_RE_gam(d=d, X="t2_ln_agp", Y="hcz_t2",  W=NULL)
preds_unadj <- predict_gam_diff(fit=res_unadj$fit, d=res_unadj$dat, quantile_diff=c(0.25,0.75), Xvar="t2_ln_agp", Yvar="hcz_t2")
preds_unadj$res
simul_plot_unadj <- gam_simul_CI(res_unadj$fit, res_unadj$dat, xlab="t2_ln_agp", ylab="hcz_t2", title="unadj") 
simul_plot_unadj$p + 
  geom_vline(xintercept = preds_unadj$res$q1) +
  geom_vline(xintercept = preds_unadj$res$q3) +
  ggtitle(paste0(round(preds_unadj$res$point.diff,2)," (",
                 round(preds_unadj$res$lb.diff,2), ", ",
                 round(preds_unadj$res$ub.diff,2), ")")) +
  geom_smooth(aes(x=t2_ln_agp, y=hcz_t2), method="lm",  se=F, data=d)


simul_plot_unadj$p + 
  geom_vline(xintercept = preds_unadj$res$q1) +
  geom_vline(xintercept = preds_unadj$res$q3) +
  geom_smooth(aes(x=t2_ln_agp, y=hcz_t2), method="lm",  se=F, data=d)


simul_plot_unadj$p + 
  geom_vline(xintercept = preds_unadj$res$q1) +
  geom_vline(xintercept = preds_unadj$res$q3) +
  geom_smooth(aes(x=t2_ln_agp, y=hcz_t2), method="lm",  se=F, data=d) +
  geom_smooth(aes(x=t2_ln_agp, y=hcz_t2), method="loess",  se=F, color="red", data=d) +
  geom_point(aes(x=t2_ln_agp, y=hcz_t2), alpha=0.1, data=d) +
  coord_cartesian(xlim=c(-1,1), ylim=c(-3,-1))



#---------------------------
# adjusted, categorical growth
#---------------------------


#GAM fit
res_adj <- fit_RE_gam(d=d, X="t2_ln_agp", Y="hcz_t2",  W=W2_immmune.W2_anthro)
preds_adj <- predict_gam_diff(fit=res_adj$fit, d=res_adj$dat, quantile_diff=c(0.25,0.75), Xvar="t2_ln_agp", Yvar="hcz_t2")
preds_adj$res
res_adj$covars

#glm predictions
glm.data <- d %>% select("t2_ln_agp", "hcz_t2", !!(res_adj$covars)) %>% droplevels()
res <- glm(hcz_t2~., data=glm.data)
summary(res)
preds <- predict(res,type="response")
glm.data.pred <- glm.data
for(i in 3:ncol(glm.data.pred)){
  if(class(glm.data.pred[,i])=="character"|class(glm.data.pred[,i])=="factor"){
    glm.data.pred[,i] <- Mode(glm.data.pred[,i])
  }else{
    glm.data.pred[,i] <- median(glm.data.pred[,i], na.rm=T)
  }
}
glm.data.pred <- glm.data.pred %>% arrange(t2_ln_agp)
preds <- predict(res,newdata=glm.data.pred,type="response")

predY.glm <- data.frame(x=glm.data.pred$t2_ln_agp, y=preds)

simul_plot <- gam_simul_CI(res_adj$fit, res_adj$dat, xlab="t2_ln_agp", ylab="hcz_t2", title="adj")
simul_plot$p + 
  geom_vline(xintercept = preds_adj$res$q1) +
  geom_vline(xintercept = preds_adj$res$q3) +
  geom_line(aes(x=x, y=y), data=predY.glm, color="blue")



#---------------------------
#continious, missing data
#---------------------------
dim(d)
dsub <- d %>% filter(!is.na(hcz_t1), !is.na(waz_t1), !is.na(laz_t1))
dim(dsub)


res_adj_cont <- fit_RE_gam(d=dsub, X="t2_ln_agp", Y="hcz_t2",  W=W2_immmune.W2_anthro_cont)
preds_adj_cont <- predict_gam_diff(fit=res_adj_cont$fit, d=res_adj_cont$dat, quantile_diff=c(0.25,0.75), Xvar="t2_ln_agp", Yvar="hcz_t2")
preds_adj_cont$res

#glm predictions
glm.data <- dsub %>% select("t2_ln_agp", "hcz_t2", !!(res_adj_cont$covars)) %>% droplevels()
res <- glm(hcz_t2~., data=glm.data)
summary(res)
preds <- predict(res,type="response")
glm.data.pred <- glm.data
for(i in 3:ncol(glm.data.pred)){
  if(class(glm.data.pred[,i])=="character"|class(glm.data.pred[,i])=="factor"){
    glm.data.pred[,i] <- Mode(glm.data.pred[,i])
  }else{
    glm.data.pred[,i] <- median(glm.data.pred[,i], na.rm=T)
  }
}
glm.data.pred <- glm.data.pred %>% arrange(t2_ln_agp)
preds <- predict(res,newdata=glm.data.pred,type="response")

predY.glm_cont <- data.frame(x=glm.data.pred$t2_ln_agp, y=preds)

simul_plot_cont <- gam_simul_CI(res_adj_cont$fit, res_adj_cont$dat, xlab="t2_ln_agp", ylab="hcz_t2", title="adj")
simul_plot_cont$p + 
  geom_vline(xintercept = preds_adj_cont$res$q1) +
  geom_vline(xintercept = preds_adj_cont$res$q3) +
  geom_line(aes(x=x, y=y), data=predY.glm_cont, color="blue")




#---------------------------
#continious, imputed data
#---------------------------
d$hcz_t1[is.na(d$hcz_t1)] <- median(d$hcz_t1, na.rm=T)
d$laz_t1[is.na(d$laz_t1)] <- median(d$laz_t1, na.rm=T)
d$waz_t1[is.na(d$waz_t1)] <- median(d$waz_t1, na.rm=T)



res_adj_cont_imp <- fit_RE_gam(d=d, X="t2_ln_agp", Y="hcz_t2",  W=W2_immmune.W2_anthro_cont)
preds_adj_cont_imp <- predict_gam_diff(fit=res_adj_cont_imp$fit, d=res_adj_cont_imp$dat, quantile_diff=c(0.25,0.75), Xvar="t2_ln_agp", Yvar="hcz_t2")
preds_adj_cont_imp$res

#glm predictions
glm.data <- d %>% select("t2_ln_agp", "hcz_t2", !!(res_adj_cont_imp$covars)) %>% droplevels()
res <- glm(hcz_t2~., data=glm.data)
summary(res)
preds <- predict(res,type="response")
glm.data.pred <- glm.data
for(i in 3:ncol(glm.data.pred)){
  if(class(glm.data.pred[,i])=="character"|class(glm.data.pred[,i])=="factor"){
    glm.data.pred[,i] <- Mode(glm.data.pred[,i])
  }else{
    glm.data.pred[,i] <- median(glm.data.pred[,i], na.rm=T)
  }
}
glm.data.pred <- glm.data.pred %>% arrange(t2_ln_agp)
preds <- predict(res,newdata=glm.data.pred,type="response")

predY.glm_cont_imp <- data.frame(x=glm.data.pred$t2_ln_agp, y=preds)

simul_plot_cont_imp <- gam_simul_CI(res_adj_cont_imp$fit, res_adj_cont_imp$dat, xlab="t2_ln_agp", ylab="hcz_t2", title="adj")
simul_plot_cont_imp$p + 
  geom_vline(xintercept = preds_adj_cont_imp$res$q1) +
  geom_vline(xintercept = preds_adj_cont_imp$res$q3) +
  geom_line(aes(x=x, y=y), data=predY.glm_cont_imp, color="blue")

