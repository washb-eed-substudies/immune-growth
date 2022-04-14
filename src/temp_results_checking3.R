


# velocity
X <- c("t2_ln_agp")            
Y <- c( "wei_velocity_t2_t3")

  
  df = tremm.mods
  i=Xvars = X
  j=Yvars = Y
  Wvars = W2_immune.W23_anthro
  data = d
  Vvar = "tr"
  
  
  
  dfunc <- outliers(j, data)
  Wset <- add_hcz(i, j, Wvars)
  res_adj <- fit_RE_gam(d=dfunc, X=i, Y=j,  W=Wset, V = Vvar)
  res <- data.frame(X=i, Y=j, V = Vvar, int.p =res_adj$int.p, fit=I(list(res_adj$fit)), dat=I(list(res_adj$dat)), covs = as.character(paste(res_adj$covars, collapse = ", ")))
  res
  
  res_adj$int.p
  
  res_adj2 <- fit_RE_gam(d=dfunc, X=i, Y=j,  W=Wset[Wset!="tr"], V = Vvar)
  res_adj2$int.p
  
  
  res_adj_emm <- fit_RE_gam(d=dfunc, X=i, Y=j,  W=Wset[Wset!="tr"], V = Vvar)
  
  preds_adj_emm <- predict_gam_emm(fit=res_adj_emm$fit, d=res_adj_emm$dat, quantile_diff=c(0.25,0.75), Xvar="t3_ln_agp", Yvar="wei_velocity_t2_t3")
  preds_adj_emm$res
  
  emm_simul_plot <- emm_simul_CI(res_adj_emm$fit, res_adj_emm$dat, xlab="t2_ln_agp", ylab="wei_velocity_t2_t3", title="adj")
  emm_simul_plot$p + geom_vline(xintercept=-0.2763686) + geom_vline(xintercept=0.3318642)
 
  equation <- as.formula(paste0("Y~s(X, bs = \"cr\") + V + X * V + birthord + momedu + hfiacat + 
    floor + nose7d_t2 + Nlt18 + month_bt2 + s(momage, bs = \"cr\") + 
    s(momheight, bs = \"cr\") + s(HHwealth_scaled, bs = \"cr\") + 
    s(laz_t2, bs = \"cr\") + s(waz_t2, bs = \"cr\") + 
    s(ageday_at2, bs = \"cr\") + s(ageday_at3, bs = \"cr\") + 
    s(id, bs = \"re\", by = dummy)"))
  
  dfunc2 <- dfunc %>% mutate(X=t2_ln_agp, Y=wei_velocity_t2_t3, V=tr, id = clusterid, dummy=1)
   fit <- mgcv::gam(formula = equation, data = dfunc2)
   
   
   equation_null <- as.formula(paste0("Y~s(X, bs = \"cr\") + V + birthord + momedu + hfiacat + 
    floor + nose7d_t2 + Nlt18 + month_bt2 + s(momage, bs = \"cr\") + 
    s(momheight, bs = \"cr\") + s(HHwealth_scaled, bs = \"cr\") + 
    s(laz_t2, bs = \"cr\") + s(waz_t2, bs = \"cr\") + 
    s(ageday_at2, bs = \"cr\") + s(ageday_at3, bs = \"cr\") + 
    s(id, bs = \"re\", by = dummy)"))
   
   fit_null <- mgcv::gam(formula = equation_null,  data = dfunc2)
  LRp <- lrtest(fit, fit_null)
  LRp
  
  
  
  
  
  