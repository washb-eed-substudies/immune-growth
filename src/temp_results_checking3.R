
dfull=d

# velocity
X <- c("t2_ln_agp")            
Y <- c( "wei_velocity_t2_t3")

    i=Xvars = X
  j=Yvars = Y
  Wvars = W2_immune.W23_anthro
  data = dfull
  Vvar = "tr"
  
  
  
  dfunc <- outliers(j, data)
  Wset <- add_hcz(i, j, Wvars)
  res_adj <- fit_RE_gam(d=dfunc, X=i, Y=j,  W=Wset, V = Vvar)
  res <- data.frame(X=i, Y=j, V = Vvar, int.p =res_adj$int.p, fit=I(list(res_adj$fit)), dat=I(list(res_adj$dat)), covs = as.character(paste(res_adj$covars, collapse = ", ")))
  res
  
  res_adj$int.p
  
  res_adj2 <- fit_RE_gam(d=dfunc, X=i, Y=j,  W=Wset[Wset!="tr"], V = Vvar)


  res_adj_emm <- fit_RE_gam(d=dfunc, X=i, Y=j,  W=Wset[Wset!="tr"], V = Vvar)
  res_adj_emm$int.p
  summary(res_adj_emm$fit)
  preds_adj_emm <- predict_gam_emm(fit=res_adj_emm$fit, d=res_adj_emm$dat, quantile_diff=c(0.25,0.75), Xvar="t2_ln_agp", Yvar="wei_velocity_t2_t3")
  preds_adj_emm$res

  
  emm_simul_plot <- emm_simul_CI(res_adj_emm$fit, res_adj_emm$dat, xlab="t2_ln_agp", ylab="wei_velocity_t2_t3", title="adj")
  emm_simul_plot$p + geom_vline(xintercept=-0.2743744 ) + geom_vline(xintercept=0.3344291 )
 
  
  res_adj_main <- fit_RE_gam(d=dfunc, X=i, Y=j,  W=Wset[Wset!="tr"], V = NULL)
  preds_adj_gam <- predict_gam_diff(fit=res_adj_main$fit, d=res_adj_main$dat, quantile_diff=c(0.25,0.75), Xvar="t2_ln_agp", Yvar="wei_velocity_t2_t3")
  preds_adj_gam$res
  
  gam_simul_plot <- gam_simul_CI(res_adj_main$fit, res_adj_main$dat, xlab="t2_ln_agp", ylab="wei_velocity_t2_t3", title="adj")
  gam_simul_plot$p + geom_vline(xintercept=-0.2743744 ) + geom_vline(xintercept=0.3344291 )
  
  
  
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
  
  
  
  
  d=dfunc
  X=i
  Y=j
  W=Wset[Wset!="tr"]
  V = "tr"
  id = "clusterid"
  forcedW = NULL
  family = "gaussian"
  pval = 0.2
  print = TRUE
  
  cat("\nNon-prescreened covariates: ", paste(forcedW, 
                                              sep = "", collapse = ", "), "\n")
  set.seed(12345)
  require(mgcv)
  require(dplyr)
  require(faraway)
  require(lmtest)
  
  W <- subset(d, select = W)
  
  Y <- subset(d, select = Y)
  colnames(Y) <- "Y"
  X <- subset(d, select = X)
  colnames(X) <- "X"
  id <- subset(d, select = id)
  colnames(id) <- "id"
  
  Vvar <- subset(d, select = V)
  colnames(Vvar) <- "V"
  
  gamdat <- data.frame(Y, X, id, Vvar, W)
  
  if(sum(is.na(forcedW)) != 0){
    colnamesW <- names(W)
  }else{
    if(is.null(forcedW)){
      Wnames <- names(W)
      forcedW <- c(Wnames[Wnames == "tr" | grepl("age_", 
                                                 Wnames) | grepl("agedays_", Wnames) | 
                            grepl("ageday_", Wnames)])
    }
    cat("\nNon-prescreened covariates: ", paste(forcedW, 
                                                sep = "", collapse = ", "), "\n")
    colnamesW <- names(W)[!(names(W) %in% forcedW)]
  }
  screenW <- subset(gamdat, select = colnamesW)
  
  if(!is.null(screenW)){
    if(print == TRUE){
      cat("\n-----------------------------------------\nPre-screening the adjustment covariates:\n-----------------------------------------\n")
    }
    suppressWarnings(Wscreen <- washb_prescreen(Y = gamdat$Y, 
                                                Ws = screenW, family = family, pval = pval, print = print))
    if(!is.null(forcedW)){
      Wscreen <- c(as.character(Wscreen), as.character(forcedW))
    }
    W <- subset(gamdat, select = Wscreen)
    Wdf <- W
    Wdf$constant <- rep(1, nrow(gamdat))
    for(i in 1:ncol(W)){
      tmp <- glm(constant ~ ., data = Wdf, family = family)
      todrop <- NULL
      todrop <- suppressWarnings(names(tmp$coefficients)[-1][as.vector(vif(tmp)) > 
                                                               10][1])
      if(!is.null(todrop) & !is.na(todrop)){
        collinear_vars <- c(collinear_vars, todrop)
        Wdf <- Wdf[, colnames(Wdf) != todrop]
      }
    }
    to_keep <- colnames(W)[!(colnames(W) %in% collinear_vars)]
    if(length(to_keep) != length(colnames(W))){
      cat("\nDropped for collinearity with other covariates:\n", 
          colnames(W)[!(colnames(W) %in% to_keep)])
    }
    W_processed <- W[which(colnames(W) %in% to_keep)]
    Wscreen <- colnames(W_processed)
    cat("\n\nCovariated included in model:\n", Wscreen)
  }
  
  d <- subset(gamdat, select = c("Y", "X", "id", "V", Wscreen))
  
  fullrows <- nrow(d)
  d <- d %>% filter(!is.na(Y))
  Yrows <- nrow(d)
  cat("\nRows dropped due to missing outcome: ", fullrows - 
        Yrows, "\n")
  d <- d %>% filter(!is.na(X))
  Xrows <- nrow(d)
  cat("Rows dropped due to missing exposure: ", Yrows - 
        Xrows, "\n")
  if(!is.null(W) & length(Wscreen) > 0){
    cat("Percent missingness by covariate:\n")
    print(sapply(d[, -c(1:3)], function(x) round(sum(is.na(x))/nrow(X) * 
                                                   100, 1)))
    d <- d[complete.cases(d), ]
    cat("\nRows dropped due to missing covariates: ", 
        Xrows - nrow(d), "\n")
  }
  cat("Final sample size: ", nrow(d), "\n")
  d$dummy <- 1
  
  
  Ws <- subset(gamdat, select = c(Wscreen))
  W_factors <- colnames(Ws)[(grepl("factor", sapply(Ws, 
                                                    class)) | grepl("character", sapply(Ws, class)))]
  W_numeric <- colnames(Ws)[(grepl("integer", sapply(Ws, 
                                                     class)) | grepl("numeric", sapply(Ws, class)))]
  indicator_vec <- rep(TRUE, length(W_numeric))
  for (i in 1:length(W_numeric)){
    N_unique <- length(unique(Ws[, W_numeric[i]]))
    if(N_unique > 20){
      indicator_vec[i] <- FALSE
    }
  }
  
  W_indicator <- W_numeric[indicator_vec]
  W_continious <- W_numeric[!indicator_vec]
  if(length(W_continious) > 0){
    eq_num <- paste0("s(", W_continious, ", bs=\"cr\")", 
                     collapse = " + ")
  }else{
    eq_num = NULL
  }
  if(length(W_factors) + length(W_indicator) > 0){
    eq_fact <- paste0(" + ", paste0(c(W_factors, 
                                      W_indicator), collapse = " + "))
  }else{
    eq_fact = NULL
  }
  
  
  form <- paste0("Y~s(X, bs=\"cr\")+ V + X*V+", 
                 eq_fact, " +", eq_num, "+ s(id,bs=\"re\",by=dummy)")
  form <- gsub("+ +", "+", form, fixed = TRUE)
  equation <- as.formula(form)
  form_null <- paste0("Y~s(X, bs=\"cr\")+ V + ", 
                      eq_fact, " +", eq_num, "+ s(id,bs=\"re\",by=dummy)")
  form_null <- gsub("+ +", "+", form_null, fixed = TRUE)
  form_null <- gsub("+  +", "+", form_null, fixed = TRUE)
  equation_null <- as.formula(form_null)
  
  
  fit <- mgcv::gam(formula = equation, data = d)
  fit_null <- mgcv::gam(formula = equation_null, data = d)
  LRp <- lrtest(fit, fit_null)[2, 5]
  summary(fit)
  summary(fit_null)
  
  AIC(fit, fit_null)
  LRp
  
  glm_formula <- gsub("\"","",form)
  glm_formula <- gsub(", bs=cr","",glm_formula)
  glm_formula <- gsub(",bs=re,by=dummy","",glm_formula)
  glm_formula <- gsub(")","",glm_formula)
  glm_formula <- gsub("s\\(","",glm_formula)
  glm_equation = as.formula(glm_formula)
  
  glm_formula_null <- gsub("\"","",form_null)
  glm_formula_null <- gsub(", bs=cr","",glm_formula_null)
  glm_formula_null <- gsub(",bs=re,by=dummy","",glm_formula_null)
  glm_formula_null <- gsub(")","",glm_formula_null)
  glm_formula_null <- gsub("s\\(","",glm_formula_null)
  glm_equation_null = as.formula(glm_formula_null)
  
  glm.fit <- glm(formula = glm_equation, data = d)
  glm.fit_null <- glm(formula = glm_equation_null, data = d)
  glm.LRp <- lrtest(glm.fit, glm.fit_null)[2, 5]
  glm.LRp
  
  summary(glm.fit)
  summary(glm.fit_null)
  
  
  