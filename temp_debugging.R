
data.full<-d



#main effect analysis
X <- c("t2_ln_agp", "t2_ln_crp", "newcrp_t2")            
Y <- c("delta_laz_t2_t3", "delta_whz_t2_t3", "delta_hcz_t2_t3", "delta_waz_t2_t3")
adj_models <- gam.adj(df = adj_models,  Xvars = "t2_ln_agp", Yvars = "delta_waz_t2_t3", 
                      Wvars = W2_immune.W23_anthro, data = d)

main_preds <- predict_gam_diff(fit=adj_models$fit[[1]], 
                          d=adj_models$dat[[1]], quantile_diff=c(0.25,0.75), Xvar="t2_ln_agp", Yvar="delta_waz_t2_t3")

main_preds$res



#look at the predicted results
adj.res <- adj_res %>% filter(Y=="delta_waz_t2_t3", X=="t2_ln_agp")

adj.res


tremm <- tremm.res %>% filter(Y=="delta_waz_t2_t3", X=="t2_ln_agp")
tremm

tremm.mod <- tremm.mods %>% filter(Y=="delta_waz_t2_t3", X=="t2_ln_agp")


summary(tremm.mod$fit[[1]])


#preds <- predict_gam_emm(fit=tremm.mods$fit[i][[1]], d=tremm.mods$dat[i][[1]], quantile_diff=c(0.25,0.75), Xvar=tremm.mods$X[i], Yvar=tremm.mods$Y[i])
main_effect_emm_mod <- predict_gam_diff(fit=tremm.mod$fit[[1]], d=tremm.mod$dat[[1]], quantile_diff=c(0.25,0.75), Xvar=tremm.mod$X, Yvar=tremm.mod$Y)
main_effect_emm_mod$res
main_preds$res

emm_preds_emm_mod <- predict_gam_emm(fit=tremm.mod$fit[[1]], d=tremm.mod$dat[[1]], quantile_diff=c(0.25,0.75), Xvar=tremm.mod$X, Yvar=tremm.mod$Y)
emm_preds_emm_mod$res


# fit model without tr being included twice
X <- c("t2_ln_agp")            
Y <- c("delta_waz_t2_t3")
# tremm.mods2 <- gam.emm(df = tremm.mods,  Xvars = X, Yvars = Y, 
#                       Wvars = W2_immune.W23_anthro[-13], data = d, Vvar = "tr")

dfunc <- outliers(Y, d)
Wset <- add_hcz(X, Y, W2_immune.W23_anthro[-13])
res_adj <- fit_RE_gam(d=dfunc, X=X, Y=Y,  W=Wset, V = "tr", forcedW = "sex")
tremm.mods2 <- data.frame(X=X, Y=Y, V = "tr", int.p =res_adj$int.p,
                  fit=I(list(res_adj$fit)), dat=I(list(res_adj$dat)), 
                  covs = as.character(paste(res_adj$covars, collapse = ", ")))

emm_preds_emm_mod2 <- predict_gam_emm(fit=tremm.mods2$fit[[1]], d=tremm.mods2$dat[[1]], quantile_diff=c(0.25,0.75), Xvar=tremm.mod$X, Yvar=tremm.mod$Y)
emm_preds_emm_mod2$res
main_effect_emm_mod2 <- predict_gam_diff(fit=tremm.mods2$fit[[1]], d=tremm.mods2$dat[[1]], quantile_diff=c(0.25,0.75), Xvar=tremm.mod$X, Yvar=tremm.mod$Y)
main_effect_emm_mod2$res

summary(tremm.mods2$fit[[1]])

#step by step
fit=tremm.mods2$fit[[1]]
d=tremm.mods2$dat[[1]]
quantile_diff=c(0.25,0.75)
Xvar=tremm.mod$X
Yvar=tremm.mod$Y
binaryX = FALSE
  
  set.seed(12345)
  require(mgcv)
  require(dplyr)
  d$dummy <- 0
  
  Wvars <- colnames(d)[!(colnames(d) %in% c("Y", "X", "V",
                                            "id", "dummy"))]
  for (i in Wvars) {
    if (class(d[, i]) == "character" | class(d[, i]) ==
        "factor") {
      d[, i] <- Mode(d[, i])
    }
    else {
      d[, i] <- median(d[, i])
    }
  }
  d <- d[order(d$X), ]
  
  
  if (binaryX == F) {
    q1 <- unname(quantile(d$X, quantile_diff[1]))
    q3 <- unname(quantile(d$X, quantile_diff[2]))
    q1_pos <- which(abs(d$X - q1) == min(abs(d$X - q1)))[1]
    q3_pos <- which(abs(d$X - q3) == min(abs(d$X - q3)))[1]
    d$X[q1_pos] <- q1
    d$X[q3_pos] <- q3
  }
  if (binaryX == T) {
    q1 <- min(d$X)
    q3 <- max(d$X)
    q1_pos <- 1
    q3_pos <- nrow(d)
    d$X[q1_pos] <- min(d$X)
    d$X[q3_pos] <- max(d$X)
  }
  
  #Add specific rows for the quartile predictions
  Nrows <- nrow(d)
  d <- bind_rows(d, d[1:2,])
  d$X[c(Nrows+1,Nrows+2)] <- c(q1, q3)
  
  dfull <- d
  plotdf <- res <- NULL
  
  #Detect if modifier is
  if(class(dfull$V)!="numeric"){
    for(i in unique(dfull$V)){
      d <- dfull
      d$V <- i
      preds <- predict(fit, newdata = d, type = "response")
      Xp <- predict(fit, newdata = d, type = "lpmatrix")
      #Xp <- Xp[order(d$X), ] #XXXX this is the issue!
      #Xp <- Xp[Nrows+1,] #prediction at the 25th percentile
      #diff <- t(apply(Xp, 1, function(x) x - Xp[Nrows+1, ]))
      diff <- t(apply(Xp, 1, function(x) x - Xp[Nrows+1, ]))
      point.diff <- diff %*% coef(fit)
      se.diff <- sqrt(diag(diff %*% vcov(fit) %*% t(diff)))
      lb.diff <- point.diff - 1.96 * se.diff
      ub.diff <- point.diff + 1.96 * se.diff
      Zval <- abs(point.diff/se.diff)
      Pval <- exp(-0.717 * Zval - 0.416 * Zval^2)
      resdf <- data.frame(Y = Yvar, X = Xvar, Vlevel=i, N = nrow(d), q1 = d$X[q1_pos],
                          q3 = d$X[q3_pos], pred.q1 = preds[q1_pos], pred.q3 = preds[q3_pos],
                          point.diff, lb.diff = lb.diff, ub.diff = ub.diff, Pval = Pval)
      
      if(binaryX==T){
        temp_res <- resdf[1, ]
      }else{
        temp_res <- resdf[nrow(resdf), ]
      }
      
      
      plotdf <- bind_rows(plotdf, resdf)
      
      res <- bind_rows(res, temp_res)
    }
    
  }else{
    
    quartiles <- as.numeric(summary(dfull$V)[c(2,5)])
    
    for(i in unique(quartiles)){
      d <- dfull
      d$V <- i
      preds <- predict(fit, newdata = d, type = "response")
      Xp <- predict(fit, newdata = d, type = "lpmatrix")
      Xp <- Xp[order(d$X), ]
      diff <- t(apply(Xp, 1, function(x) x - Xp[Nrows+1, ]))
      point.diff <- diff %*% coef(fit)
      se.diff <- sqrt(diag(diff %*% vcov(fit) %*% t(diff)))
      lb.diff <- point.diff - 1.96 * se.diff
      ub.diff <- point.diff + 1.96 * se.diff
      Zval <- abs(point.diff/se.diff)
      Pval <- exp(-0.717 * Zval - 0.416 * Zval^2)
      resdf <- data.frame(Y = Yvar, X = Xvar, Vlevel=i, N = nrow(d), q1 = d$X[q1_pos],
                          q3 = d$X[q3_pos], pred.q1 = preds[q1_pos], pred.q3 = preds[q3_pos],
                          point.diff, lb.diff = lb.diff, ub.diff = ub.diff, Pval = Pval)
      
      #if(binaryX == T){
      temp_res <- resdf[nrow(resdf), ]
      # }else{
      #   temp_res <- resdf[q3_pos, ]
      # }
      
      plotdf <- bind_rows(plotdf, resdf[1:Nrows,])
      
      res <- bind_rows(res, temp_res)
    }
    
  }






