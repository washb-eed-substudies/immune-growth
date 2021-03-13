rm(list=ls())
source(here::here("0-config.R"))

######################
###Load in data
######################

#Load in enrollment data,blinded tr data, stool data for adjusted analysis. Use read.dta() to read the .dta files, or read.csv() to 
#read .csv files. Use stringAsFactors=TRUE so that any character-based variable will be read in as a factor.
d <- read.csv(paste0(dropboxDir,"Data/Cleaned/Audrie/bangladesh-dm-ee-immune-growth-covariates-immunelab-anthro.csv"), stringsAsFactors = TRUE)

summary(d$t2_ln_il2)
summary(d$t3_ln_il2)
summary(d$t3_ln_il2 - d$t2_ln_il2)
summary(d$il1_t3 - d$il1_t2)
summary(log(d$il1_t3) - log(d$il1_t2))
summary(log(d$il1_t3 - d$il1_t2))
summary(d$d23_il2)


#drop Z-score, sd, and ratio measures
d <- d[,!(grepl("^(z_)",colnames(d)) | grepl("^(sd_)",colnames(d)))]
d <- d %>% subset(., select = -c(t2_ratio_pro_il10,
                                 t2_ratio_th1_il10,
                                 t2_ratio_th2_il10,
                                 t2_ratio_th17_il10,
                                 t2_ratio_th1_th2,
                                 t2_ratio_th1_th17,
                                 t3_ratio_pro_il10,
                                 t3_ratio_th1_il10,
                                 t3_ratio_th2_il10,
                                 t3_ratio_th17_il10,
                                 t3_ratio_th1_th2,
                                 t3_ratio_th1_th17,
                                 d23_ratio_pro_il10, 
                                 d23_ratio_th1_il10,  
                                 d23_ratio_th2_il10,  
                                 d23_ratio_th17_il10,
                                 d23_ratio_th1_th2,
                                 d23_ratio_th1_th17))

x=c("il1_t2", "il6_t2", "tnfa_t2")[1]
summary(as.vector(scale(d[,x], center = FALSE, scale = apply(as.matrix(d[,x]), 2, sd, na.rm = TRUE))))
x=c("il1_t2", "il6_t2", "tnfa_t2")[2]
summary(as.vector(scale(d[,x], center = FALSE, scale = apply(as.matrix(d[,x]), 2, sd, na.rm = TRUE))))
x=c("il1_t2", "il6_t2", "tnfa_t2")[3]
summary(as.vector(scale(d[,x], center = FALSE, scale = apply(as.matrix(d[,x]), 2, sd, na.rm = TRUE))))


#function to create composite score
create_score <- function(d, numerator_vars=c("il1_t2", "il6_t2", "tnfa_t2"), denominator_vars="il10_t2", varname="t2_ratio_pro_il10"){
  for(i in numerator_vars){
    if(i==numerator_vars[1]){
      x = as.vector(scale(d[,i], center = FALSE, scale = apply(as.matrix(d[,i]), 2, sd, na.rm = TRUE)))
    }else{
      x = x + as.vector(scale(d[,i], center = FALSE, scale = apply(as.matrix(d[,i]), 2, sd, na.rm = TRUE)))
    }
  }
  summary(x)
  
  
  for(i in denominator_vars){
    if(i==denominator_vars[1]){
      y = as.vector(scale(d[,i], center = FALSE, scale = apply(as.matrix(d[,i]), 2, sd, na.rm = TRUE)))
    }else{
      y = y + as.vector(scale(d[,i], center = FALSE, scale = apply(as.matrix(d[,i]), 2, sd, na.rm = TRUE)))
    }
  }
  summary(y)
  
  score=log(x/y)
  summary(score)
  d$score <- score
  colnames(d)[ncol(d)] <- varname
  return(d)
}





# *Pro-inflammatory cytokines / IL-10
# *(IL-1 + IL-6 + TNF-a) / IL-10
d <- create_score(d, numerator_vars=c("il1_t2", "il6_t2", "tnfa_t2"), denominator_vars="il10_t2", varname="t2_ratio_pro_il10")
summary(d$t2_ratio_pro_il10)
ggplot(d, aes(x=t2_ratio_pro_il10)) + geom_density()

# *Th1 / IL-10
# *(IL-12 + IFN) / IL-10
# gen t2_ratio_th1_il10 = (il12_t2 + ifng_t2) / il10_t2
d <- create_score(d, numerator_vars=c("il12_t2", "ifng_t2"), denominator_vars="il10_t2", varname="t2_ratio_th1_il10")
summary(d$t2_ratio_th1_il10)
ggplot(d, aes(x=t2_ratio_th1_il10)) + geom_density()

# *Th2 / IL-10 
# *(IL-4 + IL-5 + IL-13) / IL-10
# gen t2_ratio_th2_il10 = (il4_t2 + il5_t2 + il13_t2) / il10_t2
d <- create_score(d, numerator_vars=c("il4_t2", "il5_t2", "il13_t2"), denominator_vars="il10_t2", varname="t2_ratio_th2_il10")
summary(d$t2_ratio_th2_il10)
ggplot(d, aes(x=t2_ratio_th2_il10)) + geom_density()


# *Th17 / IL-10
# *(IL-17A + IL-21) / IL-10
# gen t2_ratio_th17_il10 = (il17_t2 + il21_t2) / il10_t2
d <- create_score(d, numerator_vars=c("il17_t2", "il21_t2"), denominator_vars="il10_t2", varname="t2_ratio_th17_il10")
summary(d$t2_ratio_th17_il10)
ggplot(d, aes(x=t2_ratio_th17_il10)) + geom_density()


# *Th1 / Th2
# *(IL-12 + IFN) / (IL-4 + IL-5 + IL-13)
# gen t2_ratio_th1_th2 = (il12_t2 + ifng_t2) / (il4_t2 + il5_t2 + il13_t2)
d <- create_score(d, numerator_vars=c("il12_t2", "ifng_t2"), denominator_vars=c("il4_t2", "il5_t2", "il13_t2"), varname="t2_ratio_th1_th2")
summary(d$t2_ratio_th1_th2)
ggplot(d, aes(x=t2_ratio_th1_th2)) + geom_density()


# *Th1 / Th17
# *(IL-12+IFN) / (IL-17A + IL-21)
# gen t2_ratio_th1_th17 = (il12_t2 + ifng_t2) / (il17_t2 + il21_t2)
d <- create_score(d, numerator_vars=c("il12_t2", "ifng_t2"), denominator_vars=c("il17_t2", "il21_t2"), varname="t2_ratio_th1_th17")
summary(d$t2_ratio_th1_th17)
ggplot(d, aes(x=t2_ratio_th1_th17)) + geom_density()


# 
# *Pro-inflammatory cytokines / IL-10
# *(IL-1 + IL-6 + TNF-a) / IL-10
# gen t3_ratio_pro_il10 = (il1_t3 + il6_t3 + tnf_t3) / il10_t3
d <- create_score(d, numerator_vars=c("il1_t3", "il6_t3", "tnfa_t3"), denominator_vars="il10_t3", varname="t3_ratio_pro_il10")
summary(d$t3_ratio_pro_il10)
ggplot(d, aes(x=t3_ratio_pro_il10)) + geom_density()


# *Th1 / IL-10
# *(IL-12 + IFN) / IL-10
# gen t3_ratio_th1_il10 = (il12_t3 + ifn_t3) / il10_t3
d <- create_score(d, numerator_vars=c("il12_t3", "ifng_t3"), denominator_vars="il10_t3", varname="t3_ratio_th1_il10")
summary(d$t3_ratio_th1_il10)
ggplot(d, aes(x=t3_ratio_th1_il10)) + geom_density()


# *Th2 / IL-10 
# *(IL-4 + IL-5 + IL-13) / IL-10
# gen t3_ratio_th2_il10 = (il4_t3 + il5_t3 + il13_t3) / il10_t3
d <- create_score(d, numerator_vars=c("il4_t3", "il5_t3", "il13_t3"), denominator_vars="il10_t3", varname="t3_ratio_th2_il10")
summary(d$t3_ratio_th2_il10)
ggplot(d, aes(x=t3_ratio_th2_il10)) + geom_density()


# *Th17 / IL-10
# *(IL-17A + IL-21) / IL-10
# gen t3_ratio_th17_il10 = (il17_t3 + il21_t3) / il10_t3
d <- create_score(d, numerator_vars=c("il17_t3", "il21_t3"), denominator_vars="il10_t3", varname="t3_ratio_th17_il10")
summary(d$t3_ratio_th17_il10)
ggplot(d, aes(x=t3_ratio_th17_il10)) + geom_density()


# *Th1 / Th2
# *(IL-12 + IFN) / (IL-4 + IL-5 + IL-13)
# gen t3_ratio_th1_th2 = (il12_t3 + ifn_t3) / (il4_t3 + il5_t3 + il13_t3)
d <- create_score(d, numerator_vars=c("il12_t3", "ifng_t3"), denominator_vars=c("il4_t3", "il5_t3", "il13_t3"), varname="t3_ratio_th1_th2")
summary(d$t3_ratio_th1_th2)
ggplot(d, aes(x=t3_ratio_th1_th2)) + geom_density()


# *Th1 / Th17
# *(IL-12+IFN) / (IL-17A + IL-21)
# gen t3_ratio_th1_th17 = (il12_t3 + ifn_t3) / (il17_t3 + il21_t3)
d <- create_score(d, numerator_vars=c("il12_t3", "ifng_t3"), denominator_vars=c("il17_t3", "il21_t3"), varname="t3_ratio_th1_th17")
summary(d$t3_ratio_th1_th17)
ggplot(d, aes(x=t3_ratio_th1_th17)) + geom_density()

d <- d %>% mutate(
  d23_ratio_pro_il10 = t3_ratio_pro_il10 - t2_ratio_pro_il10,
  d23_ratio_th1_il10 = t3_ratio_th1_il10 - t2_ratio_th1_il10,
  d23_ratio_th2_il10 = t3_ratio_th2_il10 - t2_ratio_th2_il10,
  d23_ratio_th17_il10 = t3_ratio_th17_il10 - t2_ratio_th17_il10,
  d23_ratio_th1_th2 = t3_ratio_th1_th2 - t2_ratio_th1_th2, 
  d23_ratio_th1_th17 = t3_ratio_th1_th17 - t2_ratio_th1_th17)

ggplot(d, aes(x=d23_ratio_pro_il10)) + geom_density()
ggplot(d, aes(x=d23_ratio_th1_il10)) + geom_density()
ggplot(d, aes(x=d23_ratio_th2_il10)) + geom_density()
ggplot(d, aes(x=d23_ratio_th17_il10)) + geom_density()
ggplot(d, aes(x=d23_ratio_th1_th2)) + geom_density()
ggplot(d, aes(x=d23_ratio_th1_th17)) + geom_density()


# add sum score
sum_score <- read.csv(here('results/child sum score/child immune sum scores.csv')) %>% select(-X)
sum_score_d <- left_join(d, sum_score, by='childid')


# add hhwealth
d_hhwealth <- read.csv("C:/Users/Sophia/Documents/ee-secondary/sophia scripts/hhwealth.csv")
dfull <- left_join(sum_score_d, d_hhwealth, by="dataid")


# check covariate missingness
Wvars<-c("sex","birthord", "momage","momheight","momedu", 
                  "hfiacat", "Nlt18","Ncomp", "watmin", "walls", 
                  "floor", "HHwealth", "tr", "cesd_sum_t2", 
                  "ari7d_t2", "diar7d_t2", "nose7d_t2", "life_viol_any_t3")

#Add in time varying covariates:
Wvars2<-c("ageday_bt2", "ageday_at2",  "month_bt2", "month_at2", "laz_t1", "waz_t1") 
Wvars3<-c("ageday_bt3", "ageday_at3", "month_bt3", "month_at3", 
          "laz_t2", "waz_t2", "cesd_sum_ee_t3", "pss_sum_mom_t3", 
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

add_hcz <- function(j, W){
  if (j=="hcz_t3"){Wset=c(W, "hcz_t2")}
  else if (j=="hcz_t2"){Wset=c(W, "hcz_t1")}
  else {Wset=W}
  return(Wset)
}

generate_miss_tbl <- function(Wvars, d){
  W <- d %>% select(all_of(Wvars))  
  miss <- data.frame(name = names(W), missing = colSums(is.na(W))/nrow(W), row.names = c(1:ncol(W)))
  for (i in 1:nrow(miss)) {
    miss$class[i] <- class(W[,which(colnames(W) == miss[i, 1])])
  }
  miss 
}

generate_miss_tbl(Wvars, dfull)

# add missingness category to IPV covariate
dfull$life_viol_any_t3<-as.factor(dfull$life_viol_any_t3)
summary(dfull$life_viol_any_t3)
dfull$life_viol_any_t3<-addNA(dfull$life_viol_any_t3)
levels(dfull$life_viol_any_t3)[length(levels(dfull$life_viol_any_t3))]<-"Missing"
summary(dfull$life_viol_any_t3)

# add missingness category to caregiver report covariates
summary(dfull$diar7d_t2)
dfull$diar7d_t2<-as.factor(dfull$diar7d_t2)
dfull$diar7d_t2<-addNA(dfull$diar7d_t2)
levels(dfull$diar7d_t2)[length(levels(dfull$diar7d_t2))]<-"Missing"
summary(dfull$diar7d_t2)

summary(dfull$ari7d_t2)
dfull$ari7d_t2<-as.factor(dfull$ari7d_t2)
dfull$ari7d_t2<-addNA(dfull$ari7d_t2)
levels(dfull$ari7d_t2)[length(levels(dfull$ari7d_t2))]<-"Missing"
summary(dfull$ari7d_t2)

summary(dfull$nose7d_t2)
dfull$nose7d_t2<-as.factor(dfull$nose7d_t2)
dfull$nose7d_t2<-addNA(dfull$nose7d_t2)
levels(dfull$nose7d_t2)[length(levels(dfull$nose7d_t2))]<-"Missing"
summary(dfull$nose7d_t2)

generate_miss_tbl(Wvars, dfull)

Xvars <- c("t2_ratio_pro_il10", "t2_ln_agp", "t2_ln_crp", "sumscore_t2_Z", "t2_ln_ifn")         
Yvars <- c("laz_t2", "waz_t2", "whz_t2" ,"hcz_t2")

for (i in Xvars){
  for (j in Yvars){
    print(i)
    print(j)
    Wvars = add_hcz(j, W2_immmune.W2_anthro)
    d_sub <- subset(dfull, !is.na(dfull[,i]) & !is.na(dfull[,j]))
    print(generate_miss_tbl(Wvars, d_sub))
  }
}

Yvars <- c("laz_t3", "waz_t3", "whz_t3", "hcz_t3")

for (i in Xvars){
  for (j in Yvars){
    print(i)
    print(j)
    Wvars = add_hcz(j, W2_immune.W3_anthro)
    d_sub <- subset(dfull, !is.na(dfull[,i]) & !is.na(dfull[,j]))
    print(generate_miss_tbl(Wvars, d_sub))
  }
}

Yvars <- c("len_velocity_t2_t3", "wei_velocity_t2_t3", "hc_velocity_t2_t3")

for (i in Xvars){
  for (j in Yvars){
    print(i)
    print(j)
    Wvars = W2_immune.W23_anthro
    d_sub <- subset(dfull, !is.na(dfull[,i]) & !is.na(dfull[,j]))
    print(generate_miss_tbl(Wvars, d_sub))
  }
}

Xvars <- c("t3_ratio_pro_il10", "sumscore_t3_Z", "t3_ln_ifn", "t3_ln_agp", "t3_ln_crp")            
Yvars <- c("laz_t3", "waz_t3", "whz_t3" ,"hcz_t3") 


for (i in Xvars){
  for (j in Yvars){
    print(i)
    print(j)
    Wvars = add_hcz(j, W3_immune.W3_anthro)
    d_sub <- subset(dfull, !is.na(dfull[,i]) & !is.na(dfull[,j]))
    print(generate_miss_tbl(Wvars, d_sub))
  }
}

# add missingness category to caregiver report covariates
summary(dfull$diar7d_t3)
dfull$diar7d_t3<-as.factor(dfull$diar7d_t3)
dfull$diar7d_t3<-addNA(dfull$diar7d_t3)
levels(dfull$diar7d_t3)[length(levels(dfull$diar7d_t3))]<-"Missing"
summary(dfull$diar7d_t3)

summary(dfull$ari7d_t3)
dfull$ari7d_t3<-as.factor(dfull$ari7d_t3)
dfull$ari7d_t3<-addNA(dfull$ari7d_t3)
levels(dfull$ari7d_t3)[length(levels(dfull$ari7d_t3))]<-"Missing"
summary(dfull$ari7d_t3)

summary(dfull$nose7d_t3)
dfull$nose7d_t3<-as.factor(dfull$nose7d_t3)
dfull$nose7d_t3<-addNA(dfull$nose7d_t3)
levels(dfull$nose7d_t3)[length(levels(dfull$nose7d_t3))]<-"Missing"
summary(dfull$nose7d_t3)

# create factor variable with missingness level for growth measurements at year 1 and year 2
growth.var <- c("laz_t1", "waz_t1", "hcz_t1", "laz_t2", "waz_t2", "hcz_t2")
for (i in growth.var) {
  cutpoints <- c(-3, -2, -1, -0)  
  cuts <- c(min(dfull[[i]], na.rm = T), cutpoints, max(dfull[[i]], na.rm = T))
  new_var <- paste(i, "_cat", sep="")
  dfull[[new_var]] <- cut(dfull[[i]], 
                          cuts,
                          right = FALSE,
                          include.lowest = TRUE)
  dfull[[new_var]] <- as.factor(dfull[[new_var]])
  dfull[[new_var]] <- fct_explicit_na(dfull[[new_var]], "Missing")
  dfull[[new_var]] <- factor(dfull[[new_var]], levels = levels(dfull[[new_var]]))
}

# check missingness of categorical growth covariates
generate_miss_tbl(paste(growth.var, "_cat", sep=""), dfull)

names(dfull)

saveRDS(dfull, paste0(dropboxDir,"Data/Cleaned/Audrie/bangladesh-immune-growth-analysis-dataset.rds"))








