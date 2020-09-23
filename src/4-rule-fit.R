

rm(list=ls())

source(here::here("0-config.R"))
library(RCurl)
library(xrf)
#https://github.com/marjoleinF/pre
library(pre)

# grabbing data from uci
d<-readRDS(paste0(dropboxDir,"Data/Cleaned/Audrie/bangladesh-immune-analysis-dataset.rds"))


df <- data.frame(momheight=d$momheight, d[,grep('t2_ln_', names(d))])
df <- data.frame(d[,grep('t2_ln_', names(d))])
df <- df[complete.cases(df),]


# m_xrf <- xrf(t2_ln_il13 ~ ., df,
#          xgb_control = list(nrounds = 20, max_depth = 2),
#          family = 'gaussian', deoverlap = TRUE)
# coef(m_xrf, lambda = 'lambda.1se') %>%
#   filter(coefficient_lambda.1se != 0) %>%
#   arrange(-(abs(coefficient_lambda.1se)))
# 
# coef(m_xrf) %>%
#   filter(coefficient_lambda.min != 0) %>%
#   arrange(rev(abs(coefficient_lambda.min)))



airq <- airquality[complete.cases(airquality), ]
set.seed(42)
airq.ens <- pre(Ozone ~ ., data = airq)
airq.ens
plot(airq.ens, nterms = 9, cex = .5)
coefs <- coef(airq.ens)
coefs[1:6, ]

predict(airq.ens, newdata = airq[1:6, ])
set.seed(43)
airq.cv <- cvpre(airq.ens)
imps <- importance(airq.ens, round = 4)

par(mfrow = c(1, 2))
expl <- explain(airq.ens, newdata = airq[1:2, ], cex = .8)

singleplot(airq.ens, varname = "Temp")
pairplot(airq.ens, varnames = c("Temp", "Wind"))


set.seed(44)
nullmods <- bsnullinteract(airq.ens)
int <- interact(airq.ens, nullmods = nullmods)
corplot(airq.ens)
