rm(list=ls())
source(here::here("0-config.R"))

d <- readRDS(paste0(dropboxDir,"Data/Cleaned/Audrie/bangladesh-immune-growth-analysis-dataset.rds"))

y1_exposure <- select(d, grep("t2_ln", names(d), value=T), grep("t2_ratio", names(d), value=T))
y1_exposure <- y1_exposure[,-c(39:44)]
  
y2_exposure <- select(d, grep("t3_ln", names(d), value=T), grep("t3_ratio", names(d), value=T))
y1_exposure <- y1_exposure[,-c(37:42)]

sapply(y1_exposure, sd, na.rm=T)
summary(sapply(y1_exposure, is.na))

y1_exposure_complete <- y1_exposure[complete.cases(y1_exposure),]
y1_pca <- prcomp(y1_exposure_complete)
y1_pca$rotation[,1]
summary(y1_pca)

library(corrplot)
corrplot(cor(y1_exposure_complete), method="ellipse")
