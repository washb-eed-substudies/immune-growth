


rm(list=ls())

source(here::here("0-config.R"))
source(here::here("src/0-gam-functions.R"))


d<-readRDS(paste0(dropboxDir,"Data/Cleaned/Audrie/bangladesh-immune-analysis-dataset.rds"))


saveRDS(dfull, paste0(dropboxDir,"Data/Cleaned/Andrew/immune_growth_data.RDS"))


