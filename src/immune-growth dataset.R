rm(list=ls())

source(here::here("0-config.R"))

d <- readRDS(paste0(dropboxDir,"Data/Cleaned/Audrie/bangladesh-immune-growth-analysis-dataset.rds"))

#sum score
sum_score <- read.csv(here('results/child sum score/child immune sum scores.csv')) %>% select(-X)
sum_score_d <- left_join(d, sum_score, by='childid')

#hhwealth
d_hhwealth <- read.csv("C:/Users/Sophia/Documents/ee-secondary/sophia scripts/hhwealth.csv")
d_full <- left_join(sum_score_d, d_hhwealth, by="dataid")

saveRDS(d_full, paste0(dropboxDir,"Data/Cleaned/Audrie/bangladesh-immune-growth-analysis-dataset.rds"))
