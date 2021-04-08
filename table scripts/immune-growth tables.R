rm(list=ls())

source(here::here("0-config.R"))
source(here::here("table-functions.R"))

# load enrollment characteristics and results
d <- readRDS(paste0(dropboxDir,"Data/Cleaned/Audrie/bangladesh-immune-growth-analysis-dataset.rds"))
H1 <- readRDS(here('results/unadjusted/H1_res.RDS'))
H2 <- readRDS(here('results/unadjusted/H2_res.RDS'))
H3 <- readRDS(here('results/unadjusted/H3_res.RDS'))
delta_growth <- readRDS(here('results/unadjusted/delta_growth_res.RDS'))
H1adj <- readRDS(here('results/adjusted/H1_adj_nofever_res.RDS'))
H2adj <- readRDS(here('results/adjusted/H2_adj_nofever_res.RDS'))
H3adj <- readRDS(here('results/adjusted/H3_adj_nofever_res.RDS'))
delta_growth_adj <- readRDS(here('results/adjusted/delta_growth_adj_nofever_res.RDS'))

full_res <- rbind(H1, H2, H3, delta_growth)
full_adj_res <- rbind(H1adj, H2adj, H3adj, delta_growth_adj)
#### MAIN TABLES ####
#### Table 1 ####
# Characteristics of participants
# nperc <- function(vector){
#   n <- sum(vector==1, na.rm=T)
#   perc <- round(n/sum(!is.na(vector))*100)
#   paste(n, " (", perc, "%)", sep="")
# }
# 
# mediqr <- function(vector){
#   quantiles <- round(quantile(vector, na.rm=T), 2)
#   paste(quantiles[3], " (", quantiles[2], ", ", quantiles[4], ")", sep="")
# }
# 
# n_med_col <- c(nperc(d$sex), mediqr(d$t2_f2_8ip), mediqr(d$t2_f2_23d), mediqr(d$t2_f2_VI), mediqr(d$t2_f2_12i),
#                mediqr(d$t3_cort_slope), mediqr(d$t3_residual_cort), mediqr(d$t3_saa_slope), mediqr(d$t3_residual_saa),
#                mediqr(d$t3_map), mediqr(d$t3_hr_mean), mediqr(d$t3_gcr_mean), mediqr(d$t3_gcr_cpg12),
#                mediqr(d$laz_t2), mediqr(d$waz_t2), mediqr(d$whz_t2), mediqr(d$hcz_t2),
#                mediqr(d$laz_t3), mediqr(d$waz_t3), mediqr(d$whz_t3), mediqr(d$hcz_t3),
#                nperc(d$diar7d_t2), nperc(d$diar7d_t3), mediqr(d$momage), mediqr(d$momheight), 
#                mediqr(d$momeduy), mediqr(d$cesd_sum_t2), mediqr(d$cesd_sum_ee_t3), mediqr(d$pss_sum_mom_t3), 
#                nperc(d$life_viol_any_t3))
# 
# 
# #immune markers, anthropometry (Y1, Y2)
# tbl1 <- data.table(" " = c("Child","","","","","","","","","","","","","","","","","","","","","","","Mother","","","","","",""),
#                    " " = c("", "Urinary F2-isoprostanes (Year 1)","","","", "Salivary cortisol reactivity (Year 2)","", "sAA reactivity (Year 2)","",
#                            "SAM biomarkers (Year 2)","", "Glucocorticoid receptor","", "Anthropometry (14 months, Year 1)","","","",
#                            "Anthropometry (28 months, Year 2)","","","", "Diarrhea (14 months, Year 1)", "Diarrhea (28 months, Year 2)","",
#                            "Anthropometry at enrollment", "Education", "Depression at Year 1", "Depression at Year 2", "Perceived stress at Year 2", 
#                            "Intimate partner violence"),
#                    " " = c("Female", "iPF(2a)-III", "2,3-dinor-iPF(2a)-III", "iPF(2a-VI", "8,12-iso-iPF(2a)-VI", 
#                            "Change in slope between pre- and post-stressor cortisol", "Cortisol residualized gain score", 
#                            "Change in slope between pre- and post-stressor sAA change", "sAA residualized gain score",
#                            "Mean arterial pressure", "Resting heart rate", "NR3C1 exon 1F promoter methylation", "NGFI-A transcription factor binding site methylation",
#                            "Length-for-age Z score", "Weight-for-age Z score", "Weight-for-length Z score", "Head circumference-for-age Z score",
#                            "Length-for-age Z score", "Weight-for-age Z score", "Weight-for-length Z score", "Head circumference-for-age Z score",
#                            "Caregiver-reported 7-day recall", "Caregiver-reported 7-day recall", "Age (years)", "Height (cm)", "Schooling completed (years)",
#                            "CES-D score", "CES-D score", "Perceived Stress Scale score", "Any lifetime exposure"),
#                    "n (%) or median (IQR)" = n_med_col)


#### Table 2 ####

exposure <- c("t2_ln_agp", "t2_ln_crp") 
outcome <- c("laz_t2", "waz_t2", "whz_t2", "hcz_t2", "laz_t3", "waz_t3", "whz_t3", "hcz_t3",
             "len_velocity_t2_t3", "wei_velocity_t2_t3", "hc_velocity_t2_t3",
             "delta_laz_t2_t3", "delta_waz_t2_t3", "delta_whz_t2_t3", "delta_hcz_t2_t3")
expo_var <- c("Ln AGP", "Ln CRP") 
out_var <- c("LAZ Year 1", "WAZ Year 1", "WLZ Year 1", "HCZ Year 1",
             "LAZ Year 2", "WAZ Year 2", "WLZ Year 2", "HCZ Year 2",
             "Length velocity Year 1 to Year 2",
             "Weight velocity (kg/month) Year 1 to Year 2",
             "Head circumference velocity (cm/month) Year 1 to Year 2",
             "Change in LAZ Year 1 to Year 2", "Change in WAZ from Year 1 to Year 2",
             "Change in WLZ from Year 1 to Year 2", "Change in HCZ from Year 1 to Year 2")

tbl2 <- growth_tbl("CRP and AGP Year 1", expo_var, out_var, exposure, outcome, full_res, full_adj_res, T)
tbl2flex <- growth_tbl_flex("CRP and AGP Year 1", expo_var, out_var, exposure, outcome, full_res, full_adj_res, T)
tbl1supp <- growth_tbl("CRP and AGP Year 1", expo_var, out_var, exposure, outcome, full_res, full_adj_res,)
tbl1flexsupp <- growth_tbl_flex("CRP and AGP Year 1", expo_var, out_var, exposure, outcome, full_res, full_adj_res,)


#### crp/agp year 2 ####
exposure <- c("t3_ln_agp", "t3_ln_crp") 
outcome <- c("laz_t3", "waz_t3", "whz_t3", "hcz_t3")
expo_var <- c("Ln AGP", "Ln CRP") 
out_var <- c("LAZ Year 2", "WAZ Year 2", "WLZ Year 2", "HCZ Year 2")

tbl3 <- growth_tbl("CRP and AGP Year 2", expo_var, out_var, exposure, outcome, full_res, full_adj_res, T)
tbl3flex <- growth_tbl_flex("CRP and AGP Year 2", expo_var, out_var, exposure, outcome, full_res, full_adj_res, T)
tbl2supp <- growth_tbl("CRP and AGP Year 2", expo_var, out_var, exposure, outcome, full_res, full_adj_res,)
tbl2flexsupp <- growth_tbl_flex("CRP and AGP Year 2", expo_var, out_var, exposure, outcome, full_res, full_adj_res,)


#### Table 4 ####

exposure <- c("t2_ratio_pro_il10", "t2_ratio_il2_il10", "t2_ratio_gmc_il10", "t2_ratio_th1_il10", "t2_ratio_th2_il10",     
              "t2_ratio_th17_il10", "t2_ratio_th1_th2", "t2_ratio_th1_th17", "t2_ln_ifn", "sumscore_t2_Z") 
outcome <- c("laz_t2", "waz_t2", "whz_t2", "hcz_t2")
expo_var <- c("Ln Pro-inflammatory cytokines/IL-10", "Ln IL-2/IL-10", "Ln GM-CSF/IL-10", "Ln Th1/IL-10", "Ln Th2/IL-10",     
              "Ln Th17/IL-10", "Ln Th1/Th2", "Ln Th1/Th17", "Ln IFN-y", "Sum score of 13 cytokines") 
out_var <- c("LAZ Year 1", "WAZ Year 1", "WLZ Year 1", "HCZ Year 1")

tbl4 <- growth_tbl("Immune Status Year 1", expo_var, out_var, exposure, outcome, H1, H1adj, T)
tbl4flex <- growth_tbl_flex("Immune Status Year 1", expo_var, out_var, exposure, outcome, H1, H1adj, T)
tbl3supp <- growth_tbl("Immune Status Year 1", expo_var, out_var, exposure, outcome, H1, H1adj)
tbl3flexsupp <- growth_tbl_flex("Immune Status Year 1", expo_var, out_var, exposure, outcome, H1, H1adj)


#### Table 5 ####

exposure <- c("t3_ratio_pro_il10", "t3_ratio_il2_il10", "t3_ratio_gmc_il10", "t3_ratio_th1_il10", "t3_ratio_th2_il10",     
              "t3_ratio_th17_il10", "t3_ratio_th1_th2", "t3_ratio_th1_th17", "t3_ln_ifn", "sumscore_t3_Z")   
outcome <- c("laz_t3", "waz_t3", "whz_t3", "hcz_t3")
expo_var <- c("Ln Pro-inflammatory cytokines/IL-10", "Ln IL-2/IL-10", "Ln GM-CSF/IL-10", "Ln Th1/IL-10", "Ln Th2/IL-10",     
              "Ln Th17/IL-10", "Ln Th1/Th2", "Ln Th1/Th17", "Ln IFN-y", "Sum score of 13 cytokines") 
out_var <- c("LAZ Year 2", "WAZ Year 2", "WLZ Year 2", "HCZ Year 2")

tbl5 <- growth_tbl("Immune Status Year 2", expo_var, out_var, exposure, outcome, H1, H1adj, T)
tbl5flex <- growth_tbl_flex("Immune Status Year 2", expo_var, out_var, exposure, outcome, H1, H1adj, T)
tbl4supp <- growth_tbl("Immune Status Year 2", expo_var, out_var, exposure, outcome, H1, H1adj)
tbl4flexsupp <- growth_tbl_flex("Immune Status Year 2", expo_var, out_var, exposure, outcome, H1, H1adj)


#### Table 6 ####

exposure <- c("t2_ratio_pro_il10", "t2_ratio_il2_il10", "t2_ratio_gmc_il10", "t2_ratio_th1_il10", "t2_ratio_th2_il10",     
              "t2_ratio_th17_il10", "t2_ratio_th1_th2", "t2_ratio_th1_th17", "t2_ln_ifn", "sumscore_t2_Z") 
outcome <- c("laz_t3", "waz_t3", "whz_t3", "hcz_t3")
expo_var <- c("Ln Pro-inflammatory cytokines/IL-10", "Ln IL-2/IL-10", "Ln GM-CSF/IL-10", "Ln Th1/IL-10", "Ln Th2/IL-10",     
              "Ln Th17/IL-10", "Ln Th1/Th2", "Ln Th1/Th17","Ln IFN-y", "Sum score of 13 cytokines") 
out_var <- c("LAZ Year 2", "WAZ Year 2", "WLZ Year 2", "HCZ Year 2")

tbl6 <- growth_tbl("Immune Status Year 1", expo_var, out_var, exposure, outcome, H2, H2adj, T)
tbl6flex <- growth_tbl_flex("Immune Status Year 1", expo_var, out_var, exposure, outcome, H2, H2adj, T)
tbl5supp <- growth_tbl("Immune Status Year 1", expo_var, out_var, exposure, outcome, H2, H2adj)
tbl5flexsupp <- growth_tbl_flex("Immune Status Year 1", expo_var, out_var, exposure, outcome, H2, H2adj)


#### Table 7 ####

exposure <- c("t2_ratio_pro_il10", "t2_ratio_il2_il10", "t2_ratio_gmc_il10", "t2_ratio_th1_il10", "t2_ratio_th2_il10",     
              "t2_ratio_th17_il10", "t2_ratio_th1_th2", "t2_ratio_th1_th17", "t2_ln_ifn", "sumscore_t2_Z") 
outcome <- c("len_velocity_t2_t3", "wei_velocity_t2_t3", "hc_velocity_t2_t3")
expo_var <- c("Ln Pro-inflammatory cytokines/IL-10", "Ln IL-2/IL-10", "Ln GM-CSF/IL-10", "Ln Th1/IL-10", "Ln Th2/IL-10",     
              "Ln Th17/IL-10", "Ln Th1/Th2", "Ln Th1/Th17", "Ln IFN-y", "Sum score of 13 cytokines") 
out_var <- c("Length velocity (cm/month) Year 1 to Year 2",
             "Weight velocity (kg/month) Year 1 to Year 2",
             "Head circumference velocity (cm/month) Year 1 to Year 2")

tbl7 <- growth_tbl("Immune Status Year 1", expo_var, out_var, exposure, outcome, H3, H3adj, T)
tbl7flex <- growth_tbl_flex("Immune Status Year 1", expo_var, out_var, exposure, outcome, H3, H3adj, T)
tbl6supp <- growth_tbl("Immune Status Year 1", expo_var, out_var, exposure, outcome, H3, H3adj)
tbl6flexsupp <- growth_tbl_flex("Immune Status Year 1", expo_var, out_var, exposure, outcome, H3, H3adj)


#### Table 8 ####

exposure <- c("t2_ratio_pro_il10", "t2_ratio_il2_il10", "t2_ratio_gmc_il10", "t2_ratio_th1_il10", "t2_ratio_th2_il10",     
              "t2_ratio_th17_il10", "t2_ratio_th1_th2", "t2_ratio_th1_th17", "t2_ln_ifn", "sumscore_t2_Z") 
outcome <- c("delta_laz_t2_t3", "delta_waz_t2_t3", "delta_whz_t2_t3", "delta_hcz_t2_t3")
expo_var <- c("Ln Pro-inflammatory cytokines/IL-10", "Ln IL-2/IL-10", "Ln GM-CSF/IL-10", "Ln Th1/IL-10", "Ln Th2/IL-10",     
              "Ln Th17/IL-10", "Ln Th1/Th2", "Ln Th1/Th17", "Ln IFN-y", "Sum score of 13 cytokines") 
out_var <- c("Change in LAZ Year 1 to Year 2", "Change in WAZ Year 1 to Year 2",
             "Change in WLZ Year 1 to Year 2", "Change in HCZ Year 1 to Year 2")

tbl8 <- growth_tbl("Immune Status Year 1", expo_var, out_var, exposure, outcome, delta_growth, delta_growth_adj, T)
tbl8flex <- growth_tbl_flex("Immune Status Year 1", expo_var, out_var, exposure, outcome, delta_growth, delta_growth_adj, T)
tbl7supp <- growth_tbl("Immune Status Year 1", expo_var, out_var, exposure, outcome, delta_growth, delta_growth_adj)
tbl7flexsupp <- growth_tbl_flex("Immune Status Year 1", expo_var, out_var, exposure, outcome, delta_growth, delta_growth_adj)


#### SAVE TABLES ####

write.csv(tbl2, here('tables/main/immune-growth-table1.csv'))
write.csv(tbl3, here('tables/main/immune-growth-table2.csv'))
write.csv(tbl4, here('tables/main/immune-growth-table3.csv'))
write.csv(tbl5, here('tables/main/immune-growth-table4.csv'))
write.csv(tbl6, here('tables/main/immune-growth-table5.csv'))
write.csv(tbl7, here('tables/main/immune-growth-table6.csv'))
write.csv(tbl8, here('tables/main/immune-growth-table8.csv'))

write.csv(tbl1supp, here('tables/supplementary/immune-growth-supptable1.csv'))
write.csv(tbl2supp, here('tables/supplementary/immune-growth-supptable2.csv'))
write.csv(tbl3supp, here('tables/supplementary/immune-growth-supptable3.csv'))
write.csv(tbl4supp, here('tables/supplementary/immune-growth-supptable4.csv'))
write.csv(tbl5supp, here('tables/supplementary/immune-growth-supptable5.csv'))
write.csv(tbl6supp, here('tables/supplementary/immune-growth-supptable6.csv'))
write.csv(tbl7supp, here('tables/supplementary/immune-growth-supptable7.csv'))

save_as_docx("Table 1: Association between CRP and AGP at Year 1 and Growth" = tbl2flex, 
             "Table 2: Association between CRP and AGP at Year 2 and Growth" = tbl3flex, 
             "Table 3: Association between Immune Status and Growth at Year 1" = tbl4flex, 
             "Table 4: Association between Immune Status and Growth at Year 2" = tbl5flex, 
             "Table 5: Association between Immune Status at Year 1 and Growth at Year 2" = tbl6flex, 
             "Table 6: Association between Immune Status and Growth Velocity Between Year 1 and Year 2" = tbl7flex, 
             "Table 7: Association between Immune Status and Change in Growth Between Year 1 and Year 2" = tbl8flex, 
             path=here('C:/Users/Sophia/Documents/WASH/WASH Immune and Growth/immune-growth main v4.docx'))

save_as_docx("Table S1: Association between CRP and AGP at Year 1 and Growth" = tbl1flexsupp, 
             "Table S2: Association between CRP and AGP at Year 1 and Growth" = tbl2flexsupp, 
             "Table S3: Association between Immune Status and Growth at Year 1" = tbl3flexsupp, 
             "Table S4: Association between Immune Status and Growth at Year 2" = tbl4flexsupp, 
             "Table S5: Association between Immune Status at Year 1 and Growth at Year 2" = tbl5flexsupp, 
             "Table S6: Association between Immune Status and Growth Velocity Between Year 1 and Year 2" = tbl6flexsupp, 
             "Table S7: Association between Immune Status and Change in Growth Between Year 1 and Year 2" = tbl7flexsupp, 
             path=here('C:/Users/Sophia/Documents/WASH/WASH Immune and Growth/immune-growth-supplementary v4.docx'))

