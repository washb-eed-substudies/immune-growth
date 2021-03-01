rm(list=ls())
source(here::here("0-config.R"))
library(cowplot)
library(patchwork)
theme_set(theme_ki())

#load spline data
H1_spline <- readRDS(here("figure-data/H1_adj_nofever_spline_data.RDS"))
H2_spline <- readRDS(here("figure-data/H2_adj_nofever_spline_data.RDS"))
H3_spline <- readRDS(here("figure-data/H3_adj_nofever_spline_data.RDS"))
dgrowth_spline <- readRDS(here("figure-data/delta_growth_adj_nofever_spline_data.RDS"))

#load results for quartiles
H1_quartiles <- readRDS(here("results/adjusted/H1_adj_nofever_res.RDS"))
H2_quartiles <- readRDS(here("results/adjusted/H2_adj_nofever_res.RDS"))
H3_quartiles <- readRDS(here("results/adjusted/H3_adj_nofever_res.RDS"))
dgrowth_quartiles <- readRDS(here("results/adjusted/delta_growth_adj_nofever_res.RDS"))

d_for_plot <- function(x_name, y_name, x_var, y_var, spline, quart){
    d <- NULL
  for (i in 1:length(y_var)){
    new <- data.frame(x=x_name, y=y_name[i], spline%>%filter(Xvar==x_var, Yvar==y_var[i]), quart%>%filter(X==x_var, Y==y_var[i])%>%select(q1, q3))
    d <- rbind(d, new)
  }
  d
}

d1 <- d_for_plot("Pro-inflammatory/IL-10 \nYear 1", c("LAZ Year 1", "WAZ Year 1", "WLZ Year 1", "HCZ Year 1"), "t2_ratio_pro_il10", c("laz_t2", "waz_t2", "whz_t2", "hcz_t2"), H1_spline, H1_quartiles)
d2 <- d_for_plot("IL-2/IL-10 Year 1", c("LAZ Year 1", "WAZ Year 1", "WLZ Year 1", "HCZ Year 1"), "t2_ratio_il2_il10", c("laz_t2", "waz_t2", "whz_t2", "hcz_t2"), H1_spline, H1_quartiles)
d3 <- d_for_plot("GM-CSF/IL-10 Year 1", c("LAZ Year 1", "WAZ Year 1", "WLZ Year 1", "HCZ Year 1"), "t2_ratio_gmc_il10", c("laz_t2", "waz_t2", "whz_t2", "hcz_t2"), H1_spline, H1_quartiles)
d4 <- d_for_plot("Th1/IL-10 Year 1", c("LAZ Year 1", "WAZ Year 1", "WLZ Year 1", "HCZ Year 1"), "t2_ratio_th1_il10", c("laz_t2", "waz_t2", "whz_t2", "hcz_t2"), H1_spline, H1_quartiles)
d5 <- d_for_plot("Th2/IL-10 Year 1", c("LAZ Year 1", "WAZ Year 1", "WLZ Year 1", "HCZ Year 1"), "t2_ratio_th2_il10", c("laz_t2", "waz_t2", "whz_t2", "hcz_t2"), H1_spline, H1_quartiles)
d6 <- d_for_plot("Th17/IL-10 Year 1", c("LAZ Year 1", "WAZ Year 1", "WLZ Year 1", "HCZ Year 1"), "t2_ratio_th17_il10", c("laz_t2", "waz_t2", "whz_t2", "hcz_t2"), H1_spline, H1_quartiles)
d7 <- d_for_plot("Th1/Th2 Year 1", c("LAZ Year 1", "WAZ Year 1", "WLZ Year 1", "HCZ Year 1"), "t2_ratio_th1_th2", c("laz_t2", "waz_t2", "whz_t2", "hcz_t2"), H1_spline, H1_quartiles)
d8 <- d_for_plot("Th1/Th17 Year 1", c("LAZ Year 1", "WAZ Year 1", "WLZ Year 1", "HCZ Year 1"), "t2_ratio_th1_th17", c("laz_t2", "waz_t2", "whz_t2", "hcz_t2"), H1_spline, H1_quartiles)
d9 <- d_for_plot("AGP Year 1", c("LAZ Year 1", "WAZ Year 1", "WLZ Year 1", "HCZ Year 1"), "t2_ln_agp", c("laz_t2", "waz_t2", "whz_t2", "hcz_t2"), H1_spline, H1_quartiles)
d10<- d_for_plot("CRP Year 1", c("LAZ Year 1", "WAZ Year 1", "WLZ Year 1", "HCZ Year 1"), "t2_ln_crp", c("laz_t2", "waz_t2", "whz_t2", "hcz_t2"), H1_spline, H1_quartiles)

d11 <- d_for_plot("Pro-inflammatory/IL-10 Year 2", c("LAZ Year 1", "WAZ Year 2", "WLZ Year 2", "HCZ Year 2"), "t3_ratio_pro_il10", c("laz_t3", "waz_t3", "whz_t3", "hcz_t3"), H1_spline, H1_quartiles)
d12 <- d_for_plot("IL-2/IL-10 Year 2", c("LAZ Year 2", "WAZ Year 2", "WLZ Year 2", "HCZ Year 2"), "t3_ratio_il2_il10", c("laz_t3", "waz_t3", "whz_t3", "hcz_t3"), H1_spline, H1_quartiles)
d13 <- d_for_plot("GM-CSF/IL-10 Year 2", c("LAZ Year 2", "WAZ Year 2", "WLZ Year 2", "HCZ Year 2"), "t3_ratio_gmc_il10", c("laz_t3", "waz_t3", "whz_t3", "hcz_t3"), H1_spline, H1_quartiles)
d14 <- d_for_plot("Th1/IL-10 Year 2", c("LAZ Year 2", "WAZ Year 2", "WLZ Year 2", "HCZ Year 2"), "t3_ratio_th1_il10", c("laz_t3", "waz_t3", "whz_t3", "hcz_t3"), H1_spline, H1_quartiles)
d15 <- d_for_plot("Th2/IL-10 Year 2", c("LAZ Year 2", "WAZ Year 2", "WLZ Year 2", "HCZ Year 2"), "t3_ratio_th2_il10", c("laz_t3", "waz_t3", "whz_t3", "hcz_t3"), H1_spline, H1_quartiles)
d16 <- d_for_plot("Th17/IL-10 Year 2", c("LAZ Year 2", "WAZ Year 2", "WLZ Year 2", "HCZ Year 2"), "t3_ratio_th17_il10", c("laz_t3", "waz_t3", "whz_t3", "hcz_t3"), H1_spline, H1_quartiles)
d17 <- d_for_plot("Th1/Th2 Year 2", c("LAZ Year 2", "WAZ Year 2", "WLZ Year 2", "HCZ Year 2"), "t3_ratio_th1_th2", c("laz_t3", "waz_t3", "whz_t3", "hcz_t3"), H1_spline, H1_quartiles)
d18 <- d_for_plot("Th1/Th17 Year 2", c("LAZ Year 2", "WAZ Year 2", "WLZ Year 2", "HCZ Year 2"), "t3_ratio_th1_th17", c("laz_t3", "waz_t3", "whz_t3", "hcz_t3"), H1_spline, H1_quartiles)

d19 <- d_for_plot("Pro-inflammatory/IL-10 \nYear 1", c("LAZ Year 2", "WAZ Year 2", "WLZ Year 2", "HCZ Year 2"), "t2_ratio_pro_il10", c("laz_t3", "waz_t3", "whz_t3", "hcz_t3"), H2_spline, H2_quartiles)
d20 <- d_for_plot("IL-2/IL-10 Year 1", c("LAZ Year 2", "WAZ Year 2", "WLZ Year 2", "HCZ Year 2"), "t2_ratio_il2_il10", c("laz_t3", "waz_t3", "whz_t3", "hcz_t3"), H2_spline, H2_quartiles)
d21 <- d_for_plot("GM-CSF/IL-10 Year 1", c("LAZ Year 2", "WAZ Year 2", "WLZ Year 2", "HCZ Year 2"), "t2_ratio_gmc_il10", c("laz_t3", "waz_t3", "whz_t3", "hcz_t3"), H2_spline, H2_quartiles)
d22 <- d_for_plot("Th1/IL-10 Year 1", c("LAZ Year 2", "WAZ Year 2", "WLZ Year 2", "HCZ Year 2"), "t2_ratio_th1_il10", c("laz_t3", "waz_t3", "whz_t3", "hcz_t3"), H2_spline, H2_quartiles)
d23 <- d_for_plot("Th2/IL-10 Year 1", c("LAZ Year 2", "WAZ Year 2", "WLZ Year 2", "HCZ Year 2"), "t2_ratio_th2_il10", c("laz_t3", "waz_t3", "whz_t3", "hcz_t3"), H2_spline, H2_quartiles)
d24 <- d_for_plot("Th17/IL-10 Year 1", c("LAZ Year 2", "WAZ Year 2", "WLZ Year 2", "HCZ Year 2"), "t2_ratio_th17_il10", c("laz_t3", "waz_t3", "whz_t3", "hcz_t3"), H2_spline, H2_quartiles)
d25 <- d_for_plot("Th1/Th2 Year 1", c("LAZ Year 2", "WAZ Year 2", "WLZ Year 2", "HCZ Year 2"), "t2_ratio_th1_th2", c("laz_t3", "waz_t3", "whz_t3", "hcz_t3"), H2_spline, H2_quartiles)
d26 <- d_for_plot("Th1/Th17 Year 1", c("LAZ Year 2", "WAZ Year 2", "WLZ Year 2", "HCZ Year 2"), "t2_ratio_th1_th17", c("laz_t3", "waz_t3", "whz_t3", "hcz_t3"), H2_spline, H2_quartiles)
d27 <- d_for_plot("AGP Year 1", c("LAZ Year 2", "WAZ Year 2", "WLZ Year 2", "HCZ Year 2"), "t2_ln_agp", c("laz_t3", "waz_t3", "whz_t3", "hcz_t3"), H2_spline, H2_quartiles)
d28 <- d_for_plot("CRP Year 1", c("LAZ Year 2", "WAZ Year 2", "WLZ Year 2", "HCZ Year 2"), "t2_ln_crp", c("laz_t3", "waz_t3", "whz_t3", "hcz_t3"), H2_spline, H2_quartiles)

d29 <- d_for_plot("Pro-inflammatory/IL-10 \nYear 1", c("Length velocity", "Weight velocity", "Head \n circumference \n velocity"), "t2_ratio_pro_il10", c("len_velocity_t2_t3", "wei_velocity_t2_t3", "hc_velocity_t2_t3"), H3_spline, H3_quartiles)
d30 <- d_for_plot("IL-2/IL-10 Year 1", c("Length velocity", "Weight velocity", "Head \n circumference \n velocity"), "t2_ratio_il2_il10", c("len_velocity_t2_t3", "wei_velocity_t2_t3", "hc_velocity_t2_t3"), H3_spline, H3_quartiles)
d31 <- d_for_plot("GM-CSF/IL-10 Year 1", c("Length velocity", "Weight velocity", "Head \n circumference \n velocity"), "t2_ratio_gmc_il10", c("len_velocity_t2_t3", "wei_velocity_t2_t3", "hc_velocity_t2_t3"), H3_spline, H3_quartiles)
d32 <- d_for_plot("Th1/IL-10 Year 1", c("Length velocity", "Weight velocity", "Head \n circumference \n velocity"), "t2_ratio_th1_il10", c("len_velocity_t2_t3", "wei_velocity_t2_t3", "hc_velocity_t2_t3"), H3_spline, H3_quartiles)
d33 <- d_for_plot("Th2/IL-10 Year 1", c("Length velocity", "Weight velocity", "Head \n circumference \n velocity"), "t2_ratio_th2_il10", c("len_velocity_t2_t3", "wei_velocity_t2_t3", "hc_velocity_t2_t3"), H3_spline, H3_quartiles)
d34 <- d_for_plot("Th17/IL-10 Year 1", c("Length velocity", "Weight velocity", "Head \n circumference \n velocity"), "t2_ratio_th17_il10", c("len_velocity_t2_t3", "wei_velocity_t2_t3", "hc_velocity_t2_t3"), H3_spline, H3_quartiles)
d35 <- d_for_plot("Th1/Th2 Year 1", c("Length velocity", "Weight velocity", "Head \n circumference \n velocity"), "t2_ratio_th1_th2", c("len_velocity_t2_t3", "wei_velocity_t2_t3", "hc_velocity_t2_t3"), H3_spline, H3_quartiles)
d36 <- d_for_plot("Th1/Th17 Year 1", c("Length velocity", "Weight velocity", "Head \n circumference \n velocity"), "t2_ratio_th1_th17", c("len_velocity_t2_t3", "wei_velocity_t2_t3", "hc_velocity_t2_t3"), H3_spline, H3_quartiles)
d37 <- d_for_plot("AGP Year 1", c("Length velocity", "Weight velocity", "Head \n circumference \n velocity"), "t2_ln_agp", c("len_velocity_t2_t3", "wei_velocity_t2_t3", "hc_velocity_t2_t3"), H3_spline, H3_quartiles)
d38 <- d_for_plot("CRP Year 1", c("Length velocity", "Weight velocity", "Head \n circumference \n velocity"), "t2_ln_crp", c("len_velocity_t2_t3", "wei_velocity_t2_t3", "hc_velocity_t2_t3"), H3_spline, H3_quartiles)

d39 <- d_for_plot("Pro-inflammatory/IL-10 \nYear 1", c("Change in LAZ", "Change in WAZ", "Change in WLZ", "Change in HCZ"), "t2_ratio_pro_il10", c("delta_laz_t2_t3", "delta_waz_t2_t3", "delta_whz_t2_t3", "delta_hcz_t2_t3"), dgrowth_spline, dgrowth_quartiles)
d40 <- d_for_plot("IL-2/IL-10 Year 1", c("Change in LAZ", "Change in WAZ", "Change in WLZ", "Change in HCZ"), "t2_ratio_il2_il10", c("delta_laz_t2_t3", "delta_waz_t2_t3", "delta_whz_t2_t3", "delta_hcz_t2_t3"), dgrowth_spline, dgrowth_quartiles)
d41 <- d_for_plot("GM-CSF/IL-10 Year 1", c("Change in LAZ", "Change in WAZ", "Change in WLZ", "Change in HCZ"), "t2_ratio_gmc_il10", c("delta_laz_t2_t3", "delta_waz_t2_t3", "delta_whz_t2_t3", "delta_hcz_t2_t3"), dgrowth_spline, dgrowth_quartiles)
d42 <- d_for_plot("Th1/IL-10 Year 1", c("Change in LAZ", "Change in WAZ", "Change in WLZ", "Change in HCZ"), "t2_ratio_th1_il10", c("delta_laz_t2_t3", "delta_waz_t2_t3", "delta_whz_t2_t3", "delta_hcz_t2_t3"), dgrowth_spline, dgrowth_quartiles)
d43 <- d_for_plot("Th2/IL-10 Year 1", c("Change in LAZ", "Change in WAZ", "Change in WLZ", "Change in HCZ"), "t2_ratio_th2_il10", c("delta_laz_t2_t3", "delta_waz_t2_t3", "delta_whz_t2_t3", "delta_hcz_t2_t3"), dgrowth_spline, dgrowth_quartiles)
d44 <- d_for_plot("Th17/IL-10 Year 1", c("Change in LAZ", "Change in WAZ", "Change in WLZ", "Change in HCZ"), "t2_ratio_th17_il10", c("delta_laz_t2_t3", "delta_waz_t2_t3", "delta_whz_t2_t3", "delta_hcz_t2_t3"), dgrowth_spline, dgrowth_quartiles)
d45 <- d_for_plot("Th1/Th2 Year 1", c("Change in LAZ", "Change in WAZ", "Change in WLZ", "Change in HCZ"), "t2_ratio_th1_th2", c("delta_laz_t2_t3", "delta_waz_t2_t3", "delta_whz_t2_t3", "delta_hcz_t2_t3"), dgrowth_spline, dgrowth_quartiles)
d46 <- d_for_plot("Th1/Th17 Year 1", c("Change in LAZ", "Change in WAZ", "Change in WLZ", "Change in HCZ"), "t2_ratio_th1_th17", c("delta_laz_t2_t3", "delta_waz_t2_t3", "delta_whz_t2_t3", "delta_hcz_t2_t3"), dgrowth_spline, dgrowth_quartiles)
d47 <- d_for_plot("AGP Year 1", c("Change in LAZ", "Change in WAZ", "Change in WLZ", "Change in HCZ"), "t2_ln_agp", c("delta_laz_t2_t3", "delta_waz_t2_t3", "delta_whz_t2_t3", "delta_hcz_t2_t3"), dgrowth_spline, dgrowth_quartiles)
d48 <- d_for_plot("CRP Year 1", c("Change in LAZ", "Change in WAZ", "Change in WLZ", "Change in HCZ"), "t2_ln_crp", c("delta_laz_t2_t3", "delta_waz_t2_t3", "delta_whz_t2_t3", "delta_hcz_t2_t3"), dgrowth_spline, dgrowth_quartiles)

d1$y <- factor(d1$y)
d2$y <- factor(d2$y)
d3$y <- factor(d3$y)
d4$y <- factor(d4$y)
d5$y <- factor(d5$y)
d6$y <- factor(d6$y)
d7$y <- factor(d7$y)
d8$y <- factor(d8$y)
d9$y <- factor(d9$y)
d10$y <- factor(d10$y)
d11$y <- factor(d11$y)
d12$y <- factor(d12$y)
d13$y <- factor(d13$y)
d14$y <- factor(d14$y)
d15$y <- factor(d15$y)
d16$y <- factor(d16$y)
d17$y <- factor(d17$y)
d18$y <- factor(d18$y)
d19$y <- factor(d19$y)
d20$y <- factor(d20$y)
d21$y <- factor(d21$y)
d22$y <- factor(d22$y)
d23$y <- factor(d23$y)
d24$y <- factor(d24$y)
d25$y <- factor(d25$y)
d26$y <- factor(d26$y)
d27$y <- factor(d27$y)
d28$y <- factor(d28$y)
d29$y <- factor(d29$y)
d30$y <- factor(d30$y)
d31$y <- factor(d31$y)
d32$y <- factor(d32$y)
d33$y <- factor(d33$y)
d34$y <- factor(d34$y)
d35$y <- factor(d35$y)
d36$y <- factor(d36$y)
d37$y <- factor(d37$y)
d38$y <- factor(d38$y)
d39$y <- factor(d39$y)
d40$y <- factor(d40$y)
d41$y <- factor(d41$y)
d42$y <- factor(d42$y)
d43$y <- factor(d43$y)
d44$y <- factor(d44$y)
d45$y <- factor(d45$y)
d46$y <- factor(d46$y)
d47$y <- factor(d47$y)
d48$y <- factor(d48$y)


#spline plot function
spline_plot_functions <- function(d){
  
  color_levels = c("LAZ Year 1", "WAZ Year 1", "WLZ Year 1", "HCZ Year 1",
                   "LAZ Year 2", "WAZ Year 2", "WLZ Year 2", "HCZ Year 2",
                   "Change in LAZ", "Change in WAZ", "Change in WLZ", "Change in HCZ",
                   "Length velocity", "Weight velocity", "Head \n circumference \n velocity")
  
  nlevels <- length(levels(d$y))
  
  quantiles <- d %>% group_by(y) %>%
    summarize(
      x.lb=as.numeric(quantile(X, probs = seq(0, 1, 0.05))[2]),
      x.ub=as.numeric(quantile(X, probs = seq(0, 1, 0.05))[20]),
      y.lb=as.numeric(quantile(Y, probs = seq(0, 1, 0.05))[2]),
      y.ub=as.numeric(quantile(Y, probs = seq(0, 1, 0.05))[20])
    )
  
  d <- left_join(d, quantiles, by="y")
  p1 <- d[d$y==levels(d$y)[1],] %>% {ggplot(.,aes(x = X)) +
      geom_smooth(aes(y = fit, color=y), se = F) +
      geom_rug(aes(x=q1), sides="b", length = unit(0.15, "npc"), size=1, color="grey30") +
      geom_rug(aes(x=q3), sides="b", length = unit(0.15, "npc"), size=1, color="grey30") +
      geom_point(aes(y=Y), alpha=0.5) +
      geom_ribbon(aes(ymin=lwrS, ymax=uprS, fill=y, color=y), alpha=0.5) +
      coord_cartesian(xlim = c(.$x.lb[1], .$x.ub[1]), ylim = c(.$y.lb[1], .$y.ub[1])) +
      scale_colour_manual(values=tableau10[c(1:4,1:4,1:4,5:7)], drop=TRUE, limits=color_levels) + 
      scale_fill_manual(values=tableau10[c(1:4,1:4,1:4,5:7)], drop=TRUE, limits=color_levels) + 
      xlab(.$x[1]) + ylab(.$y[1])
  }
  if(nlevels>1){
    p2 <- d[d$y==levels(d$y)[2],] %>% {ggplot(.,aes(x = X)) +
        geom_smooth(aes(y = fit, color=y), se = F) +
        geom_ribbon(aes(ymin=lwrS, ymax=uprS, fill=y, color=y), alpha=0.5) +
        geom_rug(aes(x=q1), sides="b", length = unit(0.15, "npc"), size=1, color="grey30") +
        geom_rug(aes(x=q3), sides="b", length = unit(0.15, "npc"), size=1, color="grey30") +
        geom_point(aes(y=Y), alpha=0.5) +
        coord_cartesian(xlim = c(.$x.lb[1], .$x.ub[1]), ylim = c(.$y.lb[1], .$y.ub[1])) +
        scale_colour_manual(values=tableau10[c(1:4,1:4,1:4,5:7)], drop=TRUE, limits=color_levels) +
        scale_fill_manual(values=tableau10[c(1:4,1:4,1:4,5:7)], drop=TRUE, limits=color_levels) +
        xlab(.$x[1]) + ylab(.$y[1])
    }
  }
  if(nlevels>2){
    p3 <- d[d$y==levels(d$y)[3],] %>% {ggplot(.,aes(x = X)) +
        geom_smooth(aes(y = fit, color=y), se = F) +
        geom_rug(aes(x=q1), sides="b", length = unit(0.15, "npc"), size=1, color="grey30") +
        geom_rug(aes(x=q3), sides="b", length = unit(0.15, "npc"), size=1, color="grey30") +
        geom_ribbon(aes(ymin=lwrS, ymax=uprS, fill=y, color=y), alpha=0.5) +
        geom_point(aes(y=Y), alpha=0.5) +
        coord_cartesian(xlim = c(.$x.lb[1], .$x.ub[1]), ylim = c(.$y.lb[1], .$y.ub[1])) +
        scale_colour_manual(values=tableau10[c(1:4,1:4,1:4,5:7)], drop=TRUE, limits=color_levels) +
        scale_fill_manual(values=tableau10[c(1:4,1:4,1:4,5:7)], drop=TRUE, limits=color_levels) +
        xlab(.$x[1]) + ylab(.$y[1])
    }
    if(nlevels==4){
      p4 <- d[d$y==levels(d$y)[4],] %>% {ggplot(.,aes(x = X)) +
          geom_smooth(aes(y = fit, color=y), se = F) +
          geom_point(aes(y=Y), alpha=0.5) +
          geom_rug(aes(x=q1), sides="b", length = unit(0.15, "npc"), size=1, color="grey30") +
          geom_rug(aes(x=q3), sides="b", length = unit(0.15, "npc"), size=1, color="grey30") +
          geom_ribbon(aes(ymin=lwrS, ymax=uprS, fill=y, color=y), alpha=0.5) +
          coord_cartesian(xlim = c(.$x.lb[1], .$x.ub[1]), ylim = c(.$y.lb[1], .$y.ub[1])) +
          scale_colour_manual(values=tableau10[c(1:4,1:4,1:4,5:7)], drop=TRUE, limits=color_levels) +
          scale_fill_manual(values=tableau10[c(1:4,1:4,1:4,5:7)], drop=TRUE, limits=color_levels) +
          xlab(.$x[1]) + ylab(.$y[1])
      }
    }
  }
  
  if(nlevels==2){
    return(list(p1, p2))
  } else if (nlevels==3){
    return(list(p1, p2, p3))
  } else if (nlevels==4) {
    return (list(p1, p2, p3, p4))
  }
  return (list(p1))
}

plist1 <- spline_plot_functions(d1)
plist2 <- spline_plot_functions(d2)
plist3 <- spline_plot_functions(d3)
plist4 <- spline_plot_functions(d4)
plist5 <- spline_plot_functions(d5)
plist6 <- spline_plot_functions(d6)
plist7 <- spline_plot_functions(d7)
plist8 <- spline_plot_functions(d8)
plist9 <- spline_plot_functions(d9)
plist10 <- spline_plot_functions(d10)
plist11 <- spline_plot_functions(d11)
plist12 <- spline_plot_functions(d12)
plist13 <- spline_plot_functions(d13)
plist14 <- spline_plot_functions(d14)
plist15 <- spline_plot_functions(d15)
plist16 <- spline_plot_functions(d16)
plist17 <- spline_plot_functions(d17)
plist18 <- spline_plot_functions(d18)
plist19 <- spline_plot_functions(d19)
plist20 <- spline_plot_functions(d20)
plist21 <- spline_plot_functions(d21)
plist22 <- spline_plot_functions(d22)
plist23 <- spline_plot_functions(d23)
plist24 <- spline_plot_functions(d24)
plist25 <- spline_plot_functions(d25)
plist26 <- spline_plot_functions(d26)
plist27 <- spline_plot_functions(d27)
plist28 <- spline_plot_functions(d28)
plist29 <- spline_plot_functions(d29)
plist30 <- spline_plot_functions(d30)
plist31 <- spline_plot_functions(d31)
plist32 <- spline_plot_functions(d32)
plist33 <- spline_plot_functions(d33)
plist34 <- spline_plot_functions(d34)
plist35 <- spline_plot_functions(d35)
plist36 <- spline_plot_functions(d36)
plist37 <- spline_plot_functions(d37)
plist38 <- spline_plot_functions(d38)
plist39 <- spline_plot_functions(d39)
plist40 <- spline_plot_functions(d40)
plist41 <- spline_plot_functions(d41)
plist42 <- spline_plot_functions(d42)
plist43 <- spline_plot_functions(d43)
plist44 <- spline_plot_functions(d44)
plist45 <- spline_plot_functions(d45)
plist46 <- spline_plot_functions(d46)
plist47 <- spline_plot_functions(d47)
plist48 <- spline_plot_functions(d48)


p1 <- plot_grid(plist1[[1]], plist1[[2]], plist1[[3]], plist1[[4]], ncol=4, labels=c("",""))
p2 <- plot_grid(plist2[[1]], plist2[[2]], plist2[[3]], plist2[[4]], ncol=4, labels=c("",""))
p3 <- plot_grid(plist3[[1]], plist3[[2]], plist3[[3]], plist3[[4]], ncol=4, labels=c("",""))
p4 <- plot_grid(plist4[[1]], plist4[[2]], plist4[[3]], plist4[[4]], ncol=4, labels=c("",""))
p5 <- plot_grid(plist5[[1]], plist5[[2]], plist5[[3]], plist5[[4]], ncol=4, labels=c("",""))
p6 <- plot_grid(plist6[[1]], plist6[[2]], plist6[[3]], plist6[[4]], ncol=4, labels=c("",""))
p7 <- plot_grid(plist7[[1]], plist7[[2]], plist7[[3]], plist7[[4]], ncol=4, labels=c("",""))
p8 <- plot_grid(plist8[[1]], plist8[[2]], plist8[[3]], plist8[[4]], ncol=4, labels=c("",""))
p9 <- plot_grid(plist9[[1]], plist9[[2]], plist9[[3]], plist9[[4]], ncol=4, labels=c("",""))
p10 <- plot_grid(plist10[[1]], plist10[[2]], plist10[[3]], plist10[[4]], ncol=4, labels=c("",""))

p11 <- plot_grid(plist11[[1]], plist11[[2]], plist11[[3]], plist11[[4]], ncol=4, labels=c("",""))
p12 <- plot_grid(plist12[[1]], plist12[[2]], plist12[[3]], plist12[[4]], ncol=4, labels=c("",""))
p13 <- plot_grid(plist13[[1]], plist13[[2]], plist13[[3]], plist13[[4]], ncol=4, labels=c("",""))
p14 <- plot_grid(plist14[[1]], plist14[[2]], plist14[[3]], plist14[[4]], ncol=4, labels=c("",""))
p15 <- plot_grid(plist15[[1]], plist15[[2]], plist15[[3]], plist15[[4]], ncol=4, labels=c("",""))
p16 <- plot_grid(plist16[[1]], plist16[[2]], plist16[[3]], plist16[[4]], ncol=4, labels=c("",""))
p17 <- plot_grid(plist17[[1]], plist17[[2]], plist17[[3]], plist17[[4]], ncol=4, labels=c("",""))
p18 <- plot_grid(plist18[[1]], plist18[[2]], plist18[[3]], plist18[[4]], ncol=4, labels=c("",""))

p19 <- plot_grid(plist19[[1]], plist19[[2]], plist19[[3]], plist19[[4]], ncol=4, labels=c("",""))
p20 <- plot_grid(plist20[[1]], plist20[[2]], plist20[[3]], plist20[[4]], ncol=4, labels=c("",""))
p21 <- plot_grid(plist21[[1]], plist21[[2]], plist21[[3]], plist21[[4]], ncol=4, labels=c("",""))
p22 <- plot_grid(plist22[[1]], plist22[[2]], plist22[[3]], plist22[[4]], ncol=4, labels=c("",""))
p23 <- plot_grid(plist23[[1]], plist23[[2]], plist23[[3]], plist23[[4]], ncol=4, labels=c("",""))
p24 <- plot_grid(plist24[[1]], plist24[[2]], plist24[[3]], plist24[[4]], ncol=4, labels=c("",""))
p25 <- plot_grid(plist25[[1]], plist25[[2]], plist25[[3]], plist25[[4]], ncol=4, labels=c("",""))
p26 <- plot_grid(plist26[[1]], plist26[[2]], plist26[[3]], plist26[[4]], ncol=4, labels=c("",""))
p27 <- plot_grid(plist27[[1]], plist27[[2]], plist27[[3]], plist27[[4]], ncol=4, labels=c("",""))
p28 <- plot_grid(plist28[[1]], plist28[[2]], plist28[[3]], plist28[[4]], ncol=4, labels=c("",""))

p29 <- plot_grid(plist29[[1]], plist29[[2]], plist29[[3]], ncol=3, labels=c("",""))
p30 <- plot_grid(plist30[[1]], plist30[[2]], plist30[[3]], ncol=3, labels=c("",""))
p31 <- plot_grid(plist31[[1]], plist31[[2]], plist31[[3]], ncol=3, labels=c("",""))
p32 <- plot_grid(plist32[[1]], plist32[[2]], plist32[[3]], ncol=3, labels=c("",""))
p33 <- plot_grid(plist33[[1]], plist33[[2]], plist33[[3]], ncol=3, labels=c("",""))
p34 <- plot_grid(plist34[[1]], plist34[[2]], plist34[[3]], ncol=3, labels=c("",""))
p35 <- plot_grid(plist35[[1]], plist35[[2]], plist35[[3]], ncol=3, labels=c("",""))
p36 <- plot_grid(plist36[[1]], plist36[[2]], plist36[[3]], ncol=3, labels=c("",""))
p37 <- plot_grid(plist37[[1]], plist37[[2]], plist37[[3]], ncol=3, labels=c("",""))
p38 <- plot_grid(plist38[[1]], plist38[[2]], plist38[[3]], ncol=3, labels=c("",""))

p39 <- plot_grid(plist39[[1]], plist39[[2]], plist39[[3]], plist39[[4]], ncol=4, labels=c("",""))
p40 <- plot_grid(plist40[[1]], plist40[[2]], plist40[[3]], plist40[[4]], ncol=4, labels=c("",""))
p41 <- plot_grid(plist41[[1]], plist41[[2]], plist41[[3]], plist41[[4]], ncol=4, labels=c("",""))
p42 <- plot_grid(plist42[[1]], plist42[[2]], plist42[[3]], plist42[[4]], ncol=4, labels=c("",""))
p43 <- plot_grid(plist43[[1]], plist43[[2]], plist43[[3]], plist43[[4]], ncol=4, labels=c("",""))
p44 <- plot_grid(plist44[[1]], plist44[[2]], plist44[[3]], plist44[[4]], ncol=4, labels=c("",""))
p45 <- plot_grid(plist45[[1]], plist45[[2]], plist45[[3]], plist45[[4]], ncol=4, labels=c("",""))
p46 <- plot_grid(plist46[[1]], plist46[[2]], plist46[[3]], plist46[[4]], ncol=4, labels=c("",""))
p47 <- plot_grid(plist47[[1]], plist47[[2]], plist47[[3]], plist47[[4]], ncol=4, labels=c("",""))
p48 <- plot_grid(plist48[[1]], plist48[[2]], plist48[[3]], plist48[[4]], ncol=4, labels=c("",""))

#adjusted differences for analyses of concurrent growth at year 1
pcomb1 <- plot_grid(p1,p2,p3,p4,p5,p6,p7,p8,
                    nrow=8,
                    labels = c("","",""),
                    hjust=0.5, vjust=0.5,
                    rel_heights = c(1, 1, 1))

#adjusted differences for analyses of concurrent growth at year 2
pcomb2 <- plot_grid(p11,p12,p13,p14,p15,p16,p17,p18,
                    nrow=8,
                    labels = c("","","",""),
                    hjust=0.5,vjust=0.5,
                    rel_heights = c(1, 1, 1, 1))
#adjusted differences for immmune status and subsequent growth
pcomb3 <- plot_grid(p19,p20,p21,p22,p23,p24,p25,p26,
                    nrow=8,
                    labels = c("",""),
                    hjust=0.5,vjust=0.5,
                    rel_heights = c(1, 1, 1, 1))
#adjusted differences for immune status and growth velocity
pcomb4 <- plot_grid(p29,p30,p31,p32,p33,p34,p35,p36,
                    nrow=8,
                    labels = c("",""),
                    hjust=0.5,vjust=0.5,
                    rel_heights = c(1, 1, 1, 1))
pcomb5 <- plot_grid(p39,p40,p41,p42,p43,p44,p45,p46,
                    nrow=8,
                    labels = c("",""),
                    hjust=0.5,vjust=0.5,
                    rel_heights = c(1, 1, 1, 1))

# agp only
agp <- plot_grid(p9,p27,p47,p37,
                    nrow=4,
                    labels = c("",""),
                    hjust=0.5,vjust=0.5,
                    rel_heights = c(1, 1, 1, 1))
# crp only
crp <- plot_grid(p10,p28,p48,p38,
                    nrow=4,
                    labels = c("",""),
                    hjust=0.5,vjust=0.5,
                    rel_heights = c(1, 1, 1, 1))

ggsave(pcomb1, file = here("figures/concurrent-growth-y1.tiff"), height=16, width=10)
ggsave(pcomb2, file = here("figures/concurrent-growth-y2.tiff"), height=16, width=10)
ggsave(pcomb3, file = here("figures/subsequent-growth.tiff"), height=16, width=10)
ggsave(pcomb4, file = here("figures/growth-velocity.tiff"), height=16, width=10)
ggsave(pcomb5, file = here("figures/change-in-growth.tiff"), height=16, width=10)

ggsave(agp, file=here("figures/agp.tiff"), height=14, width=10)
ggsave(crp, file=here("figures/crp.tiff"), height=14, width=10)

# ggplot(d,aes(x = X)) +
#   geom_smooth(aes(y = fit, color=y), se = F) +
#   geom_ribbon(aes(ymin=lwrS, ymax=uprS, fill=y, color=y), alpha=0.5) +
#   geom_rug(aes(y=Y)) +
#   facet_wrap(~y)


