rm(list=ls())
source(here::here("0-config.R"))
library(cowplot)
library(patchwork)
library(ggrepel)
theme_set(theme_ki())

H1adj <- readRDS(here('results/adjusted/H1_adj_nofever_res.RDS'))


plotdf <- H1adj %>% mutate(contrast=paste0(Y,"-",X)) %>%
arrange(Y, -point.diff) %>% 
  mutate(#contrast=factor(contrast, levels=unique(contrast)),
         order= row_number(),
         Pval_cat = case_when(
           corrected.Pval < 0.001 ~ "***",
           corrected.Pval < 0.01 ~ "**",
           corrected.Pval < 0.05 ~ "*",
           corrected.Pval > 0.05 ~ ""
         ))
table(plotdf$Pval_cat)

ggplot(plotdf, aes(x=order, y=point.diff)) + geom_point() + 
  geom_linerange(aes(ymin=lb.diff, ymax=ub.diff)) + 
  #geom_text_repel(aes(label=Pval_cat)) +
  geom_text(aes(label=Pval_cat), position = position_nudge(y = 0.05, x=0.25)) +
  
  facet_wrap(~Y, scales="free") +
  #scale_x_discrete(labels= X) +
  scale_x_continuous(
    breaks = plotdf$order,
    labels = plotdf$X#,
    #expand = c(0,0)
  ) +
  geom_hline(yintercept = 0) +
  theme(axis.text.x=element_text(size=8,angle=45,hjust=1, vjust=1))

