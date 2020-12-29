rm(list=ls())
source(here::here("0-config.R"))
library(viridis)
library(reshape2)
library(caret)
library(RANN)
library(dplyr)
library(ggplot2)
library(GGally)

#read in data
d <- readRDS(paste0(dropboxDir,"Data/Cleaned/Audrie/bangladesh-immune-growth-analysis-dataset.rds"))

#-------------------------Year 1-----------------------------#
#select immune variables
y1.var <- select(d, childid, grep("t2_ln", names(d), value=T))
y1.var <- y1.var[,c(1, 5:17)]

### deal with missingness
#remove observations with all missing data except childid
y1.var <- y1.var[rowSums(is.na(y1.var)) < 13,]
nrow(y1.var) #remaining obs

#check missingness with children that have at least one measurement
sum(is.na(y1.var)) #number total obs
colSums(is.na(y1.var)) #number by cytokine

missing <- y1.var %>% 
  mutate(missing = rowSums(is.na(y1.var)),
         flag.y1 = ifelse(missing>0, 1, 0)) #create missing flag

length(which(missing$missing >0)) #number of children with missing data
mean(missing$missing[missing$missing >0]) #average missing per child w/ missing data

#store ids and flag
y1.id <- missing[,c(1,16)]
#remove ids from PCA dataset
y1.var <- y1.var[,-1] 

#check correlation prior to imputation
corrplot <- y1.var
colnames(corrplot) <- c("GMC", "IFN-g", "IL-10", "IL-12",
                        "IL-13", "IL-17", "IL-1", "IL-2", "IL-21",
                         "IL-4", "IL-5", "IL-6", "TNF-a")
ggcorr(corrplot, label = TRUE, label_round = 2, label_size = 3)

#cormat <- round(cor(y1.var, use="pairwise.complete.obs"),2)
#cormat[lower.tri(cormat)]<- NA
#all positively correlated

#save means and SD for back transformation
y1.means <- sapply(y1.var, mean, na.rm = TRUE)
y1.sd <- sapply(y1.var, sd, na.rm = TRUE)

#impute using kNN
set.seed(12345)
missing.model = preProcess(y1.var, "knnImpute")
y1.impute.Z = predict(missing.model, y1.var); sum(is.na(y1.impute.Z))
  
### run pca
y1.pca <- prcomp(y1.impute.Z)
summary(y1.pca) 

#diagnostic plots
screeplot(y1.pca, npcs = 10, type = "lines")
biplot(y1.pca)

#predict out
y1.pred <- as.data.frame(predict(y1.pca, newdata=y1.impute.Z))
colnames(y1.pred) <- paste(colnames(y1.pred), "t2", sep = "_")

### create sum score
y1.sumscore <- y1.impute.Z %>%
  mutate(sumscore_t2_Z = scale(rowSums(y1.impute.Z), center = TRUE, scale = TRUE)) %>%
  select(sumscore_t2_Z)

#test correlation with PC1
cor.test(y1.sumscore$sumscore_t2_Z, y1.pred$PC1_t2)

### create imputed ratios
#backtransform imputed variables from Z-scores to log
y1.impute <- sweep(y1.impute.Z, MARGIN=2, y1.sd, `*`)
y1.impute <- sweep(y1.impute, MARGIN=2, y1.means, `+`)

#take antilog
y1.cluster_imp <- exp(y1.impute)
names(y1.cluster_imp) = gsub(pattern = "_ln", replacement = "", x = names(y1.cluster_imp))

#scale
y1.cluster_imp <- as.data.frame(scale(y1.cluster_imp, center = FALSE, 
                                      scale = apply(y1.cluster_imp, 2, sd, na.rm = TRUE)))

#combine ratios and take log
y1.cluster_imp <- y1.cluster_imp %>%
  mutate(pro = t2_il1 + t2_il6 + t2_tnf,
         th1 = t2_il12 + t2_ifn,
         th2 = t2_il4 + t2_il5 + t2_il13,
         th17 = t2_il17 + t2_il21) %>%
  mutate(t2_ratio_pro_il10_imp = log(pro/t2_il10),
         t2_ratio_th1_th2_imp = log(th1/th2),
         t2_ratio_th1_th17_imp = log(th1/th17),
         t2_ratio_th1_il10_imp = log(th1/t2_il10),
         t2_ratio_th2_il10_imp = log(th2/t2_il10),
         t2_ratio_th17_il10_imp = log(th17/t2_il10),
         t2_ratio_gmc_il10_imp = log(t2_gmc/t2_il10),
         t2_ratio_il2_il10_imp = log(t2_il2/t2_il10)) %>%
  select(t2_ratio_pro_il10_imp, t2_ratio_th1_th2_imp, t2_ratio_th1_th17_imp, t2_ratio_th1_il10_imp, 
         t2_ratio_th2_il10_imp, t2_ratio_th17_il10_imp, t2_ratio_gmc_il10_imp, t2_ratio_il2_il10_imp)

#bind all back to ids
y1.pc.ids <- as.data.frame(cbind(y1.id, y1.cluster_imp, y1.sumscore, y1.pred[,1:10]))

#-------------------------Year 2-----------------------------#
#select PCA variables
y2.var <- select(d, childid, grep("t3_ln", names(d), value=T))
y2.var <- y2.var[,c(1, 3:15)]

### deal with missingness
#remove observations with all missing data except childid
y2.var <- y2.var[rowSums(is.na(y2.var)) < 13,]
nrow(y2.var) #remaining obs

#check missingness with children that have at least one measurement
sum(is.na(y2.var)) #number total obs
colSums(is.na(y2.var)) #number by cytokine
missing <- y2.var %>% 
  mutate(missing = rowSums(is.na(y2.var)),
         flag.y2 = ifelse(missing>0, 1, 0))
length(which(missing$missing >0)) #number of children with missing data
mean(missing$missing[missing$missing >0]) #average missing per child

#store ids and flag
y2.id <- missing[,c(1,16)]
#remove ids from PCA dataset
y2.var <- y2.var[,-1] 

#check correlation prior to imputation
corrplot <- y2.var
colnames(corrplot) <- c("GMC", "IFN-g", "IL-10", "IL-12",
                        "IL-13", "IL-17", "IL-1", "IL-2", "IL-21",
                        "IL-4", "IL-5", "IL-6", "TNF-a")
ggcorr(corrplot, label = TRUE, label_round = 2, label_size = 3)

#cormat <- round(cor(y2.var, use="pairwise.complete.obs"),2)
#cormat[lower.tri(cormat)]<- NA
#all positively correlated

#save means and SD for back transformation
y2.means <- sapply(y2.var, mean, na.rm = TRUE)
y2.sd <- sapply(y2.var, sd, na.rm = TRUE)

#impute using kNN
set.seed(12345)
missing.model = preProcess(y2.var, "knnImpute")
y2.impute.Z = predict(missing.model, y2.var); sum(is.na(y2.impute.Z))

#run pca
y2.pca <- prcomp(y2.impute.Z)
summary(y2.pca) 

#diagnostic plots
screeplot(y2.pca, npcs = 10, type = "lines")
biplot(y2.pca)

#predict out
y2.pred <- as.data.frame(predict(y2.pca, newdata=y2.impute.Z))
colnames(y2.pred) <- paste(colnames(y2.pred), "t3", sep = "_")

### create sum score
y2.sumscore <- y2.impute.Z %>%
  mutate(sumscore_t3_Z = scale(rowSums(y2.impute.Z), center = TRUE, scale = TRUE)) %>%
  select(sumscore_t3_Z)

#test correlation with PC1
cor(y2.sumscore$sumscore_t3_Z, y2.pred$PC1_t3)

### create imputed ratios
#backtransform imputed variables from Z-scores
y2.impute <- sweep(y2.impute.Z, MARGIN=2, y2.sd, `*`)
y2.impute <- sweep(y2.impute, MARGIN=2, y2.means, `+`)

#take antilog
y2.cluster_imp <- exp(y2.impute)
names(y2.cluster_imp) = gsub(pattern = "_ln", replacement = "", x = names(y2.cluster_imp))

#scale
y2.cluster_imp <- as.data.frame(scale(y2.cluster_imp, center = FALSE, 
                                      scale = apply(y2.cluster_imp, 2, sd, na.rm = TRUE)))

#combine and take log
y2.cluster_imp <- y2.cluster_imp %>%
  mutate(pro = t3_il1 + t3_il6 + t3_tnf,
         th1 = t3_il12 + t3_ifn,
         th2 = t3_il4 + t3_il5 + t3_il13,
         th17 = t3_il17 + t3_il21) %>%
  mutate(t3_ratio_pro_il10_imp = log(pro/t3_il10),
         t3_ratio_th1_th2_imp = log(th1/th2),
         t3_ratio_th1_th17_imp = log(th1/th17),
         t3_ratio_th1_il10_imp = log(th1/t3_il10),
         t3_ratio_th2_il10_imp = log(th2/t3_il10),
         t3_ratio_th17_il10_imp = log(th17/t3_il10),
         t3_ratio_gmc_il10_imp = log(t3_gmc/t3_il10),
         t3_ratio_il2_il10_imp = log(t3_il2/t3_il10)) %>%
  select(t3_ratio_pro_il10_imp, t3_ratio_th1_th2_imp, t3_ratio_th1_th17_imp, t3_ratio_th1_il10_imp, 
         t3_ratio_th2_il10_imp, t3_ratio_th17_il10_imp, t3_ratio_gmc_il10_imp, t3_ratio_il2_il10_imp)

#bind back to ids
y2.pc.ids <- as.data.frame(cbind(y2.id, y2.cluster_imp, y2.sumscore, y2.pred[,1:10]))

##### merge year 1 and year 2
pca.results <- merge(y1.pc.ids, y2.pc.ids, all = TRUE)
pca.results$flag.y1[is.na(pca.results$flag.y1)] <- 0
pca.results$flag.y2[is.na(pca.results$flag.y2)] <- 0

#childid with no immune data (n=15)
d$childid[!(d$childid %in% pca.results$childid)]

##### export results
write.csv(pca.results,
          file = "~/Documents/immune-growth/results/clustering pca/PCA results.csv")

##### crosscheck imputed values
d2 <- readRDS(paste0(dropboxDir,"Data/Cleaned/Audrie/bangladesh-immune-growth-analysis-dataset.rds"))

Xvars <- c("t3_ratio_pro_il10", "t3_ratio_il2_il10", "t3_ratio_gmc_il10", "t3_ratio_th1_il10", "t3_ratio_th2_il10",     
           "t3_ratio_th17_il10", "t3_ratio_th1_th2", "t3_ratio_th1_th17") 

d2 <- d2 %>% select(childid, all_of(Xvars))

d3 <- merge(d2, pca.results, by = "childid", all = TRUE)

d3 <- d3 %>% 
  select(childid, flag.y1, starts_with("t3_ratio_gmc"))

#-------------------------------- Visualizations --------------------------#
#Year 1
#PC loadings
y1.loadings <- as.matrix(y1.pca$rotation[,c(1:10)])
y1.loadings_long <- as.data.frame(cbind(cytokine = rownames(y1.loadings), y1.loadings[,1:10]))

barchart <- y1.loadings_long
barchart$PC1 <- as.numeric(barchart$PC1)
y1.PC1.barchart <- ggplot (data = barchart, aes(x = cytokine, y = PC1, label=PC1)) +
  geom_col() +
  geom_text(aes(label = round(PC1, 2)), hjust = -0.1)+
  ylim(0,1)+
  coord_flip()

y1.loadings_long <- pivot_longer(y1.loadings_long, cols = starts_with("PC"), 
                                 names_to = "PC", values_to = "value")

y1.loadings_long$value <- as.numeric(y1.loadings_long$value)
y1.loadings_long$PC <- factor(y1.loadings_long$PC, levels = c("PC1","PC2","PC3","PC4","PC5",
                                                              "PC6","PC7","PC8","PC9","PC10"))

cytokine.levels.y1 <- c("t2_ln_il10", "t2_ln_il21", "t2_ln_il17",  
                        "t2_ln_il13", "t2_ln_il5","t2_ln_il4","t2_ln_ifn","t2_ln_il12",
                        "t2_ln_il2", "t2_ln_gmc","t2_ln_tnf", "t2_ln_il6", "t2_ln_il1")
y1.loadings_long$cytokine <- factor(y1.loadings_long$cytokine, 
                                    levels = cytokine.levels.y1,
                                    labels = c("IL-10", "IL-21","IL-17",   
                                               "IL-13", "IL-5","IL-4","IFNg","IL-12",
                                               "IL-2","GMCSF","TNFa","IL-6","IL-1"))

y1.loadings.plot <- ggplot(data = y1.loadings_long) + 
  geom_tile(aes(x = PC, y = cytokine, fill = value)) +
  theme(legend.position = "right", axis.title.y = element_blank()) +
  scale_fill_viridis(name = "Value") +
  labs(x = "Principal component", subtitle = "Year 1 PCA loadings")

#density plots of distribution of principal components
long <- pivot_longer(as.data.frame(y1.pc.ids[,c(1,10:19)]), cols = starts_with("PC"), 
                     names_to = "PC", values_to = "value")

long$PC <- factor(long$PC, levels = c("PC1_t2","PC2_t2","PC3_t2","PC4_t2","PC5_t2",
                                      "PC6_t2","PC7_t2","PC8_t2","PC9_t2","PC10_t2"))

ggplot(data = long) +
  geom_density(aes(x = value))+
  facet_wrap(~PC)

#plot sum score
ggplot(data = y1.pc.ids) +
  geom_density(aes(x = sumscore_t2_Z))+
  labs(x = "Sum Score (Year 1)")

ggplot(data = y1.pc.ids, aes(x = sumscore_t2_Z, y = PC1_t2)) +
  geom_point()+
  labs(x = "Sum Score (Year 1)", y = "PC1 (Year 1)")+
  ylim(-9,6.25)

#Year 2
#visualize PC loadings
y2.loadings <- as.matrix(y2.pca$rotation[,c(1:10)])
y2.loadings_long <- as.data.frame(cbind(cytokine = rownames(y2.loadings), y2.loadings[,1:10]))

barchart <- y2.loadings_long
barchart$PC1 <- as.numeric(barchart$PC1)
ggplot (data = barchart, aes(x = cytokine, y = PC1, label=PC1)) +
  geom_col() +
  geom_text(aes(label = round(PC1, 2)), hjust = -0.1)+
  ylim(0,1)+
  coord_flip()

y2.loadings_long <- pivot_longer(y2.loadings_long, cols = starts_with("PC"), 
                                 names_to = "PC", values_to = "value")

y2.loadings_long$value <- as.numeric(y2.loadings_long$value)
y2.loadings_long$PC <- factor(y2.loadings_long$PC, levels = c("PC1","PC2","PC3","PC4","PC5",
                                                              "PC6","PC7","PC8","PC9","PC10"))

cytokine.levels.y2 <- c("t3_ln_il10", "t3_ln_il21", "t3_ln_il17",  
                        "t3_ln_il13", "t3_ln_il5","t3_ln_il4","t3_ln_ifn","t3_ln_il12",
                        "t3_ln_il2", "t3_ln_gmc","t3_ln_tnf", "t3_ln_il6", "t3_ln_il1")
y2.loadings_long$cytokine <- factor(y2.loadings_long$cytokine, 
                                    levels = cytokine.levels.y2,
                                    labels = c("IL-10", "IL-21","IL-17",   
                                               "IL-13", "IL-5","IL-4","IFNg","IL-12",
                                               "IL-2","GMCSF","TNFa","IL-6","IL-1"))

y2.loadings.plot <- ggplot(data = y2.loadings_long) + 
  geom_tile(aes(x = PC, y = cytokine, fill = value)) +
  theme(legend.position = "right", axis.title.y = element_blank()) +
  scale_fill_viridis(name = "Value") +
  labs(x = "Principal component", subtitle = "Year 2 PCA loadings")

#density plots of distribution of principal components
long <- pivot_longer(as.data.frame(y2.pc.ids[,c(1,10:19)]), cols = starts_with("PC"), 
                     names_to = "PC", values_to = "value")

long$PC <- factor(long$PC, levels = c("PC1_t3","PC2_t3","PC3_t3","PC4_t3","PC5_t3",
                                      "PC6_t3","PC7_t3","PC8_t3","PC9_t3","PC10_t3"))

ggplot(data = long) +
  geom_density(aes(x = value))+
  facet_wrap(~PC)

#plot sum score
ggplot(data = y2.pc.ids) +
  geom_density(aes(x = sumscore_t3_Z))+
  labs(x = "Sum Score (Year 2)")

ggplot(data = y2.pc.ids, aes(x = sumscore_t3_Z, y = PC1_t3)) +
  geom_point()+
  labs(x = "Sum Score (Year 2)", y = "PC1 (Year 2)")+
  ylim(-9,6.25)

