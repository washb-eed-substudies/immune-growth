rm(list=ls())
source(here::here("0-config.R"))
library(viridis)
library(corrplot)
library(reshape2)
library(caret)
library(RANN)
library(dplyr)

#read in data
d <- readRDS(paste0(dropboxDir,"Data/Cleaned/Audrie/bangladesh-immune-growth-analysis-dataset.rds"))

#####  Year 1
#select immune variables
y1.exposure <- select(d, childid, grep("t2_ln", names(d), value=T))
y1.exposure <- y1.exposure[,c(1, 3:17)]

#remove observations with all missing data except childid
y1.exposure <- y1.exposure[rowSums(is.na(y1.exposure)) < 13,]
nrow(y1.exposure) #remaining obs

#check missingness with children that have at least one measurement
sum(is.na(y1.exposure)) #number total obs
colSums(is.na(y1.exposure)) #number by cytokine
missing <- y1.exposure %>% 
  mutate(missing = rowSums(is.na(y1.exposure)),
         flag.y1 = ifelse(missing>0, 1, 0))
length(which(missing$missing >0)) #number of children with missing data
mean(missing$missing[missing$missing >0]) #average missing per child

#store ids and flag
y1.id <- missing[,c(1,18)]
#remove ids from PCA dataset
y1.exposure <- y1.exposure[,-1] 

#check correlation prior to imputation
cormat <- round(cor(y1.exposure, use="pairwise.complete.obs"),2)
cormat[lower.tri(cormat)]<- NA
melted_cormat <- melt(cormat, na.rm = TRUE)
ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "red", high = "green", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1),
        axis.title = element_blank())+
  coord_fixed()
#all positively correlated

#save means and SD for back transformation
y1.means <- sapply(y1.exposure, mean, na.rm = TRUE)
y1.sd <- sapply(y1.exposure, sd, na.rm = TRUE)

#impute using kNN
set.seed(12345)
missing.model = preProcess(y1.exposure, "knnImpute")
y1.impute.Z = predict(missing.model, y1.exposure); sum(is.na(y1.impute.Z))
  
#run pca
y1.pca <- prcomp(y1.impute.Z)
summary(y1.pca) 

#diagnostic plots
screeplot(y1.pca, npcs = 10, type = "lines")
biplot(y1.pca)

#visualize PC loadings
y1.loadings <- as.matrix(y1.pca$rotation[,c(1:10)])
y1.loadings_long <- as.data.frame(cbind(cytokine = rownames(y1.loadings), y1.loadings[,1:10]))
y1.loadings_long <- pivot_longer(y1.loadings_long, cols = starts_with("PC"), 
                                 names_to = "PC", values_to = "value")

y1.loadings_long$value <- as.numeric(y1.loadings_long$value)
y1.loadings_long$PC <- factor(y1.loadings_long$PC, levels = c("PC1","PC2","PC3","PC4","PC5",
                                                              "PC6","PC7","PC8","PC9","PC10"))

cytokine.levels.y1 <- c("t2_ln_agp", "t2_ln_crp","t2_ln_il10", "t2_ln_il21", "t2_ln_il17",  
                        "t2_ln_il13", "t2_ln_il5","t2_ln_il4","t2_ln_ifn","t2_ln_il12",
                        "t2_ln_il2", "t2_ln_gmc","t2_ln_tnf", "t2_ln_il6", "t2_ln_il1")
y1.loadings_long$cytokine <- factor(y1.loadings_long$cytokine, 
                                    levels = cytokine.levels.y1,
                                    labels = c("AGP","CRP","IL-10", "IL-21","IL-17",   
                                               "IL-13", "IL-5","IL-4","IFNg","IL-12",
                                               "IL-2","GMCSF","TNFa","IL-6","IL-1"))

y1.loadings.plot <- ggplot(data = y1.loadings_long) + 
  geom_tile(aes(x = PC, y = cytokine, fill = value)) +
  theme(legend.position = "right", axis.title.y = element_blank()) +
  scale_fill_viridis(name = "Value") +
  labs(x = "Principal component", subtitle = "Year 1 PCA loadings")

#backtransform imputed variables from Z-scores
y1.impute <- sweep(y1.impute.Z, MARGIN=2, y1.sd, `*`)
y1.impute <- sweep(y1.impute, MARGIN=2, y1.means, `+`)

#create imputed ratios
y1.cluster_imp <- exp(y1.impute)
names(y1.cluster_imp) = gsub(pattern = "_ln", replacement = "", x = names(y1.cluster_imp))
y1.cluster_imp <- as.data.frame(scale(y1.cluster_imp, center = FALSE))
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

#create sum score
y1.sumscore <- y1.impute.Z %>%
  mutate(sumscore_t2_Z = scale(rowSums(y1.impute.Z), center = TRUE, scale = TRUE)) %>%
  select(sumscore_t2_Z)

#bind back to ids
pred <- predict(y1.pca, newdata=y1.impute.Z)
colnames(pred) <- paste(colnames(pred), "t2", sep = "_")
y1.pc.ids <- as.data.frame(cbind(y1.id, y1.cluster_imp, y1.sumscore, pred[,1:10]))

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

ggplot(data = y1.pc.ids, aes(x = sumscore_t2_Z, y = PC3_t2)) +
  geom_point()+
  labs(x = "Sum Score (Year 1)", y = "PC1 (Year 1)")+
  ylim(-9,6.25)


#-------------------------Year 2-----------------------------#
#select PCA variables
y2.exposure <- select(d, childid, grep("t3_ln", names(d), value=T))
y2.exposure <- y2.exposure[,c(1, 3:15)]

#remove observations with all missing data except childid
y2.exposure <- y2.exposure[rowSums(is.na(y2.exposure)) < 13,]
nrow(y2.exposure) #remaining obs

#check missingness with children that have at least one measurement
sum(is.na(y2.exposure)) #number total obs
colSums(is.na(y2.exposure)) #number by cytokine
missing <- y2.exposure %>% 
  mutate(missing = rowSums(is.na(y2.exposure)),
         flag.y2 = ifelse(missing>0, 1, 0))
length(which(missing$missing >0)) #number of children with missing data
mean(missing$missing[missing$missing >0]) #average missing per child

#store ids and flag
y2.id <- missing[,c(1,16)]
#remove ids from PCA dataset
y2.exposure <- y2.exposure[,-1] 

#check correlation prior to imputation
cormat <- round(cor(y2.exposure, use="pairwise.complete.obs"),2)
cormat[lower.tri(cormat)]<- NA
melted_cormat <- melt(cormat, na.rm = TRUE)
ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "red", high = "green", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1),
        axis.title = element_blank())+
  coord_fixed()
#all positively correlated

#save means and SD for back transformation
y2.means <- sapply(y2.exposure, mean, na.rm = TRUE)
y2.sd <- sapply(y2.exposure, sd, na.rm = TRUE)

#impute using kNN
set.seed(12345)
missing.model = preProcess(y2.exposure, "knnImpute")
y2.impute.Z = predict(missing.model, y2.exposure); sum(is.na(y2.impute.Z))

#run pca
y2.pca <- prcomp(y2.impute.Z)
summary(y2.pca) 

#diagnostic plots
screeplot(y2.pca, npcs = 10, type = "lines")
biplot(y2.pca)

#visualize PC loadings
y2.loadings <- as.matrix(y2.pca$rotation[,c(1:10)])
y2.loadings_long <- as.data.frame(cbind(cytokine = rownames(y2.loadings), y2.loadings[,1:10]))
y2.loadings_long <- pivot_longer(y2.loadings_long, cols = starts_with("PC"), 
                                 names_to = "PC", values_to = "value")

y2.loadings_long$value <- as.numeric(y2.loadings_long$value)
y2.loadings_long$PC <- factor(y2.loadings_long$PC, levels = c("PC1","PC2","PC3","PC4","PC5",
                                                              "PC6","PC7","PC8","PC9","PC10"))

cytokine.levels.y2 <- c("t3_ln_agp", "t3_ln_crp","t3_ln_il10", "t3_ln_il21", "t3_ln_il17",  
                        "t3_ln_il13", "t3_ln_il5","t3_ln_il4","t3_ln_ifn","t3_ln_il12",
                        "t3_ln_il2", "t3_ln_gmc","t3_ln_tnf", "t3_ln_il6", "t3_ln_il1")
y2.loadings_long$cytokine <- factor(y2.loadings_long$cytokine, 
                                    levels = cytokine.levels.y2,
                                    labels = c("AGP","CRP","IL-10", "IL-21","IL-17",   
                                               "IL-13", "IL-5","IL-4","IFNg","IL-12",
                                               "IL-2","GMCSF","TNFa","IL-6","IL-1"))

y2.loadings.plot <- ggplot(data = y2.loadings_long) + 
  geom_tile(aes(x = PC, y = cytokine, fill = value)) +
  theme(legend.position = "right", axis.title.y = element_blank()) +
  scale_fill_viridis(name = "Value") +
  labs(x = "Principal component", subtitle = "Year 1 PCA loadings")

#backtransform imputed variables from Z-scores
y2.impute <- sweep(y2.impute.Z, MARGIN=2, y2.sd, `*`)
y2.impute <- sweep(y2.impute, MARGIN=2, y2.means, `+`)

#create imputed ratios
y2.cluster_imp <- exp(y2.impute)
names(y2.cluster_imp) = gsub(pattern = "_ln", replacement = "", x = names(y2.cluster_imp))
y2.cluster_imp <- as.data.frame(scale(y2.cluster_imp, center = FALSE))
y2.cluster_imp <- y1.cluster_imp %>%
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

#create sum score
y2.sumscore <- y2.impute.Z %>%
  mutate(sumscore_t3_Z = scale(rowSums(y2.impute.Z), center = TRUE, scale = TRUE)) %>%
  select(sumscore_t3_Z)

#bind back to ids
pred <- predict(y2.pca, newdata=y2.impute.Z)
colnames(pred) <- paste(colnames(pred), "t3", sep = "_")
y2.pc.ids <- as.data.frame(cbind(y2.id, y2.cluster_imp, y2.sumscore, pred[,1:10]))

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
  labs(x = "Sum Score (Year 1)")

ggplot(data = y2.pc.ids, aes(x = sumscore_t3_Z, y = PC1_t3)) +
  geom_point()+
  labs(x = "Sum Score (Year 1)", y = "PC1 (Year 1)")+
  ylim(-9,6.25)

##### merge year 1 and year 2
pca.results <- merge(y1.pc.ids, y2.pc.ids, all = TRUE)
pca.results$flag.y1[is.na(pca.results$flag.y1)] <- 0
pca.results$flag.y2[is.na(pca.results$flag.y2)] <- 0

##### crosscheck imputed values
d2 <- readRDS(paste0(dropboxDir,"Data/Cleaned/Audrie/bangladesh-immune-growth-analysis-dataset.rds"))

Xvars <- c("t2_ratio_pro_il10", "t2_ratio_il2_il10", "t2_ratio_gmc_il10", "t2_ratio_th1_il10", "t2_ratio_th2_il10",     
           "t2_ratio_th17_il10", "t2_ratio_th1_th2", "t2_ratio_th1_th17", "t2_ln_agp", "t2_ln_crp") 

d2<-d2 %>% select(childid, Xvars)

d3 <- merge(d2, pca.results, by = "childid", all = TRUE)

d3 <- d3 %>% 
  select(childid, flag.y1, starts_with("t2_ratio_gmc"))

##### export results
write.csv(pca.results,
          file = "~/Documents/immune-growth/results/clustering pca.csv")

####### sensitivity analysis - PCA with complete cases
#year 1
y1.complete <- y1.exposure[complete.cases(y1.exposure),]
nrow(y1.complete)

#scale and perform PCA
y1.complete_scale <- scale(y1.complete, center = FALSE)
y1.complete.pca <- prcomp(y1.complete_scale)
y1.complete.loadings <- as.matrix(y1.complete.pca$rotation[,c(1:10)])

loadings.diff.y1 = y1.complete.loadings - y1.loadings
melted <- melt(loadings.diff.y1)

loaddiff.plot.y1 <- ggplot(data = melted, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Difference") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 10, hjust = 1),
        axis.title = element_blank())

### year 2
y2.complete <- y2_exposure[complete.cases(y2_exposure),]
nrow(y2.complete)

#scale and perform PCA
y2.complete_scale <- scale(y2.complete, center = FALSE)
y2.complete.pca <- prcomp(y2.complete_scale)
y2.complete.loadings <- as.matrix(y2.complete.pca$rotation[,c(1:10)])

loadings.diff.y2 = y2.complete.loadings - y2_loadings
melted <- melt(loadings.diff.y2)

loaddiff.plot.y2 <- ggplot(data = melted, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Difference") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 10, hjust = 1),
        axis.title = element_blank())
