rm(list=ls())
source(here::here("0-config.R"))
library(viridis)
library(corrplot)
library(reshape2)
library(caret)
library(RANN)

#read in data
d <- readRDS(paste0(dropboxDir,"Data/Cleaned/Audrie/bangladesh-immune-growth-analysis-dataset.rds"))

#####  Year 1
#select immune variables
y1_exposure <- select(d, childid, grep("t2_ln", names(d), value=T))
y1_exposure <- y1_exposure[,c(1, 3:17)]

#remove observations with all missing data except childid
y1_exposure <- y1_exposure[rowSums(is.na(y1_exposure)) != ncol(y1_exposure)-1,]
nrow(y1_exposure) #remaining obs
y1_id <- y1_exposure[,1] #store ids
y1_exposure <- y1_exposure[,-1] #remove ids from PCA dataset

#check correlation prior to imputation
cormat <- round(cor(y1_exposure, use="pairwise.complete.obs"),2)
cormat[lower.tri(cormat)]<- NA
melted_cormat <- melt(cormat, na.rm = TRUE)
y1.corr <- ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
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

#check missingness
sum(is.na(y1_exposure)) #number total obs
colSums(is.na(y1_exposure)) #number by cytokine
missing <- y1_exposure %>% mutate(missing = rowSums(is.na(y1_exposure))) %>% select(missing)
length(which(missing$missing >0)) #number of children with missing data
mean(missing$missing[missing$missing >0]) #average missing per child

#save means and SD for back transformation
y1.means <- sapply(y1_exposure, mean, na.rm = TRUE)
y1.sd <- sapply(y1_exposure, sd, na.rm = TRUE)

#impute using kNN
set.seed(12345)
missing.model = preProcess(y1_exposure, "knnImpute")
y1_impute.Z = predict(missing.model, y1_exposure); sum(is.na(y1_impute.Z))

#check distributions
densityplot(y1_exposure$t2_ln_tnf)
densityplot(y1_impute.Z$t2_ln_tnf)
  
#run pca
y1_pca <- prcomp(y1_impute.Z)
summary(y1_pca) 
#diagnostic plots
screeplot(y1_pca, npcs = 10, type = "lines")
biplot(y1_pca)

#backtransform imputed variables from Z-scores
y1_impute <- sweep(y1_impute.Z, MARGIN=2, y1.sd, `*`)
y1_impute <- sweep(y1_impute, MARGIN=2, y1.means, `+`)

#bind back to ids
pred <- predict(y1_pca, newdata=y1_impute.Z)
colnames(pred) <- paste(colnames(pred), "t2", sep = "_")
colnames(y1_impute) <- paste(colnames(y1_impute), "imp", sep = "_")
colnames(y1_impute.Z) <- paste(colnames(y1_impute.Z), "imp_z", sep = "_")
y1.pc.ids <- as.data.frame(cbind(childid = y1_id, y1_impute, y1_impute.Z, pred[,1:10]))

#create sum score
y1.pc.ids$sumscore_t2 <- rowSums(y1.pc.ids[,c(17:31)])
y1.pc.ids$sumscore_t2_Z <- scale(y1.pc.ids$sumscore_t2, center = TRUE, scale = TRUE)

#plot sum score
ggplot(data = y1.pc.ids) +
  geom_density(aes(x = sumscore_t2))+
  labs(x = "Sum Score (Year 1)")

ggplot(data = y1.pc.ids, aes(x = sumscore_t2, y = PC1_t2)) +
  geom_point()+
  labs(x = "Sum Score (Year 1)", y = "PC1 (Year 1)")+
  ylim(-9,6.25)

#visualize PC loadings
y1_loadings <- as.matrix(y1_pca$rotation[,c(1:10)])
y1_loadings_long <- as.data.frame(cbind(cytokine = rownames(y1_loadings), y1_loadings[,1:10]))
y1_loadings_long <- pivot_longer(y1_loadings_long, cols = starts_with("PC"), 
                              names_to = "PC", values_to = "value")

y1_loadings_long$value <- as.numeric(y1_loadings_long$value)
y1_loadings_long$PC <- factor(y1_loadings_long$PC, levels = c("PC1","PC2","PC3","PC4","PC5",
                                                        "PC6","PC7","PC8","PC9","PC10"))

cytokine.levels.y1 <- c("t2_ln_agp", "t2_ln_crp","t2_ln_il10", "t2_ln_il21", "t2_ln_il17",  
                     "t2_ln_il13", "t2_ln_il5","t2_ln_il4","t2_ln_ifn","t2_ln_il12",
                     "t2_ln_il2", "t2_ln_gmc","t2_ln_tnf", "t2_ln_il6", "t2_ln_il1")
y1_loadings_long$cytokine <- factor(y1_loadings_long$cytokine, 
                                 levels = cytokine.levels.y1,
                                 labels = c("AGP","CRP","IL-10", "IL-21","IL-17",   
                                            "IL-13", "IL-5","IL-4","IFNg","IL-12",
                                            "IL-2","GMCSF","TNFa","IL-6","IL-1"))

y1.loadings.plot <- ggplot(data = y1_loadings_long) + 
  geom_tile(aes(x = PC, y = cytokine, fill = value)) +
  theme(legend.position = "right", axis.title.y = element_blank()) +
  scale_fill_viridis(name = "Value") +
  labs(x = "Principal component", subtitle = "Year 1 PCA loadings")

#----------Year 2-------------#
#select PCA variables
y2_exposure <- select(d, childid, grep("t3_ln", names(d), value=T))
y2_exposure <- y2_exposure[,c(1, 3:15)]

#remove observations with all missing data except childid
y2_exposure <- y2_exposure[rowSums(is.na(y2_exposure)) != ncol(y2_exposure)-1,]
nrow(y2_exposure) #remaining obs
y2_id <- y2_exposure[,1] #store ids
y2_exposure <- y2_exposure[,-1] #remove ids from PCA dataset

#check correlation prior to imputation
cormat <- round(cor(y2_exposure, use="pairwise.complete.obs"),2)
cormat[lower.tri(cormat)]<- NA
melted_cormat <- melt(cormat, na.rm = TRUE)
y2.corr <- ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
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

#check missingness
sum(is.na(y2_exposure)) #number total obs
colSums(is.na(y2_exposure)) #number by cytokine
missing <- y2_exposure %>% mutate(missing = rowSums(is.na(y2_exposure))) %>% select(missing)
length(which(missing$missing >0)) #number of children with missing data
mean(missing$missing[missing$missing >0]) #average missing per child

#save means and SD for back transformation
y2.means <- sapply(y2_exposure, mean, na.rm = TRUE)
y2.sd <- sapply(y2_exposure, sd, na.rm = TRUE)

#impute using kNN
set.seed(12345)
missing.model = preProcess(y2_exposure, "knnImpute")
y2_impute.Z = predict(missing.model, y2_exposure); sum(is.na(y2_impute.Z))

#check distributions
densityplot(y2_exposure$t3_ln_tnf)
densityplot(y2_impute.Z$t3_ln_tnf)

#run pca
y2_pca <- prcomp(y2_impute.Z)
summary(y2_pca) 
#diagnostic plots
screeplot(y2_pca, npcs = 10, type = "lines")
biplot(y2_pca)

#backtransform imputed variables from Z-scores
y2_impute <- sweep(y2_impute.Z, MARGIN=2, y2.sd, `*`)
y2_impute <- sweep(y2_impute, MARGIN=2, y2.means, `+`)

#bind back to ids
pred <- predict(y2_pca, newdata=y2_impute.Z)
colnames(pred) <- paste(colnames(pred), "t3", sep = "_")
colnames(y2_impute) <- paste(colnames(y2_impute), "imp", sep = "_")
colnames(y2_impute.Z) <- paste(colnames(y2_impute.Z), "imp_z", sep = "_")
y2.pc.ids <- as.data.frame(cbind(childid = y2_id, y2_impute, y2_impute.Z, pred[,1:10]))

#create sum score
y2.pc.ids$sumscore_t3 <- rowSums(y2.pc.ids[,c(17:31)])
y2.pc.ids$sumscore_t3_Z <- scale(y2.pc.ids$sumscore_t3, center = TRUE, scale = TRUE)

#plot sum score
ggplot(data = y2.pc.ids) +
  geom_density(aes(x = sumscore_t3))+
  labs(x = "Sum Score (Year 1)")

ggplot(data = y2.pc.ids, aes(x = sumscore_t3_Z, y = PC1_t3)) +
  geom_point()+
  labs(x = "Sum Score (Year 1)", y = "PC1 (Year 1)")+
  ylim(-9,6.25)

#visualize PC loadings
y2_loadings <- as.matrix(y2_pca$rotation[,c(1:10)])
y2_loadings_long <- as.data.frame(cbind(cytokine = rownames(y2_loadings), y2_loadings[,1:10]))
y2_loadings_long <- pivot_longer(y2_loadings_long, cols = starts_with("PC"), 
                                 names_to = "PC", values_to = "value")

y2_loadings_long$value <- as.numeric(y2_loadings_long$value)
y2_loadings_long$PC <- factor(y2_loadings_long$PC, levels = c("PC1","PC2","PC3","PC4","PC5",
                                                              "PC6","PC7","PC8","PC9","PC10"))

cytokine.levels.y2 <- c("t3_ln_il10", "t3_ln_il21", "t3_ln_il17",  
                        "t3_ln_il13", "t3_ln_il5","t3_ln_il4","t3_ln_ifn","t3_ln_il12",
                        "t3_ln_il2", "t3_ln_gmc","t3_ln_tnf", "t3_ln_il6", "t3_ln_il1")
y2_loadings_long$cytokine <- factor(y2_loadings_long$cytokine, 
                                    levels = cytokine.levels.y2,
                                    labels = c("IL-10", "IL-21","IL-17",   
                                               "IL-13", "IL-5","IL-4","IFNg","IL-12",
                                               "IL-2","GMCSF","TNFa","IL-6","IL-1"))

y2.loadings.plot <- ggplot(data = y2_loadings_long) + 
  geom_tile(aes(x = PC, y = cytokine, fill = value)) +
  theme(legend.position = "right", axis.title.y = element_blank()) +
  scale_fill_viridis(name = "Value") +
  labs(x = "Principal component", subtitle = "Year 2 PCA loadings")

#export results
pca.results <- merge(y1.pc.ids, y2.pc.ids, all = TRUE)
write.csv(pca.results,
          file = "/Volumes/External 1TB Drive 2/immune-growth/results/clustering pca/PCA results.csv")

#density plots of distribution of principal components
#y1
long <- pivot_longer(as.data.frame(y1.pc.ids), cols = starts_with("PC"), 
                     names_to = "PC", values_to = "value")

long$PC <- factor(long$PC, levels = c("PC1_y1","PC2_y1","PC3_y1","PC4_y1","PC5_y1",
                                      "PC6_y1","PC7_y1","PC8_y1","PC9_y1","PC10_y1"))

y1.score.density.plots <- ggplot(data = long) +
  geom_density(aes(x = value))+
  facet_wrap(~PC)

#y2
long <- pivot_longer(as.data.frame(y2.pc.ids), cols = starts_with("PC"), 
                     names_to = "PC", values_to = "value")

long$PC <- factor(long$PC, levels = c("PC1_y2","PC2_y2","PC3_y2","PC4_y2","PC5_y2",
                                      "PC6_y2","PC7_y2","PC8_y2","PC9_y2","PC10_y2"))

y2.score.density.plots <- ggplot(data = long) +
  geom_density(aes(x = value))+
  facet_wrap(~PC)

####### sensitivity analysis - PCA with complete cases
#year 1
y1.complete <- y1_exposure[complete.cases(y1_exposure),]
nrow(y1.complete)

#scale and perform PCA
y1.complete_scale <- scale(y1.complete, center = FALSE)
y1.complete.pca <- prcomp(y1.complete_scale)
y1.complete.loadings <- as.matrix(y1.complete.pca$rotation[,c(1:10)])

loadings.diff.y1 = y1.complete.loadings - y1_loadings
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
