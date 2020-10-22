rm(list=ls())
source(here::here("0-config.R"))
library(viridis)
library(corrplot)

#read in data
d <- readRDS(paste0(dropboxDir,"Data/Cleaned/Audrie/bangladesh-immune-growth-analysis-dataset.rds"))

#select PCA variables
y1_exposure <- select(d, childid, grep("_t2", names(d), value=T))
y1_exposure <- y1_exposure[,c(1, 4:18)]

#remove observations with all missing data except childid
y1_exposure <- y1_exposure[rowSums(is.na(y1_exposure)) != ncol(y1_exposure)-1,]
nrow(y1_exposure) #remaining obs
summary(sapply(y1_exposure, is.na)) #check missingness
y1_id <- y1_exposure[,1] #store ids
y1_exposure <- y1_exposure[,-1] #remove ids

#impute median
y1_impute <- y1_exposure
for(i in 1:ncol(y1_impute)){
  y1_impute[is.na(y1_impute[,i]), i] <- median(y1_impute[,i], na.rm = TRUE)
}

#check correlations
corrplot(cor(y1_impute), method="ellipse") 
#all positively correlated, would have preferred some negative correlation

#scale
y1_impute_scale <- scale(y1_impute, center = FALSE)
  
#run pca
y1_pca <- prcomp(y1_impute_scale)
summary(y1_pca) #first PC accounts for 25.2% of variance

#bind back to dataset
y1.pc.ids <- y1_pca$x[,c(1:10)]
colnames(y1.pc.ids) <- paste(colnames(y1.pc.ids), "y1", sep = "_")
y1.pc.ids <- cbind(childid = y1_id, y1.pc.ids)

#visualize
y1_loadings <- as.matrix(y1_pca$rotation[,c(1:10)])
y1_loadings_long <- as.data.frame(cbind(cytokine = rownames(y1_loadings), y1_loadings[,1:10]))
y1_loadings_long <- pivot_longer(y1_loadings_long, cols = starts_with("PC"), 
                              names_to = "PC", values_to = "value")

y1_loadings_long$value <- as.numeric(y1_loadings_long$value)
y1_loadings_long$PC <- factor(y1_loadings_long$PC, levels = c("PC1","PC2","PC3","PC4","PC5",
                                                        "PC6","PC7","PC8","PC9","PC10"))

cytokine <- unique(y1_loadings_long$cytokine)
y1_loadings_long$cytokine <- factor(y1_loadings_long$cytokine, 
                                 levels = cytokine,
                                 labels = c("CRP", "AGP", "GMCSF", "IFN-g", "IL-10", "IL-12",
                                            "IL-13", "IL-17", "IL-1", "IL-2", "IL-21", "IL-4",
                                            "IL-5", "IL-6", "TNF-a"))

y1.loadings.plot <- ggplot(data = y1_loadings_long) + 
  geom_tile(aes(x = PC, y = cytokine, fill = value)) +
  theme(legend.position = "right") +
  scale_fill_viridis() +
  labs(legend = "Loading", y = "Cytokine", x = "Principal component")

#----------Year 2-------------#
#select PCA variables
y2_exposure <- select(d, childid, grep("_t3", names(d), value=T))
y2_exposure <- y2_exposure[,c(1, 4:16)]

#remove observations with all missing data except childid
y2_exposure <- y2_exposure[rowSums(is.na(y2_exposure)) != ncol(y2_exposure)-1,]
nrow(y2_exposure) #remaining obs
summary(sapply(y2_exposure, is.na)) #check missingness
y2_id <- y2_exposure[,1] #store ids
y2_exposure <- y2_exposure[,-1] #remove ids

#impute median
y2_impute <- y2_exposure
for(i in 1:ncol(y2_impute)){
  y2_impute[is.na(y2_impute[,i]), i] <- median(y2_impute[,i], na.rm = TRUE)
}

#check correlations
corrplot(cor(y2_impute), method="ellipse") 
#all positively correlated, moreso than year 1

#scale
y2_impute_scale <- scale(y2_impute, center = FALSE)

#run pca
y2_pca <- prcomp(y2_impute_scale)
summary(y2_pca) #first PC accounts for 26.5% of variance

#bind back to dataset
y2.pc.ids <- y2_pca$x[,c(1:10)] 
colnames(y2.pc.ids) <- paste(colnames(y2.pc.ids), "y2", sep = "_")
y2.pc.ids <- cbind(childid = y2_id, y2.pc.ids)

#visualize 
y2_loadings <- as.matrix(y2_pca$rotation[,c(1:10)])
y2_loadings_long <- as.data.frame(cbind(cytokine = rownames(y2_loadings), y2_loadings[,1:10]))
y2_loadings_long <- pivot_longer(y2_loadings_long, cols = starts_with("PC"), 
                              names_to = "PC", values_to = "value")

y2_loadings_long$value <- as.numeric(y2_loadings_long$value)
y2_loadings_long$PC <- factor(y2_loadings_long$PC, levels = c("PC1","PC2","PC3","PC4","PC5",
                                                        "PC6","PC7","PC8","PC9","PC10"))

cytokine <- unique(y2_loadings_long$cytokine)
y2_loadings_long$cytokine <- factor(y2_loadings_long$cytokine, 
                                    levels = cytokine,
                                    labels = c("GMCSF", "IFN-g", "IL-10", "IL-12",
                                            "IL-13", "IL-17", "IL-1", "IL-2", "IL-21", "IL-4",
                                            "IL-5", "IL-6", "TNF-a"))

y2.loadings.plot <- ggplot(data = y2_loadings_long) + 
  geom_tile(aes(x = PC, y = cytokine, fill = value)) +
  theme(legend.position = "right") +
  scale_fill_viridis() +
  labs(legend = "Loading", y = "Cytokine", x = "Principal component")

pca.results <- merge(y1.pc.ids, y2.pc.ids, by = "childid")
write.csv(pca.results, file = "~/Documents/immune-growth/results/clustering pca/PCA results.csv")
