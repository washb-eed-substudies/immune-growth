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
y1_id <- y1_exposure[,1] #store ids
y1_exposure <- y1_exposure[,-1] #remove ids

#check missingness
sum(is.na(y1_exposure)) #number total obs
colSums(is.na(y1_exposure)) #number by cytokine
missing <- y1_exposure %>% mutate(missing = rowSums(is.na(y1_exposure))) %>% select(missing)
length(which(missing$missing >0)) #number of children with missing data
mean(missing$missing[missing$missing >0]) #average missing per child

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
screeplot(y1_pca, npcs = 10, type = "lines")
#diagnostic plots
screeplot(y1_pca, npcs = 10, type = "lines")
biplot(y1_pca)

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

cytokine.levels.y1 <- c("agp_t2", "crp_t2","il10_t2", "il21_t2", "il17_t2",  
                     "il13_t2", "il5_t2","il4_t2","ifng_t2","il12_t2",
                     "il2_t2", "gmcsf_t2","tnfa_t2", "il6_t2", "il1_t2")
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
y2_exposure <- select(d, childid, grep("_t3", names(d), value=T))
y2_exposure <- y2_exposure[,c(1, 4:16)]

#remove observations with all missing data except childid
y2_exposure <- y2_exposure[rowSums(is.na(y2_exposure)) != ncol(y2_exposure)-1,]
nrow(y2_exposure) #remaining obs
y2_id <- y2_exposure[,1] #store ids
y2_exposure <- y2_exposure[,-1] #remove ids

#check missingness
sum(is.na(y2_exposure)) #number total obs
colSums(is.na(y2_exposure)) #number by cytokine
missing <- y2_exposure %>% mutate(missing = rowSums(is.na(y2_exposure))) %>% select(missing)
length(which(missing$missing >0)) #number of children with missing data
mean(missing$missing[missing$missing >0]) #average missing per child

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
#diagnostic plots
screeplot(y2_pca, npcs = 10, type = "lines")
biplot(y2_pca)

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

cytokine.levels.y2 <- c("il10_t3", "il21_t3", "il17_t3",  
                     "il13_t3", "il5_t3","il4_t3","ifng_t3","il12_t3",
                     "il2_t3", "gmcsf_t3","tnfa_t3", "il6_t3", "il1_t3")

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

pca.results <- merge(y1.pc.ids, y2.pc.ids, by = "childid")

#density plots
#y1
long <- pivot_longer(as.data.frame(y1.pc.ids), cols = starts_with("PC"), 
                     names_to = "PC", values_to = "value")

long$PC <- factor(long$PC, levels = c("PC1_y1","PC2_y1","PC3_y1","PC4_y1","PC5_y1",
                                      "PC6_y1","PC7_y1","PC8_y1","PC9_y1","PC10_y1"))

y1.density.plots <- ggplot(data = long) +
  geom_density(aes(x = value))+
  facet_wrap(~PC)+
  xlim(-2.5,2.5)+
  ylim(0, 3)

#y2
long <- pivot_longer(as.data.frame(y2.pc.ids), cols = starts_with("PC"), 
                     names_to = "PC", values_to = "value")

long$PC <- factor(long$PC, levels = c("PC1_y2","PC2_y2","PC3_y2","PC4_y2","PC5_y2",
                                      "PC6_y2","PC7_y2","PC8_y2","PC9_y2","PC10_y2"))

y2.density.plots <- ggplot(data = long) +
  geom_density(aes(x = value))+
  facet_wrap(~PC)+
  xlim(-2.5,2.5)+
  ylim(0, 3)

#write.csv(pca.results, file = "~/Documents/immune-growth/results/clustering pca/PCA results.csv")
