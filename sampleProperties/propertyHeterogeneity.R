library(vegan)
library(DescTools)
df<-read.csv('sampleProperties.csv',row.names=1) # Table containing edaphic and climatic properties of each sample (Supplementary Table 1)
scaled<-scale(df[, c("pH_H2O","pH_CaCl2","EC","OC","CaCO3","P","N","K","monthly_precipitation","monthly_air_temperature")]) # Standardize values
d <- vegdist(scaled, method = "euclidean") # Calculate Euclidean distances between samples based on standardized properties
groups<-factor(df$Environment) # Define ecosystem types as factors
mod <- betadisper(d, groups) # Calculate beta-dispersion
df$betadispersion <- mod$distances # Extract individual distances of samples to the centroid of their ecosystem type in the Euclidean distance matrix
results<-df[,c('betadispersion','Environment')]
write.csv(results, 'betadisper_soilProp.csv') # Export for plotting

K<-kruskal.test(betadispersion~Environment,data=results) # Kruskal-Wallis test for effect of ecosystem type on beta-dispersion
D<-DunnTest(betadispersion~Environment,data=results) # Post-hoc Dunn test comparing ecosystem types pairwise
print(K)
print(D)