# Load required package
library(vegan)

# Read data
data <- read.csv("networkCharacteristics.csv", row.names = 'Network') # Obtained from running getNetworkProperties_50NetworksPerEcosystem.py or getNetworkProperties_nullModelNetworks.py
data<-data[data$"Environment" != "Null", ] # Null models to be excluded from this analysis

metadata<-read.csv('allSubsampleMetadata.csv',row.names='X')
metadata<-metadata[metadata$"Environment" != "Null", ] # Null models to be excluded from this analysis


for (org in c("Bacteria","Fungi")){
     # 1. Keep only Bacteria/Fungi at once
     print(paste0('========',org,'========'))
     data_org<-data[data$"Organism" == org, ]
     metadata_org<-metadata[grep(paste0("^",substr(org, 1, 1)), rownames(metadata)), ]


     # 2. Separate metadata (e.g., environment) from numeric variables
     # Adjust column name if needed
     env <- data_org$Environment
     vars <- data_org[, names(data_org) %in% c("Node.number","Edge.number","Mean.interaction.coefficient","Degree.heterogeneity","Network.density","Clustering.coefficient","Modularity..Q.")]

     # 3. Standardize variables (z-score)
     vars_scaled <- scale(vars)
     #vars_numeric <- as.data.frame(apply(vars_scaled, 2, as.numeric))

     # 4. Compute Euclidean distance matrix on network characteristics
     dist_matrix <- vegdist(vars_scaled, method = "euclidean")
     adonisresult <- adonis2(dist_matrix ~ Environment, data = data_org)
     #print(adonisresult)

     # 5. Perform PCoA (classical multidimensional scaling)
     pcoa <- cmdscale(dist_matrix, eig = TRUE, k = 2)
     eig_vals <- pcoa$eig
     pct_var <- 100 * eig_vals / sum(eig_vals) # get percentage of variance:
     pcoa_points <- as.data.frame(pcoa$points) # Extract coordinates
     colnames(pcoa_points) <- c("PCoA1", "PCoA2")
     pcoa_points$environment <- env # Add environment back for plotting
     write.csv(pcoa_points,paste0('pcoaNetworks.',org,'.coords.csv')) # export coordinates for plotting

     # 6. Plot PCoA
     pdf(paste0('pcoaNetworks.',org,'.pcoa.pdf'))
     plot(pcoa_points$PCoA1, pcoa_points$PCoA2,
          col = as.factor(pcoa_points$environment),
          pch = 19,
          cex = 2,
          xlab = "PCoA1",
          ylab = "PCoA2")
          legend("topright", legend = levels(as.factor(env)),
          col = 1:length(unique(env)), pch = 19)
     dev.off()

     # 7. Variance partitioning
     vp<-varpart(dist_matrix,
        ~ betadispersion+alpha_div,# microbial diversity/community properties
        ~ pH_H2O+pH_CaCl2+EC+OC+CaCO3+P+N+K,# edaphic properties
        ~ monthly_precipitation+monthly_air_temperature,# climate
        ~ centroid_x+centroid_y+convex_hull_area+mean_pairwise_distance,# geography
        data = metadata_org)
     print(str(summary(vp)))
     pdf(paste('pcoaNetworks',org,'venn','pdf',sep='.'))
     plot(vp)
     dev.off()

     # Extract the Individual fractions table as a data frame
     ind_fractions <- as.data.frame(vp$part$indfract)
     # Save as CSV
     write.csv(ind_fractions, file = paste('pcoaNetworks',org,'individual_fractions','csv',sep='.'))


     # 8. Fitting Environmental Vectors
     fit <- envfit(pcoa ~ betadispersion+alpha_div+pH_H2O+pH_CaCl2+EC+OC+CaCO3+P+N+K+monthly_precipitation+monthly_air_temperature+centroid_x+centroid_y+convex_hull_area+mean_pairwise_distance,
              data = metadata_org,
              permutations = 999)
     pdf(paste0('pcoaNetworks.',org,'.pcoaWithEnvFit.pdf'))
     plot(pcoa_points$PCoA1, pcoa_points$PCoA2,
          col = as.factor(pcoa_points$environment),
          cex = 2,
          pch = 19,
          xlab = "PCoA1",
          ylab = "PCoA2")
          legend("topright", legend = levels(as.factor(env)),
          col = 1:length(unique(env)), pch = 19)
     plot(fit)
     dev.off()
     

     #9. plot ordination ellipses
     env_factor <- as.factor(pcoa_points$environment)
     unique_env <- levels(env_factor)
     centroids <- data.frame(
          PCoA1 = tapply(pcoa_points$PCoA1, env_factor, mean),
          PCoA2 = tapply(pcoa_points$PCoA2, env_factor, mean)
     )
     rownames(centroids) <- unique_env
     # Empty plot with same limits
     pdf(paste0('pcoaNetworks.',org,'.pcoaWithEnvFit_ellipses_sd.pdf'))
     plot(pcoa_points$PCoA1, pcoa_points$PCoA2, type = "n",
          xlab = "PCoA1", ylab = "PCoA2",
          xlim = range(pcoa_points$PCoA1), ylim = range(pcoa_points$PCoA2))

     # Add ellipses (95% CI; adjust conf=0.95)
     cols <- 1:length(unique_env)
     ordiellipse(pcoa, env_factor, conf = 0.95, col = cols, draw = "polygon", kind = "sd")  # or "ehull"

     # Add centroid points (large, matching legend)
     points(centroids$PCoA1, centroids$PCoA2, pch = 19, cex = 3, col = cols)
     plot(fit)
     # Legend
     legend("topright", legend = unique_env, col = cols, pch = 19, pt.cex = 1.5)
     dev.off

}
