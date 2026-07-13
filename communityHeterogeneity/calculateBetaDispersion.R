library(vegan)


writeDist <- function(x, file = "", format = "phylip", ...) { ## This function extracts a file to be used as an input to the HOMOVA test
  format <- match.arg(format, c("phylip", "nexus"))
  if (format == "phylip") {
    x <- as.matrix(x)
    # maybe x <- format(x, digits = digits, justify = "none")
    cat(ncol(x), "\n", file = file)
    write.table(x, file, append = TRUE, quote = FALSE, col.names = FALSE)
  }
  else write.nexus.dist(x, file = file, ...)
}


for (org in c('Bacteria','Fungi')){
    otu<-read.csv(paste0(org,'.OTU_table.csv'),row.names=1) # OTU table including all samples across the three ecosystem types
    env<-read.csv(paste0(org,'.design'),row.names=1,sep='\t') # TSV file linking Sample name to ecosystem types
    env<-env$treatment
    dist_obj <- vegdist(otu, method = "bray") # Bray-Curtis distance matrix between samples
    dist_mat <- as.matrix(dist_obj)
    rownames(dist_mat) <- colnames(dist_mat) <- rownames(otu)
    #writeDist(dist_mat, file = paste0(org,'.brayCurtisMatrix.phy'), format = "phylip") # input for HOMOVA test

    pcoa_res <- wcmdscale(dist_obj, eig = TRUE, k = 2) # Calculate a PCoA based on the Bray-Curtis distances
    groups <- factor(env)
    mod <- betadisper(dist_obj, groups) # Calculate beta dispersion
    dists <- data.frame(mod$distances) # Extract individual distances of samples to the centroid of their ecosystem type

    # Export
    write.csv(dists, paste0(org, ".betadispersion.csv"))
}