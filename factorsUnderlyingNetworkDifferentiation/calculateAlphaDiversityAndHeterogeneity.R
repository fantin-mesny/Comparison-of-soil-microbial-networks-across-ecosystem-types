library(vegan)

## Parse the complete dataset
otu_bac<- read.csv('Bacteria.OTU_table.csv',row.names=1) # OTU table including all samples across the three ecosystem types
otu_bac<- t(otu_bac)
otu_fun<- read.csv('Fungi.OTU_table.csv',row.names=1) # OTU table including all samples across the three ecosystem types
otu_fun<-t(otu_fun)


for (org in c('Bacteria','Fungi')){
    print(org)
    ifelse(org=='Bacteria', otu0<-otu_bac, otu0<-otu_fun)

    bc_dist <- vegdist(as.matrix(otu0), method = "bray") ## Calculate a bray-curtis distance matrix

    for (env in c('Cropland','Grassland','Woodland')){
        for (rep in c(0:49)){
            otu<-read.csv(paste0(env,'_',org,'_',as.character(rep),'.tsv'),sep='\t',row.names='X') # Parse OTU table of the subset of samples only
            otu<-as.matrix(otu)

            # calculate alpha diversity (Shannon Index) in subset of samples
            alpha_div<-diversity(otu, index = "shannon")
            alpha_div<-data.frame(alpha_div)
            
            # calculate betadispersion (heterogeneity) of subset
            groups<- factor(ifelse(rownames(otu0) %in% rownames(prop), "In_Subset", "Out_Subset")) # Use the full-dataset OTU table but distinguishing samples in subset and not in subset
            mod <- betadisper(bc_dist, groups)
            betadisper_res<-data.frame(mod$distance) # Distances to the centroid of the sample subset in the complete dataset bray-curtis matrix
            
            write.csv(betadisper_res,paste0(env,'_',org,'_',as.character(rep),'.betadisp.csv'))
            write.csv(alpha_div,paste0(env,'_',org,'_',as.character(rep),'.shannon.csv'))

        }
    }
}