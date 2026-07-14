# Comparison of soil microbial networks across ecosystem types

This repository includes scripts used in the study *The complexity and robustness of soil microbial networks depend on ecosystem type* (Mesny *et al.*)

### Microbial community homogeneity analysis
- Beta-dispersion (heterogeneity) of microbial communities
- HOMOVA test

&rarr; See [communityHeterogeneity](https://github.com/fantin-mesny/Comparison-of-soil-microbial-networks-across-ecosystem-types/tree/main/communityHeterogeneity)

### Analyses of soil sample properties 
- PCA of soil sample edaphic and climatic metadata
- Beta-dispersion analysis testing for heterogeneity of environmental properties
- Moran's I and join-counts tests for spatial autocorrelation

&rarr; See [sampleProperties](https://github.com/fantin-mesny/Comparison-of-soil-microbial-networks-across-ecosystem-types/tree/main/sampleProperties)

### Indicator genus identification
- OTU data transformation (genus level, presence/absence)
- Indicator taxa analysis

&rarr; See [indicatorTaxa](https://github.com/fantin-mesny/Comparison-of-soil-microbial-networks-across-ecosystem-types/tree/main/indicatorTaxa)

### Network property calculation
- Network property extraction from networks reconstructed from complete sets of samples
- Subsampling the dataset to reconstruct networks from 50 subsets of 50 samples per ecosystem types
- Network property extraction and analysis from networks reconstructed from 50 subsets of 50 samples per ecosystem type
- Shuffling the dataset to reconstruct null-model networks
- Comparison of null-model networks to networks reconstructed from 50 subsets of 50 samples per ecosystem type
- Export node degree, closeness and betweenness from networks for further analyses (*e.g.* keystone identification)

&rarr; See [networkProperties](https://github.com/fantin-mesny/Comparison-of-soil-microbial-networks-across-ecosystem-types/tree/main/networkProperties)

### Analysis of factors associated to network topology differentiation
- Calculation of community alpha-diversity and sample subset beta-dispersion
- Calculation of geographical descriptors from latitudes and longitudes
- Variance partitioning analysis

&rarr; See [factorsUnderlyingNetworkDifferentiation](https://github.com/fantin-mesny/Comparison-of-soil-microbial-networks-across-ecosystem-types/tree/main/factorsUnderlyingNetworkDifferentiation)

### Node removal analyses
- Removal of high-closeness nodes from networks reconstructed from complete sets of samples and plotting
- Removal of high-closeness nodes from networks reconstructed from 50 subsets of 50 samples per ecosystem type, AUC calculation and plotting

&rarr; See [nodeRemoval](https://github.com/fantin-mesny/Comparison-of-soil-microbial-networks-across-ecosystem-types/tree/main/nodeRemoval)

### Analysis of bacterial functional groups and fungal lifestyles
- Calculating proportions of nodes from each functional group in multikingdom networks
- Analysing the degree and closeness centrality of nodes from each functional group in multikingdom networks

&rarr; See [microbialFunctions](https://github.com/fantin-mesny/Comparison-of-soil-microbial-networks-across-ecosystem-types/tree/main/microbialFunctions)

### Keystone taxa identification and analyses
- Identification of microbial keystones in multikingdom networks and analysis of their specificity (Venn diagrams, closeness/degree in each ecosystem type)
- Analysis of the prevalence and relative abundance of each keystone in each ecosystem type-specific multikingdom network
- Regression analysis testing for association between keystone richness and environmental variables

&rarr; See [microbialKeystones](https://github.com/fantin-mesny/Comparison-of-soil-microbial-networks-across-ecosystem-types/tree/main/microbialKeystones)


---

### Data availability:
Networks can be downloaded from the following ESDAC repository: [Bacterial and fungal co-occurrence networks in European croplands, grasslands and woodlands](https://esdac.jrc.ec.europa.eu/content/bacterial-and-fungal-co-occurrence-networks-european-croplands-grasslands-and-woodlands)

