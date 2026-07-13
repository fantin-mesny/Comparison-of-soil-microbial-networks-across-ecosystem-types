from sklearn.decomposition import PCA
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.preprocessing import StandardScaler


## Parse sample metadata
woodland=pd.read_csv('./metadata/Woodland.soilProp.tsv',sep='\t')
woodland['Environment']='Woodland'
cropland=pd.read_csv('./metadata/Cropland.soilProp.tsv',sep='\t')
cropland['Environment']='Cropland'
grassland=pd.read_csv('./metadata/Grassland.soilProp.tsv',sep='\t')
grassland['Environment']='Grassland'
all=pd.concat([cropland,grassland,woodland]).set_index('SampleID')
all_forPCA=all.drop(columns='Environment')[['pH_H2O','pH_CaCl2','EC','OC','CaCO3','P','N','K','monthly_precipitation','monthly_air_temperature']]

## Standardize data
Xstd = StandardScaler().fit_transform(all_forPCA)

## Compute a PCA on standardized sample properties
n_components = 4
pca = PCA(n_components=n_components)
reduced = pca.fit_transform(Xstd)

## Append the principle components to the 'all' dataframe
for i in range(0, n_components):
    pc=round(pca.explained_variance_ratio_[i]*100,2)
    all['PC' + str(i + 1)+' ('+str(pc)+'%)'] = reduced[:, i]


## Prepare figure
(fig, ax) = plt.subplots(1,2, figsize=(13, 5),width_ratios=(1,0.7))
palette={'Cropland':'#8da0cb','Grassland':'#66c2a5','Woodland':'#fc8d62'}

## Plot the PCA
g = sns.scatterplot(ax=ax[0],x=[col for col in all.columns if col.startswith('PC1')][0],y=[col for col in all.columns if col.startswith('PC2')][0],hue='Environment',data=all.sample(frac = 1),palette=palette,edgecolors='black',alpha=1,sizes=500)

# Calculate and plot a variable factor map (sample properties as vectors) for the first PCs
for i in range(0, pca.components_.shape[1]):
    ax[1].arrow(0,
             0,  # Start the arrow at the origin
             pca.components_[0, i],  #0 for PC1
             pca.components_[1, i],  #1 for PC2
             head_width=0.02,
             head_length=0.02,
             color='black')
    plt.text(pca.components_[0, i] + 0.05,
             pca.components_[1, i] + 0.05,
             all.columns.values[i])
an = np.linspace(0, 2 * np.pi, 100)
plt.plot(np.cos(an), np.sin(an),color='black')  # Add a unit circle for scale
ax[1].axvline(x=0, ymin=-1, ymax=1,ls='--',color='black')
ax[1].axhline(y=0, xmin=-1, xmax=1,ls='--',color='black')
plt.axis('equal')
ax[1].set_title('Variable factor map')
plt.savefig('PCA.pdf') # Save PCA as PDF


