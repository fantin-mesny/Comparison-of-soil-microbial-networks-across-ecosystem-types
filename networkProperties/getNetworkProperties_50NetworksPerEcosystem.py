import networkx as nx
import pandas as pd
import os
import seaborn as sns
import matplotlib.pyplot as plt
import networkx.algorithms.community as nx_comm
from cdlib import algorithms
from networkx_robustness import networkx_robustness

def getModularity(G):
    Louvain_partition = algorithms.louvain(G,weight='Correlation coefficient')
    Q=nx_comm.modularity(G, Louvain_partition.communities, weight='Correlation coefficient')
    return Q

networkDirectory='LUCAS_networks_2026/networksFromRandomSampleSubsets/' # to be downloaded here: https://esdac.jrc.ec.europa.eu/content/bacterial-and-fungal-co-occurrence-networks-european-croplands-grasslands-and-woodlands
networks=sorted([File for File in os.listdir('networkDirectory') if File.endswith('.edgelist')]) # each network (50 per ecosystem type)
env_vars=['pH_H2O','pH_CaCl2','EC','OC','CaCO3','P','N','K','monthly_precipitation','monthly_air_temperature'] # sample metadata used in FlashWeave


data=[]
for network in networks: # to process each of the network independently
    org=network.split('_')[-1].replace('.edgelist','') # reads in file name if Bacteria|Fungi|multikingdom
    df0=pd.read_csv(network,sep='\t',comment="#",header=None) # parse .edgelist
    
    df=df0.rename(columns={0:'Node 1',1:'Node 2',2:'Correlation coefficient'})
    df=df[~(df['Node 1'].isin(env_vars))] # exclude metadata nodes that are incorporated by FlashWeave
    df=df[~(df['Node 2'].isin(env_vars))] # exclude metadata nodes that are incorporated by FlashWeave

    df_signif=df[df['Correlation coefficient']>0] # Only consider positive co-occurrence association
    G=nx.from_pandas_edgelist(df_signif, source='Node 1', target='Node 2', edge_attr='Correlation coefficient') # Reconstruct network with Networkx

    ## Calculate network topological properties:
    data.append({
        'Network':network,
        'Organism':org,
        'Repeat': network.split('_')[1],
        'Environment':network.split('_')[0],
        'Node number':G.number_of_nodes(),
        'Edge number': G.number_of_edges(),
        'Mean interaction coefficient': df_signif['Correlation coefficient'].mean(),
        'Degree heterogeneity':pd.DataFrame(G.degree())[1].var(),
        'Network density':nx.density(G),
        'Clustering coefficient':nx.average_clustering(G),
        #'Edge/node ratio': G.number_of_edges()/G.number_of_nodes(),
        'Modularity (Q)': getModularity(G)
        #'Robustness (Molloy-Reed)': networkx_robustness.molloy_reed(G)
        })

data=pd.DataFrame(data)
data.to_csv('networkCharacteristics.csv')

## Plot and export properties as boxplots:
flierprops = dict(marker='o', markerfacecolor='black', markersize=1,linestyle='none')
palette={'Cropland':'#ffffb3','Grassland':'#ccebc5','Woodland':'#fdc086'}
order=['Cropland','Grassland','Woodland']
fig,ax=plt.subplots(3,6,figsize=(15,15),sharex=True)

sns.boxplot(x='Environment',y='Node number',data=data[data['Organism']=='Bacteria'],ax=ax[0][0],order=order, palette=palette,flierprops=flierprops,linecolor='black')
sns.boxplot(x='Environment',y='Edge number',data=data[data['Organism']=='Bacteria'],ax=ax[0][1],order=order, palette=palette,flierprops=flierprops,linecolor='black')
sns.boxplot(x='Environment',y='Average degree',data=data[data['Organism']=='Bacteria'],ax=ax[0][2],order=order, palette=palette,flierprops=flierprops,linecolor='black')
sns.boxplot(x='Environment',y='Network density',data=data[data['Organism']=='Bacteria'],ax=ax[0][3],order=order, palette=palette,flierprops=flierprops,linecolor='black')
sns.boxplot(x='Environment',y='Clustering coefficient',data=data[data['Organism']=='Bacteria'],ax=ax[0][4],order=order, palette=palette,flierprops=flierprops,linecolor='black')
sns.boxplot(x='Environment',y='Modularity (Q)',data=data[data['Organism']=='Bacteria'],ax=ax[0][5],order=order, palette=palette,flierprops=flierprops,linecolor='black')
ax[0][0].set_title('Bacteria', loc='left')

sns.boxplot(x='Environment',y='Node number',data=data[data['Organism']=='Fungi'],ax=ax[1][0],order=order, palette=palette,flierprops=flierprops,linecolor='black')
sns.boxplot(x='Environment',y='Edge number',data=data[data['Organism']=='Fungi'],ax=ax[1][1],order=order, palette=palette,flierprops=flierprops,linecolor='black')
sns.boxplot(x='Environment',y='Average degree',data=data[data['Organism']=='Fungi'],ax=ax[1][2],order=order, palette=palette,flierprops=flierprops,linecolor='black')
sns.boxplot(x='Environment',y='Network density',data=data[data['Organism']=='Fungi'],ax=ax[1][3],order=order, palette=palette,flierprops=flierprops,linecolor='black')
sns.boxplot(x='Environment',y='Clustering coefficient',data=data[data['Organism']=='Fungi'],ax=ax[1][4],order=order, palette=palette,flierprops=flierprops,linecolor='black')
sns.boxplot(x='Environment',y='Modularity (Q)',data=data[data['Organism']=='Fungi'],ax=ax[1][5],order=order, palette=palette,flierprops=flierprops,linecolor='black')
ax[1][0].set_title('Fungi', loc='left')

sns.boxplot(x='Environment',y='Node number',data=data[data['Organism']=='multikingdom'],ax=ax[2][0],order=order, palette=palette,flierprops=flierprops,linecolor='black')
sns.boxplot(x='Environment',y='Edge number',data=data[data['Organism']=='multikingdom'],ax=ax[2][1],order=order, palette=palette,flierprops=flierprops,linecolor='black')
sns.boxplot(x='Environment',y='Average degree',data=data[data['Organism']=='multikingdom'],ax=ax[2][2],order=order, palette=palette,flierprops=flierprops,linecolor='black')
sns.boxplot(x='Environment',y='Network density',data=data[data['Organism']=='multikingdom'],ax=ax[2][3],order=order, palette=palette,flierprops=flierprops,linecolor='black')
sns.boxplot(x='Environment',y='Clustering coefficient',data=data[data['Organism']=='multikingdom'],ax=ax[2][4],order=order, palette=palette,flierprops=flierprops,linecolor='black')
sns.boxplot(x='Environment',y='Modularity (Q)',data=data[data['Organism']=='multikingdom'],ax=ax[2][5],order=order, palette=palette,flierprops=flierprops,linecolor='black')
ax[2][0].set_title('multikingdom Bacteria+Fungi', loc='left')

for i in [0,1,2,3,4,5]:
    for i2 in [0,1,2]:
        bottom, top=ax[i2][i].get_ylim()
        ax[i2][i].set_ylim(bottom, top+((top-bottom)/10))

plt.tight_layout()
plt.savefig('boxplots_propertiesNetworks.pdf') #Fig.2
plt.close()