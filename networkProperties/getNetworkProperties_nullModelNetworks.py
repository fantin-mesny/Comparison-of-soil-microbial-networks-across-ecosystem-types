import networkx as nx
import pandas as pd
import os
import seaborn as sns
import matplotlib.pyplot as plt
import networkx.algorithms.community as nx_comm
from cdlib import algorithms
from networkx_robustness import networkx_robustness
import itertools as it
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed

def getModularity(G): ## Function that calculates a modularity score Q from a Louvain partition of the network
    Louvain_partition = algorithms.louvain(G,weight='Correlation coefficient')
    Q=nx_comm.modularity(G, Louvain_partition.communities, weight='Correlation coefficient')
    return Q

def process_network(network0): # Function to process every network independently, to be run in parallel for computational time reasons
    network=network0.split('/')[-1]
    if 'Bacteria' in network or 'Fungi' in network: #not to consider multikingdom networks here
        org=network.split('_')[-1].replace('.edgelist','')
        env=network.split('_')[0]
        rep=network.split('_')[1]

        df0=pd.read_csv(network0,sep='\t',comment="#",header=None) # Parse network file
        df=df0.rename(columns={0:'Node 1',1:'Node 2',2:'Correlation coefficient'})
        df=df[~(df['Node 1'].isin(env_vars))] # exclude metadata nodes that are incorporated by FlashWeave
        df=df[~(df['Node 2'].isin(env_vars))] # exclude metadata nodes that are incorporated by FlashWeave
        df_signif=df[df['Correlation coefficient']>0] # To only consider positive co-occurrence associations

        G=nx.from_pandas_edgelist(df_signif, source='Node 1', target='Node 2', edge_attr='Correlation coefficient') # Reconstruct network with Networkx
        return {
            'Network':network,
            'Organism':org,
            'Repeat': rep,
            'Environment':env,
            'Node number':G.number_of_nodes(),
            'Edge number': G.number_of_edges(),
            'Mean interaction coefficient': df_signif['Correlation coefficient'].mean(),
            'Degree heterogeneity':pd.DataFrame(G.degree())[1].var(),
            'Network density':nx.density(G),
            'Clustering coefficient':nx.average_clustering(G),
            #'Edge/node ratio': G.number_of_edges()/G.number_of_nodes(),
            'Modularity (Q)': getModularity(G)
            #'Robustness (Molloy-Reed)': networkx_robustness.molloy_reed(G)
        }


## Directory with networks reconstructed from  50 subsets of 50 samples
networkDirectory='LUCAS_networks_2026/networksFromRandomSampleSubsets/' # to be downloaded here: https://esdac.jrc.ec.europa.eu/content/bacterial-and-fungal-co-occurrence-networks-european-croplands-grasslands-and-woodlands
networks=sorted([File for File in os.listdir('.') if File.endswith('.edgelist')])

## Directory with null model networks
nullNetworkDirectory='./nullNetworks'
networks=networks+sorted([nullNetworkDirectory+'/'+File for File in os.listdir(nullNetworkDirectory) if File.endswith('.edgelist')]) # both together
networks=networks[::-1]

env_vars=['pH_H2O','pH_CaCl2','EC','OC','CaCO3','P','N','K','monthly_precipitation','monthly_air_temperature']
orgs={'B':'Bacteria','F':'Fungi'}
multiking=True
data=[]

## Parse and analyse every network in parallel to save time
max_workers = os.cpu_count() or 1 
with ProcessPoolExecutor(max_workers=max_workers) as executor:
    futures = [executor.submit(process_network, network0) for network0 in networks]
    for fut in as_completed(futures):
        result = fut.result()
        if result is not None:
            data.append(result)
data=pd.DataFrame(data)
data.to_csv('networkCharacteristics.csv')

## Make a boxplots figure:
flierprops = dict(marker='o', markerfacecolor='black', markersize=1,linestyle='none')
palette={'Null':'white','Cropland':'#ffffb3','Grassland':'#ccebc5','Woodland':'#fdc086'}
order=['Null','Cropland','Grassland','Woodland']
fig,ax=plt.subplots(2,5,figsize=(15,15),sharex=True)
sns.boxplot(x='Environment',y='Node number',data=data[data['Organism']=='Bacteria'],ax=ax[0][0],order=order, palette=palette,flierprops=flierprops,linecolor='black')
sns.boxplot(x='Environment',y='Edge number',data=data[data['Organism']=='Bacteria'],ax=ax[0][1],order=order, palette=palette,flierprops=flierprops,linecolor='black')
#sns.boxplot(x='Environment',y='Average degree',data=data[data['Organism']=='Bacteria'],ax=ax[0][2],order=order, palette=palette,flierprops=flierprops,linecolor='black')
sns.boxplot(x='Environment',y='Network density',data=data[data['Organism']=='Bacteria'],ax=ax[0][2],order=order, palette=palette,flierprops=flierprops,linecolor='black')
sns.boxplot(x='Environment',y='Clustering coefficient',data=data[data['Organism']=='Bacteria'],ax=ax[0][3],order=order, palette=palette,flierprops=flierprops,linecolor='black')
#sns.boxplot(x='Environment',y='Edge/node ratio',data=data[data['Organism']=='Bacteria'],ax=ax[0][5],order=order, color='grey')
sns.boxplot(x='Environment',y='Modularity (Q)',data=data[data['Organism']=='Bacteria'],ax=ax[0][4],order=order, palette=palette,flierprops=flierprops,linecolor='black')
#sns.boxplot(x='Environment',y='Robustness (Molloy-Reed)',data=data[data['Organism']=='Bacteria'],ax=ax[0][6],order=order, color='grey')
#sns.violinplot(x='Environment',y='Mean interaction coefficient',data=data[data['Organism']=='Bacteria'],ax=ax[0][7],order=order, color='grey')
ax[0][0].set_title('Bacteria', loc='left')

sns.boxplot(x='Environment',y='Node number',data=data[data['Organism']=='Fungi'],ax=ax[1][0],order=order, palette=palette,flierprops=flierprops,linecolor='black')
sns.boxplot(x='Environment',y='Edge number',data=data[data['Organism']=='Fungi'],ax=ax[1][1],order=order, palette=palette,flierprops=flierprops,linecolor='black')
#sns.boxplot(x='Environment',y='Average degree',data=data[data['Organism']=='Fungi'],ax=ax[1][2],order=order, palette=palette,flierprops=flierprops,linecolor='black')
sns.boxplot(x='Environment',y='Network density',data=data[data['Organism']=='Fungi'],ax=ax[1][2],order=order, palette=palette,flierprops=flierprops,linecolor='black')
sns.boxplot(x='Environment',y='Clustering coefficient',data=data[data['Organism']=='Fungi'],ax=ax[1][3],order=order, palette=palette,flierprops=flierprops,linecolor='black')
#sns.boxplot(x='Environment',y='Edge/node ratio',data=data[data['Organism']=='Fungi'],ax=ax[1][5],order=order, color='grey')
sns.boxplot(x='Environment',y='Modularity (Q)',data=data[data['Organism']=='Fungi'],ax=ax[1][4],order=order, palette=palette,flierprops=flierprops,linecolor='black')
#sns.boxplot(x='Environment',y='Robustness (Molloy-Reed)',data=data[data['Organism']=='Fungi'],ax=ax[1][6],order=order, color='grey')
#sns.violinplot(x='Environment',y='Mean interaction coefficient',data=data[data['Organism']=='Fungi'],ax=ax[1][7],order=order, color='grey')
ax[1][0].set_title('Fungi', loc='left')

plt.tight_layout()
plt.savefig('boxplots_withNullNetworks.pdf')
plt.close()