import pandas as pd
import os
import networkx as nx
import networkx.algorithms.community as nx_comm
from cdlib import algorithms
from networkx_robustness import networkx_robustness



def getModularity(G): ## Function that calculates a modularity score Q from a Louvain partition of the network
    Louvain_partition = algorithms.louvain(G,weight='Correlation coefficient')
    Q=nx_comm.modularity(G, Louvain_partition.communities, weight='Correlation coefficient')
    return Q



networkDirectory='LUCAS_networks_2026/fullSizeNetworks/' # to be downloaded here: https://esdac.jrc.ec.europa.eu/content/bacterial-and-fungal-co-occurrence-networks-european-croplands-grasslands-and-woodlands

nInteractionTypes=[] # to store number of associations per group (B-B, F-F, B-F)
Network_stats=[] # to store network properties

for netw in os.listdir(networkDirectory):
    print(netw)
    env=netw.split('_')[0]
    org=netw.split('_')[1]

    envVars=['pH_H2O','pH_CaCl2','EC','OC','CaCO3','P','N','K','monthly_precipitation','monthly_air_temperature'] # sample metadata used in network inferrence with FlashWeave
    if org=='multikingdom':
        fungalOTUs=list(pd.read_csv('./OTU_tables/%s_Fungi.tsv' % (env),sep='\t').set_index('Unnamed: 0').columns) # OTU table of fungi parsed to get a list of fungal OTUs
        bacterialOTUs=list(pd.read_csv('./OTU_tables/%s_Bacteria.tsv' % (env),sep='\t').set_index('Unnamed: 0').columns) # OTU table of bacteria parsed to get a list of bacterial OTUs
        OTUs=bacterialOTUs+fungalOTUs # all OTUs in a list
    else:
        OTUs=list(pd.read_csv('./OTU_tables/%s_%s.tsv' % (env,org),sep='\t').set_index('Unnamed: 0').columns) # OTU table of bacteria|fungi parsed to get a list of OTUs


    ## Create a directory where each OTU is linked to its kingdom        
    otu_map={}
    if org=='multikingdom':
        for o in fungalOTUs:
            otu_map[o]='Fungi'
        for o in bacterialOTUs:
            otu_map[o]='Bacteria'
    else:
        otu_map={o:org for o in OTUs}
    for o in envVars:
        otu_map[o]='Environmental variable'

    ## Parse .edgelist file
    df0=pd.read_csv(networkDirectory+netw,sep='\t',comment="#",header=None)
    df=df0.rename(columns={0:'Node 1',1:'Node 2',2:'Correlation coefficient'})
    df_signif=df[df['Correlation coefficient']!=0]
    df_signif=df_signif[~(df_signif['Node 1'].isin(envVars))] # exclude metadata nodes that are incorporated by FlashWeave
    df_signif=df_signif[~(df_signif['Node 2'].isin(envVars))] # exclude metadata nodes that are incorporated by FlashWeave
    Npos=len(df_signif[df_signif['Correlation coefficient']>0]) # number of positive co-occurrence associations
    Nneg=len(df_signif[df_signif['Correlation coefficient']<0]) # number of negative co-occurrence associations


    ## Look at each association in the dataframe and classify as B-B, F-F or B-F
    for ind in df_signif.index:
        if otu_map[df_signif.loc[ind,'Node 1']]==otu_map[df_signif.loc[ind,'Node 2']] and otu_map[df_signif.loc[ind,'Node 1']]=='Fungi':
            df_signif.loc[ind,'Interaction type']='F-F'
        elif otu_map[df_signif.loc[ind,'Node 1']]==otu_map[df_signif.loc[ind,'Node 2']] and otu_map[df_signif.loc[ind,'Node 1']]=='Bacteria':
            df_signif.loc[ind,'Interaction type']='B-B'
        elif otu_map[df_signif.loc[ind,'Node 1']]!=otu_map[df_signif.loc[ind,'Node 2']]:
            df_signif.loc[ind,'Interaction type']='B-F'
        else:
            df_signif.loc[ind,'Interaction type']="?????"
    nInteractionTypes.append(df_signif['Interaction type'].value_counts().to_dict())
    nInteractionTypes[-1]['Network']=netw
    nInteractionTypes[-1]['Environment']=env

    ## Only keep positive co-occurrence association for future analyses:
    df_signif=df_signif[df_signif['Correlation coefficient']>float(0)].sort_values(by='Correlation coefficient',ascending=True)

    ## Create a metwork in NetworkX
    G=nx.from_pandas_edgelist(df_signif, source='Node 1', target='Node 2', edge_attr='Correlation coefficient')

    ## Extract topological properties of each network
    Network_stats.append({
        'Network':netw,
        'Organism':org,
        'Environment':env,
        'Number of predicted associations':Npos+Nneg,
        'Number of predicted positive associations':Npos,
        'Number of predicted negative associations':Nneg,
        'Proportion of negative associations':Nneg/(Nneg+Npos),
        'Node number':G.number_of_nodes(),
        'Edge number': G.number_of_edges(),
        'Mean coefficient': df_signif['Correlation coefficient'].mean(),
        'Average degree':pd.DataFrame(G.degree())[1].mean(),
        'Network density':nx.density(G),
        'Clustering coefficient':nx.average_clustering(G),
        'Edge/node ratio': G.number_of_edges()/G.number_of_nodes(),
        'Modularity (Q)': getModularity(G),
        'Robustness (Molloy-Reed)': networkx_robustness.molloy_reed(G)
        })

## Export tables
nInteractionTypes=pd.DataFrame(nInteractionTypes)
nInteractionTypes.to_csv('numberOfAssociationsPerType.csv')
Network_stats=pd.DataFrame(Network_stats)
Network_stats.to_csv('networkProperties.csv')
