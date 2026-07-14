import networkx as nx
import pandas as pd
import os
from multiprocessing import Pool


def simple_closeness_attack(G, closeness, numberOfNodesToAttack): #removes the top-closeness nodes in network
    Gnew = G.copy()
    sortedNodes=list(closeness['Node'])
    attacked=[]
    for i in range(0, numberOfNodesToAttack):
        Gnew.remove_node(sortedNodes[i]) # remove the node with the highest degree
        attacked.append(sortedNodes[i])
    return Gnew,attacked

def attackNetwork(G,closeness,maxNodesToRemove,step): # removes increasing numbers of nodes (with the highest closeness) in network 
    Range=range(0,maxNodesToRemove,step)
    disruptionDf=pd.DataFrame(index=range(0,len(Range),1),columns=['Density'])
    Nind=0
    for numberOfNodesToAttack in Range:
        attackedG, attacked = simple_closeness_attack(G, closeness, numberOfNodesToAttack)#refined: 10000
        disruptionDf.loc[Nind,'Density']=nx.density(attackedG)
        disruptionDf.loc[Nind,'Attacked nodes']='-'.join(attacked)
        disruptionDf.loc[Nind,'Number of attacked nodes']=numberOfNodesToAttack
        Nind+=1
    return disruptionDf

def parseAttackAnalyse(network): # node removal for one network
    org=network.split('_')[2].replace('.edgelist','')

    # Parse network:
    df0=pd.read_csv(network,sep='\t',comment="#",header=None) # Parse network
    df=df0.rename(columns={0:'Node 1',1:'Node 2',2:'Correlation coefficient'})
    df=df[~(df['Node 1'].isin(env_vars))] # exclude metadata node
    df=df[~(df['Node 2'].isin(env_vars))] # exclude metadata node
    df_signif=df[df['Correlation coefficient']>0] # only keep positive co-occurrence associations
    G=nx.from_pandas_edgelist(df_signif, source='Node 1', target='Node 2', edge_attr='Correlation coefficient') # create network in NetworkX
    initialDensity=nx.density(G)

    # Calculate node closeness:
    closeness=pd.DataFrame(nx.closeness_centrality(G),index=['Closeness']).T.reset_index(drop=False).rename(columns={'index':'Node'}).sort_values(by='Closeness',ascending=False)
    Nnodes=G.number_of_nodes()
    if org=='Bacteria':
        attack=attackNetwork(G,closeness,3853,1) #because mean bacterial network size is 38525.066666666666, so 10% of it
    elif org=='Fungi':
        attack=attackNetwork(G,closeness,44,1) #because mean bacterial network size is 437.38, so 10% of it
    attack['Reduction in density (in %)']=((attack['Density']-initialDensity)/initialDensity)*100
    attack.to_csv('nodeRemoval_'+network.replace('.edgelist','.csv'))
    return Nnodes




networkDirectory='./LUCAS_networks_2026/networksFromRandomSampleSubsets'
Bnetworks=sorted([File for File in os.listdir(networkDirectory) if File.endswith('_Bacteria.edgelist')])
Fnetworks=sorted([File for File in os.listdir(networkDirectory) if File.endswith('_Fungi.edgelist')]) 
env_vars=['pH_H2O','pH_CaCl2','EC','OC','CaCO3','P','N','K','monthly_precipitation','monthly_air_temperature']



## Run analysis in parallel (pools of 50 CPUs)
data=[]
with Pool() as pool:
    Fdata = pool.map(parseAttackAnalyse, Fnetworks)
print('fungi finished')

with Pool(50) as bpool1:
    Bdata1 = bpool1.map(parseAttackAnalyse, Bnetworks[:50])
print('bacterial batch 1 finished')
with Pool(50) as bpool2:
    Bdata2 = bpool2.map(parseAttackAnalyse, Bnetworks[50:100])
print('bacterial batch 2 finished')
with Pool(50) as bpool3:
    Bdata3 = bpool3.map(parseAttackAnalyse, Bnetworks[100:])
print('bacterial batch 3 finished')

