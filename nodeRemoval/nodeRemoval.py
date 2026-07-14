import pandas as pd
import os
import networkx as nx
import networkx.algorithms.community as nx_comm


def simple_closeness_attack(G,closeness, attack_fraction=0.1): #removes the top-closeness nodes in network
    Gnew = G.copy()
    G_nodes = G.number_of_nodes()
    sortedNodes=list(closeness['OTU'])
    attacked=[]
    for i in range(0, int(G_nodes * attack_fraction)):
        # remove the node with the highest degree
        Gnew.remove_node(sortedNodes[i])
        attacked.append(sortedNodes[i])
    return Gnew,attacked

def attackNetwork(G,closeness): # removes increasing proportions of nodes (with the highest closeness) in network 
    propRange=range(0,60,1) # from 0 to 6% with increment of 1 (divided by 1000 lower)
    disruptionDf=pd.DataFrame(index=range(0,len(propRange),1),columns=['Molloy-Reed'])
    Nind=0
    for prop in propRange:
        attackedG, attacked = simple_closeness_attack(G, closeness, attack_fraction=prop/1000)
        disruptionDf.loc[Nind,'Proportion of attacked nodes']=prop/1000 
        disruptionDf.loc[Nind,'Density']=nx.density(attackedG)
        disruptionDf.loc[Nind,'Attacked nodes']='-'.join(attacked)
        print(prop/1000,'done')
        Nind+=1
    return disruptionDf


if __name__ == '__main__':

    ## Parse data obtained with networkProperties/exportNodeDegreeClosenessBetweenness.py
    degree0=pd.read_csv('degrees.csv')
    betweenness0=pd.read_csv('betweenness.csv')
    closeness0=pd.read_csv('closeness.csv')

    Nax=0
    envs=['Cropland','Grassland','Woodland']
    envVars=['pH_H2O','pH_CaCl2','EC','OC','CaCO3','P','N','K','monthly_precipitation','monthly_air_temperature'] # sample metadata used in FlashWeave
    networkDirectory='./LUCAS_networks_2026/fullSizeNetworks'
    netws=[net for net in os.listdir(networkDirectory)] # to be downloaded here: https://esdac.jrc.ec.europa.eu/content/bacterial-and-fungal-co-occurrence-networks-european-croplands-grasslands-and-woodlands
    Network_stats=pd.DataFrame(index=envs,columns=['Environment','Percentage of negative interactions','Network density','Clustering coefficient'])

    attacked_Df=[]
    for env in envs:
        for netw in [netEnv for netEnv in netws if env in netEnv]: # Parse networks of one ecosystem-type
            print(env,netw)

            ## Subset dataframes to keep only the considered network
            degree=degree0[degree0["Network"]==netw]
            betweenness=betweenness0[betweenness0["Network"]==netw]
            closeness=closeness0[closeness0["Network"]==netw]
            all_dat=degree.merge(betweenness,on="OTU")
            all_dat=all_dat.merge(closeness,on="OTU")[['OTU','Degree','Betweenness','Closeness']]
            all_dat['Group']=all_dat['OTU'].str.startswith('zot')
            all_dat['Group']=all_dat['Group'].map({True:'Bacteria',False:'Fungi'})

            ## Parse network:
            fileName=networkDirectory+netw
            df0=pd.read_csv(fileName,sep='\t',comment="#",header=None) # parse .edgelist file
            df=df0.rename(columns={0:'Node 1',1:'Node 2',2:'Correlation coefficient'})
            df_signif=df[df['Correlation coefficient']>0] # Only positive co-occurrence associations
            df_signif=df_signif[~(df_signif['Node 1'].isin(envVars))] # Exclude metadata nodes
            df_signif=df_signif[~(df_signif['Node 2'].isin(envVars))] # Exclude metadata nodes
            G=nx.from_pandas_edgelist(df_signif, source='Node 1', target='Node 2', edge_attr='Correlation coefficient')

            # Sort node properties by decreasing closeness
            closeness = all_dat.sort_values(by='Closeness',ascending=False)
            removal=attackNetwork(G,closeness)
            removal['Environment']=env
            removal['Network']=netw
            attacked_Df.append(removal)

    df=pd.concat(attacked_Df)
    df.to_csv('nodeRemoval.csv')

























