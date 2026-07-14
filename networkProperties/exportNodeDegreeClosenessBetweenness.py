import pandas as pd
import sys
import os
import networkx as nx
from concurrent.futures import ProcessPoolExecutor, as_completed


def runAnalysis(netw):
    envVars=['pH_H2O','pH_CaCl2','EC','OC','CaCO3','P','N','K','monthly_precipitation','monthly_air_temperature']   

    df0=pd.read_csv(netw,sep='\t',comment="#",header=None) # parse .edgelist network file
    df_signif=df0.rename(columns={0:'Node 1',1:'Node 2',2:'Correlation coefficient'})
    df_signif=df_signif[~(df_signif['Node 1'].isin(envVars))] # exclude metadata nodes
    df_signif=df_signif[~(df_signif['Node 2'].isin(envVars))] # exclude metadata nodes
    df_signif=df_signif[df_signif['Correlation coefficient']>0].sort_values(by='Correlation coefficient',ascending=True) # Keep only positive co-occurrence associations

    G=nx.from_pandas_edgelist(df_signif, source='Node 1', target='Node 2', edge_attr='Correlation coefficient') # create network in NetworkX

    # Calculate node betweenness, closeness and degree
    print(netw,'betweenness...')
    bet=pd.DataFrame(nx.betweenness_centrality(G),index=[netw.split('/')[-1]]).T
    print(netw,'closeness...')
    clo=pd.DataFrame(nx.closeness_centrality(G),index=[netw.split('/')[-1]]).T
    print(netw,'degrees...')
    deg=pd.DataFrame(nx.degree_centrality(G),index=[netw.split('/')[-1]]).T

    return bet, clo, deg



if __name__ == '__main__':
    a = get_params(sys.argv[1:])

    Node_betweenness=[]
    Node_degree=[]
    Node_closeness=[]
    max_workers = os.cpu_count()  # available CPUs on the machine
    networkDirectory='./LUCAS_networks_2026/fullSizeNetworks/'
    nets = os.listdir(networkDirectory)

    ## Parse every network and get their betweenness, closeness and degree - runs in parallel on multiple cores
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(runAnalysis, netw): networkDirectory+netw for netw in nets}
        for future in as_completed(futures):
            bet, clo, deg = future.result()
            Node_degree.append(deg)
            Node_betweenness.append(bet)
            Node_closeness.append(clo)

    
    ## Concatenate all outputs and save individual tables
    print('    degree')
    Node_degree=pd.concat(Node_degree,axis=1).stack().reset_index(drop=False).rename(columns={'level_0':'OTU','level_1':'Network',0:'Degree'}).sort_values(by='Network')
    Node_degree.to_csv('degrees.csv',index=False)

    print('    betweenness')
    Node_betweenness=pd.concat(Node_betweenness,axis=1).stack().reset_index(drop=False).rename(columns={'level_0':'OTU','level_1':'Network',0:'Betweenness'}).sort_values(by='Network')
    Node_betweenness.to_csv('betweenness.csv',index=False)

    print('    closeness')
    Node_closeness=pd.concat(Node_closeness,axis=1).stack().reset_index(drop=False).rename(columns={'level_0':'OTU','level_1':'Network',0:'Closeness'}).sort_values(by='Network')
    Node_closeness.to_csv('closeness.csv',index=False)
