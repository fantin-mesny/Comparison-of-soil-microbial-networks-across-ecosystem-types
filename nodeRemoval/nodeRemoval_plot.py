import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from itertools import combinations
import networkx.algorithms.community as nx_comm
from cdlib import algorithms
from networkx_robustness import networkx_robustness
import random


if __name__ == '__main__':
    envs=['Cropland','Grassland','Woodland']
    netws=['B_','F_']

    # Parse both the properties of complete networks and those of networks that underwent node removal
    stats=pd.read_csv('networkProperties.csv') #obtained by running networkProperties/getNetworkProperties.py
    df=pd.read_csv('nodeRemoval.csv') # obtained by running nodeRemoval.py
    df['Environment']=df['Network'].str.split('_').str[1]
    df['Number of attacked nodes']=df['Attacked nodes'].str.split('-').str.len()
    df['Number of attacked nodes']=df['Number of attacked nodes'].fillna(0)
    df['Network initial density']=df['Network'].map(stats.set_index('Network')['Network density'].to_dict()) # get the density before node removal
    df['Reduction in density (in %)']=((df['Density']-df['Network initial density'])/df['Network initial density'])*100 # value to be plotted


    # Plot:
    fig,ax=plt.subplots(1,2,figsize=(12.5,3.5))
    palette={'Cropland':'#ffd237','Grassland':'#73dc5c','Woodland':'#f88a22'}
    for netw in netws:
        sns.lineplot(data=df[(df['Network'].str.startswith(netw))], marker="o",x='Number of attacked nodes',y='Reduction in density (in %)',hue='Environment',ax=ax[netws.index(netw)],palette=palette)
        ax[netws.index(netw)].set_title(netw)
    plt.savefig('nodeRemoval.pdf')
    plt.close()

























