import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from itertools import combinations
from statsmodels.stats import multitest


def returnLifestylesAboveThreshold(df,colName,threshold): # filter lifestyles that are very rare
    return df[df['Lifestyle'].isin(set(df[df[colName]>threshold]['Lifestyle']))]


def fisherTests(df,ls_list): # compute pairwise Fisher's 
    outputs=[]
    for l in ls_list:
        output_df=[]
        df_l=df[df['Lifestyle']==l]
        df_l.index=df_l['Lifestyle']+' '+df_l['Environment']
        pairs = list(combinations(['Cropland','Grassland','Woodland'], 2))
        for pair in pairs:
            table=[[df_l.loc[l+' '+pair[0],'inNetwork'],df_l.loc[l+' '+pair[0],'Num in network']-df_l.loc[l+' '+pair[0],'inNetwork']],
                   [df_l.loc[l+' '+pair[1],'inNetwork'],df_l.loc[l+' '+pair[1],'Num in network']-df_l.loc[l+' '+pair[1],'inNetwork']]]
            st,p=stats.fisher_exact(table)
            output_df.append({'Lifestyle':l,'Comparison':pair[0]+'-'+pair[1],'Stat':st,'P-value':p})
        output_df=pd.DataFrame(output_df)
        output_df['FDR']=multitest.multipletests(output_df['P-value'])[1]
        outputs.append(output_df)
    return pd.concat(outputs)



if __name__ == '__main__':

    networkDirectory='./LUCAS_networks_2026/fullSizeNetworks/'
    envs=['Cropland','Grassland','Woodland']

    ## Parse functional annotations of fungi and bacteria
    annot=pd.read_excel('./functionalAnnotations/fungi.annot.xlsx')
    annot_dict=annot.set_index('OTU')['Primary lifestyle'].to_dict()
    annot_bacteria=pd.read_csv('./functionalAnnotations/bacteria.annot.csv').set_index('Unnamed: 0')


    nodes_df=[]
    bacNodes_df=[]
    for env in envs: # iterate over multikingdom networks
        netw=env+'_multikingdom.edgelist'

        envVars=['pH_H2O','pH_CaCl2','EC','OC','CaCO3','P','N','K','monthly_precipitation','monthly_air_temperature'] # sample metadata used in FlashWeave
        fungalOTUs=list(pd.read_csv('./OTU_tables/%s_Fungi.tsv' % (env),sep='\t').set_index('Unnamed: 0').columns) # OTU table of fungi parsed to get a list of fungal OTUs
        bacterialOTUs=list(pd.read_csv('./OTU_tables/%s_Bacteria.tsv' % (env),sep='\t').set_index('Unnamed: 0').columns) # OTU table of bacteria parsed to get a list of bacterial OTUs
        OTUs=bacterialOTUs+fungalOTUs

        otu_map={}
        for o in fungalOTUs:
            otu_map[o]='Fungi'
        for o in bacterialOTUs:
            otu_map[o]='Bacteria'

        df0=pd.read_csv(networkDirectory+netw,sep='\t',comment="#",header=None) # Parse network (.edgelist)
        df=df0.rename(columns={0:'Node 1',1:'Node 2',2:'Correlation coefficient'})
        df_signif=df[df['Correlation coefficient']!=0]
        df_signif=df_signif[~(df_signif['Node 1'].isin(envVars))] # exclude metadata nodes
        df_signif=df_signif[~(df_signif['Node 2'].isin(envVars))] # exclude metadata nodes

        df_signif=df_signif[df_signif['Correlation coefficient']>0].sort_values(by='Correlation coefficient',ascending=True) # keep only positive cooccurrence associations
        df_signif['Node 1 kingdom']=df_signif['Node 1'].map(otu_map) # add node kingdom in table
        df_signif['Node 2 kingdom']=df_signif['Node 2'].map(otu_map) # add node kingdom in table


        # Analyse the fungi in network 
        fungiInNetwork=list(set(list(df_signif[(df_signif['Node 1 kingdom']=='Fungi')]['Node 1'])+list(df_signif[(df_signif['Node 2 kingdom']=='Fungi')]['Node 2']))) # list of all fungal nodes in network
        fungiInNetwork_df=pd.DataFrame(index=fungiInNetwork,columns=['Lifestyle']) # added to a new dataframe
        fungiInNetwork_df['Lifestyle']=fungiInNetwork_df.index.map(annot_dict) # lifestyle added
        counts=pd.DataFrame(fungiInNetwork_df['Lifestyle'].value_counts()).rename(columns={'count':'inNetwork'}) # count lifestyle nodes
        counts['Percentage of fungal nodes']=counts['inNetwork']/len(fungiInNetwork)
        counts['Num in network']=len(fungiInNetwork)
        counts['Environment']=env
        nodes_df.append(counts)

        ## Bacteria
        bacteriaInNetwork=list(set(list(df_signif[(df_signif['Node 1 kingdom']=='Bacteria')]['Node 1'])+list(df_signif[(df_signif['Node 2 kingdom']=='Bacteria')]['Node 2']))) # list of all fungal nodes in network
        bacteriaInNetwork_df=annot_bacteria[annot_bacteria.index.isin(bacteriaInNetwork)] # added to a new dataframe
        counts=pd.DataFrame(bacteriaInNetwork_df.sum(axis=0)).rename(columns={0:'inNetwork'}) # function added added
        counts['Percentage of bacterial nodes']=counts['inNetwork']/len(bacteriaInNetwork) # count function nodes
        counts['Num in network']=len(bacteriaInNetwork)
        counts['Environment']=env
        bacNodes_df.append(counts)


## Plot fungi
nodes_df=pd.concat(nodes_df).reset_index(drop=False)
palette={'Cropland':'#ffd237','Grassland':'#73dc5c','Woodland':'#f88a22'}
top_ls=['soil_saprotroph','litter_saprotroph','wood_saprotroph','dung_saprotroph','pollen_saprotroph','unspecified_saprotroph','arbuscular_mycorrhizal','ectomycorrhizal','plant_pathogen','root_endophyte','mycoparasite','animal_parasite']
fig,ax=plt.subplots(1,2,sharey=True,figsize=(10,15))
nodes_df_4plot=returnLifestylesAboveThreshold(nodes_df,'Percentage of fungal nodes',0.01)
sns.barplot(x='Percentage of fungal nodes',y='Lifestyle',hue='Environment',ax=ax[0],data=nodes_df_4plot,palette=palette,order=top_ls)
plt.close()
## Fisher test for fungi:
print(fisherTests(nodes_df_4plot,top_ls))

## Plot bacteria
bacNodes_df=pd.concat(bacNodes_df).reset_index(drop=False).rename(columns={'index':'Lifestyle'})
palette={'Cropland':'#ffd237','Grassland':'#73dc5c','Woodland':'#f88a22'}
top_ls=['chemoheterotrophy','aerobic_chemoheterotrophy','anaerobic_chemoheterotrophy','fermentation','nitrate_reduction','ureolysis']
bacNodes_df_4plot=returnLifestylesAboveThreshold(bacNodes_df,'Percentage of bacterial nodes',0.01)
fig,ax=plt.subplots(1,2,sharey=True,figsize=(10,15))
sns.barplot(x='Percentage of bacterial nodes',y='Lifestyle',hue='Environment',ax=ax[0],data=bacNodes_df_4plot,palette=palette,order=top_ls)
plt.savefig('bacterialLifestyles.pdf')
plt.close()
## Fisher test for bacteria 
print(fisherTests(bacNodes_df_4plot,top_ls))