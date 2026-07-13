import pandas as pd
from statsmodels.stats import multitest

## This script parses the R output of script runIndicSpecies_genuslevel.R to generate a clean table with additional information 

envs=pd.read_csv('allEnvs.env.csv').set_index('group') # every sample associated to its ecosystem type
ranks=['Phylum','Class','Order','Family','Genus']
dfs=[]
for org in ['bacteria','fungi']:
    df=pd.read_csv(org+'_indicGenus.csv').set_index('Unnamed: 0') # parse the R output
    df=df[["index","stat","p.value"]].dropna()
    df['Kingdom']=org
    df['Ecosystem type']=df['index'].map({1:'Cropland',2:'Grassland',3:'Woodland'})
    for rank in ranks:
        df[rank]=df.index.str.split('__').str[ranks.index(rank)]
    df['FDR']=multitest.multipletests(df['p.value'], method='fdr_bh')[1] # Benjamini-hochberg correction of P-values
    df=df.reset_index(drop=False)[['Unnamed: 0','Ecosystem type','Kingdom']+ranks+["stat","p.value",'FDR']].rename(columns={'stat':'Phi','p.value':'P-value'})
    dfs.append(df)
df=pd.concat(dfs).sort_values(by=['Ecosystem type','Kingdom'])

indicators=df[df['FDR']<0.01].reset_index()
nonDupIndex=list(indicators[['Unnamed: 0']].drop_duplicates().index)
indicators=indicators[indicators.index.isin(nonDupIndex)]

# Let's calculate the presence/absence ratios in each ecosystem type
for org in ['bacteria','fungi']:
    otu=pd.read_csv('allEnvs.Genus.pa.'+org+'.csv').set_index('Unnamed: 0') # table used as input in R, 
    otu.columns=otu.columns.str.replace('-','.')
    otu.columns=otu.columns.str.replace('/','.')
    for et in ['Cropland','Grassland','Woodland']:
        otu_env=otu[otu.index.isin(list(envs[envs['treatment']==et].index))]
        otu_notEnv=otu[~(otu.index.isin(list(envs[envs['treatment']==et].index)))]
        otu_env_num=otu_env.sum(axis=0).to_dict()
        otu_notEnv_num=otu_notEnv.sum(axis=0).to_dict()
        for ind in indicators.index:
            if (indicators.loc[ind,'Ecosystem type']==et) & (indicators.loc[ind,'Kingdom']==org):
                indicators.loc[ind,'Ratio in ecosystem type']=str(otu_env_num[indicators.loc[ind,'Unnamed: 0']])+'/'+str(len(otu_env))
                indicators.loc[ind,'Ratio in other ecosystem types']=str(otu_notEnv_num[indicators.loc[ind,'Unnamed: 0']])+'/'+str(len(otu_notEnv))

indicators.drop(columns=['index','Unnamed: 0']).to_csv('FinalTable.csv',index=None)