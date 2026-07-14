import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib_venn import venn3


## Parse data obtained with networkProperties/exportNodeDegreeClosenessBetweenness.py
degree0=pd.read_csv('degrees.csv')
betweenness0=pd.read_csv('betweenness.csv')
closeness0=pd.read_csv('closeness.csv')

## Parse taxonomy data
funTax=pd.read_excel('Taxonomy.xlsx',sheet_name=1).set_index('OTU_ID')
for c in funTax.columns:
    funTax[c]=funTax[c].str.split('__').str[1]
bacTax=pd.read_excel('Taxonomy.xlsx',sheet_name=0).set_index('zOTU')
bacTax=bacTax.rename(columns={'phylum':'Phylum','class':'Class','order':'Order','family':'Family','genus':'Genus','species':'Species'})
bacTax = bacTax.replace("ukn", "")


K_env={}
pK_env={}
allData_env={}
envs=['Cropland','Grassland','Woodland']

for env in envs: # to iterate over multikingdom networks
    netw=env+'_multikingdom.edgelist'

    ## Subset dataframes to only keep the network currently considered then merge them in a single dataframe
    degree=degree0[degree0["Network"]==netw]
    betweenness=betweenness0[betweenness0["Network"]==netw]
    closeness=closeness0[closeness0["Network"]==netw]
    all_dat=degree.merge(betweenness,on='OTU')
    all_dat=all_dat.merge(closeness,on='OTU')[['OTU','Degree','Betweenness','Closeness']]
    allData_env[env]=all_dat.set_index('OTU')

    ## Identify keystones
    thresh=0.99 # percentile threshold for keystone identification
    keystones=all_dat[(all_dat['Degree']>=all_dat['Degree'].quantile(thresh)) & (all_dat['Closeness']>=all_dat['Closeness'].quantile(thresh))] # subsetting based on threshold
    keystones['Organism']=keystones['OTU'].str.startswith('zot')
    keystones['Organism']=keystones['Organism'].map({True:'Bacteria',False:'Fungi'}) # identify if bacterium or fungus
    keystones=keystones.sort_values(by='Organism')

    # New dataframe to add taxonomy information
    keystones_withTax=[]
    keystones_withTax.append(keystones[keystones['Organism']=='Bacteria'].set_index('OTU'))
    keystones_withTax[-1]=keystones_withTax[-1].merge(bacTax,left_index=True,right_index=True,how='left')
    keystones_withTax.append(keystones[keystones['Organism']=='Fungi'].set_index('OTU'))
    keystones_withTax[-1]=keystones_withTax[-1].merge(funTax,left_index=True,right_index=True,how='left')
    keystones_withTax=pd.concat(keystones_withTax)
    K_env[env]=keystones_withTax
    keystones_withTax.to_csv(env+'.keystoneTable.csv') # final table

    ## Identify prime keystones
    prime=keystones_withTax[(keystones_withTax['Betweenness']<=keystones_withTax['Betweenness'].quantile(0.1))] #keystones + low Betweenness (<10th percentile among keystones)
    pK_env[env]=prime
    prime.to_csv(env+'.primeKeystoneTable.csv') # final table

    ## Draw plot showing degree and closeness distribution and keystones highlighted
    j=sns.jointplot(x='Degree',y='Closeness',data=all_dat,color='lightgrey',marginal_kws={'bins':80,'color':"grey"})
    sns.scatterplot(x='Degree',y='Closeness',data=keystones,color='black')
    plt.savefig(netw.replace('.edgelist','.keystones.png'),dpi=300)
    plt.close()

## Prepare Venn diagram showing overlap of keystone sets
palette=['#ffd237','#73dc5c','#f88a22']
envs=['Cropland','Grassland','Woodland']
venn3([set(K_env[env].index) for env in envs], envs,set_colors=palette)
plt.savefig('venn.keystones.pdf')
# print('intersection C-G:',set(K_env['Cropland'].index).intersection(set(K_env['Grassland'].index))) # to print the overlap between cropland and grassland

## Prepare Venn diagram showing overlap of prime keystone sets
venn3([set(pK_env[env].index) for env in envs], envs,set_colors=palette)
plt.savefig('primeKeystones.venn.pdf')
plt.close()

## Get the degree, betweenness and closeness of keystones in every ecosystem types:
vars=['Degree','Betweenness','Closeness']
keystonesAcrossEnv=[]
primeKeystonesAcrossEnv=[]
for env in envs:
    keystones=K_env[env][vars].rename(columns={c:c+'_'+env for c in vars})
    for env2 in envs:
        if env2!=env:
            keystones=keystones.merge(allData_env[env2][vars].rename(columns={c:c+'_'+env2 for c in vars}),left_index=True,right_index=True,how='left')
    primek=keystones[keystones.index.isin(list(pK_env[env].index))]

    keystonesAcrossEnv.append(pd.DataFrame(keystones.stack()).reset_index(drop=False))
    keystonesAcrossEnv[-1]['Var']=keystonesAcrossEnv[-1]['level_1'].str.split('_').str[0]
    keystonesAcrossEnv[-1]['Value in']=keystonesAcrossEnv[-1]['level_1'].str.split('_').str[1]
    keystonesAcrossEnv[-1]['Keystone in']=env

    primeKeystonesAcrossEnv.append(pd.DataFrame(primek.stack()).reset_index(drop=False))
    primeKeystonesAcrossEnv[-1]['Var']=primeKeystonesAcrossEnv[-1]['level_1'].str.split('_').str[0]
    primeKeystonesAcrossEnv[-1]['Value in']=primeKeystonesAcrossEnv[-1]['level_1'].str.split('_').str[1]
    primeKeystonesAcrossEnv[-1]['Prime keystone in']=env
keystonesAcrossEnv=pd.concat(keystonesAcrossEnv)
keystonesAcrossEnv_bact=keystonesAcrossEnv[keystonesAcrossEnv['OTU'].str.startswith('zot')]
keystonesAcrossEnv_fung=keystonesAcrossEnv[~(keystonesAcrossEnv['OTU'].str.startswith('zot'))]
primeKeystonesAcrossEnv=pd.concat(primeKeystonesAcrossEnv)
primeKeystonesAcrossEnv_bact=primeKeystonesAcrossEnv[primeKeystonesAcrossEnv['OTU'].str.startswith('zot')]
primeKeystonesAcrossEnv_fung=primeKeystonesAcrossEnv[~(primeKeystonesAcrossEnv['OTU'].str.startswith('zot'))]


## Plot the degree and closeness of keystones and prime keystone in every ecosystem types on boxplots
palette={'Cropland':'#ffffb3','Grassland':'#ccebc5','Woodland':'#fdc086'}
flierprops = dict(marker='o', markerfacecolor='black', markersize=0.5,linestyle='none')
fig,ax=plt.subplots(2,4,figsize=(12,6),sharey=True)
for var in ['Degree','Closeness']:
    sns.boxplot(data=keystonesAcrossEnv_bact[keystonesAcrossEnv_bact['Var']==var],x=0,y='Keystone in',hue='Value in',palette=palette,hue_order=envs,ax=ax[0][['Degree','Closeness'].index(var)],flierprops=flierprops,linecolor='black',legend=None)
    ax[0][['Degree','Closeness'].index(var)].set_xlabel(var)
    keystonesAcrossEnv_bact[keystonesAcrossEnv_bact['Var']==var].to_csv('data4stats/keystones.%s.Bacteria.csv' % (var))
    sns.boxplot(data=primeKeystonesAcrossEnv_bact[primeKeystonesAcrossEnv_bact['Var']==var],x=0,y='Prime keystone in',hue='Value in',palette=palette,hue_order=envs,ax=ax[0][['Degree','Closeness'].index(var)+2],flierprops=flierprops,linecolor='black',legend=None)
    ax[0][['Degree','Closeness'].index(var)+2].set_xlabel(var)
    primeKeystonesAcrossEnv_bact[primeKeystonesAcrossEnv_bact['Var']==var].to_csv('data4stats/prime.%s.Bacteria.csv' % (var))

    sns.boxplot(data=keystonesAcrossEnv_fung[keystonesAcrossEnv_fung['Var']==var],x=0,y='Keystone in',hue='Value in',palette=palette,hue_order=envs,ax=ax[1][['Degree','Closeness'].index(var)],flierprops=flierprops,linecolor='black',legend=None)
    ax[1][['Degree','Closeness'].index(var)].set_xlabel(var)
    keystonesAcrossEnv_fung[keystonesAcrossEnv_fung['Var']==var].to_csv('data4stats/keystones.%s.Fungi.csv' % (var))
    sns.boxplot(data=primeKeystonesAcrossEnv_fung[primeKeystonesAcrossEnv_fung['Var']==var],x=0,y='Prime keystone in',hue='Value in',palette=palette,hue_order=envs,ax=ax[1][['Degree','Closeness'].index(var)+2],flierprops=flierprops,linecolor='black',legend=None)
    ax[1][['Degree','Closeness'].index(var)+2].set_xlabel(var)
    primeKeystonesAcrossEnv_fung[primeKeystonesAcrossEnv_fung['Var']==var].to_csv('data4stats/prime.%s.Fungi.csv' % (var))
plt.tight_layout()
plt.savefig('keystone.centrality.pdf')
plt.close()