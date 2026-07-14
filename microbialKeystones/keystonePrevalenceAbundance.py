import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def getPrevalence(OTUtable,environment): #calculates prevalence of taxa in an OTU table
    funPrevalence=OTUtable>0
    TotalSample=len(funPrevalence)
    funPrevalence=funPrevalence.T
    funPrevalence['Nsample']=funPrevalence.sum(axis=1)
    funPrevalence['Prevalence']=funPrevalence['Nsample']/TotalSample
    funPrevalence=funPrevalence[['Prevalence']].rename(columns={c:c+'_'+environment for c in funPrevalence.columns})
    return funPrevalence.reset_index(drop=False).rename(columns={'index':'OTU'})

def getRA(OTUtable,environment): #calculates relative abundances of taxa in an OTU table
    funAbundance=OTUtable.T
    funAbundance=funAbundance.div(funAbundance.sum())
    funAbundance['Mean relative abundance']=funAbundance.sum(axis=1)/len(funAbundance.columns)
    funAbundance=funAbundance[['Mean relative abundance']].rename(columns={c:c+'_'+environment for c in funAbundance.columns})
    return funAbundance.reset_index(drop=False).rename(columns={'index':'OTU'})

### Parse OTU tables and calculate OTU abundance/prevalence data
envs=['Cropland','Grassland','Woodland']
orgs=['Bacteria','Fungi']
otuTables={}
abundances={}
prevalences={}
for env in envs:
    for org in orgs:
        otuTables['%s_%s' % (env,org)]=pd.read_csv('./OTU_tables/%s_%s.tsv' % (env,org),sep='\t').set_index('Unnamed: 0')
        prevalences['%s_%s' % (env,org)]=getPrevalence(otuTables['%s_%s' % (env,org)],env)
        abundances['%s_%s' % (env,org)]=getRA(otuTables['%s_%s' % (env,org)],env)


## Load Keystones and get their abundance/prevalence across ecosystem types
keystones_out={org:[] for org in orgs}
prime_out={org:[] for org in orgs}
for env in envs:
    keystones=pd.read_csv(env+'.keystoneTable.csv') # tables generated with keystoneIdentification.py
    prime=pd.read_csv(env+'.primeKeystoneTable.csv') # tables generated with keystoneIdentification.py
    for org in orgs:
        keystones_org=keystones[keystones['Organism']==org]
        prime_org=prime[prime['Organism']==org]
        keystones_org=keystones[keystones['Organism']==org]
        prime_org=prime[prime['Organism']==org]
        for env2 in envs:
            keystones_org=keystones_org.merge(prevalences['%s_%s' % (env2,org)],on='OTU',how='left').merge(abundances['%s_%s' % (env2,org)],on='OTU',how='left')
            prime_org=prime_org.merge(prevalences['%s_%s' % (env2,org)],on='OTU',how='left').merge(abundances['%s_%s' % (env2,org)],on='OTU',how='left')
        keystones_org['Keystones in']=env
        keystones_org['Kingdom']=org
        prime_org['Keystones in']=env
        prime_org['Kingdom']=org
        keystones_out[org].append(keystones_org)
        prime_out[org].append(prime_org)


## Plot in boxplots
flierprops = dict(marker='o', markerfacecolor='black', markersize=1,linestyle='none')
palette={'Cropland':'#ffffb3','Grassland':'#ccebc5','Woodland':'#fdc086'}

fig,ax=plt.subplots(2,4,figsize=(12,6),sharey=True)
Nax=0
for org in orgs:
    keystones=pd.concat(keystones_out[org])
    keystones=pd.DataFrame(keystones[['OTU','Keystones in']+[c for c in keystones.columns if 'abundance' in c or 'Prevalence' in c]].set_index(['OTU','Keystones in']).stack()).reset_index(drop=False)
    keystones['Var']=keystones['level_2'].str.split('_').str[0]
    keystones['Samples considered']=keystones['level_2'].str.split('_').str[1]
    keystones=keystones.rename(columns={0:'Value'})

    prime=pd.concat(prime_out[org])
    prime=pd.DataFrame(prime[['OTU','Keystones in']+[c for c in prime.columns if 'abundance' in c or 'Prevalence' in c]].set_index(['OTU','Keystones in']).stack()).reset_index(drop=False)
    prime['Var']=prime['level_2'].str.split('_').str[0]
    prime['Samples considered']=prime['level_2'].str.split('_').str[1]
    prime=prime.rename(columns={0:'Value'})

    sns.boxplot(x='Value',y='Keystones in',hue='Samples considered',data=keystones[keystones['Var']=='Mean relative abundance'],ax=ax[Nax][0],legend=None,flierprops=flierprops,palette=palette,linecolor='black',hue_order=['Cropland','Grassland','Woodland'])
    ax[Nax][0].set_xlabel('Mean relative abundance')
    sns.boxplot(x='Value',y='Keystones in',hue='Samples considered',data=keystones[keystones['Var']=='Prevalence'],ax=ax[Nax][1],legend=None,flierprops=flierprops,palette=palette,linecolor='black',hue_order=['Cropland','Grassland','Woodland'])
    ax[Nax][1].set_xlabel('Prevalence')
    sns.boxplot(x='Value',y='Keystones in',hue='Samples considered',data=prime[prime['Var']=='Mean relative abundance'],ax=ax[Nax][2],legend=None,flierprops=flierprops,palette=palette,linecolor='black',hue_order=['Cropland','Grassland','Woodland'])
    ax[Nax][2].set_xlabel('Mean relative abundance')
    sns.boxplot(x='Value',y='Keystones in',hue='Samples considered',data=prime[prime['Var']=='Prevalence'],ax=ax[Nax][3],legend=None,flierprops=flierprops,palette=palette,linecolor='black',hue_order=['Cropland','Grassland','Woodland'])
    ax[Nax][3].set_xlabel('Prevalence')
    ax[Nax][0].set_xscale('symlog', linthresh=1e-6)
    ax[Nax][0].set_xlim(-1e-7,1e-2)
    ax[Nax][2].set_xscale('symlog', linthresh=1e-6)
    ax[Nax][2].set_xlim(-1e-7,1e-3)
    Nax+=1

plt.savefig('keystone.relativeAbundanceAndPrevalence.pdf')
plt.close()

