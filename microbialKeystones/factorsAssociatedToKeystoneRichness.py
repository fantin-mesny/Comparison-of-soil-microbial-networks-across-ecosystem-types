import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np 
from statsmodels.stats import multitest
import statsmodels.api as sm
from sklearn.metrics import r2_score


def getRA(OTUtable): #calculates relative abundances of taxa in an OTU table
    funAbundance=OTUtable.T
    funAbundance=funAbundance.div(funAbundance.sum())
    return funAbundance

## Parse data obtained with networkProperties/exportNodeDegreeClosenessBetweenness.py
degree0=pd.read_csv('degrees.csv')
betweenness0=pd.read_csv('betweenness.csv')
closeness0=pd.read_csv('closeness.csv')

## Parse taxonomy of OTUs
funTax=pd.read_excel('Taxonomy.xlsx',sheet_name=1).set_index('OTU_ID')
for c in funTax.columns:
    funTax[c]=funTax[c].str.split('__').str[1]
bacTax=pd.read_excel('Taxonomy.xlsx',sheet_name=0).set_index('zOTU')
bacTax=bacTax.rename(columns={'phylum':'Phylum','class':'Class','order':'Order','family':'Family','genus':'Genus','species':'Species'})


hubs_output=[]
primeHubs_output=[]

statCor_df=[]
all_pvalues=[]

envs=['Cropland','Grassland','Woodland']
for env in envs: # to iterate over multikingdom networks
    netw=env+'_multikingdom.edgelist'

    ## Subset and merge dataframes for the considered network
    degree=degree0[degree0["Network"]==netw]
    betweenness=betweenness0[betweenness0["Network"]==netw]
    closeness=closeness0[closeness0["Network"]==netw]
    all_dat=degree.merge(betweenness,on='OTU')
    all_dat=all_dat.merge(closeness,on='OTU')[['OTU','Degree','Betweenness','Closeness']]

    ## identify keystones and prime keystones
    thresh=0.99
    hubs=all_dat[(all_dat['Degree']>=all_dat['Degree'].quantile(thresh)) & (all_dat['Closeness']>=all_dat['Closeness'].quantile(thresh))].set_index('OTU')
    primeHubs=all_dat[(all_dat['Degree']>=all_dat['Degree'].quantile(thresh)) & (all_dat['Closeness']>=all_dat['Closeness'].quantile(thresh)) & (all_dat['Betweenness']<=hubs['Betweenness'].quantile(0.1))].set_index('OTU')
    hubs['Threshold']=thresh
    hubs['Network']=netw
    hubs['Environment']=env
    hubs['Organism']=hubs.index.str.startswith('zot')
    hubs['Organism']=hubs['Organism'].map({True:'Bacteria',False:'Fungi'})
    hubs_output.append(hubs)
    primeHubs['Threshold']=thresh
    primeHubs['Network']=netw
    primeHubs['Environment']=env
    primeHubs_output.append(primeHubs)

    ## Parse environmental variables and OTU tables for the ecosystem type
    envVars=pd.read_csv('metadata/%s.soilProp.tsv' % (env),sep='\t').set_index('SampleID')
    envVars_vars=list(envVars.columns)
    botu=pd.read_csv('OTU_tables/%s_Bacteria.tsv' % (env),sep='\t').set_index('Unnamed: 0')
    fotu=pd.read_csv('OTU_tables/%s_Fungi.tsv' % (env),sep='\t').set_index('Unnamed: 0')

    ## Calculate relative abundances and get presence absence of keystones in each sample
    #Bacteria:
    botu=getRA(botu)
    botu_hubs=botu[botu.index.isin(hubs.index)].T
    bpa=(botu_hubs>0).astype(int) # Presence/Absence table: if>0 then present, otherwise absent
    bpa_dic=bpa.sum(axis=1).to_dict()
    envVars['Nbac_keystones']=envVars.index.map(bpa_dic) # add a keystone richness column in the same dataframe as environmental variables
    #Fungi:
    fotu=getRA(fotu)
    fotu_hubs=fotu[fotu.index.isin(hubs.index)].T
    fpa=(fotu_hubs>0).astype(int)# Presence/Absence table: if>0 then present, otherwise absent
    fpa_dic=fpa.sum(axis=1).to_dict()
    envVars['Nfun_keystones']=envVars.index.map(fpa_dic) # add a keystone richness column in the same dataframe as environmental variables

    ## Calculate polynomial regression between the numnber of keystones and each environmental variable
    statCor=pd.DataFrame(columns=envVars_vars)
    for rich in ['Nbac_keystones','Nfun_keystones']:
        for var in envVars_vars:
            x=envVars[var]
            y=envVars[rich]
            X = np.column_stack([x**i for i in range(4)])  # Specify the polynomial regression and its degree
            X = sm.add_constant(X)  # Optional: adds a constant term if not already included
            model = sm.OLS(y, X).fit() # Fit model
            statCor.loc['r2_%s_%s' % (env,rich[1:4]),var]=model.rsquared # R2
            statCor.loc['p_%s_%s' % (env,rich[1:4]),var]=model.f_pvalue # P-value of the model
            all_pvalues.append(model.f_pvalue)
    statCor_df.append(statCor)

    ## Null models: random sets of OTUs are picked and their richness tested for association to environmental variables
    null_R2={'Bacteria':[],'Fungi':[]}
    null_F={'Bacteria':[],'Fungi':[]}
    for org in ['Bacteria','Fungi']:
        for nul in range(9999): # 9999 sets of OTUs (same number as bacterial|fungal keystones)
            if org=='Bacteria':
                otu_nullSample=botu.sample(n=len(hubs[hubs['Organism']==org]),random_state=nul).T
            elif org=='Fungi':
                otu_nullSample=fotu.sample(n=len(hubs[hubs['Organism']==org]),random_state=nul).T
            pa=(otu_nullSample>0).astype(int) # presence/absence of OTUs of the random set
            pa_dic=pa.sum(axis=1).to_dict()
            envVars['Nnull']=envVars.index.map(pa_dic)
            envVars=envVars.sort_values(by='pH_CaCl2')
            x=envVars['pH_CaCl2']
            y=envVars['Nnull']
            X = np.column_stack([x**i for i in range(4)])  #  Specify the polynomial regression and its degree
            X = sm.add_constant(X)  # Optional: adds a constant term if not already included
            model = sm.OLS(y, X).fit()# Fit model
            null_R2[org].append(model.rsquared) # R2
            null_F[org].append(model.fvalue) # P-value of the model

    ## For pH, fit the model and plot the polynomial regression
    #Bacteria: 
    envVars=envVars.sort_values(by='pH_CaCl2')
    x=envVars['pH_CaCl2']
    y_bac=envVars['Nbac_keystones']
    X = np.column_stack([x**i for i in range(4)])   #  Specify the polynomial regression and its degree
    X = sm.add_constant(X)  # Optional: adds a constant term if not already included
    model_bac = sm.OLS(y_bac, X).fit() # Fit model
    R2_bac=model_bac.rsquared
    fvalue_bac=model_bac.fvalue
    #Fungi:
    y_fun=envVars['Nfun_keystones']
    model_fun = sm.OLS(y_fun, X).fit() # Fit model
    R2_fun=model_fun.rsquared
    fvalue_fun=model_fun.fvalue
    #print(model.summary()) # Show summary including coefficients, standard errors, t-values, and p-values
    x_smooth = np.linspace(x.min(), x.max(), 200) # Make a smooth range for plotting the curve
    X_smooth = np.column_stack([x_smooth**i for i in range(4)])
    X_smooth = sm.add_constant(X_smooth)
    y_pred_bac = model_bac.predict(X_smooth)
    y_pred_fun = model_fun.predict(X_smooth)


    ## Plot on a figure
    fig,ax=plt.subplots(2,2,figsize=(10,6))
    # Data and polynomial regressions
    sns.scatterplot(x='pH_CaCl2',y='Nbac_keystones',data=envVars,color='grey',ax=ax[0][0]) # bacteria
    sns.scatterplot(x='pH_CaCl2',y='Nfun_keystones',data=envVars,color='grey',ax=ax[1][0]) # fungi
    sns.lineplot(x=x_smooth, y=y_pred_bac, color='blue', linewidth=2.5,ax=ax[0][0]) # regression for bacteria
    sns.lineplot(x=x_smooth, y=y_pred_fun, color='blue', linewidth=2.5,ax=ax[1][0]) # regression for fungi
    ax[0][0].set_xlim(2.5,8.1)
    ax[1][0].set_xlim(2.5,8.1)
    # Histograms showing the result of the null-model
    sns.histplot(x=null_R2['Bacteria'],ax=ax[0][1],bins=100,color='lightgrey')
    ax[0][1].axvline(R2_bac,color='red') # red vertical bar showing the observed R2 for bacteria
    ax[0][1].axvline(pd.Series(null_R2['Bacteria']).quantile(0.95),color='black') # black bar showing the 95th percentile of the null model distribution
    sns.histplot(x=null_R2['Fungi'],ax=ax[1][1],bins=100,color='lightgrey')
    ax[1][1].axvline(R2_fun,color='red') # red vertical bar showing the observed R2 for fungi
    ax[1][1].axvline(pd.Series(null_R2['Fungi']).quantile(0.95),color='black') # black bar showing the 95th percentile of the null model distribution
    ax[0][1].set_xlim(0,0.75)
    ax[1][1].set_xlim(0,0.75)
    plt.savefig(env+'.keystone_richness.pdf')
    plt.close()



## Plot a general figure showing all models linking keystone richness to individual environmental variables
fdr=multitest.multipletests(all_pvalues)[1] # Benjamini-Hochberg correction for multiple testing
fdr_mapper = dict(zip(all_pvalues, fdr))
statCor_df=pd.concat(statCor_df).T.reset_index(drop=False)
for col in statCor_df.columns:
    if col.startswith('p_'):
        statCor_df[col.replace('p_','fdr_')]=statCor_df[col].map(fdr_mapper)
norm = plt.Normalize(0, 1) # setting up a colormap
sm = plt.cm.ScalarMappable(cmap="Reds", norm=norm)
sm.set_array([])
r_cols=sorted([col for col in statCor_df.columns if col.startswith('r2_')])
fig,ax=plt.subplots(1,1,figsize=(20,9))
N=0
for col in r_cols:
    sns.scatterplot(ax=ax,data=statCor_df,y=-N,x='index',size=col,hue=col,palette='seismic',sizes=(10,1500),hue_norm=(-1,1))
    N+=1
ax.get_legend().remove()
ax.figure.colorbar(sm,ax=ax)
ax.set_ylim(-6,1)
plt.savefig('regressionSignificantce.keystone_richness.pdf')
plt.close()
