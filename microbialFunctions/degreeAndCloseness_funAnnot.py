import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
from scipy import stats
if __name__ == '__main__':

    # Parse bacterial and fungal functional annotations
    annot=pd.read_excel('./functionalAnnotations/fungi.annot.xlsx')
    annot_dict=annot.set_index('OTU')['Primary lifestyle'].to_dict()
    annot_b=pd.read_csv('./functionalAnnotations/bacteria.annot.csv').set_index('Unnamed: 0')

    # Parse data obtained with networkProperties/exportNodeDegreeClosenessBetweenness.py
    degree0=pd.read_csv('degrees.csv')
    betweenness0=pd.read_csv('betweenness.csv')
    closeness0=pd.read_csv('closeness.csv')

    envs=['Cropland','Grassland','Woodland']

    Network_stats=pd.DataFrame(index=envs,columns=['Environment','Percentage of negative interactions','Network density','Clustering coefficient'])
    output=[]
    for env in envs: # iterate over multikingdom networks
        netw=env+'_multikingdom.edgelist'

        ## subset dataframes to only keep the considered network
        degree=degree0[degree0["Network"]==netw]
        betweenness=betweenness0[betweenness0["Network"]==netw]
        closeness=closeness0[closeness0["Network"]==netw]
        all_dat=degree.merge(betweenness,on="OTU")
        all_dat=all_dat.merge(closeness,on="OTU")[['OTU','Degree','Betweenness','Closeness']]
        all_dat['Group']=all_dat['OTU'].str.startswith('zot')
        all_dat['Group']=all_dat['Group'].map({True:'Bacteria',False:'Fungi'})
        all_dat['Environment']=env
        all_dat['Degree pct'] = all_dat['Degree'].rank(pct=True) # percentile rank transformation
        all_dat['Closeness pct'] = all_dat['Closeness'].rank(pct=True) # percentile rank transformation
        output.append(all_dat)
    output=pd.concat(output)


    palette={'Cropland':'#ffffb3','Grassland':'#ccebc5','Woodland':'#fdc086'}
    flierprops = dict(marker='o', markerfacecolor='black', markersize=0.5,linestyle='none')

    ## Plot fungal data:
    output_fungi=output[(output['Group']=='Fungi')]
    output_fungi['Lifestyle']=output_fungi['OTU'].map(annot_dict) # assign fungal nodes to their primary lifestyle
    fig,ax=plt.subplots(1,2,figsize=(7,10),sharey=True)

    top_ls=['soil_saprotroph','litter_saprotroph','wood_saprotroph','dung_saprotroph','pollen_saprotroph','unspecified_saprotroph','arbuscular_mycorrhizal','ectomycorrhizal','plant_pathogen','root_endophyte','mycoparasite','animal_parasite']
    output_fungi=output_fungi[output_fungi['Lifestyle'].isin(top_ls)]
    output_fungi['order']=output_fungi['Lifestyle'].map({l:top_ls.index(l) for l in top_ls})
    output_fungi=output_fungi.sort_values(by='order',ascending=True)

    sns.boxplot(hue='Environment',x='Degree pct',y='Lifestyle', data=output_fungi.dropna(),ax=ax[0],hue_order=envs,flierprops=flierprops,palette=palette,linecolor='black')
    sns.boxplot(hue='Environment',x='Closeness pct',y='Lifestyle', data=output_fungi.dropna(),ax=ax[1],hue_order=envs,flierprops=flierprops,palette=palette,linecolor='black')
    ax[0].set_xlim(-0.025,1.025)
    ax[1].set_xlim(-0.025,1.025)
    plt.tight_layout()
    plt.savefig('fungalLifestyles_centrality.pdf')
    plt.close()


    ## Plot bacterial data:
    output_bac=[]
    top_ls=['chemoheterotrophy','aerobic_chemoheterotrophy','anaerobic_chemoheterotrophy','fermentation','nitrate_reduction','ureolysis']
    for func in top_ls:
        output_bac.append(output[output['OTU'].isin(annot_b[annot_b[func]>0].index)]) ## assigns bacterial OTUs to their functional group
        output_bac[-1]['Lifestyle']=func
    output_bac=pd.concat(output_bac)

    flierprops = dict(marker='o', markerfacecolor='black', markersize=0.5,linestyle='none')
    fig,ax=plt.subplots(1,2,figsize=(7,5),sharey=True)
    output_bac=output_bac[output_bac['Lifestyle'].isin(top_ls)]
    sns.boxplot(hue='Environment',x='Degree pct',y='Lifestyle', data=output_bac,ax=ax[0],hue_order=envs,flierprops=flierprops,palette=palette,linecolor='black')
    sns.boxplot(hue='Environment',x='Closeness pct',y='Lifestyle', data=output_bac,ax=ax[1],hue_order=envs,flierprops=flierprops,palette=palette,linecolor='black')
    ax[0].set_xlim(-0.025,1.025)
    ax[1].set_xlim(-0.025,1.025)
    plt.tight_layout()
    plt.savefig('bacterialLifestyles_centrality.pdf')
    plt.close()

