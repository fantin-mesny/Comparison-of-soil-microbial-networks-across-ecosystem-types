import pandas as pd
import argparse
import sys
import numpy as np
import os
import math


filterOutRareOTUs=True
samplesPerRandom={}
soilProp=pd.read_csv('sampleProperties.tsv',sep='\t').set_index('SampleID') # Sample edaphic and climatic variables used in input of FlashWeave (Supplementary Table 1)

for org in ['Bacteria','Fungi']:
    # parse OTU table:
    otuFile=org+'.OTU_table.csv'
    otu=pd.read_csv(otuFile) 
    otu=otu.rename(columns={'Unnamed: 0':'#OTU ID'}).set_index('#OTU ID')
    otu=otu.dropna()
    otu=otu.reindex(sorted(otu.columns), axis=1)
    otu=otu[list(soilProp.index)]
        

    ## Select 50 random samples across all ecosystem types
    if len(samplesPerRandom)==0: # if org==Fungi, then this has already been done for bacteria, don't redo, use the same sets of samples
        for i in range(0,50): # 50 random selection...
            samplesPerRandom[i]=list(pd.Series(otu.columns).sample(n=50,random_state=i)) # ... of 50 samples

    
    for i in range(0,50):
        otu_sub50=otu[samplesPerRandom[i]] #Take each of the random subsets of 50 samples
        for col in otu_sub50.columns: # for each sample...
            otu_sub50[col]=otu_sub50[col].sample(frac=1) # ...shuffle the read counts across OTUs within each column (i.e. sample)

        if filterOutRareOTUs: # prevalence filter: drop OTUs present in less than 10% (5) samples
            print('Before filtering:', len(otu_sub50), 'OTUs in dataframe')
            threshold=5 # 
            toFilterOut=[o for o in otu_sub50.index if sum(otu_sub50.loc[o]>0)<threshold]
            otu_sub50=otu_sub50[~(otu_sub50.index.isin(toFilterOut))]
            #print('Removing '+str(len(toFilterOut))+' OTUs from table... Left: '+str(len(otu_sub50)))

        otu_sub50.T.to_csv('Null_'+org+'_'+str(i)+'.tsv',sep='\t') # Export the generated null OTU tables...
        soilProp.reindex(index=samplesPerRandom[i]).to_csv('null_new/Null_'+org+'_'+str(i)+'.soilProp.tsv',sep='\t') #... and a matching metadata sample




