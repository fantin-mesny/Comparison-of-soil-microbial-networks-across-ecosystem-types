import pandas as pd
import argparse
import sys
import numpy as np
import os
import math


filterOutRareOTUs=True
samplesPerRandom={}

for env in ['Cropland','Grassland','Woodland']:
    soilProp=pd.read_csv(env+'.soilProp.tsv',sep='\t').set_index('SampleID')
    for org in ['Bacteria','Fungi']:

        # Parse the OTU table of the ecosystem type:
        otu=pd.read_csv('%s_%s.tsv' % (env,org),sep='\t')
        otu=otu.rename(columns={'Unnamed: 0':'#OTU ID'}).set_index('#OTU ID')
        otu=otu.dropna()
        otu=otu.reindex(sorted(otu.columns), axis=1)
        otu=otu[list(soilProp.index)] # samples in the same order as in the metadata table
        
        ## Select 50 random samples across all ecosystem types
        if len(samplesPerRandom)==0: # if org==Fungi, then this has already been done for bacteria, don't redo, use the same sets of samples
            for i in range(0,50): # 50 random selection...
                samplesPerRandom[i]=list(pd.Series(otu.columns).sample(n=50,random_state=i)) # ... of 50 samples

    

        for i in range(0,50): # for each of the 50 subset
            otu_sub50=otu[samplesPerRandom[i]] # subsample the OTU table to only keep the 50 samples
            
            if filterOutRareOTUs: # prevalence filter: drop OTUs present in less than 10% (5) samples
                threshold=5
                #print('Before filtering:', len(otu_sub50), 'OTUs in dataframe')
                toFilterOut=[o for o in otu_sub50.index if sum(otu_sub50.loc[o]>0)<threshold]
                otu_sub50=otu_sub50[~(otu_sub50.index.isin(toFilterOut))]
                #print('Removing '+str(len(toFilterOut))+' OTUs from table... Left: '+str(len(otu_sub50)))
            otu_sub50.T.to_csv(env+'_'+org+'_'+str(i)+'.tsv',sep='\t') # Export the subsampled OTU table...
            soilProp.reindex(index=samplesPerRandom[i]).to_csv(env+'_'+org+'_'+str(i)+'.soilProp.tsv',sep='\t') #... and a matching metadata sample




