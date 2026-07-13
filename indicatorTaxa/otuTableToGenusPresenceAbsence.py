import pandas as pd

# Parse data
otu={'fungi':pd.read_csv('Fungi.OTU_table.csv').set_index('Unnamed: 0').T, 'bacteria':pd.read_csv('Bacteria.OTU_table.csv').set_index('Unnamed: 0').T} # OTU tables containing all samples
env_metadata=pd.read_csv('../Bacteria.design',sep='\t').set_index('group') # TSV file linking each sample to its ecosystem type
otu['fungi']=otu['fungi'][list(env_metadata.index)] # Put fungal data in the same order as bacterial data and metadata
otu['bacteria']=otu['bacteria'][list(env_metadata.index)] # Put bacterial data in the same order as fungal data and metadata

# Parse OTU taxonomy
taxo={'bacteria':pd.read_excel('../../Taxonomy.xlsx',sheet_name=0).set_index('zOTU'),'fungi':pd.read_excel('../../Taxonomy.xlsx',sheet_name=1).set_index('OTU_ID')}
taxo['bacteria']=taxo['bacteria'].rename(columns={'phylum':'Phylum','class':'Class','order':'Order','family':'Family','genus':'Genus','species':'Species'})
for c in taxo['fungi'].columns:
    taxo['fungi'][c]=taxo['fungi'][c].str.split('__').str[1]
taxo_speciesKnown={'bacteria':taxo['bacteria'][taxo['bacteria']['Genus']!='ukn'],'fungi':taxo['fungi'][(taxo['fungi']['Genus']!='.') & (taxo['fungi']['Genus']!='')]}

# Concat
for org in ['bacteria','fungi']:
    taxo_speciesKnown[org]['downToGenus']=taxo_speciesKnown[org]['Phylum']+'__'+taxo_speciesKnown[org]['Class']+'__'+taxo_speciesKnown[org]['Order']+'__'+taxo_speciesKnown[org]['Family']+'__'+taxo_speciesKnown[org]['Genus'] 
    otu_here=otu[org].merge(taxo_speciesKnown[org][['downToGenus']],left_index=True,right_index=True) # Each OTU assigned to its genus
    sumBySpecies=otu_here.groupby('downToGenus').sum() # Sum to the genus level
    sumBySpecies[sumBySpecies>1]=1 # If superior to 1 read, then presence (1), otherwise absence (0)
    sumBySpecies.T.to_csv('genus.presenceAbsence.'+org+'.csv') # Save file for indicator analysis in R