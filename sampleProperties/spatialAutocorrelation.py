import pandas as pd
import numpy as np
from esda.join_counts import Join_Counts
import geopandas as gpd
from libpysal.weights import KNN
from esda.moran import Moran

## Parse latitudes and longitudes:
df=pd.read_csv('sampleLocations.csv').set_index('SampleID')
gdf = gpd.GeoDataFrame(
    df,
    geometry=gpd.points_from_xy(df.TH_LONG, df.TH_LAT),
    crs="EPSG:4326"
)
gdf = gdf.to_crs(epsg=3857)
w = KNN.from_dataframe(gdf, k=5)
w.transform = "r"  # row-standardized
y = gdf["Ecosystem type"].map({'Cropland':3,'Grassland':2,'Woodland':1}) # associate ecosystem types a number reflecting disturbance (needed for Moran's test)

## Calculate Moran test:
moran = Moran(y, w)
print("Moran's I:", moran.I)
print("Expected I:",moran.EI)
print('Standard deviation of I:',moran.seI_norm)
print("p-value:", moran.p_sim)

## Calculate individual Join-count tests (better for categorical labels):
gdf["Cropland"] = (gdf["Ecosystem type"] == 'Cropland').astype(int)
gdf["Grassland"] = (gdf["Ecosystem type"] == 'Grassland').astype(int)
gdf["Woodland"] = (gdf["Ecosystem type"] == 'Woodland').astype(int)
for col in ["Cropland", "Grassland", "Woodland"]:
    jc = Join_Counts(gdf[col].values, w)
    print(f"{col}: BB={jc.bb}, chi2={jc.chi2}, p={jc.p_sim_bb}")
