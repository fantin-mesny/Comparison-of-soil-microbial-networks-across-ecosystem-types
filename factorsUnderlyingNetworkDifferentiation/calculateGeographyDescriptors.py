import numpy as np
from scipy.spatial.distance import pdist
from shapely.geometry import MultiPoint
import pandas as pd
import geopandas as gpd

def spatial_descriptors(coords): # Function that calculate geography indicators for a set of locations (latitude/longitude coordinates)
    centroid = coords.mean(axis=0)
    area = MultiPoint(coords).convex_hull.area
    mean_dist = pdist(coords).mean()
    return {
        "centroid_x": centroid[0],
        "centroid_y": centroid[1],
        "convex_hull_area": area,
        "mean_pairwise_distance": mean_dist
    }


df=pd.read_csv('sampleLocations.csv').set_index('SampleID')

data={}
for env in ['Cropland','Grassland','Woodland']:
    for org in ['Bacteria','Fungi']:
        for rep in range(0,50):
            samplesInSubset=list(pd.read_csv('%s_%s_%s.soilProp.tsv' % (env,org,str(rep)),sep='\t').set_index('SampleID').index)
            df_subset=df[df.index.isin(samplesInSubset)]
            gdf = gpd.GeoDataFrame(
                df_subset,
                geometry=gpd.points_from_xy(df_subset.TH_LONG, df_subset.TH_LAT),
                crs="EPSG:4326"
            )
            gdf = gdf.to_crs(epsg=3857)  # Web Mercator
            coords = np.column_stack([gdf.geometry.x, gdf.geometry.y])
            data['%s_%s_%s' % (org[0],str(rep),env)]=spatial_descriptors(coords)

data=pd.DataFrame(data)
data.T.to_csv('geography.csv')