'''
***********************************************************************************************************************
*                                 create raster module main script                                                    *
*                                         Created by P.Nezval                                                         *
*                                             version 0.1                                                             *
*                                               May 2021                                                              *
***********************************************************************************************************************
                                                License:
                                    Mozilla Public License Version 2.0
                                    ==================================
'''

import os
import rasterio
from rasterio.plot import show
from create_raster import functions as crfun
from pre_processor import functions as pfun
from shapely.geometry import Polygon, Point, LineString
import geopandas as gpd
import matplotlib.pyplot as plt

shapefile_foot = gpd.read_file("data/foot_print_simplified.shp")
img = rasterio.open('data/Energyyearroof.tif')
ini_rast = crfun.InitialRaster(img)
#show(img)
path = 'data/Energyyearwall.txt'
js_path = 'data/study_area.json'

data = pfun.read_json(js_path)
obj_id = pfun.get_obj_ids(data)

wallFile = crfun.WallDataObject(path)
wallFile.processInput()

voxels = []
i = 0
for voxel_entry in wallFile.file_data:
    voxel = crfun.VoxelObject(voxel_entry['column'], voxel_entry['row'], voxel_entry['data'], i)
    voxel.defineCoordinates(ini_rast)
    i += 1
    voxels.append(voxel)

build_obj_list = []
i = 0
for key in obj_id:
    x = crfun.extract_footprints(key, shapefile_foot)
    building = crfun.BuildingObject(key, i, x, shapefile_foot.height[i])
    building.getWallList()
    building.assignBuildingIndex(voxels)
    building.createCornerAreas()
    building.assignWallIndex()
    building.setVoxelCorners()
    build_obj_list.append(building)
    i += 1



'''

for building in build_obj_list:
    building.assignWallIndex(voxels)

'''

import pandas
import geopandas

outfp = 'data/voxel_wall_id_all7.shp'

y = []
x = []
id = []
bid = []
vid = []
for voxel in voxels:
    if voxel.wall_start_corner or voxel.wall_end_corner:
        x.append(voxel.coord_x)
        y.append(voxel.coord_y)
        id.append(voxel.wall_index)
        bid.append(voxel.building_index)
        vid.append(voxel.voxel_id)


df = pandas.DataFrame({'id':id, 'bid':bid, 'vid':vid, 'x':x, 'y':y})
crs = {'init': 'epsg:3008'}
gdf = geopandas.GeoDataFrame(df,geometry=geopandas.points_from_xy(df.x, df.y))
gdf.to_file(outfp)

outfp = 'data/corner_areas2.shp'

y = []
x = []
id = []
bid = []
vid = []
for building in build_obj_list:
    i = 0
    for geom in building.corner_area.geoms:
        id.append(i)
        bid.append(building.building_id)
        x, y = geom.exterior.xy
        vid.append(Polygon(zip(x,y)))
        i += 1

df = pandas.DataFrame({'id':id, 'bid':bid})
gdf = geopandas.GeoDataFrame(df, geometry=vid)
gdf.to_file(outfp)

