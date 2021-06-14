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
    x = crfun.extract_footprints(data, key)
    c_obj = data['CityObjects'][key]
    ver_list = pfun.get_vertex_list(c_obj)
    ver_heights = pfun.vertex_heights(data)
    build_height = pfun.get_build_height(ver_list, ver_heights)
    i += 1
    building = crfun.BuildingObject(key, i, x, build_height)
    building.getWallList()
    building.assignBuildingIndex(voxels)
    building.assignWallIndex(voxels)
    build_obj_list.append(building)



'''
budova = build_obj_list[4]
voxel_list1 = []
for voxel in voxels:
    if voxel.building_index == budova.building_id:
        voxel_list1.append(voxel)

stena = budova.wall_list[1]
line = LineString([Point(stena['wall_start'][0], stena['wall_start'][1]), Point(stena['wall_end'][0], stena['wall_end'][1])])
line = line.buffer(1.5)
vox_close_wall = []
for voxel in voxel_list1:
    p1 = Point(voxel.coord_x, voxel.coord_y)
    if p1.within(line):
        vox_close_wall.append(voxel)

start = [stena['wall_start'][0], stena['wall_start'][1]]
end = [stena['wall_end'][0], stena['wall_end'][1]]

correct_vox = []
for voxel in vox_close_wall:
    point = [voxel.coord_x, voxel.coord_y]
    if budova.correctSide(start, end, point):
        correct_vox.append(voxel)
        voxel.wall_index = stena['wall_id']

import pandas
import geopandas

outfp = 'data/voxel_wall_id6.shp'

y = []
x = []
id = []
for voxel in voxel_list1:
    #if voxel.building_index == 'none':
    x.append(voxel.coord_x)
    y.append(voxel.coord_y)
    id.append(voxel.wall_index)

df = pandas.DataFrame({'id':id, 'x':x, 'y':y})
gdf = geopandas.GeoDataFrame(df, geometry=geopandas.points_from_xy(df.x, df.y))
gdf.to_file(outfp)


'''

for building in build_obj_list:
    building.assignWallIndex(voxels)


import pandas
import geopandas

outfp = 'data/voxel_wall_id_all.shp'

y = []
x = []
id = []
bid = []
for voxel in voxels:
    if voxel.wall_index != 999:
        x.append(voxel.coord_x)
        y.append(voxel.coord_y)
        id.append(voxel.wall_index)
        bid.append(voxel.building_index)


df = pandas.DataFrame({'id':id, 'bid':bid, 'x':x, 'y':y})
gdf = geopandas.GeoDataFrame(df, geometry=geopandas.points_from_xy(df.x, df.y))
gdf.to_file(outfp)