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
from create_raster import functions as crfun
from pre_processor import functions as pfun
from extract_window import win_functions as wfun
from shapely.geometry import Polygon, Point, LineString
import geopandas as gpd
import matplotlib.pyplot as plt

shapefile_foot = gpd.read_file("data2/ft_height.shp")
img = rasterio.open('data2/Energyyearroof.tif')
ini_rast = crfun.InitialRaster(img)
#show(img)
path = 'data2/Energyyearwall.txt'
js_path = 'data2/win_lund3008_red.json'

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
wall_rasters = []
i = 0
windows = []
for key in obj_id:
    x = crfun.extract_footprints(key, shapefile_foot)
    building = crfun.BuildingObject(key, i, x, int(shapefile_foot.height[i]))
    building.getWallList()
    building.assignBuildingIndex(voxels)
    building.createCornerAreas()
    building.assignWallIndex()
    building.setVoxelCorners()
    build_obj_list.append(building)
    win_geom = wfun.get_win_from_json(data['CityObjects'][key], data)
    for wall in building.wall_list:
        wallRaster = crfun.WallRasterObject(building.building_id, wall['wall_id'], wall['wall_start'], wall['wall_end'], building.building_height)
        wallRaster.getRasterDimension()
        wallRaster.getWallVoxelList(building)
        wallRaster.calculateLineCoords()
        wallRaster.populateGrid()
        wallRaster.interpolate()
        wall_rasters.append(wallRaster)
    j = 0
    for window in win_geom:
        win_obj = wfun.Window(window, key, j)
        win_obj.set_extent()
        for wall in building.wall_list:
            wall_win_buff = crfun.create_buffer(wall['wall_start'], wall['wall_end'], 1.5)
            win_obj.assign_wall(wall, wall_win_buff)
            building.windows_geometry.append(win_obj)
            windows.append(win_obj)
        j += 1
    i += 1



import pandas
import geopandas

outfp = 'data2/wall_window_id.shp'

y = []
x = []
id = []
bid = []
vid = []

for window in windows:
    if window.valid :
        x.append(window.edge1[0])
        x.append(window.edge2[0])
        y.append(window.edge1[1])
        y.append(window.edge2[1])
        id.append(window.id)
        id.append(window.id)
        bid.append(window.bid)
        bid.append(window.bid)
        vid.append(window.wall_id)
        vid.append(window.wall_id)


'''
for voxel in voxels:
    #if voxel.wall_start_corner or voxel.wall_end_corner:
        x.append(voxel.coord_x)
        y.append(voxel.coord_y)
        id.append(voxel.wall_index)
        bid.append(voxel.building_index)
        vid.append(voxel.voxel_id)
'''

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

