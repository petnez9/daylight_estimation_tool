'''
***********************************************************************************************************************
*                                Main process file daylight estimation tool                                           *
*                                         Created by P.Nezval                                                         *
*                                             version 0.1                                                             *
*                                         May 2021 - August 2021                                                      *
***********************************************************************************************************************
                                                License:
                                               MIT License
                                    ==================================
'''


import rasterio
import geopandas as gpd
from det_tool import functions as fun
from det_tool import classes as cls
from pre_processor import functions as pfun
from config import cityModel_path, foot_height_shp, roof_raster, facade_solar_energy, updated_JSON, report_path


def process(model_path, shp_path, raster_path, facade_file):

    print("Reading config file ...")
    print("Reading input data ...")
    print("Loading input data ...")

    shapefile_foot = gpd.read_file(shp_path)
    data = pfun.read_json(model_path)
    obj_id = pfun.get_obj_ids(data)

    print("3D City model contains {} city objects ...\n".format(len(obj_id)))

    print("Initializing InitialRaster Object Class ...")
    img = rasterio.open(raster_path)
    ini_rast = cls.InitialRaster(img)

    print("Initializing and processing FacadeFile Object Class ...")
    facadeFile = cls.FacadeData(facade_file)
    facadeFile.processInput()

    print("Initializing and processing Voxel Object Class ...")
    voxels = []
    i = 0
    for voxel_entry in facadeFile.file_data:
        voxel = cls.Voxel(voxel_entry['column'], voxel_entry['row'], voxel_entry['data'], i)
        voxel.defineCoordinates(ini_rast)
        i += 1
        voxels.append(voxel)

    print("Initializing and processing Building Object Classes ...")
    build_obj_list = []
    i = 0
    for key in obj_id:
        print("Building {} in being processed ...\n".format(key))
        building = cls.Building(key, i, fun.extract_footprints(key, shapefile_foot), int(shapefile_foot.height[i]), fun.get_surfaces_info(data, key))
        fun.process_building(building, voxels)
        for facade in building.facade_list:
            facadeRaster = cls.FacadeRaster(facade.feature_id, facade.height, facade.start, facade.end, facade.facade_index)
            fun.process_fac_raster(facadeRaster, building, facade)
        win_geom, wind_bound_id = fun.extract_window_geometry(data['CityObjects'][key], data)
        for n in range(len(win_geom)):
            win_obj = cls.Window(win_geom[n], key, n, wind_bound_id[n])
            for facade in building.facade_list:
                fun.process_window(win_obj, facade, building)
                facade.getFacadeIrradiance()
        roof = cls.Roof(key)
        fun.process_roof(roof, data['CityObjects'][key], data, img)
        building.roof = roof
        build_obj_list.append(building)
        i += 1
    build_irr_values = []
    for building in build_obj_list:
        for window in building.windows_geometry:
            build_irr_values.append(window.total_irradiance)
    return build_obj_list, build_irr_values


def main():
    build_obj_list, irr_values = process(cityModel_path, foot_height_shp, roof_raster, facade_solar_energy)
    fun.update_CityJSON(build_obj_list, cityModel_path, updated_JSON, irr_values)
    fun.generate_report(build_obj_list, report_path)
    print("Done!")