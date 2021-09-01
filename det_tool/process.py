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

import datetime as dt
import rasterio
import geopandas as gpd
from det_tool import functions as fun
from det_tool import classes as cls
from pre_processor import functions as pfun
from config import cityModel_path, foot_height_shp, roof_raster, facade_solar_energy, updated_JSON, report_path


def process(model_path, shp_path, raster_path, facade_file):

    print("Reading config file ...")
    print("Reading and loading input data ...")

    shapefile_foot = gpd.read_file(shp_path)
    data = pfun.read_json(model_path)
    obj_id = pfun.get_obj_ids(data)

    print("3D City model contains {} city objects ...\n".format(len(obj_id)))

    print("Initializing InitialRaster Object Class ...")
    img = rasterio.open(raster_path)
    ini_rast = cls.InitialRaster(img)

    print("Initializing and processing FacadeFile Object Class ...")
    facadeFile = cls.FacadeData(facade_file)
    facadeFile.process_input()

    print("Initializing and processing Voxel Object Class ...")
    voxels = []
    i = 0
    for voxel_entry in facadeFile.file_data:
        voxel = cls.Voxel(voxel_entry['column'], voxel_entry['row'], voxel_entry['data'], i)
        voxel.define_coordinates(ini_rast)
        i += 1
        voxels.append(voxel)

    print("Initializing and processing Building Object Classes ...\n")
    build_obj_list = []
    build_irr_values = []
    stat_data = {}
    i = 0
    for key in obj_id:
        b_start = dt.datetime.now()
        print("Building {} in being processed ...".format(key))

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
                #build_irr_values.append(win_obj.total_irradiance)
        stat_data['{}'.format(key)] = {'win_assigned': len(building.windows_geometry)}

        roof = cls.Roof(key)
        fun.process_roof(roof, data['CityObjects'][key], data, ini_rast)
        building.roof = roof
        end_time = dt.datetime.now()
        dif = end_time - b_start
        print("The processing time of building {a} is: {b}\n".format(a=key, b=end_time - b_start))
        stat_data['{}'.format(key)]['proc_time'] = float(dif.total_seconds())
        build_obj_list.append(building)

        uni_ver = []
        bound = data['CityObjects'][key]['geometry'][0]['boundaries'][0]
        for ver_list in bound:
            for element in ver_list:
                if element not in uni_ver:
                    uni_ver.append(element)
        stat_data['{}'.format(key)]['n_vertices'] = len(uni_ver)
        i += 1

    for building in build_obj_list:
        for window in building.windows_geometry:
            build_irr_values.append(window.total_irradiance)
    return build_obj_list, build_irr_values


def main():
    start = dt.datetime.now()
    build_obj_list, irr_values = process(cityModel_path, foot_height_shp, roof_raster, facade_solar_energy)
    process_end = dt.datetime.now()
    fun.update_CityModel(build_obj_list, cityModel_path, updated_JSON, irr_values)
    fun.generate_report(build_obj_list, report_path)
    print("Done!\nProcessing time of the whole dataset: {a}\n"
          "The running time of the tool, including the user interaction: {b}"
          .format(a=process_end-start, b=dt.datetime.now() - start))