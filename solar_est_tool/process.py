'''
***********************************************************************************************************************
*                                Main process file solar estimation tool                                              *
*                                         Created by P.Nezval                                                         *
*                                             version 0.1                                                             *
*                                         May 2021 - August 2021                                                      *
***********************************************************************************************************************
                                                License:
                                               MIT License
                                    ==================================
'''

import json
import rasterio
from solar_est_tool import functions as fun
from solar_est_tool import classes as cls
from pre_processor import functions as pfun
from CityJSON_gen_upd import classes as ct_cls
import geopandas as gpd
from config import cityModel_path, foot_height_shp, roof_raster, facade_solar_energy, updated_JSON


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
        building = cls.Building(key, i, fun.extract_footprints(key, shapefile_foot), int(shapefile_foot.height[i]))
        fun.process_building(building, voxels)
        roof = cls.Roof(key)
        fun.process_roof(roof, data['CityObjects'][key], data, img)
        building.roof = roof
        for facade in building.facade_list:
            facadeRaster = cls.FacadeRaster(building.fid, building.building_height, facade.start, facade.end, facade.facade_index)
            fun.process_fac_raster(facadeRaster, building, facade)
        win_geom, wind_bound_id = fun.get_win_from_json(data['CityObjects'][key], data)
        for n in range(len(win_geom)):
            win_obj = cls.Window(win_geom[n], key, n, wind_bound_id[n])
            for facade in building.facade_list:
                win_obj.set_extent(facade)
                facade_win_buff = fun.create_buffer(facade.start, facade.end, 1)
                win_obj.assign_facade(facade, facade_win_buff)
                if win_obj not in building.windows_geometry:
                    building.windows_geometry.append(win_obj)
                facade.getFacadeWindows(building)
                facade.getWindowIrradiance()
        build_obj_list.append(building)
        i += 1
    return build_obj_list


def update_CityJSON(build_obj_list, model_path, output_path):
    print("Data are being written to CityJSON ...")
    data = pfun.read_json(model_path)
    nm = ['irradiation1', 'irradiation2', 'irradiation3', 'irradiation4']
    cl = [[0, 0, 0.5], [0.13, 0.55, 0.13], [0.93, 0.93, 0], [0.98, 0.11, 0.18]]
    mat = ct_cls.Material(nm, cl)
    mat.printEachMaterial()
    data['appearance'] = {'materials': mat.materials}
    for building in build_obj_list:
        matVal = ct_cls.MaterialValues()
        matVal.populateMaterialValues(nm, building, data['CityObjects'][building.fid]['geometry'][0]['boundaries'][0])
        data['CityObjects'][building.fid]['geometry'][0]['material'] = matVal.material
    with open(output_path, 'w') as jsonFile:
        json.dump(data, jsonFile)


def main():
    build_obj_list = process(cityModel_path, foot_height_shp, roof_raster, facade_solar_energy)
    update_CityJSON(build_obj_list, cityModel_path, updated_JSON)
    print("Done!\n")


if __name__ == "__main__":
    main()



