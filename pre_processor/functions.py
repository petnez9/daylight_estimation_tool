'''
***********************************************************************************************************************
*                                           Functions file                                                            *
*                                         Created by P.Nezval                                                         *
*                                             version 0.1                                                             *
*                                               May 2021                                                              *
***********************************************************************************************************************
                                                License:
                                    Mozilla Public License Version 2.0
                                    ==================================
'''

import json
import geopandas
import pandas
from shapely.geometry import Polygon


######################################################################################################################

def read_json(input_data):
    with open(input_data) as file:
        data = json.load(file)
    return data


def get_obj_ids(input_data):
    city_obj_id = []
    for i in input_data['CityObjects'].keys():
        city_obj_id.append(i)
    return city_obj_id


def vertex_heights(input_data):
    heights = {}
    for i in range(len(input_data['vertices'])):
        x = {i: input_data['vertices'][i][2]}
        heights.update(x)
    return heights


def get_vertex_list(obj_build):
    vertex = []
    for i in range(len(obj_build['geometry'][0]['boundaries'][0])):
        for j in range(len(obj_build['geometry'][0]['boundaries'][0][i][0])):
            x = obj_build['geometry'][0]['boundaries'][0][i][0][j]
            if x not in vertex:
                vertex.append(x)
    return vertex


def get_build_height(vert_list, heights):
    height_min = 10000
    height_max = 0
    for x in vert_list:
        y = heights[x]
        if y > height_max:
            height_max = y
        if y < height_min:
            height_min = y
    height = height_max - height_min
    return height


def get_footprint_vertex_list(building_bound, input_data):
    footprints_ver = []
    min_sum = 9999999
    for polygon in building_bound:
        sum_height = 0
        n = 0
        for vertex in polygon[0]:
            n += 1
            sum_height += input_data['vertices'][vertex][2]
        sum_height = sum_height/n
        if sum_height < min_sum:
            min_sum = sum_height
            footprints_ver = polygon[0]
    return footprints_ver


def extract_footprints_own(data, key):
    building = data['CityObjects'][key]
    bound_list = building['geometry'][0]['boundaries'][0]
    foot_vertex_list = get_footprint_vertex_list(bound_list, data)
    coordinates = []
    for vertex in foot_vertex_list:
        coordinates.append(data['vertices'][vertex][0:2])
    return coordinates


def create_foot_SHP(input_data, oid, out_shp):
    out = out_shp
    fid = []
    heights = []
    vid = []
    for i in oid:
        fid.append(i)
        building = input_data['CityObjects'][i]
        ver_list = get_vertex_list(building)
        ver_heights = vertex_heights(input_data)
        heights.append(get_build_height(ver_list, ver_heights))
        coords = extract_footprints_own(input_data, i)
        x, y = [],[]
        for pair in coords:
            coor_x, coor_y = pair[0],pair[1]
            x.append(coor_x)
            y.append(coor_y)
        vid.append(Polygon(zip(x, y)))
    df = pandas.DataFrame({'id': fid, 'heights': heights})
    gdf = geopandas.GeoDataFrame(df, geometry=vid)
    gdf.to_file(out)


def main(input_json, output_shp):
    data = read_json(input_json)
    obj_id = get_obj_ids(data)
    create_foot_SHP(data, obj_id, output_shp)