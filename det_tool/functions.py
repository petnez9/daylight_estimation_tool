'''
***********************************************************************************************************************
*                                Functions for daylight estimation tool                                               *
*                                         Created by P.Nezval                                                         *
*                                             version 0.1                                                             *
*                                         May 2021 - August 2021                                                      *
***********************************************************************************************************************
                                                License:
                                               MIT License
                                    ==================================
'''

from CityJSON_gen_upd import classes as ct_cls
from matplotlib import pyplot as plt
from shapely.geometry import Polygon, Point, LineString
import json
import numpy as np
from pre_processor import functions as pfun


def handleZeroDivision(x,y):
    try:
        return x/y
    except ZeroDivisionError:
        return 0


def extract_footprints(obj_id, shapefile):
    i = 0
    coords = []
    while len(coords) == 0:
        if shapefile.fid[i] == obj_id:
            x, y = shapefile.exterior[i].xy
            j = 0
            for j in range(len(x)):
                coords.append([x[j], y[j]])
        i += 1
    return coords


def get_roof_id(semantics):
    for i in range(len(semantics['surfaces'])):
        if semantics['surfaces'][i]['type'] == 'RoofSurface':
            return i


def get_roof_bound_id(id, val):
    ids = []
    for i in range(len(val)):
        if val[i] == id:
            ids.append(i)
    return ids


def get_roof_geom(ids, bound):
    geom = []
    for i in range(len(bound)):
        if i in ids:
            geom.append(bound[i][0])
    return geom


def get_all_vertices_list(data):
    vertex = []
    for i in range(len(data['vertices'])):
        x = {'ver_id':i, 'ver_coords':data['vertices'][i]}
        vertex.append(x)
    return vertex


def assign_roof_ver_coords(roof_geom, ver_list):
    val_list = []
    for dict in ver_list:
        x = list(dict.values())[0]
        val_list.append(x)
    win_geom_list = []
    for obj in roof_geom:
        win_coord_list = []
        for vertex in obj:
            position = val_list.index(vertex)
            win_coord_list.append(ver_list[position])
        win_geom_list.append(win_coord_list)
    return win_geom_list


def create_buffer(point1, point2, buff_len):
    line = LineString([Point(point1), Point(point2)])
    left = line.parallel_offset(buff_len / 2, 'left')
    right = line.parallel_offset(buff_len / 2, 'right')
    p1 = left.boundary[1]
    p2 = right.boundary[0]
    p3 = right.boundary[1]
    p4 = left.boundary[0]
    return Polygon([p1, p2, p3, p4, p1])


def my_fun(e):
    return e['distance']


def get_vertex_list(obj_build, data):
    vertex = []
    for i in range(len(obj_build['geometry'][0]['boundaries'][0])):
        for j in range(len(obj_build['geometry'][0]['boundaries'][0][i][0])):
            x = obj_build['geometry'][0]['boundaries'][0][i][0][j]
            x = {'ver_id':x, 'ver_coords':data['vertices'][x]}
            if x not in vertex:
                vertex.append(x)
    return vertex


def extract_footprints_own(data, key):
    building = data['CityObjects'][key]
    bound_list = building['geometry'][0]['boundaries'][0]
    foot_vertex_list = get_footprint_vertex_list(bound_list, data)
    coordinates = []
    for vertex in foot_vertex_list:
        coordinates.append(data['vertices'][vertex][0:2])
    x = {'id': key, 'coordinates': coordinates}
    return x


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


def get_surfaces_info(data, key):
    surf = data['CityObjects'][key]['geometry'][0]['semantics']['surfaces']
    surface = {}
    i = 0
    for dict in surf:
        surface[dict['type']] = {'index': i, 'json_type': dict}
        i += 1
    return surface


def get_win_id(semantics):
    for i in range(len(semantics['surfaces'])):
        if semantics['surfaces'][i]['type'] == 'Window':
            return i


def get_win_bound_id(id, val):
    ids = []
    for i in range(len(val)):
        if val[i] == id:
            ids.append(i)
    return ids


def get_win_geom(ids, bound):
    geom = []
    for i in range(len(bound)):
        if i in ids:
            geom.append(bound[i][0])
    return geom


def assign_ver_coords(win_geom, ver_list):
    val_list = []
    for dict in ver_list:
        x = list(dict.values())[0]
        val_list.append(x)
    win_geom_list = []
    for obj in win_geom:
        win_coord_list = []
        for vertex in obj:
            position = val_list.index(vertex)
            win_coord_list.append(ver_list[position])
        win_geom_list.append(win_coord_list)
    return win_geom_list


def extract_window_geometry(building, data):
    semantics = building['geometry'][0]['semantics']
    id = get_win_id(semantics)
    win_b_ids = get_win_bound_id(id, semantics['values'][0])
    geom = get_win_geom(win_b_ids, building['geometry'][0]['boundaries'][0])
    ver_list = get_all_vertices_list(data)
    win_geom_list = assign_ver_coords(geom, ver_list)
    return win_geom_list, win_b_ids


def process_building(building, voxels):
    building.getFacadeList()
    building.assignBuildingIndex(voxels)
    building.createCornerAreas()
    building.assignFacadeIndex()


def process_roof(roof, feature, data, img):
    roof.getRoofGeometry(feature, data)
    roof.getRoofPolygon()
    roof.extractRoofIrradiance(img)


def process_fac_raster(facadeRaster, building, facade):
    facadeRaster.getRasterDimension()
    facadeRaster.getFacadeVoxelList(building)
    facade.facade_voxel_list = facadeRaster.voxelList
    facadeRaster.calculateLineCoords()
    facadeRaster.populateGrid()
    facadeRaster.interpolate()
    facade.facade_raster = facadeRaster.grid


def process_window(win_obj, facade, building):
    win_obj.set_extent(facade)
    facade_win_buff = create_buffer(facade.start, facade.end, 1)
    win_obj.assign_facade(facade, facade_win_buff)
    if win_obj not in building.windows_geometry:
        building.windows_geometry.append(win_obj)
    facade.getFacadeWindows(building)
    facade.getWindowIrradiance()


def process_output(building, nm, data, matVal, semAttr):
    matVal.populateMaterialValues(nm, building, data['CityObjects'][building.fid]['geometry'][0]['boundaries'][0])
    data['CityObjects'][building.fid]['geometry'][0]['material'] = matVal.material
    val = data['CityObjects'][building.fid]['geometry'][0]['semantics']['values'][0]
    semAttr.populateSemanticValues(val)
    semAttr.populateSemanticSurface(val)
    data['CityObjects'][building.fid]['geometry'][0]['semantics']['values'] = semAttr.values
    data['CityObjects'][building.fid]['geometry'][0]['semantics']['surfaces'] = semAttr.surfaces
    data['version'] = '1.0'
    data['metadata']['referenceSystem'] = 'urn:ogc:def:crs:EPSG::3008'


def update_CityModel(build_obj_list, model_path, output_path, irr_values):
    data = pfun.read_json(model_path)
    nm = ['irradiation1', 'irradiation2', 'irradiation3', 'irradiation4']
    cl = [[0, 0, 0.5], [0.13, 0.55, 0.13], [0.93, 0.93, 0], [0.98, 0.11, 0.18]]
    mat = ct_cls.Material(nm, cl)
    mat.printEachMaterial()
    data['appearance'] = {'materials': mat.materials}
    irr_values = [x for x in irr_values if x > 0]
    createHist(irr_values)
    stat = makeStat(irr_values)
    threshold = getInput(stat)
    print("Data are being written to CityJSON ...")
    for building in build_obj_list:
        matVal = ct_cls.MaterialValues(threshold)
        semAttr = ct_cls.SemanticAttributes(building)
        process_output(building, nm, data, matVal, semAttr)
    with open(output_path, 'w') as jsonFile:
        json.dump(data, jsonFile)


def generate_report(build_obj_list, path):
    report = ct_cls.ReportGenerator(build_obj_list)
    report.print_report()
    with open(path, 'w') as jsonFile:
        json.dump(report.report, jsonFile)


def getInput(stat):
    print("Based on basic statistics: \nmax - {a}, min - {b}, avg - {c}\n"
          "and histogram, choose your threshold values for windows irradiance suitability".format(a=stat[0],b=stat[1], c=stat[2]))
    val1 = int(input("Threshold 1 - Not suitable (values less then this one are considered as not suitable)\n"
                 "Enter the value:"))
    val2 = int(input("Threshold 2 - Average suitability (values less then this one are considered as fairly suitable)\n"
                 "Enter the value:"))
    val3 = int(input("Threshold 3 - Good suitability (values less then this one are considered as well suitable),"
                 "all values higher than this one are considered as highly suitable\n"
                 "Enter the value:"))
    return [val1, val2, val3]


def createHist(irr_list):
    irr_list = [x for x in irr_list if x > 0]
    irr = np.array(irr_list)
    plt.style.use('ggplot')
    plt.title("Histogram of windows irradiance values;\nmax - {a}, min - {b}, average - {c}"
              .format(a=format(irr.max(),'.2f'),b=format(irr.min(), '.2f'), c=format(irr.mean(), '.2f')))
    plt.hist(irr_list, bins=round(irr.max()/irr.min()))
    plt.show()


def makeStat(irr_list):
    irr_list = [x for x in irr_list if x > 0]
    irr = np.array(irr_list)
    return [format(irr.max(), '.2f'), format(irr.min(), '.2f'), format(irr.mean(), '.2f')]

