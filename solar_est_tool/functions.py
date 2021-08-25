'''
***********************************************************************************************************************
*                                Functions for solar estimation tool                                                  *
*                                         Created by P.Nezval                                                         *
*                                             version 0.1                                                             *
*                                         May 2021 - August 2021                                                      *
***********************************************************************************************************************
                                                License:
                                               MIT License
                                    ==================================
'''


from shapely.geometry import Polygon, Point, LineString


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


def get_win_from_json(building, data):
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
    roof.getRoofRawGeometry(feature, data)
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