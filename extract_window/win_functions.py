'''
***********************************************************************************************************************
*                                Functions file for extract windows module                                            *
*                                         Created by P.Nezval                                                         *
*                                             version 0.1                                                             *
*                                         May 2021 - July 2021                                                        *
***********************************************************************************************************************
                                                License:
                                    Mozilla Public License Version 2.0
                                    ==================================
'''

from shapely.geometry import Point
import math


def get_vertex_list(obj_build, data):
    vertex = []
    for i in range(len(obj_build['geometry'][0]['boundaries'][0])):
        for j in range(len(obj_build['geometry'][0]['boundaries'][0][i][0])):
            x = obj_build['geometry'][0]['boundaries'][0][i][0][j]
            x = {'ver_id':x, 'ver_coords':data['vertices'][x]}
            if x not in vertex:
                vertex.append(x)
    return vertex


def get_all_vertices_list(data):
    vertex = []
    for i in range(len(data['vertices'])):
        x = {'ver_id':i, 'ver_coords':data['vertices'][i]}
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
    return win_geom_list

class Window:

    def __init__(self, geom, bid, id):
        self.geometry = geom
        self.bid = bid
        self.id = id
        self.width = 0
        self.height = 0
        self.win_start_xy = [0,0]
        self.win_end_xy = 0
        self.win_height_top = 0
        self.facade_id = None
        self.valid = False

    def set_extent(self, facade):
        i = 0
        z_min = 999
        z_max = 0
        self.win_end_xy = facade.start
        for geo_dict in self.geometry:
            x,y,z = geo_dict['ver_coords']
            if Point(facade.start).distance(Point([x,y])) < Point(facade.start).distance(Point(self.win_start_xy)):
                self.win_start_xy = [x,y]
            if Point(facade.start).distance(Point([x,y])) > Point(facade.start).distance(Point(self.win_end_xy)):
                self.win_end_xy = [x,y]
            if z > z_max:
                z_max = z
            if z < z_min:
                z_min = z
        self.width = Point(self.win_start_xy).distance(Point(self.win_end_xy))
        self.height = z_max - z_min
        self.win_height_top = z_max
        self.valid = True

    def assign_facade(self, facade, buffer):
        if self.valid:
            if Point(self.win_start_xy).within(buffer) and Point(self.win_end_xy).within(buffer):
                self.facade_id = facade.facade_index