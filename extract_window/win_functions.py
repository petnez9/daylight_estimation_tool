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

def extract_footprints(data, key):
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
        self.edge1 = 0
        self.edge2 = 0
        self.wall_id = None
        self.valid = False

    def check_win_geom_validity(self):
        x = len(self.geometry)
        if x != 4 and x != 2:
            print('Window {a} is invalid'.format(a=self.id))
            print(x)
        else:
            self.valid = True

    def set_extent(self):
        self.check_win_geom_validity()
        if self.valid:
            if len(self.geometry) == 4:
                x1 = self.geometry[0]['ver_coords'][0]
                x2 = self.geometry[2]['ver_coords'][0]
                y1 = self.geometry[0]['ver_coords'][1]
                y2 = self.geometry[2]['ver_coords'][1]
                z1 = self.geometry[0]['ver_coords'][2]
                z2 = self.geometry[2]['ver_coords'][2]
                p1 = Point([x1, y1])
                p2 = Point([x2, y2])
                a = round(p1.distance(p2))
                b = abs(z1 - z2)
            elif len(self.geometry) == 2:
                x1 = self.geometry[0]['ver_coords'][0]
                x2 = self.geometry[1]['ver_coords'][0]
                y1 = self.geometry[0]['ver_coords'][1]
                y2 = self.geometry[1]['ver_coords'][1]
                z1 = self.geometry[0]['ver_coords'][2]
                z2 = self.geometry[1]['ver_coords'][2]
                p1 = Point([x1, y1])
                p2 = Point([x2, y2])
                a = round(p1.distance(p2))
                b = abs(z1 - z2)
            if a*b > 0:
                self.width = round(p1.distance(p2))
                self.height = abs(z1 - z2)
                self.edge1 = [x1, y1]
                self.edge2 = [x2, y2]
            else:
                self.valid = False

    def assign_wall(self, wall, buffer):
        if self.valid:
            if Point(self.edge1).within(buffer) and Point(self.edge2).within(buffer):
                self.wall_id = wall['wall_id']
                print('Window {a} is has been assigned to wall id {b}'.format(a=self.id, b=self.wall_id))
