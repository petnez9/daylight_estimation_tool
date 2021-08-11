'''
***********************************************************************************************************************
*                                Functions file for create raster module                                              *
*                                         Created by P.Nezval                                                         *
*                                             version 0.1                                                             *
*                                         May 2021 - July 2021                                                        *
***********************************************************************************************************************
                                                License:
                                    Mozilla Public License Version 2.0
                                    ==================================
'''
import math
import os
from itertools import chain
from rasterio import mask
from shapely.geometry import Polygon, Point, LineString, MultiPolygon
import numpy as np

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

def assign_ver_coords(roof_geom, ver_list):
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


class InitialRaster:

    def __init__(self, input_raster):
        self.top = input_raster.bounds.top
        self.left = input_raster.bounds.left
        self.resolution = input_raster.res
        self.rows = input_raster.shape[0]
        self.columns = input_raster.shape[1]
        self.CRS = 'EPSG:3008'


class Voxel:

    def __init__(self, col_index, row_index, voxel_data, id):
        self.voxel_id = id
        self.col_index = col_index
        self.row_index = row_index
        self.voxel_data = voxel_data
        self.building_index = 'none'
        self.facade_index = 999
        self.facade_start_corner = False
        self.facade_end_corner = False

    def defineCoordinates(self, init_raster):
        self.coord_x = init_raster.left + self.col_index * init_raster.resolution[0] - init_raster.resolution[0] / 2
        self.coord_y = init_raster.top - self.row_index * init_raster.resolution[1] + init_raster.resolution[1] / 2


class FacadeData:

    def __init__(self, path):
        self.path = path
        self.file_size = os.stat(path).st_size / 1000
        self.file_data = []

    def processInput(self):
        with open(self.path) as f:
            facadeData = f.readlines()[1:]
            for line in facadeData:
                line = line.split(' ')
                line = [x for x in line if x]
                irr_data = [float(x) for x in line[2:-1] if float(x) > 0]
                self.file_data.append({'row': int(line[0]), 'column': int(line[1]), 'data': irr_data})

class Facade:

    def __init__(self, fid, height, start, end, fac_id):
        self.feature_id = fid
        self.facade_index = fac_id
        self.start = start
        self.end = end
        self.facade_buffer = create_buffer(start, end, 3)
        self.height = height
        self.lineString = LineString([Point(start), Point(end)])
        self.facade_voxel_list = []
        self.facade_windows = []
        self.facade_raster = []
        self.facade_windows_raster = []

    def getFacadeWindows(self, building):
        for window in building.windows_geometry:
            if window.facade_id == self.facade_index:
                self.facade_windows.append(window)

    def createFacadeWindowRaster(self):
        raster_window = np.zeros(list(self.facade_raster.shape))
        for window in self.facade_windows:
            p1 = Point(self.start)
            p2 = Point(window.win_start_xy)
            j = round(p1.distance(p2))
            i = self.height - window.win_height_top
            n = round(window.width)
            raster_window[i:i + window.height, j:j + n] = 1
        self.facade_windows_raster = raster_window


class Building:

    def __init__(self, fid, bid, foot, height):
        self.fid = fid
        self.building_id = bid
        self.footprints = foot
        self.building_height = height
        self.facade_list = []
        self.build_voxel_list = []
        self.corner_area = 0
        self.corner_list = []
        self.facade_lines = []
        self.windows_geometry = []
        self.raster_list = []

    def assignBuildingIndex(self, voxels):
        pointList = []
        for pair in self.footprints:
            point = Point(pair[0], pair[1])
            pointList.append(point)
        poly = Polygon(pointList)
        poly = poly.buffer(2)
        for voxel in voxels:
            p1 = Point(voxel.coord_x, voxel.coord_y)
            if p1.within(poly):
                voxel.building_index = self.building_id
                self.build_voxel_list.append(voxel)

    def getFacadeList(self):
        facade_list = []
        for i in range(len(self.footprints)-1):
            facade_start = self.footprints[i]
            facade_end = self.footprints[i + 1]
            facade_id = i
            facade_list.append(Facade(self.fid, self.building_height, facade_start, facade_end, facade_id))
            self.corner_list.append(Point(facade_start[0], facade_start[1]))
            self.facade_lines.append(LineString([Point(facade_start), Point(facade_end)]))
        self.facade_list = facade_list

    def checkSide(self, x, y, z):
        # first two points are coords of a line and third point is actually the evaluated point
        area = (x[1] * (y[0] - z[0]) + y[1] * (z[0] - x[0]) + z[1] * (x[0] - y[0])) / 2
        if area > 0:  # if the area is positive the points is on one side
            return True  # if the area is negative the points is on the other side
        elif area < 0:
            return False
        else:
            return True

    def createConvexCornerArea(self, poly1, poly2, corner):
        extract = poly1.union(poly2)
        corner_area = corner.buffer(1.25).difference(extract)
        if corner_area.geom_type == 'MultiPolygon':
            polygons = list(corner_area)
            max_area = 0
            for poly in polygons:
                if poly.area > max_area:
                    corner_area = poly
        return corner_area

    def createConcaveCornerArea(self, poly1, poly2, corner):
        extract = poly1.intersection(poly2)
        corner_area = corner.buffer(1.25).intersection(extract)
        if corner_area.geom_type == 'MultiPolygon':
            polygons = list(corner_area)
            max_area = 0
            for poly in polygons:
                if poly.area > max_area:
                    corner_area = poly
        return corner_area

    def createCornerAreas(self):
        corner_list = []
        for i in range(len(self.facade_list)):
            if i < len(self.facade_list) - 1:
                if self.checkSide(self.facade_list[i].start, self.facade_list[i].end, self.facade_list[i + 1].end):
                    corner = Point(self.facade_list[i].end[0], self.facade_list[i].end[1])
                    corner_area = self.createConvexCornerArea(self.facade_list[i].facade_buffer, self.facade_list[i + 1].facade_buffer, corner)
                    corner_list.append(corner_area)
                else:
                    corner = Point(self.facade_list[i].end[0], self.facade_list[i].end[1])
                    corner_area = self.createConcaveCornerArea(self.facade_list[i].facade_buffer, self.facade_list[i + 1].facade_buffer, corner)
                    corner_list.append(corner_area)
            else:
                if self.checkSide(self.facade_list[i].start, self.facade_list[i].end, self.facade_list[0].end):
                    corner = Point(self.facade_list[i].end[0], self.facade_list[i].end[1])
                    corner_area = self.createConvexCornerArea(self.facade_list[i].facade_buffer, self.facade_list[0].facade_buffer, corner)
                    corner_list.append(corner_area)
                else:
                    corner = Point(self.facade_list[i].end[0], self.facade_list[i].end[1])
                    corner_area = self.createConcaveCornerArea(self.facade_list[i].facade_buffer, self.facade_list[0].facade_buffer, corner)
                    corner_list.append(corner_area)
        try:
            corner_area_merged = MultiPolygon(corner_list)
            if corner_area_merged.is_valid:
                self.corner_area = corner_area_merged
            else:
                self.corner_area = corner_list
                print("An exception occurred - corner area is represented by list of polygons instead of multipolygon")

        except:
            self.corner_area = corner_list
            print("An exception occurred - corner area is represented by list of polygons instead of multipolygon")

    def assignFacadeIndex(self):
        corner_points = []
        for facade in self.facade_list:
            for voxel in self.build_voxel_list:
                p1 = Point(voxel.coord_x, voxel.coord_y)
                if isinstance(self.corner_area, list):
                    for poly_corner in self.corner_area:
                        if p1.within(poly_corner):
                            if voxel not in corner_points:
                                corner_points.append(voxel)
                        else:
                            if p1.within(facade.facade_buffer):
                                voxel.facade_index = facade.facade_index
                else:
                    if p1.within(self.corner_area):
                        if voxel not in corner_points:
                            corner_points.append(voxel)
                    else:
                        if p1.within(facade.facade_buffer):
                            voxel.facade_index = facade.facade_index
        self.resolveCorners(corner_points)

    def resolveCorners(self, list):
        for voxel in list:
            p1 = Point(voxel.coord_x, voxel.coord_y)
            i = 0
            min_dist = 99999
            assigned_corner = 0
            for corner in self.corner_list:
                dist = p1.distance(corner)
                if dist < min_dist:
                    min_dist = dist
                    assigned_corner = i
                i += 1
            if assigned_corner == 0:
                dist1 = p1.distance(self.facade_lines[assigned_corner])
                dist2 = p1.distance(self.facade_lines[-1])
                if dist1 < dist2:
                    voxel.facade_index = assigned_corner
                else:
                    voxel.facade_index = len(self.facade_lines) - 1
            else:
                dist1 = p1.distance(self.facade_lines[assigned_corner])
                dist2 = p1.distance(self.facade_lines[assigned_corner - 1])
                if dist1 < dist2:
                    voxel.facade_index = assigned_corner
                else:
                    voxel.facade_index = assigned_corner - 1

    def getWallVoxelList(self, i):
        selected = []
        for voxel in self.build_voxel_list:
            if voxel.facade_index == i:
                selected.append(voxel)
        return selected

    def setVoxelCorners(self):
        for i in range(len(self.facade_list)):
            facade_voxels = self.getWallVoxelList(i)
            start_v = 0
            end_v = 0
            start_dis = 999
            end_dis = 999
            for voxel in facade_voxels:
                p1 = Point(voxel.coord_x, voxel.coord_y)
                dist1 = p1.distance(Point(self.facade_list[i].start))
                dist2 = p1.distance(Point(self.facade_list[i].end))
                if dist1 < start_dis:
                    start_dis = dist1
                    start_v = voxel.voxel_id
                if dist2 < end_dis:
                    end_dis = dist2
                    end_v = voxel.voxel_id
            for voxel in facade_voxels:
                if voxel.voxel_id == start_v:
                    voxel.facade_start_corner = True
                elif voxel.voxel_id == end_v:
                    voxel.facade_end_corner = True

class Roof:
    def __init__(self, fid):
        self.feature_index = fid
        self.roof_geometry = Polygon()
        self.total_irradiance = 0
        self.roof_type = 0
        self.roof_raw_geometry = 0

    def getRoofRawGeometry(self, building, city_model):
        semantics = building['geometry'][0]['semantics']
        id = get_roof_id(semantics)
        roof_b_ids = get_roof_bound_id(id, semantics['values'][0])
        geom = get_roof_geom(roof_b_ids, building['geometry'][0]['boundaries'][0])
        ver_list = get_all_vertices_list(city_model)
        roof_geom_list = assign_ver_coords(geom, ver_list)
        self.roof_raw_geometry = roof_geom_list

    def getRoofPolygon(self):
        if len(self.roof_raw_geometry) == 1:
            roof_g_list = self.roof_raw_geometry[0]
            poly_list = []
            for edge in roof_g_list:
                poly_list.append(edge['ver_coords'])
            poly_list.append(poly_list[0])
            self.roof_geometry = [Polygon(poly_list)]
        else:
            poly = []
            for roof_feature in self.roof_raw_geometry:
                poly_list = []
                for edge in roof_feature:
                    poly_list.append(edge['ver_coords'])
                poly_list.append(poly_list[0])
                poly.append(Polygon(poly_list))
            self.roof_geometry = poly

    def extractRoofIrradiance(self, raster):
        raster_mask = mask.mask(raster, self.roof_geometry, crop=True)
        raster_values = []
        for i in range(len(raster_mask[0][0])):
            raster_values.append([x for x in raster_mask[0][0][i] if x > 0])
        raster_values = [x for x in raster_values if x]
        raster_values = list(chain(*raster_values))
        self.total_irradiance = sum(raster_values)


class FacadeRaster(Facade):

    def __init__(self, fid, height, start, end, fac_id):
        Facade.__init__(self, fid, height, start, end, fac_id)
        self.voxelList = []
        self.spatial_reference = {}
        self.columns = 0
        self.rows = 0
        self.resolution = 1
        self.grid = 0

    def getRasterDimension(self):
        start = Point(self.start)
        end = Point(self.end)
        dim = start.distance(end)
        self.columns = round(dim)
        self.rows = round(self.height)
        self.grid = np.zeros([self.rows, self.columns])

    def getFacadeVoxelList(self, building):
        selected = []
        for voxel in building.build_voxel_list:
            if voxel.facade_index == self.facade_index:
                selected.append(voxel)
        self.voxelList = selected

    def checkQuadrant(self):
        x_diff = self.end[0] - self.start[0]
        y_diff = self.end[1] - self.start[1]
        if x_diff > 0 and y_diff > 0:
            return "Q1"
        elif x_diff < 0 and y_diff > 0:
            return "Q2"
        elif x_diff < 0 and y_diff < 0:
            return "Q3"
        elif x_diff > 0 and y_diff < 0:
            return "Q4"

    def calculateLineCoords(self):
        for voxel in self.voxelList:
            quadrant = self.checkQuadrant()
            if quadrant == "Q1":
                main_dir = math.atan((self.end[1]-self.start[1])/(self.end[0]-self.start[0]))
                alpha = math.atan((voxel.coord_y - self.start[1]) / (voxel.coord_x - self.start[0])) - main_dir
                start = Point(self.start)
                p1 = Point([voxel.coord_x, voxel.coord_y])
                distance = start.distance(p1)
                new_distance = distance * math.cos(alpha)
                delta_xp1 = new_distance * math.cos(main_dir)
                delta_yp1 = new_distance * math.sin(main_dir)
                voxel.coord_x = self.start[0] + delta_xp1
                voxel.coord_y = self.start[1] + delta_yp1
            elif quadrant == "Q2":
                main_dir = math.pi - math.atan((self.end[1]-self.start[1])/(self.end[0]-self.start[0]))
                alpha = main_dir - (math.pi - math.atan((voxel.coord_y - self.start[1]) / (voxel.coord_x - self.start[0])))
                start = Point(self.start)
                p1 = Point([voxel.coord_x, voxel.coord_y])
                distance = start.distance(p1)
                new_distance = distance * math.cos(alpha)
                delta_xp1 = new_distance * math.cos(main_dir)
                delta_yp1 = new_distance * math.sin(main_dir)
                voxel.coord_x = self.start[0] - abs(delta_xp1)
                voxel.coord_y = self.start[1] + abs(delta_yp1)
            elif quadrant == "Q3":
                main_dir = math.atan((self.end[1]-self.start[1])/(self.end[0]-self.start[0]))
                alpha = main_dir - math.atan((voxel.coord_y - self.start[1]) / (voxel.coord_x - self.start[0]))
                start = Point(self.start)
                p1 = Point([voxel.coord_x, voxel.coord_y])
                distance = start.distance(p1)
                new_distance = distance * math.cos(alpha)
                delta_xp1 = new_distance * math.cos(main_dir)
                delta_yp1 = new_distance * math.sin(main_dir)
                voxel.coord_x = self.start[0] - delta_xp1
                voxel.coord_y = self.start[1] - delta_yp1
            elif quadrant == "Q4":
                main_dir = 2*math.pi - math.atan((self.end[1] - self.start[1]) / (self.end[0] - self.start[0]))
                alpha = main_dir - (2 * math.pi - math.atan((voxel.coord_y - self.start[1]) / (voxel.coord_x - self.start[0])))
                start = Point(self.start)
                p1 = Point([voxel.coord_x, voxel.coord_y])
                distance = start.distance(p1)
                new_distance = distance * math.cos(alpha)
                delta_xp1 = new_distance * math.cos(main_dir)
                delta_yp1 = new_distance * math.sin(main_dir)
                voxel.coord_x = self.start[0] + delta_xp1
                voxel.coord_y = self.start[1] - delta_yp1

    def find_free_spot_up(self, i, list):
        while i in list:
            i = i + 1
        return i

    def sort_dist_list(self):
        dist_list = []
        start = Point(self.start)
        for voxel in self.voxelList:
            voxel_loc = Point([voxel.coord_x, voxel.coord_y])
            d = start.distance(voxel_loc)
            x = {'distance': d, 'data': voxel.voxel_data}
            dist_list.append(x)
        dist_list.sort(key=my_fun)
        return dist_list

    def adjust_position(self, dist_list):
        check_list = []
        new_dist_list = []
        for i in range(len(dist_list)):
            d = round(dist_list[i]['distance'])
            if d in check_list:
                d = self.find_free_spot_up(d, check_list)
            if d <= self.columns - 1:
                x = {'distance': d, 'data': dist_list[i]['data']}
                new_dist_list.append(x)
                check_list.append(d)
        return new_dist_list

    def populate_vertical_values(self, voxel_data, d):
        if len(voxel_data) == self.rows:
            for i in range(len(voxel_data)):
                self.grid[len(voxel_data) - 1 - i][d] = voxel_data[i]
        elif len(voxel_data) < self.rows:
            diff = self.rows - len(voxel_data)
            for i in range(len(voxel_data)):
                self.grid[len(voxel_data) - diff - i][d] = voxel_data[i]
        elif len(voxel_data) > self.rows:
            diff = len(voxel_data) - self.rows
            voxel_data = voxel_data[diff:]
            for i in range(len(voxel_data)):
                self.grid[len(voxel_data) - 1 - i][d] = voxel_data[i]


    def populateGrid(self):
        if len(self.voxelList) == self.columns:
            dist_list = self.sort_dist_list()
            for i in range(len(dist_list)):
                self.populate_vertical_values(dist_list[i]['data'], i)
        elif len(self.voxelList) < self.columns:
            dist_list = self.sort_dist_list()
            dist_list = self.adjust_position(dist_list)
            for i in range(len(dist_list)):
                self.populate_vertical_values(dist_list[i]['data'], dist_list[i]['distance'])


    def create_foating_window(self, i, j, size):
        n = round(size/2)
        if i < n and j < n:
            slice = np.array(self.grid[:(i + n), :(j + n)])
            return slice
        elif i >= self.rows - n and j >= self.columns - n:
            slice = np.array(self.grid[(i - int(n / 2)):, (j - int(n / 2)):])
            return slice
        elif i < n and j >= n:
            slice = np.array(self.grid[:(i + n), (j - int(n/2)):(j + n)])
            return slice
        elif i >= n and j < n:
            slice = np.array(self.grid[(i - int(n/2)):(i + n), :(j + n)])
            return slice
        elif i >= self.rows - n and j < self.columns - n:
            slice = np.array(self.grid[(i - int(n / 2)):, (j - int(n / 2)):(j + n)])
            return slice
        elif i < self.rows - n and j >= self.columns - n:
            slice = np.array(self.grid[(i - int(n / 2)):(i + n), (j - int(n / 2)):])
            return slice
        else:
            slice = np.array(self.grid[(i - int(n/2)):(i + n), (j - int(n/2)):(j + n)])
            return slice

    def getMeanValue(self, submat):
        non_zero = []
        for x in submat:
            [non_zero.append(x) for x in x if x > 0]
        np_arr = np.array(non_zero)
        if not np_arr.any():
            return 0
        else:
            return np_arr.mean()

    def interpolate(self):
        for i in range(self.rows):
            for j in range(self.columns):
                if self.grid[i][j] == 0:
                    submat = self.create_foating_window(i, j, 3)
                    self.grid[i][j] = self.getMeanValue(submat)

    def storeRasterObject(self):
        pass

class OutputCityJSON:

    def __init__(self, data):
        self.metadata = {"geographicalExtent": [], "referenceSystem": 0}
        self.type = "CityJSON"
        self.version = "1.0"
        self.cityObjects = {}

    def printOut(self):
        json_out = {"type": self.type, "version": self.version, "metadata": self.metadata}


