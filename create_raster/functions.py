'''
***********************************************************************************************************************
*                                Functions file for create raster module                                              *
*                                         Created by P.Nezval                                                         *
*                                             version 0.1                                                             *
*                                               May 2021                                                              *
***********************************************************************************************************************
                                                License:
                                    Mozilla Public License Version 2.0
                                    ==================================
'''

import os
from shapely.geometry import Polygon, Point, LineString, MultiPolygon


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

def create_buffer(point1, point2):
    line = LineString([Point(point1[0], point1[1]), Point(point2[0], point2[1])])
    buff_length = 2.4
    left = line.parallel_offset(buff_length / 2, 'left')
    right = line.parallel_offset(buff_length / 2, 'right')
    p1 = left.boundary[1]
    p2 = right.boundary[0]
    p3 = right.boundary[1]
    p4 = left.boundary[0]
    return Polygon([p1, p2, p3, p4, p1])


class InitialRaster():

    def __init__(self, input_raster):
        self.top = input_raster.bounds.top
        self.left = input_raster.bounds.left
        self.resolution = input_raster.res
        self.rows = input_raster.shape[0]
        self.columns = input_raster.shape[1]
        self.CRS = 'EPSG:3008'


class VoxelObject():

    def __init__(self, col_index, row_index, voxel_data, id):
        self.voxel_id = id
        self.col_index = col_index
        self.row_index = row_index
        self.voxel_data = voxel_data
        self.building_index = 'none'
        self.wall_index = 999
        self.wall_start_corner = False
        self.wall_end_corner = False

    def defineCoordinates(self, init_raster):
        self.coord_x = init_raster.left + self.col_index * init_raster.resolution[0] - init_raster.resolution[0] / 2
        self.coord_y = init_raster.top - self.row_index * init_raster.resolution[1] + init_raster.resolution[1] / 2


class WallDataObject():

    def __init__(self, path):
        self.path = path
        self.file_size = os.stat(path).st_size / 1000
        self.file_data = []

    def processInput(self):
        with open(self.path) as f:
            wallData = f.readlines()[1:]
            for line in wallData:
                line = line.split(' ')
                line = [x for x in line if x]
                irr_data = [float(x) for x in line[2:-1] if float(x) > 0]
                self.file_data.append({'row': int(line[0]), 'column': int(line[1]), 'data': irr_data})


class BuildingObject():

    def __init__(self, fid, bid, foot, height):
        self.fid = fid
        self.building_id = bid
        self.footprints = foot
        self.building_height = height
        self.wall_list = []
        self.build_voxel_list = []
        self.corner_area = 0
        self.corner_list = []
        self.wall_lines = []

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

    def getWallList(self):
        wall_list = []
        for i in range(len(self.footprints)-1):
            wall_start = self.footprints[i]
            wall_end = self.footprints[i + 1]
            wall_id = i
            wall_list.append({'wall_id': wall_id, 'wall_start': wall_start, 'wall_end': wall_end, 'wall_buffer': create_buffer(wall_start, wall_end)})
            self.corner_list.append(Point(wall_start[0], wall_start[1]))
            self.wall_lines.append(LineString([Point(wall_start[0], wall_start[1]), Point(wall_end[0], wall_end[1])]))
        self.wall_list = wall_list

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
        for i in range(len(self.wall_list)):
            if i < len(self.wall_list) - 1:
                if self.checkSide(self.wall_list[i]['wall_start'], self.wall_list[i]['wall_end'], self.wall_list[i+1]['wall_end']):
                    corner = Point(self.wall_list[i]['wall_end'][0], self.wall_list[i]['wall_end'][1])
                    corner_area = self.createConvexCornerArea(self.wall_list[i]['wall_buffer'], self.wall_list[i+1]['wall_buffer'], corner)
                    corner_list.append(corner_area)
                else:
                    corner = Point(self.wall_list[i]['wall_end'][0], self.wall_list[i]['wall_end'][1])
                    corner_area = self.createConcaveCornerArea(self.wall_list[i]['wall_buffer'], self.wall_list[i+1]['wall_buffer'], corner)
                    corner_list.append(corner_area)
            else:
                if self.checkSide(self.wall_list[i]['wall_start'], self.wall_list[i]['wall_end'], self.wall_list[0]['wall_end']):
                    corner = Point(self.wall_list[i]['wall_end'][0], self.wall_list[i]['wall_end'][1])
                    corner_area = self.createConvexCornerArea(self.wall_list[i]['wall_buffer'], self.wall_list[0]['wall_buffer'], corner)
                    corner_list.append(corner_area)
                else:
                    corner = Point(self.wall_list[i]['wall_end'][0], self.wall_list[i]['wall_end'][1])
                    corner_area = self.createConcaveCornerArea(self.wall_list[i]['wall_buffer'], self.wall_list[0]['wall_buffer'], corner)
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

    def assignWallIndex(self):
        corner_points = []
        for wall in self.wall_list:
            for voxel in self.build_voxel_list:
                p1 = Point(voxel.coord_x, voxel.coord_y)
                if isinstance(self.corner_area, list):
                    for poly_corner in self.corner_area:
                        if p1.within(poly_corner):
                            if voxel not in corner_points:
                                corner_points.append(voxel)
                        else:
                            if p1.within(wall['wall_buffer']):
                                voxel.wall_index = wall['wall_id']
                else:
                    if p1.within(self.corner_area):
                        if voxel not in corner_points:
                            corner_points.append(voxel)
                    else:
                        if p1.within(wall['wall_buffer']):
                            voxel.wall_index = wall['wall_id']
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
                dist1 = p1.distance(self.wall_lines[assigned_corner])
                dist2 = p1.distance(self.wall_lines[-1])
                if dist1 < dist2:
                    voxel.wall_index = assigned_corner
                else:
                    voxel.wall_index = len(self.wall_lines) - 1
            else:
                dist1 = p1.distance(self.wall_lines[assigned_corner])
                dist2 = p1.distance(self.wall_lines[assigned_corner - 1])
                if dist1 < dist2:
                    voxel.wall_index = assigned_corner
                else:
                    voxel.wall_index = assigned_corner - 1

    def getWallVoxelList(self, i):
        selected= []
        for voxel in self.build_voxel_list:
            if voxel.wall_index == i:
                selected.append(voxel)
        return selected

    def setVoxelCorners(self):
        for i in range(len(self.wall_list)):
            wall_voxels = self.getWallVoxelList(i)
            start_v = 0
            end_v = 0
            start_dis = 999
            end_dis = 999
            for voxel in wall_voxels:
                p1 = Point(voxel.coord_x, voxel.coord_y)
                dist1 = p1.distance(Point(self.wall_list[i]['wall_start']))
                dist2 = p1.distance(Point(self.wall_list[i]['wall_end']))
                if dist1 < start_dis:
                    start_dis = dist1
                    start_v = voxel.voxel_id
                if dist2 < end_dis:
                    end_dis = dist2
                    end_v = voxel.voxel_id
            for voxel in wall_voxels:
                if voxel.voxel_id == start_v:
                    voxel.wall_start_corner = True
                elif voxel.voxel_id == end_v:
                    voxel.wall_end_corner = True


class WallRasterObject():

    def __init__(self, bid, wall_index, coord_x_sta, coord_y_sta, coord_x_end, coord_y_end, height):
        self.build_index = bid
        self.wall_index = wall_index
        self.left = coord_x_sta, coord_y_sta
        self.right = coord_x_end, coord_y_end
        self.data = []
        self.height = height
        self.spatial_reference = {}

    def getRasterDimension(self):
        a = 0

