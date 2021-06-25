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
from shapely.geometry import Polygon, Point, LineString


def extract_footprints(obj_id, shapefile):
    i = 0
    coords = []
    while len(coords) == 0:
        if shapefile.fid[i] == obj_id:
            x,y = shapefile.exterior[i].xy
            j = 0
            for j in range(len(x)):
                coords.append([x[j], y[j]])
        i += 1
    return coords


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
        self.coord_x = init_raster.left + self.col_index*init_raster.resolution[0] - init_raster.resolution[0]/2
        self.coord_y = init_raster.top - self.row_index*init_raster.resolution[1] + init_raster.resolution[1]/2


class WallObject():

    def __init__(self, wall_index, coord_x_sta, coord_y_sta, coord_x_end, coord_y_end):
        self.wall_index = wall_index
        self.coord_x_start = coord_x_sta
        self.coord_y_start = coord_y_sta
        self.coord_x_end = coord_x_end
        self.coord_y_end = coord_y_end
        self.wall_height = 0
        self.wall_dimension = 0
        self.building_id = None


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


    def assignBuildingIndex(self, voxels):
        pointList = []
        for pair in self.footprints:
            point = Point(pair[0], pair[1])
            pointList.append(point)
        poly = Polygon(pointList)
        poly = poly.buffer(1.5)
        for voxel in voxels:
            p1 = Point(voxel.coord_x, voxel.coord_y)
            if p1.within(poly):
                voxel.building_index = self.building_id
        return


    def getWallList(self):
        wall_list = []
        for i in range(len(self.footprints)):
            if i == (len(self.footprints)-1):
                wall_start = self.footprints[-1]
                wall_end = self.footprints[0]
                wall_id = i
            else:
                wall_start = self.footprints[i]
                wall_end = self.footprints[i+1]
                wall_id = i
            wall_list.append({'wall_id': wall_id, 'wall_start': wall_start, 'wall_end': wall_end})
        self.wall_list = wall_list


    def correctSide(self, x, y, z):
        # first two points are coords of a line and third point is actually the evaluated point
        area = (x[1] * (y[0] - z[0]) + y[1] * (z[0] - x[0]) + z[1] * (x[0] - y[0])) / 2
        if area > 0: #if the area is positive the points is on one side
            return False #if the area is negative the points is on the other side
        elif area < 0:
            return True
        else:
            return True


    def checkCorner(self, list, poly):
        for voxel in list:
            


    def assignWallIndex(self, voxel_list):
        build_voxel_list = []
        for voxel in voxel_list:
            if voxel.building_index == self.building_id:
                build_voxel_list.append(voxel)
        build_shape = self.footprints
        #print('The building id is: {a} and number of voxels for this building is: {b}'.format(a=self.building_id, b=len(build_voxel_list)))
        for wall in self.wall_list:
            line = LineString([Point(wall['wall_start'][0], wall['wall_start'][1]), Point(wall['wall_end'][0], wall['wall_end'][1])])
            if line.length > 0.5:
                buff_length = 2
                left = line.parallel_offset(buff_length / 2, 'left')
                right = line.parallel_offset(buff_length / 2, 'right')
                p1 = left.boundary[1]
                p2 = right.boundary[0]
                p3 = right.boundary[1]
                p4 = left.boundary[0]
                poly_buffer = Polygon([p1, p2, p3, p4, p1])
                for voxel in build_voxel_list:
                    p1 = Point(voxel.coord_x, voxel.coord_y)
                    if p1.within(poly_buffer):
                        voxel.wall_index = wall['wall_id']
                        '''start = [wall['wall_start'][0], wall['wall_start'][1]]
                        point = [voxel.coord_x, voxel.coord_y]
                        if self.correctSide(start, end, point):
                            #print('The voxel row_id is: {a}; col_id is: {b} represents the following wall: {c}'.format(
                                    #a=voxel.row_index, b=voxel.col_index, c=wall['wall_id']))
                            #voxel.wall_index = wall['wall_id']
                            if len(voxel.wall_index) > 1:
                                if p1.within(Point(wall['wall_start'][0], wall['wall_start'][1]).buffer(1.5)):
                                    voxel.wall_start_corner = True
                                elif p1.within(Point(wall['wall_end'][0], wall['wall_end'][1]).buffer(1.5)):
                                    voxel.wall_end_corner = True'''


class WallRasterObject():

    def __init__(self, wall_index, coord_x_sta, coord_y_sta, coord_x_end, coord_y_end):
        self.top = input_raster.bounds.top
        self.left = input_raster.bounds.left
        self.resolution = input_raster.res
        self.rows = input_raster.shape[0]
        self.columns = input_raster.shape[1]
        self.CRS = 'unknown'
        self.data = []


