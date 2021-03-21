'''
***********************************************************************************************************************
*                                   Shapefile update attributes script                                                *
*                                         Created by P.Nezval                                                         *
*                                             version 0.1                                                             *
*                                             March 2021                                                              *
***********************************************************************************************************************
                                                License:
                                    Mozilla Public License Version 2.0
                                    ==================================
'''

import geopandas
import pandas
import json


def create_attr_table(input_data):
    fids = []
    heights = []
    for i in range(len(input_data['building'])):
        nm = 'property_{}'.format(i)
        fids.append(input_data['building'][i][nm]['fid'])
        heights.append(input_data['building'][i][nm]['height'])
    dat = {'fid': fids, 'height': heights}
    updated = pandas.DataFrame(dat)
    return updated


def updateSHP(shape,attr):
    updated = attr.merge(shape, on='fid')
    gdf = geopandas.GeoDataFrame(updated, geometry=updated.geometry)
    return gdf


def main():
    shp = geopandas.read_file(
        '/home/petnez9/Documents/Master_Thesis/practical/data/study_area/study_area_from_cityJSON.shp')
    with open('/home/petnez9/Documents/Master_Thesis/practical/data/building_heights.json') as file:
        data = json.load(file)
        file.close()
    attr = create_attr_table(data)
    shape = updateSHP(shp, attr)
    shape.to_file('/home/petnez9/Documents/Master_Thesis/practical/data/study_area/study_area_updated.shp')


if __name__ == "__main__":
    main()