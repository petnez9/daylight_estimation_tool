'''
***********************************************************************************************************************
*                                   Building height extractor tool                                                    *
*                                         Created by P.Nezval                                                         *
*                                             version 0.1                                                             *
*                                             March 2021                                                              *
***********************************************************************************************************************
                                                License:
                                    Mozilla Public License Version 2.0
                                    ==================================
'''

import json

######################################################################################################################


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


def print_build_height(input_data, oid):
    out = {}
    out['buildings'] = []
    j = 0
    for i in oid:
        building = input_data['CityObjects'][i]
        ver_list = get_vertex_list(building)
        ver_heights = vertex_heights(input_data)
        build_height = get_build_height(ver_list, ver_heights)
        k = 'property_{}'.format(j)
        x = {k: {'fid': i, 'height': build_height}}
        out['buildings'].append(x)
        j = j + 1
    with open('/home/petnez9/Documents/Master_Thesis/practical/data/building_heights.json', 'w') as outfile:
        json.dump(out, outfile)
    return out


def main(path):
    with open(path) as file:
        data = json.load(file)
        file.close()
    obj_id = get_obj_ids(data)
    result = print_build_height(data, obj_id)


if __name__ == "__main__":
    json_file = "/home/petnez9/Documents/Master_Thesis/practical/data/study_area.json"
    main(json_file)