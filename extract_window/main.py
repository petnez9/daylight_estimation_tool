'''
***********************************************************************************************************************
*                                 extract windows module main script                                                  *
*                                         Created by P.Nezval                                                         *
*                                             version 0.1                                                             *
*                                               May 2021                                                              *
***********************************************************************************************************************
                                                License:
                                    Mozilla Public License Version 2.0
                                    ==================================
'''

from pre_processor import functions as pfun
from extract_window import win_functions as wfun

data = 'data2/win_lund3008_red.json'
data = pfun.read_json(data)
oid = pfun.get_obj_ids(data)

win_list = []
for i in oid:
    building = data['CityObjects'][i]
    win_geom = wfun.get_win_from_json(building, data)
    for window in win_geom:
        win_obj = wfun.Window(window, i)
        win_list.append(win_obj)

