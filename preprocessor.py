#!/usr/bin/python

'''
***********************************************************************************************************************
*                                               Functions                                                  *
*                                         Created by P.Nezval                                                         *
*                                             version 0.1                                                             *
*                                             May 2021                                                                *
***********************************************************************************************************************
'''

import sys
import pre_processor.functions as fun


######################################################################################################################

if __name__ == "__main__":

    arg = sys.argv[1:]

    if arg[0] == '--help':
        print('Thank you for using the urban solar energy modelling tool v0.0.1 - preprocessor\n\n'
              'Usage of the tool is as following:\n'
              'preprocessor.py <path_to_cityJSON_input_file>.json <path_to_output_shapefile>.shp')
    else:
        if len(arg) == 2:
            fun.main(arg[0], arg[1])
        else:
            print('You have entered incorrect number of arguments.\n'
                  'Please use <script_name>.py --help to see the format of input.')
