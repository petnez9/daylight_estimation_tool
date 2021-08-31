'''
***********************************************************************************************************************
*                                Configuration file for Solar estimation tool                                         *
*                                         Created by P.Nezval                                                         *
*                                             version 0.1                                                             *
*                                         May 2021 - August 2021                                                      *
***********************************************************************************************************************
                                                License:
                                               MIT License
                                    ==================================
'''

### here comes the path to 3D city model including windows features
cityModel_path = 'data/win_lund3008_red.json'

### here comes the path to shapefile with footprints and height attributes
foot_height_shp = "data/ft_height.shp"

### here comes the path to roof raster (one of the outputs of UMEP)
roof_raster = 'data/Energyyearroof.tif'

### here comes the path to the text file with energy values on the facades
facade_solar_energy = 'data/Energyyearwall.txt'

### here comes the path to the updated CityJSON
updated_JSON = 'data/city_windows.json'

### here comes the path to the report
report_path = 'report2.json'