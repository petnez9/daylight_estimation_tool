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
cityModel_path = 'data2/win_lund3008_sample.json'

### here comes the path to shapefile with footprints and height attributes
foot_height_shp = "data2/ft_height.shp"

### here comes the path to roof raster (one of the outputs of UMEP)
roof_raster = 'data2/Energyyearroof.tif'

### here comes the path to the text file with energy values on the facades
facade_solar_energy = 'data2/Energyyearwall.txt'

### here comes the path to the updated CityJSON
updated_JSON = 'data2/city_windows_irr.json'