'''
***********************************************************************************************************************
*                                Classes file for daylight estimation tool                                            *
*                                         Created by P.Nezval                                                         *
*                                             version 0.1                                                             *
*                                         May 2021 - August 2021                                                      *
***********************************************************************************************************************
                                                License:
                                               MIT License
                                    ==================================
'''

class OutputCityJSON:

    def __init__(self):
        self.metadata = {"geographicalExtent": [], "referenceSystem": "EPSG:3008"}
        self.type = "CityJSON"
        self.version = "1.0"
        self.cityObjects = {}
        self.vertices = []
        self.appearance = {}

    def populateCityObjects(self, building, mat_name, mat_color):
        city_obj = CityObject()
        geom = geometry()
        mat = Material(mat_name, mat_color)
        mat.printEachMaterial()
        geom.populateObject(building, self.vertices)
        city_obj.geometry.append(geom.printOut())
        self.cityObjects['{a}'.format(a=building.fid)] = city_obj.printOut()
        self.appearance['materials'] = mat.materials

    def getExtent(self):
        x = [x[0] for x in self.vertices]
        y = [x[1] for x in self.vertices]
        z = [x[2] for x in self.vertices]
        self.metadata["geographicalExtent"] = [
            max(x),
            6175671.146770532,
            0,
            130297.167365228,
            6176468.73581959,
            50
        ]

    def printOut(self):
        json_out = {"type": self.type,
                    "version": self.version,
                    # "metadata": self.metadata,
                    "CityObjects": self.cityObjects,
                    "vertices": self.vertices}
        # "appearance": self.appearance}
        return json_out


class CityObject:

    def __init__(self):
        self.type = "Building"
        self.attributes = {}
        self.geometry = []

    def printOut(self):
        object_out = {"type": self.type, "attributes": self.attributes, "geometry": self.geometry}
        return object_out


class geometry:

    def __init__(self):
        self.type = "Solid"
        self.boundaries = []
        self.semantics = {'surfaces': [{"type": "RoofSurface"},
                                       {"type": "GroundSurface"},
                                       {"type": "WallSurface"},
                                       {"type": "Window"}],
                          'values': []
                          }
        self.material = {}
        self.lod = 3

    def resolveFootprints(self, foot, bound, val, ver):
        curr_bound = []
        for i in range(len(foot) - 1):
            x = [foot[i][0], foot[i][1], 0]
            if x not in ver:
                ver.append(x)
            ver_id = ver.index(x)
            curr_bound.append(ver_id)
        bound.append(curr_bound)
        val.append(1)
        return bound, val

    def resolveFacades(self, fac, height, bound, val, ver):
        for facade in fac:
            curr_bound = []
            x = [facade.start[0], facade.start[1], 0]
            if x not in ver:
                ver.append(x)
            ver_id = ver.index(x)
            curr_bound.append(ver_id)
            x = [facade.end[0], facade.end[1], 0]
            if x not in ver:
                ver.append(x)
            ver_id = ver.index(x)
            curr_bound.append(ver_id)
            x = [facade.end[0], facade.end[1], height]
            if x not in ver:
                ver.append(x)
            ver_id = ver.index(x)
            curr_bound.append(ver_id)
            x = [facade.start[0], facade.start[1], height]
            if x not in ver:
                ver.append(x)
            ver_id = ver.index(x)
            curr_bound.append(ver_id)
            bound.append(curr_bound)
            val.append(2)
        return bound, val

    def resolveRoof(self, foot, height, bound, val, ver):
        curr_bound = []
        roof = foot.copy()
        for i in range(len(roof) - 1):
            x = [foot[i][0], foot[i][1], height]
            if x not in ver:
                ver.append(x)
            ver_id = ver.index(x)
            curr_bound.append(ver_id)
        bound.append(curr_bound)
        val.append(0)
        return bound, val

    def populateObject(self, building, vertices):
        boundaries = []
        values = []
        boundaries, values = self.resolveFootprints(building.footprints, boundaries, values, vertices)
        boundaries, values = self.resolveRoof(building.footprints, building.building_height, boundaries, values,
                                              vertices)
        boundaries, values = self.resolveFacades(building.facade_list, building.building_height, boundaries, values,
                                                 vertices)
        boundaries.append(boundaries.pop(1))
        values.append(values.pop(1))
        self.boundaries.append([boundaries])
        self.semantics['values'].append(values)

    def printOut(self):
        object_out = {"type": self.type,
                      "boundaries": self.boundaries,
                      "semantics": self.semantics,
                      "material": self.material,
                      "lod": self.lod}
        return object_out


class Material:

    def __init__(self, names, color_arrays):
        self.name = names
        self.ambientIntensity = 0.2000
        self.diffuseColor = color_arrays
        self.materials = []

    def printEachMaterial(self):
        for i in range(len(self.name)):
            object_out = {"name": self.name[i],
                          "ambientIntensity": self.ambientIntensity,
                          "diffuseColor": self.diffuseColor[i]}
            self.materials.append(object_out)


class MaterialValues:

    def __init__(self, threshold):
        self.material = {}
        self.m_values = []
        self.threshold = threshold

    def populateMaterialValues(self, names, building, bound):
        val1, val2, val3, val4 = [], [], [], []
        val1.extend([None] * len(bound))
        val2.extend([None] * len(bound))
        val3.extend([None] * len(bound))
        val4.extend([None] * len(bound))
        self.m_values.extend([None] * len(bound))
        for window in building.windows_geometry:
            if window.total_irradiance < self.threshold[0]:
                val1[window.bound_id] = 0
                self.m_values[window.bound_id] = 0
                window.suitability = "Not suitable"
            elif window.total_irradiance > self.threshold[0] and window.total_irradiance < self.threshold[1]:
                val2[window.bound_id] = 1
                self.m_values[window.bound_id] = 1
                window.suitability = "Fairly suitable"
            elif window.total_irradiance > self.threshold[1] and window.total_irradiance < self.threshold[2]:
                val3[window.bound_id] = 2
                self.m_values[window.bound_id] = 2
                window.suitability = "Well suitable"
            else:
                val4[window.bound_id] = 3
                #self.m_values[window.bound_id] = 3
                window.suitability = "Highly suitable"
        self.material['{}'.format(names[0])] = {'values': [val1]}
        self.material['{}'.format(names[1])] = {'values': [val2]}
        self.material['{}'.format(names[2])] = {'values': [val3]}
        self.material['{}'.format(names[3])] = {'values': [val4]}


class SemanticAttributes:

    def __init__(self, building):
        self.building = building
        self.attributes = []
        self.surfaces = []

    def populateSemanticSurface(self, values):
        for i in range(len(values)):
            if values[i] == 0:
                self.surfaces.append({"type": "RoofSurface",
                                        "total_irradiance": float(self.building.roof.total_irradiance),
                                        "total_irradiance_per_area": float(self.building.roof.total_irradiance_per_m)})

            elif values[i] == 1:
                win_not_found = True
                for window in self.building.windows_geometry:
                    if i == window.bound_id:
                        win_not_found = False
                        self.surfaces.append({"type": "Window",
                                              "total_irradiance": float(window.total_irradiance),
                                              "total_irradiance_per_area": float(window.total_irradiance_per_m)})
                if win_not_found:
                    self.surfaces.append({"type": "Window"})

            elif values[i] == 2:
                self.surfaces.append({"type": "GroundSurface"})

            elif values[i] == 3:
                self.surfaces.append({"type": "WallSurface"})

            elif values[i] == 4:
                self.surfaces.append({"type": "Door"})

'''
    def populateSemanticValues(self, values):
        for i in range(len(values)):
            if values[i] == 0:
                self.attributes.append({"attributes": { "total_irradiance": float(self.building.roof.total_irradiance),
                                                        "total_irradiance_per_area": float(self.building.roof.total_irradiance_per_m)}})
            elif values[i] == 1:
                for window in self.building.windows_geometry:
                    if i == window.bound_id:
                        self.attributes.append({"attributes": {"total_irradiance": float(window.total_irradiance),
                                                               "total_irradiance_per_area": float(window.total_irradiance_per_m)}})
            elif values[i] == 2:
                self.attributes.append(None)
            elif values[i] == 3:
                self.attributes.append(None)
            elif values[i] == 4:
                self.attributes.append(None)
'''


class ReportGenerator:

    def __init__(self, build_list):
        self.report = {}
        self.build_list = build_list

    def print_report(self):
        buil_list = []
        for building in self.build_list:
            buil_list.append(self.print_building(building))
        self.report = {'CityObjects': buil_list}

    def print_building(self, building):
        entry = {'facades': self.print_facades(building.facade_list),
                 'roof': self.print_roof(building.roof),
                 'windows': self.print_windows(building.windows_geometry)}
        return {'Building': building.fid, 'building_data': entry}

    def print_facades(self, facade_list):
        fac_list = []
        for facade in facade_list:
            entry = {'facade_id': facade.facade_index,
                     'facade_area': float(facade.area),
                     'facade_total_irradiance': float(facade.total_irradiance),
                     'facade_total_irradiance_per_area': float(facade.total_irradiance_per_m)}
            fac_list.append(entry)
        return fac_list

    def print_roof(self, roof):
        entry = {'Roof_area': float(roof.roof_area),
                 'roof_total_irradiance': float(roof.total_irradiance),
                 'roof_total_irradiance_per_area': float(roof.total_irradiance_per_m)}
        return entry

    def print_windows(self, windows):
        win_list = []
        for window in windows:
            entry = {'window_id': window.id,
                     'window_area': float(window.height * window.width),
                     'window_total_irradiance': float(window.total_irradiance),
                     'window_total_irradiance_per_area': float(window.total_irradiance_per_m),
                     'window_suitability': window.suitability}
            win_list.append(entry)
        return win_list