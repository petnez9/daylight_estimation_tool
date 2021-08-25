'''
***********************************************************************************************************************
*                                Classes file for solar estimation tool                                               *
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
                                              129835.348124142,
                                              6175671.146770532,
                                              0,
                                              130297.167365228,
                                              6176468.73581959,
                                              50
                                            ]

    def printOut(self):
        json_out = {"type": self.type,
                    "version": self.version,
                    #"metadata": self.metadata,
                    "CityObjects": self.cityObjects,
                    "vertices": self.vertices}
                    #"appearance": self.appearance}
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
        boundaries, values = self.resolveRoof(building.footprints, building.building_height, boundaries, values, vertices)
        boundaries, values = self.resolveFacades(building.facade_list, building.building_height, boundaries, values, vertices)
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
    def __init__(self):
        self.material = {}

    def populateMaterialValues(self, names, building, bound):
        val1, val2, val3, val4 = [], [], [], []
        val1.extend([None] * len(bound))
        val2.extend([None] * len(bound))
        val3.extend([None] * len(bound))
        val4.extend([None] * len(bound))
        print("Building {a} has this many boundaries: {b} and this many windows: {c}".format(a=building.fid,b=len(bound),c=len(building.windows_geometry)))
        for window in building.windows_geometry:
            if window.total_irradiance < 750:
                val1[window.bound_id] = 0
            elif window.total_irradiance > 750 and window.total_irradiance < 1500:
                val2[window.bound_id] = 1
            elif window.total_irradiance > 1500 and window.total_irradiance < 2250:
                val3[window.bound_id] = 2
            else:
                val4[window.bound_id] = 3
        ll1 = [x for x in val1 if x is not None]
        ll2 = [x for x in val2 if x is not None]
        ll3 = [x for x in val3 if x is not None]
        ll4 = [x for x in val4 if x is not None]
        x = len(ll1) + len(ll2) + len(ll3) + len(ll4)
        print("Irr1 has {a} Irr2 has {b} Irr3 has {c} Irr4 has {d}"
              " Together it is {e} features n of unassigned windows is {f}\n"
              .format(a=len(ll1), b=len(ll2),c=len(ll3),d=len(ll4),e=x,f=len(building.windows_geometry)-x))
        self.material['{}'.format(names[0])] = {'values': [val1]}
        self.material['{}'.format(names[1])] = {'values': [val2]}
        self.material['{}'.format(names[2])] = {'values': [val3]}
        self.material['{}'.format(names[3])] = {'values': [val4]}