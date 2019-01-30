'''Gives 6 lists of lists, each contains the coordinates for a face and the physical properties stored within'''
import numpy as np

def coord_charge_extract(axis_in, cuboid, axis_check):
    #atoms_coords_charge = np.column_stack((cuboid[10], cuboid[9]))
    atoms_coords_charge = cuboid[10]
    acc = []
    pts = []
    for info in atoms_coords_charge:
        acc.append(list(info))
    for inf in axis_in:
        pts.append(list(inf))

    coord_charge = []
    if axis_check.upper() == "YZ": 
        for a in pts:
            for b in acc:
                first = str(a[1])[:7]
                second = str(b[5])[:7]
                if first == second:
                    coord_charge.append(b[4:])
    if axis_check.upper() == "XZ": 
        #print("XZ")
        for a in pts:
            for b in acc:
                first = str(a[0])[:7]
                second = str(b[3])[:7]
                #print(first,second)
                if first == second:
                    coord_charge.append([b[3], b[5], b[6], b[7], b[8]])
    if axis_check.upper() == "XY": 
        for a in pts:
            for b in acc:
                first = str(a[0])[:7]
                second = str(b[3])[:7]
                if first == second:
                    coord_charge.append([b[3], b[4], b[6], b[7], b[8]])
    return coord_charge


def coord_charges(xyz_list, cuboid):
    coord_charge_list = []
    
    for i in range(2):
        YZ = np.column_stack((xyz_list[i][0], xyz_list[i][1]))
        coord_charge_list.append(coord_charge_extract(YZ, cuboid, "YZ"))
        
    for i in range(2):
        XZ = np.column_stack((xyz_list[i+2][0], xyz_list[i+2][1]))
        coord_charge_list.append(coord_charge_extract(XZ, cuboid, "XZ"))
        
    for i in range(2):
        XY = np.column_stack((xyz_list[i+4][0], xyz_list[i+4][1]))
        coord_charge_list.append(coord_charge_extract(XY, cuboid, "XY"))
        
    return coord_charge_list

class cubeFaceData:
    def __init__(self, xyz_list, cuboid):
        '''Constructor for this class.'''
        self.xyz_list = xyz_list
        self.cuboid = cuboid

    def Main(self):
        coord_charge_list = coord_charges(self.xyz_list, self.cuboid)
        return coord_charge_list