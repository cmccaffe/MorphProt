import numpy as np

def matrix_mapping(face_1, face_2, min_loc):

    print("Indentifying minimum correlation and mapping correlation back to cube faces...")
    n1 = face_1.shape[1] - 1
    n2 = face_2.shape[1] - 1
    m1 = face_1.shape[0] - 1
    m2 = face_2.shape[0] - 1

    if n1 >= n2:
        if min_loc[1] >= n1:
            matrix_1_xint = (min_loc[1] - n2,n1 + 1)
            matrix_2_xint = (0,n2 - min_loc[1] + n1 +1)
        else:
            if min_loc[1] <= n2:
                matrix_1_xint = (0, min_loc[1] + 1)
                matrix_2_xint = (n2 - min_loc[1], n2 + 1)
            else:
                matrix_1_xint = (min_loc[1] - n2, min_loc[1] + 1)
                matrix_2_xint = (0, n2 + 1)
    else:
        if min_loc[1] >= n1:
            if min_loc[1] <= n2:
                matrix_1_xint = (0, n1 + 1)
                matrix_2_xint = (n2 - min_loc[1], n2 - min_loc[1] + n1 + 1)
            else:
                matrix_1_xint = (min_loc[1] - n2,n1 + 1)
                matrix_2_xint = (0,n2 + n1 - min_loc[1] + 1)
        else:
            matrix_1_xint = (0, n1 + n2 - min_loc[1] + 1)
            matrix_2_xint = (n2 - min_loc[1], n2 + 1)

    if m1 >= m2:
        if min_loc[0] >= m1:
            matrix_1_yint = (min_loc[0] - m2,m1 + 1)
            matrix_2_yint = (0,m2 - min_loc[0] + m1 + 1)
        else:
            if min_loc[0] <= m2:
                matrix_1_yint = (0, min_loc[0] + 1)
                matrix_2_yint = (m2 - min_loc[0], m2 + 1)
            else:
                matrix_1_yint = (min_loc[0] - m2, min_loc[0] + 1)
                matrix_2_yint = (0, m2 + 1)
    else:
        if min_loc[0] >= m1:
            if min_loc[0] <= m2:
                matrix_1_yint = (0, m1 + 1)
                matrix_2_yint = (m2 - min_loc[0], m2 - min_loc[0] + m1 + 1)
            else:
                matrix_1_yint = (min_loc[0] - m2,m1 + 1)
                matrix_2_yint = (0, m2 +  m1 - min_loc[0] + 1)
        else:
            matrix_1_yint = (0, n1 + n2 - min_loc[0] + 1)
            matrix_2_yint = (m2 - min_loc[0],m2 +1)


    matrix_1 = face_1[matrix_1_yint[0]:matrix_1_yint[1], matrix_1_xint[0]:matrix_1_xint[1]]
    matrix_2 = face_2[matrix_2_yint[0]:matrix_2_yint[1], matrix_2_xint[0]:matrix_2_xint[1]]
    
    print("Matrix 1 has constraints ", "(", matrix_1_yint[0], ":", matrix_1_yint[1], ",", matrix_1_xint[0], ":", matrix_1_xint[1], ")")
    print("Matrix 2 has constraints ", "(", matrix_2_yint[0], ":", matrix_2_yint[1], ",", matrix_2_xint[0], ":", matrix_2_xint[1], ")")
    
    return matrix_1, matrix_2


def matching_atomic_coords(out_matrix, prot_face, attribute):
    
    atom_coords_list = []
    
    if attribute == str(1): 
        for a in np.nditer(out_matrix):
            if a > 0:
                a = a - 1
                for k in prot_face[0]:
                    if np.round(a, 3) == np.round(prot_face[0][k][-1][1], 3):
                        for i in prot_face[0][k][:-1]:
                            atom_coords_list.append(i)

            elif a < 0:
                a = a + 1
                for k in prot_face[0]:
                    if np.round(a, 3) == np.round(prot_face[0][k][-1][1], 3):
                        for i in prot_face[0][k][:-1]:
                            atom_coords_list.append(i)
            else:
                a = a
                for k in prot_face[0]:
                    if np.round(a, 3) == np.round(prot_face[0][k][-1][1], 3):
                        for i in prot_face[0][k][:-1]:
                            atom_coords_list.append(i)
                            
    if attribute == str(0): 
        for a in np.nditer(out_matrix):
            if a > 0:
                a = a - 1
                for k in prot_face[0]:
                    if np.round(a, 3) == np.round(prot_face[0][k][-1][0], 3):
                        for i in prot_face[0][k][:-1]:
                            atom_coords_list.append(i)

            elif a < 0:
                a = a + 1
                for k in prot_face[0]:
                    if np.round(a, 3) == np.round(prot_face[0][k][-1][0], 3):
                        for i in prot_face[0][k][:-1]:
                            atom_coords_list.append(i)
            else:
                a = a
                for k in prot_face[0]:
                    if np.round(a, 3) == np.round(prot_face[0][k][-1][0], 3):
                        for i in prot_face[0][k][:-1]:
                            atom_coords_list.append(i)
    if attribute == str(2): 
        for a in np.nditer(out_matrix):
            for k in prot_face[0]:
                if np.round(a, 3) == np.round(prot_face[0][k][-1][2], 3):
                    for i in prot_face[0][k][:-1]:
                        atom_coords_list.append(i)
                            
    print("Returning the atomic coordinate list")
    return atom_coords_list

def cube_map_to_protein(cuboid, atom_coords_list, prot_face):
    if prot_face == 4 or prot_face == 5:
        for e in atom_coords_list:
            for k in cuboid[-1]:
                if str(e[1]) == str(k[4]):
                    e.append(k[:3])
    elif prot_face == 0 or prot_face == 1 or prot_face == 2 or prot_face == 3:
        for e in atom_coords_list:
            for k in cuboid[-1]:
                if str(e[1]) == str(k[5]):
                    e.append(k[:3])
    else:
        print("The protein face cannot be identified")
    aa_id = []
    for i in atom_coords_list:
        aa_id.append(i[-1])
    return aa_id

def sphere_map_to_protein(cuboid, atom_coords_list):

    for e in atom_coords_list:
        for k in cuboid[-1]:
            if str(e[1]) == str(k[5]):
                e.append(k[:3])
    
    return atom_coords_list



class mapBack:
    def __init__(self, face1, face2, min_loc, cuboid, prot_face, attribute):
        '''Constructor for this class.'''
        self.face1 = face1
        self.face2 = face2
        self.min_loc = min_loc
        self.cuboid = cuboid
        self.out_matrix = out_matrix
        self.prot_face = prot_face
        self.attribute = attribute

    def Main(self):
        out_matrix = matrix_mapping(self.face_1, self.face_2, self.min_loc)
        matching_atomic_coords = matching_atomic_coords(out_matrix, self.prot_face, self.attribute)
        cube_map_to_protein = cube_map_to_protein(self.cuboid, matching_atomic_coords, self.prot_face)
        sphere_map_to_protein = sphere_map_to_protein(self.cuboid, matching_atomic_coords, self.prot_face)
        return cube_map_to_protein, sphere_map_to_protein
