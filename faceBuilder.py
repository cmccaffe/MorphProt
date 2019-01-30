import numpy as np

def protein_face_array_normalized(face):
    dictionary = face[0]
    length1 = len(face[1])-1
    length2 = len(face[2])-1
    protein_chrg_list = []
    protein_hydro_list = []
    protein_mutrate_list = []
    
    for k in dictionary:
        avg_hydro = float(dictionary[k][-1][0])
        if avg_hydro > 0:
            avg_hydro = avg_hydro + 1
        if avg_hydro < 0:
            avg_hydro = avg_hydro - 1
        
        avg_chrg = float(dictionary[k][-1][1])
        if avg_chrg > 0:
            avg_chrg = avg_chrg + 1
        if avg_chrg < 0:
            avg_chrg = avg_chrg - 1
            
        avg_mutrate = float(dictionary[k][-1][2])

        protein_chrg_list.append(avg_chrg)
        protein_hydro_list.append(avg_hydro)
        protein_mutrate_list.append(avg_mutrate)

    protein_chrg_array = np.array(protein_chrg_list)
    protein_chrg_array = protein_chrg_array.reshape(length1, length2)
    protein_hydro_array = np.array(protein_hydro_list)
    protein_hydro_array = protein_hydro_array.reshape(length1, length2)
    protein_mutrate_array = np.array(protein_mutrate_list)
    protein_mutrate_array = protein_mutrate_array.reshape(length1, length2)
    
    return protein_hydro_array, protein_chrg_array, protein_mutrate_array

# attribute 0 is hydropathy and attribute 1 is charge

def face_normalization(prot_face_list, attribute, shape):
    if shape == 'cube':
        prot = []
        for i in range(6):
            protein_face_normalized = np.nan_to_num(protein_face_array_normalized(prot_face_list[i])[int(attribute)])
            prot.append(protein_face_normalized)

        return prot
    
    elif shape == 'sphere':
        protein_face_normalized = np.nan_to_num(protein_face_array_normalized(prot_face_list)[int(attribute)]) 

        return protein_face_normalized

class faceBuilder:
    def __init__(self, prot_face_list, attribute, shape):
        '''Constructor for this class.'''
        self.prot_face_list = prot_face_list
        self.attribute = attribute

    def Main(self):
        prot1 = face_normalization(self.prot_face_list, self.attribute, self.shape)
        return prot1
