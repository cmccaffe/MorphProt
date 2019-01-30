import Bio
from Bio import PDB
from Bio.PDB.Polypeptide import is_aa
from Bio.PDB.AbstractPropertyMap import AbstractPropertyMap

'''Inserts the charge and hydropathy information into a list of the coordinates'''


def cuboid_dimensions(model, surface_residues, mut_rate_dict=None):
    X_coord = []
    Y_coord = []
    Z_coord = []
    charge = []
    atom_data = []
    res_hydropathy = {'ILE': 0.31, 'VAL': -0.07, 'LEU': 0.56, 'PHE': 1.13, 'CYS': 0.24, 'MET': 0.23, 'ALA': -0.17, 'GLY': -0.01, 'THR': -0.14, 'SER': -0.13, 'TRP': 1.85, 'TYR': 0.94, 'PRO': -0.45, 'HIS': -0.96, 'GLU': -2.02, 'GLN': -0.58, 'ASP': -1.23, 'ASN': -0.42, 'LYS': -0.99, 'ARG': -0.81}
    print('Parsing atomic data...')

    for chain in model:
        for residue in chain:
            y = Bio.PDB.is_aa(residue)
            if y == True:
                aa_res = residue
                res_num = aa_res.get_id()[1]
                if res_num in surface_residues:
                    for atom in aa_res:
                        XYZ = atom.get_vector()
                        X = XYZ[0]
                        Y = XYZ[1]
                        Z = XYZ[2]
                        X_coord.append(X)
                        Y_coord.append(Y)
                        Z_coord.append(Z)
                        charge.append(atom.get_occupancy())
                        if mut_rate_dict is None:
                            null_mut_rate_dict = {'ILE': 0, 'VAL': 0, 'LEU': 0, 'PHE': 0, 'CYS': 0, 'MET': 0, 'ALA': 0, 'GLY': 0, 'THR': 0, 'SER': 0, 'TRP': 0, 'TYR': 0, 'PRO': 0, 'HIS': 0, 'GLU': 0, 'GLN': 0, 'ASP': 0, 'ASN': 0, 'LYS': 0, 'ARG': 0}
                            value = [aa_res.get_resname(), res_num, atom.get_fullname(), X, Y, Z, res_hydropathy[str(aa_res.get_resname())], atom.get_occupancy(), null_mut_rate_dict[str(aa_res.get_resname())]]  
                            atom_data.append(value)
                        else:
                            try:
                                mut_rate = int(mut_rate_dict[res_num])
                            except KeyError:
                                continue
                            value = [aa_res.get_resname(), res_num, atom.get_fullname(), X, Y, Z, res_hydropathy[str(aa_res.get_resname())], atom.get_occupancy(), mut_rate]
                            atom_data.append(value)                        

    # identifies the max/min of xyz coordinates            
    max_X = max(X_coord)
    min_X = min(X_coord)
    max_Y = max(Y_coord)
    min_Y = min(Y_coord)
    max_Z = max(Z_coord)
    min_Z = min(Z_coord)
    return max_X, min_X, max_Y, min_Y, max_Z, min_Z, X_coord, Y_coord, Z_coord, charge, atom_data



class physicalProp:
    def __init__(self, model, surface_residues):
        '''Constructor for this class.'''
        self.model = model

    def Main(self):
        cuboid = physicalProp(self.model, self.surface_residues)
        return cuboid