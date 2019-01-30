import Bio
from Bio import PDB
from Bio.PDB.Polypeptide import is_aa
from Bio.PDB.AbstractPropertyMap import AbstractPropertyMap
import numpy as np



# opens PDB or PQR files and returns the model

def pdbpqr_parse(name, file):
    pdb_parser = Bio.PDB.PDBParser(QUIET = True)
    ref_structure = pdb_parser.get_structure(name, file)
    model = ref_structure[0]
    print("PDB/PQR file is being parsed...")
    return model

# calculates the Relative ASA for each residue
# any residue with a RASA > 0 is extracted and used in the charge calculation
# returns a list of residue ids that are exposed to the surface
# plan to replace dssp with msms on the new computer


def ca_depth(residue, surface):
    if not residue.has_id("CA"):
        return None
    ca = residue["CA"]
    coord = ca.get_coord() 
    return min_dist(coord, surface)

def residue_depth(residue, surface):
    atom_list = residue.get_unpacked_list() 
    length = len(atom_list)
    d = 0
    for atom in atom_list: 
        coord = atom.get_coord() 
        d = d + min_dist(coord, surface)
    return d / length

def min_dist(coord, surface):
    """Return minimum distance between coord and surface.""" 
    d = surface - coord 
    d2 = np.sum(d * d, 1)
    return np.sqrt(min(d2))

def protein_shell(model):
    depth_dict = {}
    depth_list = []
    depth_keys = []
    # get_residue 
    residue_list = PDB.Selection.unfold_entities(model, 'R')
    # make surface from PDB file using MSMS
    surface = PDB.get_surface(model, MSMS='./msms')
    # calculate rdepth for each residue
    for residue in residue_list:
        if not is_aa(residue):
            continue 
        rd = residue_depth(residue, surface)
        ca_rd = ca_depth(residue, surface)
        # Get the key
        res_id = residue.get_id() 
        chain_id = residue.get_parent().get_id() 
        depth_dict[(chain_id, res_id)] = (rd, ca_rd)
        depth_list.append((residue, (rd, ca_rd)))
        depth_keys.append((chain_id, res_id))
        # Update xtra information 
        residue.xtra['EXP_RD'] = rd
        residue.xtra['EXP_RD_CA'] = ca_rd
    depth_id = []
    for res in depth_dict:
        if depth_dict[res][0] < 5:
            depth_id.append(res[1][1])
    return depth_id


class proteinShell:
    def __init__(self, name, pqr):
        '''Constructor for this class.'''
        self.name = name
        self.pqr = pqr

    def Main(self):
        model = pdbpqr_parse(self.name, self.pqr)
        shell = protein_shell(model)
        return shell
