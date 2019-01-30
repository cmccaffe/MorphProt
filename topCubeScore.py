import numpy as np
from scipy.signal import correlate2d
from skimage.transform import rotate

# for each face perform a cross correlation analysis
def cube_corr(protein_1, protein_2, attribute):
    rotation = np.arange(0,351,35)
    score_matrix = {}
    for i, face1 in enumerate(protein_1):
        for j, face2 in enumerate(protein_2):
            for angle in rotation:
                face2_rot = rotate(face2, angle)
                score = correlate2d(face1, face2_rot)
                score_matrix[i, j, angle] = score
    return score_matrix

def top_cube_corr(score_matrix, attribute):
    scores ={} # max scores for hydropathy
    if attribute == str(1):
        for k in score_matrix:
            min_value = np.amin(score_matrix[k])
            scores[k] = [min_value, np.unravel_index(np.argmin(score_matrix[k]), score_matrix[k].shape)]
            
    elif attribute == str(0)or str(2):     
        for k in score_matrix:
            max_value = np.amax(score_matrix[k])
            scores[k] = [max_value, np.unravel_index(np.argmax(score_matrix[k]), score_matrix[k].shape)]
    else:
        print("There is an error with the score extraction.")
        
    return scores


class topCubeScore:
    def __init__(self, protein1, protein2, attribute):
        '''Constructor for this class.'''
        self.protein1 = protein1
        self.protein2 = protein2
        self.attribute = attribute

    def Main(self):
        score_matrix = cube_corr(self.protein_1, self.protein_2, self.attribute)
        top_scores = top_cube_corr(score_matrix, self.attribute)
        return top_scores