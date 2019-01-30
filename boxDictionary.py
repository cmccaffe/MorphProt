import numpy as np
import itertools
import random

# returns modular division base 5, used to create grids of size 5

def myround(x, base):
    return int(base * round(float(x)/base))

# pairs list of values into interval of values

def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = itertools.tee(iterable)
    next(b, None)
    interval = list(zip(a, b))
    return interval

def avg_charge(dictionary,k):
    pointsInBox = dictionary[str(k)]
    C = []
    D = []
    E = []
    for point in pointsInBox:
        C.append(float(point[2]))
        D.append(float(point[3]))
        E.append(float(point[4]))
    mean_C = np.mean(C)
    mean_D = np.mean(D)
    mean_E = np.mean(E)
    return str(k), mean_C, mean_D, mean_E

# reshuffles the values in each face
def shuffle2d(arr2d, rand=random):
    """Shuffes entries of 2-d array arr2d, preserving shape."""
    reshape = []
    data = []
    iend = 0
    for row in arr2d:
        data.extend(row)
        istart, iend = iend, iend+len(row)
        reshape.append((istart, iend))
    rand.shuffle(data)
    return [data[istart:iend] for (istart,iend) in reshape]

def shuffled_dict(box_dict):
    data = list(box_dict.values())
    shuffled_data = shuffle2d(data)
    for x,k in zip(shuffled_data, box_dict):
        box_dict[k][:] = x
    return box_dict

# set min and max for 2d square cubes as intervals of 5

def box_dictionary(axis, coord_charge, resolution, shuffle):
    xmin = min(axis[0])
    xmax = max(axis[0])
    ymin = min(axis[1])
    ymax = max(axis[1])

    xmin_r = myround(xmin, base=resolution)
    xmax_r = myround(xmax, base=resolution)
    ymin_r = myround(ymin, base=resolution)
    ymax_r = myround(ymax, base=resolution)

    if xmin < (xmin_r):    
        xmin = xmin_r - resolution
    else:
        xmin = xmin_r

    if xmax > xmax_r:    
        xmax = xmax_r + resolution
    else:
        xmax = xmax_r

    if ymax > ymax_r:   
        ymax = ymax_r + resolution
    else:
        ymax = ymax_r

    if ymin < ymin_r:    
        ymin = ymin_r - resolution
    else:
        ymin = ymin_r
    # get box intervals
    stepsize = resolution
    x_steps = int((xmax-xmin)/stepsize)
    y_steps = int((ymax-ymin)/stepsize)

    X_values = []
    Y_values = []

    for i in range(0,x_steps+1):
        X_values.append(xmax - i*stepsize)
    for i in range(0,y_steps+1):
        Y_values.append(ymax - i*stepsize)

    X_interval = pairwise(X_values)
    Y_interval = pairwise(Y_values)

    box_dimen = list(itertools.product(X_interval, Y_interval))

    box_name = []
    for i in range(1,len(box_dimen)+1):
        box = 'box_' + str(i)
        box_name.append(box)

    box_values = []

    for i in box_dimen:
        box =  []
        X = i[0]
        Y = i[1]
        for a in coord_charge:
            if X[1] < float(a[0]) < X[0] and Y[1] < float(a[1]) < Y[0]:
                box.append(a)
        box_values.append(box)


    box_dict = dict(zip(box_name, box_values))
    if shuffle == 'on':
        shuffled_dict(box_dict)
        for k in box_dict:
            box_dict[k].append(avg_charge(box_dict, k)[1:])
    else:
        for k in box_dict:
            box_dict[k].append(avg_charge(box_dict, k)[1:])
       
    return(box_dict, X_values, Y_values)

class boxDictionary:
    def __init__(self, axis, coord_charge, resolution, shuffle):
        '''Constructor for this class.'''
        self.axis = axis
        self.coord_charge = coord_charge
        self.shuffle = shuffle
        self.resolution = resolution

    def Main(self):
        prot_face = box_dictionary(self.axis, self.coord_charge, self.resolution, self.shuffle)
        return prot_face
