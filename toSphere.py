import math
import numpy as np
import MorphProt
from MorphProt import proteinShell
from MorphProt import physicalProp
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import style
style.use("ggplot")
from pylab import *

def myround(x, base):
    return int(base * round(float(x)/base))

def cart2sph(x,y,z):
    XsqPlusYsq = x**2 + y**2
    r = math.sqrt(XsqPlusYsq + z**2)               # r
    theta = math.atan2(z,math.sqrt(XsqPlusYsq))     # theta
    phi = math.atan2(y,x)                           # phi
    return r, theta, phi

def to_sphere(cuboid, pqr):
    x_coord = cuboid[6]
    y_coord = cuboid[7]
    z_coord = cuboid[8]
    charge = cuboid[9]
    hydropathy = []
    mutrate = []
    aa_res = []
    aa_loc = []
    atom = []
    for i in cuboid[10]:
        hydropathy.append(i[-3])
        mutrate.append(i[-1])
        aa_res.append(i[0])
        aa_loc.append(i[1])
        atom.append(i[2])
    coords_w_prop = np.column_stack((x_coord, y_coord, z_coord, hydropathy, charge, mutrate))
    prot_centroid = (np.average(x_coord), np.average(y_coord), np.average(z_coord))
    
    x_dist = cuboid[0] - cuboid[1]
    y_dist = cuboid[2] - cuboid[3]
    z_dist = cuboid[4] - cuboid[5]
    dist_list = [x_dist, y_dist, z_dist]
    radius = myround((max(dist_list)/2) + 1, 1)
    
    sphere_proj = []
    for i in coords_w_prop:
        x = prot_centroid[0] + (radius/np.abs(math.sqrt((i[0] - prot_centroid[0])**2 + (i[1] - prot_centroid[1])**2 + (i[2] - prot_centroid[2])**2))* (i[0] - prot_centroid[0]))
        y = prot_centroid[1] + (radius/np.abs(math.sqrt((i[0] - prot_centroid[0])**2 + (i[1] - prot_centroid[1])**2 + (i[2] - prot_centroid[2])**2))* (i[1] - prot_centroid[1]))
        z = prot_centroid[2] + (radius/np.abs(math.sqrt((i[0] - prot_centroid[0])**2 + (i[1] - prot_centroid[1])**2 + (i[2] - prot_centroid[2])**2))* (i[2] - prot_centroid[2]))
        sphere_data = [x, y, z, i[3], i[4], i[5]]
        sphere_proj.append(sphere_data)
        
    sphere_avg = np.mean(sphere_proj, axis=0)
    adjusted_sphere_proj = []
    for i in sphere_proj:
        x = i[0] - sphere_avg[0]
        y = i[1] - sphere_avg[1]
        z = i[2] - sphere_avg[2]
        sphere_data = [x, y, z, i[3], i[4], i[5]]
        adjusted_sphere_proj.append(sphere_data)
    
    plt.rcParams['grid.color'] = "black"    
    fig = plt.figure(figsize=(12,9))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_facecolor('xkcd:white')
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    #plt.rcParams['grid.color'] = "black"
    for i in adjusted_sphere_proj:
        ax.scatter(i[0], i[1], i[2], zdir='x')
    plt.show()    
    #plt.title(str(pqr))
    #savefig(pqr[:-4]+'.pdf', bbox_inches='tight')
    
    polar_sphere = []
    for i in adjusted_sphere_proj:
        polar_coord = cart2sph(i[0], i[1], i[2])
        polar_data = [polar_coord[1], polar_coord[2], i[3], i[4], i[5]] 
        polar_sphere.append(polar_data)
    
    plane_proj = []
    for i in polar_sphere:
        x = radius * i[1]
        y = radius * math.log(math.tan((i[0] + math.pi/2)/2))
        plane_data = [x, y, i[2], i[3], i[4]]
        plane_proj.append(plane_data)  
    plane_proj_data = list(zip(plane_proj, aa_res, aa_loc, atom))
    
    fig = plt.figure(figsize=(12,9))
    ax = fig.add_subplot(111)
    ax.set_facecolor('xkcd:white')
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    #plt.rcParams['grid.color'] = "black"
    
    for i in polar_sphere:
        ax.scatter(radius * i[1], radius * math.log(math.tan((i[0] + math.pi/2)/2)))
    plt.show()  
    return plane_proj_data
        
class faceBuilder:
    def __init__(self, cuboid, pqr):
        '''Constructor for this class.'''
        self.cuboid = cuboid
        self.pqr = pqr

    def Main(self):
        plane_proj = to_sphere(self.cuboid, self.pqr)
        return plane_proj

