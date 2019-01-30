import numpy as np
from sklearn.cluster import KMeans
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from matplotlib import style
style.use("ggplot")
from collections import Counter
from mpl_toolkits.mplot3d import Axes3D
from pylab import *

def ClusterIndicesNumpy(clustNum, labels_array): #numpy 
    return np.where(labels_array == clustNum)[0]

def protein_processing(pdb, cuboid):

    x_coord = cuboid[6]
    y_coord = cuboid[7]
    z_coord = cuboid[8]
    charge = cuboid[9]

    X = np.array(x_coord, dtype=float)
    Y = np.array(y_coord, dtype=float)
    Z = np.array(z_coord, dtype=float)

    XYZ = np.column_stack((X, Y, Z))

    maxArg = np.argmax(XYZ, axis=0)
    minArg = np.argmin(XYZ, axis=0)

    start_centroids = np.array([XYZ[maxArg[0]], XYZ[minArg[0]], XYZ[maxArg[1]], XYZ[minArg[1]], XYZ[maxArg[2]], XYZ[minArg[2]]])

    ## using bias centroids to start

    cluster_num = 6
    kmeans = KMeans(n_clusters=6, init=start_centroids, random_state=0).fit(XYZ)
    centroids = kmeans.cluster_centers_
    labels = kmeans.labels_

    print ("centroids : ")
    print (centroids)
    print ("labels : ")
    print (labels)

    colors = ["g.","r.","c.","y."]
    color = ["g", "r", "b", "c", "y", "k"]
    c = Counter(labels)
    fig = figure(figsize=(10,10))
    ax = fig.gca(projection='3d')

    for i in range(len(XYZ)):
        ax.scatter(XYZ[i][0], XYZ[i][1], XYZ[i][2], c=color[labels[i]])

    for cluster_number in range(cluster_num):
        print("Cluster {} contains {} samples".format(cluster_number, c[cluster_number]))

    ax.scatter(centroids[:, 0],centroids[:, 1], centroids[:, 2], marker = "x", s=150, linewidths = 5, zorder = 100, c=color)
    plt.title(str(pdb))
    savefig(pdb[:-4]+'.pdf', bbox_inches='tight')
    
    clust0 = XYZ[ClusterIndicesNumpy(0, labels)]
    clust1 = XYZ[ClusterIndicesNumpy(1, labels)]
    clust2 = XYZ[ClusterIndicesNumpy(2, labels)]
    clust3 = XYZ[ClusterIndicesNumpy(3, labels)]
    clust4 = XYZ[ClusterIndicesNumpy(4, labels)]
    clust5 = XYZ[ClusterIndicesNumpy(5, labels)]
    
    # diveds cluster into 2D projections x projections are yz... etc
    # need to make this into a function

    x0 = []
    y0 = []
    z0 = []
    for i in clust0:
        x0.append(i[0])
        y0.append(i[1])
        z0.append(i[2])
    yz0 = np.vstack((y0,z0))

    x1 = []
    y1 = []
    z1 = []
    for i in clust1:
        x1.append(i[0])
        y1.append(i[1])
        z1.append(i[2])
    yz1 = np.vstack((y1,z1))

    x2 = []
    y2 = []
    z2 = []
    for i in clust2:
        x2.append(i[0])
        y2.append(i[1])
        z2.append(i[2])
    xz2 = np.vstack((x2,z2))

    x3 = []
    y3 = []
    z3 = []
    for i in clust3:
        x3.append(i[0])
        y3.append(i[1])
        z3.append(i[2])
    xz3 = np.vstack((x3,z3))

    x4 = []
    y4 = []
    z4 = []
    for i in clust4:
        x4.append(i[0])
        y4.append(i[1])
        z4.append(i[2])
    xy4 = np.vstack((x4,y4))

    x5 = []
    y5 = []
    z5 = []
    for i in clust5:
        x5.append(i[0])
        y5.append(i[1])
        z5.append(i[2])
    xy5 = np.vstack((x5,y5))
    return [yz0, yz1, xz2, xz3, xy4, xy5]


class shapeProcessing:
    def __init__(self, pdb, cuboid, clustNum, labels_array):
        '''Constructor for this class.'''
        self.pdb = pdb
        self.cuboid = cuboid
        self.clustNum = clustNum
        self.labels_array= labels_array

    def Main(self):
        xyzList = protein_processing(self.pdb, self.cuboid)
        return xyzList