import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d.art3d import Line3DCollection

def plot_gnm(coords, adj_mat):
    """
    Plots the elastic network model.

    :param: coords: List of coordinates for alpha carbon atoms
    :param: adj_mat: Adjacency matrix for elastic network model
    :param: zeros: Expected zero eigenvalues
    :return: None.
    """

    i, j = np.nonzero(adj_mat)
    coords = np.asarray(coords)

    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.plot(coords[:,0], coords[:,1], coords[:,2])

    # Create edges based on adjacency matrix
    segments = np.hstack((coords[i,None,:], coords[j,None,:]))
    lc = Line3DCollection(segments, color='red', linewidths=0.1)
    ax.add_collection(lc)

    plt.show()