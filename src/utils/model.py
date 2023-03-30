import numpy as np
from Bio.PDB import *
from scipy.spatial import distance_matrix

def get_calphas(filename):
    """
    Parses PDB file and returns coordinates of alpha carbon atoms.

    :param: filename: PDB filename
    :return: calphas: List of coordinates for alpha carbon atoms
    """

    parser = PDBParser()
    io = PDBIO()

    calphas = []
    strct = parser.get_structure('Y', filename)
    for atom in strct.get_atoms():
        if atom.name == 'CA':
            calphas.append(atom.coord)
    return calphas

def get_adjacency_matrix(coords, max_dist):
    """
    Returns an adjacancy matrix representation for a list of points.

    :param: coords: List of coordinates for alpha carbon atoms
    :param: max_dist: Maximum cutoff distance for edge between nodes (default: 10 angstrom)
    :return: mat: Adjacency matrix representation
    """

    dist = distance_matrix(coords, coords)
    # Calculate distance cutoff and clear the matrix diagonal
    mat = np.where(dist > max_dist, 0, 1) - np.identity(len(coords))

    return mat

def get_normal_modes(mat, n_modes, zeros=0):
    """
    Calculates normal modes for a square matrix by finding eigenvectors /
    eigenvalues. 

    :param: mat: NxN square matrix
    :param: n_modes: Number of lowest modes to return
    :param: zeros: Expected zero eigenvalues
    :return: eigenvals, eigenvecs: Eigenvalues and eigenvectors for n_modes lowest modes.
    """
    if (n_modes <= 0) or (n_modes > mat.shape[0]-zeros):
        raise ValueError('Invalid number of nodes requested for specified ' \
                         'matrix')

    eigenvals, eigenvecs = np.linalg.eig(mat)

    # Sort eigenvalues, (lowest to highest modes)
    idx = eigenvals.argsort()
    eigenvals = eigenvals[idx]
    eigenvecs = eigenvecs[:,idx]

    return eigenvals[zeros:n_modes+zeros], eigenvecs[:,zeros:n_modes+zeros]
