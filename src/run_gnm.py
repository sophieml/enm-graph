import argparse
from utils import model, vis
from scipy.sparse.csgraph import laplacian

def main():
    parser = argparse.ArgumentParser(description='Build and visualize a ' \
                                    'Gaussian network model from a PDB file')
    parser.add_argument('filename', help='Path to PDB file')
    parser.add_argument('-d',
                        metavar='--dist',
                        help='Cutoff distance for edge distance (default: 7 ' \
                        'angstrom)',
                        type=float,
                        required=False,
                        default=7)
    args = parser.parse_args()

    calphas = model.get_calphas(args.filename)
    gph = model.get_adjacency_matrix(calphas, args.d)
    kirchoff = laplacian(gph)
    vis.plot_gnm(calphas, gph)

if __name__ == '__main__':
    main()