#! /usr/bin/env python

# Copyright (c) 2020 Martin Rosellen

# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated
# documentation files (the "Software"), to deal in the Software without restriction, including without limitation the
# rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit
# persons to whom the Software is furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the
# Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
# WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
# COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
# OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

import argparse
import matplotlib.pyplot as plt
import numpy as np
import pytraj as pt
import sys
from os.path import basename, splitext
from scipy.spatial.distance import jensenshannon


def load_trajectory(traj, top):
    traj = pt.load(traj, top)
    return traj


def calculate_angles(traj, residues=[], angles='phi psi'):
    if residues:
        data = list(pt.multidihedral(traj, dihedral_types=angles, resrange=residues))
    else:
        # get all residues numbers
        residues = [int(l[4:]) for l in list(pt.dssp(traj)[0])]
        data = pt.multidihedral(traj, dihedral_types=angles, resrange=residues)
    return data


def set_resrange(residues, traj):
    # use all residues in topology
    if not residues:
        residues = [int(l[4:]) for l in list(pt.dssp(traj)[0])]

    # range of residues
    elif '-' in residues:
        start, end = residues.split('-')
        residues = list(np.arange(int(start), int(end)))

    # list of residues
    else:
        residues = args.residues_2.split(',')


def main(args):
    parser = argparse.ArgumentParser(description=r'Calculates the Jensen-Shannon divergence between the probability '
                                                 r'distributions of backbone dihedral angles and generates .pml files'
                                                 r'that can be loaded using pymol.')
    parser.add_argument('top_1', help='topology file of first sturcture in prmtop format.')
    parser.add_argument('traj_1', help='trajectory file of first structure in netcdf format.')
    parser.add_argument('top_2', help='topology file of second structure in prmtop format.')
    parser.add_argument('traj_2', help='trajectory file of second structure in netcdf format.')
    parser.add_argument('-r1', dest='residues_1',
                        help='Dihedral angles are calculated for this range of residues for the first '
                             'structure, example: 1-156 or 1,2,3,4,... (default: all residues)', default=[])
    parser.add_argument('-r2', dest='residues_2',
                        help='Dihedral angles are calculated for this range of residues for the second structure, '
                             'example: 1-156 or 1,2,3,4,... (default: all residues)', default=[])
    parser.add_argument('-p', dest='pymol_script', help='Name of pymol script file (default: top_1 + _js.pml)')
    parser.add_argument('-s', dest='scatter', help='Name of scatter plot file (default: top_1 + .png)')

    args = parser.parse_args()

    # file name of first topology will be used for output file names
    top_1_name = splitext(basename(args.top_1))[0]

    # read in data
    traj_1 = load_trajectory(args.traj_1, args.top_1)
    traj_2 = load_trajectory(args.traj_2, args.top_2)

    # handle residue ranges
    residues_1 = set_resrange(args.residues_1, traj_1)
    residues_2 = set_resrange(args.residues_2, traj_2)

    # if residues differ in lenght only use common residues
    if len(args.residues_1) != len(list(set(residues_1) & set(residues_2))):
        residues = list(set(residues_1) & set(residues_2))
        residues_1 = residues
        residues_2 = residues

    residues_1.sort()
    residues_2.sort()

    angles_1 = calculate_angles(traj_1, residues_1)
    angles_2 = calculate_angles(traj_2, residues_2)

    # generate histograms of angle values
    js_distances = []
    for e, _ in enumerate(residues_1):
        hist_1 = list(np.histogram(angles_1[e], bins=np.arange(-180, 181, 10))[0])
        hist_2 = list(np.histogram(angles_2[e], bins=np.arange(-180, 181, 10))[0])
        js_distances.append(jensenshannon(hist_1, hist_2))

    # plot JS divergence
    plt.clf()
    plt.scatter(residues_1, js_distances)
    plt.xticks(residues_1)
    plt.xlabel('Residues')
    plt.ylabel('Jensen Shannon Divergence')
    plt.savefig(top_1_name + ".png")

    # write pdb that gets written in by pymol
    pt.write_traj(top_1_name + ".pdb", traj_1, frame_indices=[0], overwrite=True)

    # normalize values to range from zero to one for coloring
    js_distances_normalized = js_distances / max(js_distances)
    js_distances_normalized = (1. - js_distances_normalized)

    # generate output for .pymol input file
    out_pymol_file = 'load ' + top_1_name + '.pdb \n'
    out_pymol_file += 'remove solvent\n'
    out_pymol_file += 'remove name Cl-\n'
    out_pymol_file += 'hide lines\nshow cartoon\n'

    for r, j in zip(residues_1, js_distances_normalized):
        # coloring of residues according to JensenShannon divergence
        out_pymol_file += 'set_color residue_%s = [%f, %f, 1]\n' % (r, j, j)
        out_pymol_file += 'color residue_%s, resid %s\n' % (r, str(r))

    with open(top_1_name + "_js.pml", 'w') as o:
        o.write(out_pymol_file)


if __name__ == '__main__':
    main(sys.argv)
