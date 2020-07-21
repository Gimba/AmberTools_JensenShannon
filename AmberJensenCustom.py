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
import sys
from os.path import basename, splitext

import matplotlib

matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import pytraj as pt
from scipy.spatial.distance import jensenshannon


def main(args):
    parser = argparse.ArgumentParser(description=r'Calculates the Jensen-Shannon divergence between the probability '
                                                 r'distributions of backbone dihedral angles and generates .pml files'
                                                 r'that can be loaded using pymol.')
    parser.add_argument('top', help='topology file of first sturcture in prmtop format.')
    parser.add_argument('traj', help='trajectory file of first structure in netcdf format.')
    parser.add_argument('-d', dest='dihedral_angle',
                        help='A dihedral angles defined as in AmberTools, e.g. :26@CA,:27@CA,:34@CA,:50@CA\'')
    parser.add_argument('-p', dest='pymol_script', help='Name of pymol script file (default: top_1 + _js.pml)')
    parser.add_argument('-s', dest='scatter', help='Name of scatter plot file (default: top_1 + .png)')
    parser.add_argument('-f', dest='frames', help='Size of frame window', default="")

    args = parser.parse_args()

    # file name of first topology will be used for output file names
    top_1_name = splitext(basename(args.top))[0]

    traj = pt.iterload(args.traj, args.top)
    last_k = 0
    histograms = []
    x_ticks = []
    for k in range(int(args.frames), int(traj.n_frames), int(args.frames)):
        angles = pt.dihedral(traj, mask=args.dihedral_angle.replace(',', ' '), top=args.top, frame_indices=range(
            last_k, k + 1))
        hist = list(np.histogram(angles, bins=np.arange(-180, 181, 10))[0])
        histograms.append(hist)
        x_ticks.append(
            str(last_k // 1000) + "k -> " + str(k // 1000) + "k to " + str((k + 1) // 1000) + "k -> " + str((
                                                                                                             int(
                                                                                                                 args.frames)) // 1000) + "k")
        last_k = k + 1

    js_distances = []
    for c in range(len(histograms) - 1):
        js_distances.append(jensenshannon(histograms[c], histograms[c + 1]))

    print(js_distances)
    plt.scatter(x_ticks[:-1], js_distances)
    plt.xticks(x_ticks[:-1], rotation=30, ha='right')
    plt.ylabel('Jensen Shannon Divergence')
    plt.xlabel('Compared frame windows')
    plt.tight_layout()
    plt.savefig(top_1_name + ".png")


if __name__ == '__main__':
    main(sys.argv)
