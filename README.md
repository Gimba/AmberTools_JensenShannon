# AmberTools_JensenShannon
Script that calculates the Jensen-Shannon divergence between pairs of dihedral angles (default: phi, psi). If not specified this is done for all residues. The Jensen-Shannon divergence values are plotted in a scatter plot and further a pymol script is created that colors residues in shades of blue that indicate high divergence.

Requirements:
matplotlib
pytraj
scipy
numpy
```
  usage: AmberJensenShannonTool.py [-h] [-r1 RESIDUES_1] [-r2 RESIDUES_2]
                                   [-p PYMOL_SCRIPT] [-s SCATTER]
                                   top_1 traj_1 top_2 traj_2

  Calculates the Jensen-Shannon divergence between the probability distributions
  of backbone dihedral angles and generates .pml filesthat can be loaded using
  pymol.

  positional arguments:
    top_1            topology file of first sturcture in prmtop format.
    traj_1           trajectory file of first structure in netcdf format.
    top_2            topology file of second structure in prmtop format.
    traj_2           trajectory file of second structure in netcdf format.

  optional arguments:
    -h, --help       show this help message and exit
    -r1 RESIDUES_1   Dihedral angles are calculated for this range of residues
                     for the first structure, example: 1-156 or 1,2,3,4,...
                     (default: all residues)
    -r2 RESIDUES_2   Dihedral angles are calculated for this range of residues
                     for the second structure, example: 1-156 or 1,2,3,4,...
                     (default: all residues)
    -p PYMOL_SCRIPT  Name of pymol script file (default: top_1 + _js.pml)
    -s SCATTER       Name of scatter plot file (default: top_1 + .png)
```
Jensen Shannon divergence values overlaid onto structure in pymol:
![](https://github.com/Gimba/AmberTools_JensenShannon/blob/master/pymol_overlay.png)


Jensen Shannon divergence values as scatterplot:
![](https://github.com/Gimba/AmberTools_JensenShannon/blob/master/scatter_plot.png)
