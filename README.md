# AmberTools_JensenShannon
Script that calculates the Jensen-Shannon divergence between pairs of dihedral angles (default: phi, psi). If not specified this is done for all residues. The Jensen-Shannon divergence values are plotted in a scatter plot and further a pymol script is created that colors residues in shades of blue that indicate high divergence.

Requirements:
matplotlib
pytraj
scipy
numpy

Jensen Shannon divergence values overlaid onto structure in pymol:
![](https://github.com/Gimba/AmberTools_JensenShannon/blob/master/pymol_overlay.png)


Jensen Shannon divergence values as scatterplot:
![](https://github.com/Gimba/AmberTools_JensenShannon/blob/master/scatter_plot.png)
