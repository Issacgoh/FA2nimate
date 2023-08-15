# Graphimate
This is a fun package for creating iterated animations and movies of continuous graph layouts using the FA2 algorithm. The FA2 algorithm utilises Barnes Hut simulation, degree-dependent repulsive force, and local and global adaptive temperatures. 

![Simulation of the emergence of the first blood stem cells in human life from specialized endothelium (note that the mp4 outputs are smoother)]([https://raw.githubusercontent.com/YourUsername/YourRepo/main/path_to_gif/your_gif.gif](https://github.com/Issacgoh/Graphimate/blob/main/resources/Emergence_first_blood_stem_cells.mp4.gif))

## Starting notes:
This package takes as input:

  - An anndata object
  - A categorical data variable containing labels
  - A 2D array containing some XY dimensionality-reduced coordinates (PCA, UMAP, VAE, etc...)
  - A sparse 2D matrix (CSR) containing cell-cell weighted connectivities (KNN, etc...)

## Installation:
To install directly from github run below in command line

```bash
pip install git+git@github.com:Issacgoh/Graphimate.git
```

To clone and install:
```bash
git clone https://github.com/Issacgoh/Graphimate.git

cd ./graphimate

pip install .
```

## Usage notes:
Please see the example notebook under "/example_notebooks/"

## Graphing modules:
This package currently supports two custom iterative graph layout modules, 'FA' (forceatlas2) and 'FR' (Fruchterman-Reingold) modules.
More modules from the igraph collection will be added.
