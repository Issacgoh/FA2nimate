# Graphimate
This is a fun package for creating iterated animations and movies of continuous graph layouts using the FA2 algorithm. The FA2 algorithm utilises Barnes Hut simulation, degree-dependent repulsive force, and local and global adaptive temperatures. 

<video preload="none" autoplay loop muted playsinline poster="LINK.jpg" width="100%" height="100%">
    <source src="[LINK.mp4](https://github.com/Issacgoh/Graphimate/blob/main/Emergence_first_blood_stem_cells.mp4)" type="video/mp4">
</video>

## Starting notes:
This package takes as input:
  - An anndata object
  - A categorical data variable containing labels
  - A 2D array containing some XY dimensionality-reduced coordinates (PCA, UMAP, VAE, etc...)
  - A sparse 2D matrix (CSR) containing cell-cell weighted connectivities (KNN, etc...)

## Installation:
To install directly from github run below in command line

pip install git+git@github.com:Issacgoh/Graphimate.git

To clone and install:

git clone https://github.com/Issacgoh/Graphimate.git

cd ./FA2nimate



## Usage notes:
