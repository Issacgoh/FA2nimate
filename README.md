# Graphimate
This is a fun package for creating iterated animations and movies of continuous graph layouts using the FA2 algorithm. The FA2 algorithm utilises Barnes Hut simulation, degree-dependent repulsive force, and local and global adaptive temperatures. 

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
