# Graphimate

A simple and intuitive package for creating fun animations of continuous graph layouts in single cell data. 

## About
Graphimate is a fun package for creating iterated animations and movies of continuous graph layouts using a defined graph layout algorithm. The default layout uses the FA2 algorithm which utilises Barnes Hut simulation, degree-dependent repulsive force, and local and global adaptive temperatures. With these iterative simulations, we can observe global structures in our data "unfolding" in a fun and intuitive manner! 

![Simulation of the emergence of the first blood stem cells in human life from specialized endothelium (note that the mp4 outputs are smoother)](https://github.com/Issacgoh/Graphimate/blob/main/resources/Emergence_first_blood_stem_cells.mp4.gif)


## Project team
Issac Goh, Newcastle University; Sanger institute (https://haniffalab.com/team/issac-goh.html)

### Contact
Issac Goh, (ig7@sanger.ac.uk)

## Built With
[Scanpy](https://scanpy.readthedocs.io/en/stable/)
[ForceAtlas2](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0098679)  
[Python-igraph](https://python.igraph.org/en/stable/)  
[MoviePy](https://zulko.github.io/moviepy/)

## Getting Started

### Prerequisites
This package takes as input:
  - An anndata object
  - A categorical data variable containing labels
  - A 2D array containing some XY dimensionality-reduced coordinates (PCA, UMAP, VAE, etc...)
  - A sparse 2D matrix (CSR) containing cell-cell weighted connectivities (KNN, etc...)

### Installation:
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
### Running Locally
This package was designed to run within a Jupyter notebook to utilise fully the interactive display interfaces. Functions can also be run locally via a python script. 
Please see the example notebook under "/example_notebooks/"

### Production
To deploy this package for large data submitted to schedulers on HPCs or VMs, please see the example given in "/example_notebooks/". 

### Graphing modules:
This package currently supports two custom iterative graph layout modules, 'FA' (forceatlas2) and 'FR' (Fruchterman-Reingold) modules.
More modules from the igraph collection will be added.

## Roadmap

- [ ] Initial Research
- [ ] Minimum viable product
- [x] Alpha Release <-- You are Here
- [ ] Feature-Complete Release

## Acknowledgements


