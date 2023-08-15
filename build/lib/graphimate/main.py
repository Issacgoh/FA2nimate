# FA2nimate 

##### Ver:: A0_V1
##### Author(s) : Issac Goh
##### Date : 230813;YYMMDD

### Author notes

### Features to add

import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from fa2 import ForceAtlas2
from moviepy.editor import VideoClip
from moviepy.video.io.bindings import mplfig_to_npimage
from datetime import datetime
from IPython.display import HTML
from base64 import b64encode
from scipy.sparse import csr_matrix
from moviepy.video.fx import all as vfx
from moviepy.editor import VideoClip, concatenate_videoclips
import random

from .modules import CustomForceAtlas2
from .modules import FruchtermanReingold


def register_data(adata, feat_use, n_iterations, skip_iterations = 1,explosion_duration = 0.15, knn_key = 'neighbors', use_initial="X_pca", var_length=7500, **kwargs ):
    """
    General data setup module.

    Parameters:

    Returns:

    """
    #unpack kwargs
    if kwargs:
        for key, value in kwargs.items():
            globals()[key] = value
        kwargs.update(locals())
        
    # Check cooridinate param exists
    if not use_initial in adata.obsm.keys():
        print('initialisation coordinates not detected in data, proceeding to compute PCA')
        print('We are assuming your data is pre-normalised')
        try:
            sc.pp.highly_variable_genes(adata, n_top_genes=var_length)
            keep = list(adata.var['dispersions_norm'].nlargest(7500).index)
            adata = adata[:,keep].to_memory()
        except:
            print("unable to compute dispersion, is your data norm?")
            print('proceeding to compute PCA on all genes as Highly var was not possible')
        sc.pp.pca(adata,50)
        use_initial = "X_pca"
        
    positions = adata.obsm[use_initial][:, :2].copy()
    
    # Now check KNN/custom connectivites graph
    # This line had to be added for as scanpy's recent update 1.9.3 introduces a bug for pp.neighbors if diffmap is available
    if 'X_diffmap' in adata.obsm.keys():
        adata.obsm['X_diffmap_']  = adata.obsm['X_diffmap']
        del adata.obsm['X_diffmap']
    if not knn_key in adata.uns.keys():
        print("Unable to locate your KNN graph, we will look for the key in obsp")
        try:
            snn = adata.obsp[knn_key]
            print("We found your graph, we now assume that this is the actual connectivity matrix")
        except:
            print('initialisation graph not detected in data, proceeding to compute KNN')
            sc.pp.neighbors(adata,key_added=knn_key,use_rep=use_initial)
            snn = adata.obsp[adata.uns[knn_key]['connectivities_key']]
    else:
        snn = adata.obsp[adata.uns[knn_key]['connectivities_key']]
    
    return adata, positions, snn

def setup_animation(adata,feat_use,use_initial='X_pca',knn_key='neighbors',edges=False,edge_subset=0.3,desired_total_duration=10,resolution=(1080,720),dpi=150,alg ='FA2',skip_iterations= 1,explosion_duration=0.15,var_length = 7500,out_path='./',markersize=5,**kwargs):
    """
    Creates an animation of a given data set using a specified layout algorithm.

    Parameters:
    -----------
    adata : object
        The input data object containing the data information.
    feat_use : str
        The feature to use from the `adata` object.
    use_initial : str, optional (default='X_pca')
        The initial position to use for the data points.
    knn_key : str, optional (default='neighbors')
        The key for k-nearest neighbors.
    edges : bool, optional (default=False)
        Whether to include edges in the animation.
    edge_subset : float, optional (default=0.3)
        The subset of edges to consider if `edges` is True.
    desired_total_duration : int, optional (default=10)
        The total duration of the animation in seconds.
    resolution : tuple, optional (default=(1080, 720))
        The resolution of the output animation.
    dpi : int, optional (default=150)
        Dots per inch for the output animation.
    alg : str, optional (default='FA2')
        The layout algorithm to use. Options are 'FA2' for ForceAtlas2 and 'FR' for FruchtermanReingold.
    skip_iterations : int, optional (default=1)
        The number of iterations to skip in the animation.
    explosion_duration : float, optional (default=0.15)
        The duration of the explosion effect in the animation.
    var_length : int, optional (default=7500)
        Variable length.
    out_path : str, optional (default='./')
        The path to save the animation output.
    **kwargs : dict
        Additional keyword arguments for advanced customization.

    Returns:
    --------
    fpath : str
        The file path of the generated animation.
    
    Notes:
    ------
    This function creates an animated visualization of the data points using the ForceAtlas2 or 
    FruchtermanReingold layout algorithms. The animation starts with an explosion effect where 
    data points move away from the center. The main part of the animation visualizes the data 
    points' movement as they are positioned by the chosen layout algorithm. The resulting animation 
    is saved as an MP4 file.
    """
    
    if kwargs:
        for key, value in kwargs.items():
            globals()[key] = value
        kwargs.update(locals())
        now = datetime.now()
        date_time = now.strftime("%d%m%Y_%H%M")
        
    # Check KNN/custom connectivites graph
    if not knn_key in adata.uns.keys():
        print("Unable to locate your KNN graph in uns keys, we will look for the key in obsp")
        try:
            snn = adata.obsp[knn_key]
            print("We found your graph, we now assume that this is the actual connectivity matrix")
        except:
            print('initialisation graph not detected in data, did you register this dataset?')
    else:
        snn = adata.obsp[adata.uns[knn_key]['connectivities_key']]
        
    # Nested make_frame function
    def make_frame(t):
        nonlocal positions, iteration_counter, markersize
        #print(f"Time: {t}, Iteration Counter: {iteration_counter}")
        # Perform the explosion sequence
        if t < explosion_duration:
            previous_frame = explosion_frame(t)
            return previous_frame
        else:# Adjust time for FA2 iterations
            t -= explosion_duration
            #print(f"Time: {t}, FA2 Counter: {iteration_counter}")
            # Perform the ForceAtlas2 computation for every frame, but only plot and return frames
            # every skip_iterations + 1 frames
            if alg == 'FR':
                print('Using the FR approach')
                layout = FruchtermanReingold()
            else:
                layout = CustomForceAtlas2()
            if skip_iterations > 1:
                new_positions = layout.layout(G=snn, pos=positions, iterations=skip_iterations)
            else:
                new_positions = layout.layout(G=snn, pos=positions, iterations=1)
            positions = new_positions[:]

            #if iteration_counter (skip_iterations + 1) == 0:
            ax.clear()
            ax.grid(False)
            ax.set_facecolor('black')
            ax.axis('off')
            ax.set_xlim(positions[:, 0].min() * 1.3, positions[:, 0].max() * 1.3)
            ax.set_ylim(positions[:, 1].min() * 1.3, positions[:, 1].max() * 1.3)
            ax.scatter(*np.array(positions).T, c=cell_type_colors, alpha=1, s=1)

            if edges:
                for i, j in edge_positions:
                    ax.plot(*zip(positions[i], positions[j]), color='white', linewidth=0.2)

            legend_elements = [Line2D([0], [0], marker='o', color='w', label=cell_type,
                                      markerfacecolor=color, markersize=4)
                               for cell_type, color in zip(adata.obs['celltype'].cat.categories, adata.uns['celltype_colors'])]
            ax.legend(handles=legend_elements, loc='upper right', bbox_to_anchor=(1, 1.1), frameon=False, fontsize=5,
                      labelcolor='white',ncol=2)
            previous_frame = mplfig_to_npimage(fig)
            iteration_counter[0] += 1
            return previous_frame

    # Nested explosion function
    def explosion_frame(t):
        nonlocal positions, iteration_counter, markersize
        ax.clear()
        ax.grid(False)
        ax.set_facecolor('black')
        ax.axis('off')
        ax.set_xlim(positions[:, 0].min() * 1.3, positions[:, 0].max() * 1.3)
        ax.set_ylim(positions[:, 1].min() * 1.3, positions[:, 1].max() * 1.3)
        # Linear interpolation between the mean of the positions and the initial positions
        ratio = t / explosion_duration
        explode_positions = mean_position + ratio * (positions - mean_position)
        ax.scatter(*np.array(explode_positions).T, c=cell_type_colors, alpha=1, s=1)
        if edges:
            for i, j in edge_positions:
                ax.plot(*zip(explode_positions[i], explode_positions[j]), color='white', linewidth=0.1)

        legend_elements = [Line2D([0], [0], marker='o', color='w', label=cell_type,
                                  markerfacecolor=color, markersize=4)
                           for cell_type, color in zip(adata.obs['celltype'].cat.categories, adata.uns['celltype_colors'])]
        ax.legend(handles=legend_elements, loc='upper right', bbox_to_anchor=(1, 1.1), frameon=False, fontsize=5,
                  labelcolor='white',ncol=2)
        return mplfig_to_npimage(fig)

    # Center the initial positions around (0,0)
    positions = adata.obsm[use_initial][:, :2].copy() - np.mean(adata.obsm[use_initial][:, :2], axis=0)

    # Run a single iteration to find a first position set
    if alg == 'FR':
        print('Using the FR approach')
        layout = FruchtermanReingold()
    else:
        layout = CustomForceAtlas2()
    positions = layout.layout(G=snn, pos=positions, iterations=1)
    mean_position = np.array([0, 0])
    # create celltype assignment
    adata.obs['celltype'] = adata.obs[feat_use].astype('category')
    try:
        colors = list(adata.uns[feat_use+'_colors'])
    except:
        print("You do not provide a color key in .uns, we will randomly generate one")
        # create color dict
        number_of_colors = len(adata.obs['celltype'].cat.categories)
        colors = ["#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)])
                     for i in range(number_of_colors)]
        adata.uns['celltype_colors'] = colors
    col_dic = dict(zip(list(adata.obs['celltype'].cat.categories),list(colors)))
    cell_type_colors = adata.obs['celltype'].map(col_dic)
    # Calculate edge positions once outside of make_frame
    if edges:
        print("Warning, you are plotting with edges, this can be very slow, use plot_edges=False if it takes too long")
        snn_coo = snn.tocoo()
        edge_weights = snn_coo.data
        sorted_indices = np.argsort(edge_weights)
        top_10_percent = sorted_indices[-int(len(sorted_indices) * edge_subset):]
        edge_positions = []
        for idx in top_10_percent:
            i = snn_coo.row[idx]
            j = snn_coo.col[idx]
            edge_positions.append((i, j))
    fig, ax = plt.subplots(dpi=dpi)
    fig.patch.set_facecolor('black')  # Set the entire figure background to black
    ax.set_facecolor('black')
    
    # Store the previous frame to return during skipped iterations
    previous_frame = None
    # Calculate the total number of frames
    total_frames = n_iterations# // (skip_iterations + 1)
    # Define the frames per second
    fps = (n_iterations/(desired_total_duration-explosion_duration))#*1.5 #*3  # or any other value that gives a smooth animation
    # Calculate the duration of the FA2 phase
    fa2_duration = (desired_total_duration-explosion_duration)#*5
    # Set the speed factor to fit the total duration of the FA2 phase into the desired total duration
    fa2_speed_factor = fa2_duration / (desired_total_duration - explosion_duration)
    # Counter for the number of iterations
    iteration_counter = [0] #mutable list sahould be able to pass to nested functions
    
    # Create the animation from frames
    animation = VideoClip(make_frame, duration=fa2_duration + explosion_duration)
    animation.resize(width=resolution[0], height=resolution[1])
    fpath = out_path+'/'+date_time + "fa2_animation.mp4"
    animation.write_videofile(fpath, fps=int(fps))
    
    return fpath


from moviepy.editor import VideoFileClip

def video_to_gif(input_video_path, output_gif_path, fps=10):
    """
    Convert a video into gif.
    
    Parameters:
    - input_video_path: str, path to the input video file.
    - output_gif_path: str, path to save the output gif.
    - fps: int, frames per second for the gif.
    """
    with VideoFileClip(input_video_path) as clip:
        clip.write_gif(output_gif_path, fps=fps)