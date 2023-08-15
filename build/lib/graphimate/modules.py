from fa2 import ForceAtlas2
import numpy as np
import igraph as ig

class FruchtermanReingold:
    def __init__(self):
        # Any initialization or parameters can be set here if needed.
        pass

    def layout(self, G, pos, iterations=1):
        """
        Computes the Fruchterman-Reingold layout for the given graph 'G', 
        initialized with positions 'pos', for a specified number of 'iterations'.
        
        Parameters:
        - G: A scipy sparse matrix or similar adjacency matrix.
        - pos: Initial positions of the nodes.
        - iterations: Number of iterations to run the algorithm. Default is 1.
        
        Returns:
        - New positions of the nodes after specified iterations.
        """
        adjacency_list = (G.toarray() > 0).tolist()  # Convert the sparse matrix to dense and then to list
        graph = ig.Graph.Adjacency(adjacency_list)
        
        # Set the initial positions in the graph's layout attribute
        graph.vs["layout"] = pos
        
        # Now compute the layout
        layout = graph.layout_fruchterman_reingold(niter=iterations)
        
        return np.array(layout.coords)
    
class CustomForceAtlas2:
    def __init__(self, 
                 outboundAttractionDistribution=False,
                 linLogMode=False,
                 adjustSizes=False,
                 edgeWeightInfluence=1.0,
                 jitterTolerance=1.0,
                 barnesHutOptimize=True,
                 barnesHutTheta=1.2,
                 multiThreaded=False,
                 scalingRatio=2.0,
                 strongGravityMode=False,
                 gravity=1.0,
                 verbose=False):
        
        self.forceatlas2 = ForceAtlas2(
            outboundAttractionDistribution=outboundAttractionDistribution,
            linLogMode=linLogMode,
            adjustSizes=adjustSizes,
            edgeWeightInfluence=edgeWeightInfluence,
            jitterTolerance=jitterTolerance,
            barnesHutOptimize=barnesHutOptimize,
            barnesHutTheta=barnesHutTheta,
            multiThreaded=multiThreaded,
            scalingRatio=scalingRatio,
            strongGravityMode=strongGravityMode,
            gravity=gravity,
            verbose=verbose
        )

    def layout(self, G, pos, iterations=1):
        """
        Computes the layout using ForceAtlas2 for the given graph 'G', 
        initialized with positions 'pos', for a specified number of 'iterations'.
        
        Parameters:
        - G: A scipy sparse matrix or similar adjacency matrix.
        - pos: Initial positions of the nodes.
        - iterations: Number of iterations to run the algorithm. Default is 1.
        
        Returns:
        - New positions of the nodes after specified iterations.
        """
        return np.array(self.forceatlas2.forceatlas2(G=G, pos=pos, iterations=iterations))

