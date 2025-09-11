
""" 
GMW: Genomic Microbe-Wise - hybrid assembly and contamination removal tool 

Copyright (C) 2025 Wenbing Chen 
www.github.com/trainrun/gmw 

License: 
This program is free software: you can redistribute it and/or modify 
it under the terms of the GNU General Public License as published by 
the Free Software Foundation, either version 3 of the License, or 
(at your option) any later version. 

This program is distributed in the hope that it will be useful, 
but WITHOUT ANY WARRANTY; without even the implied warranty of 
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.   
See the GNU General Public License for more details. 

You should have received a copy of the GNU General Public License 
along with this program. If not, see <https://www.gnu.org/licenses/>. 
"""

"""
Function Description: 
Depth-based graph unfolding implementation.
This module provides functionality to unfold genome assembly graphs based on
sequencing depth discrepancies between connected nodes.
"""

from unfoldGraph.abstrctUnfold import AbstrctUnfolder

class DepthUnfolder(AbstrctUnfolder):
    """
    Graph unfolder based on sequencing depth discrepancy.
    
    This class implements graph unfolding by removing edges between nodes
    that have significant differences in sequencing depth, which often
    indicates separate genomic regions that should not be connected.
    """
    # def __init__(self, graph, multiple=10):
    #     self.graph = graph
    #     self.multiple = multiple
    #     self.dir_path = "/home/cwb/chen/20241006HIV_test/FLU_1/filter_assembly/mmm"
    #     self.force = True
                
    def unfold_graph(self):
        visual_in = self.out_path + "/" + self.prefix + "_depth_input.html"
        visual_out = self.out_path + "/" + self.prefix + "_depth_output.html"
        gfa_out = self.out_path + "/" + self.prefix + "_depth_output.gfa"
        self._visual_graph(visual_in)        
        
        edges_remove = []
        for node1, node2, key in self.graph.edges(keys=True):
            depth1 = self.graph.nodes[node1]['DP']
            depth2 = self.graph.nodes[node2]['DP']
            if depth1/depth2 > self.depth_discrepancy or depth2/depth1 > self.depth_discrepancy:
                edges_remove.append((node1, node2, key))
        self.graph.remove_edges_from(edges_remove)

        if self.disable_ref_unfold and self.disable_taxon_unfold:
            self.keep_unknown_components = True
            self.remove_unknown_nodes = False
        self._remove_components()
        self._merge_nodes()
        self._output_gfa(gfa_out)
        self._visual_graph(visual_out)  