
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
Graph implementation for creating GFA graphs from NetworkX structures.
This module provides functionality to convert NetworkX graph structures into
GFA format by mapping nodes and edges to segments and lines in the GFA graph.
It inherits from AbstractGraph to maintain compatibility with other GFA operations.
"""

from gfaLib.abstractions import GFALine, Orientation, reverse
from gfaLib.gfaparser import GFAParser
from gfaLib.abstractGraph import AbstractGraph

class GraphFromNetwork(AbstractGraph):
    """
    A class for creating a GFA graph from a NetworkX graph structure.
    
    This class converts an existing NetworkX graph into the GFA (Graphical Fragment Assembly)
    format by mapping nodes and edges from the network to segments and lines in the GFA graph.
    It inherits from AbstractGraph to maintain compatibility with other GFA operations.
    """
    def __init__(self, network):
        super().__init__()
        self.metadata: dict = {'version': 'unknown', 'next_node_name':  'unknown', 'with_sequence': True}
        
        # Convert network nodes to GFA segments
        for node in network.nodes:
            self.segments[node] = {**network.nodes[node]}
        
        # Convert network edges to GFA lines with orientation information
        for u, v, data in network.edges(data=True):
            #"print hehe"
            self.lines[(u,v)] = {"orientation":{(Orientation(data['label'][0]), Orientation(data['label'][2]))}, "ARG5": data['cigar']}
        