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
NetworkX representation of GFA graphs.
This module provides functionality to convert GFA graph structures into
NetworkX MultiDiGraph objects, enabling graph analysis and manipulation
using the NetworkX library's algorithms and utilities.
"""

from networkx import MultiDiGraph
from gfaLib.graphFromFile import GraphFromFile

class GFANetwork:
    """
    NetworkX representation of a GFA graph.
    This class provides methods to convert a GFA graph represented by a GraphFromFile object
    into a NetworkX MultiDiGraph, which is a directed graph with multiple edges between nodes.
    This file is modified from GFAGraphs (https://github.com/dubssieg/gfagraphs/tree/gfagraphs).
    """
    @staticmethod
    def compute_backbone(
        graph: GraphFromFile
    ) -> MultiDiGraph:
        backbone: MultiDiGraph = MultiDiGraph()

        for node_name, node_datas in graph.segments.items():
            backbone.add_nodes_from([( node_name, node_datas)])
        for (start, end, key), edge_data in graph.lines.items():
            backbone.add_edge(
                start,
                end,
                key=key,
                label=' | '.join(
                    [f'{x.value}/{y.value}' for (x, y) in edge_data["orientation"]]),
                cigar=edge_data['ARG5']
            )
        return backbone


