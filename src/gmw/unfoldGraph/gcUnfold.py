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
GC content-based graph unfolding implementation.
This module provides functionality to unfold genome assembly graphs based on
GC content discrepancies between connected nodes, which often indicates
different genomic origins.
"""

from unfoldGraph.abstrctUnfold import AbstrctUnfolder
class GCUnfolder(AbstrctUnfolder):
    # def __init__(self, graph, gc_discrepancy=0.1):
    #     self.graph = graph
    #     self.gc_discrepancy = gc_discrepancy
                
    def unfold_graph(self):
        visual_in = self.out_path + "/" + self.prefix + "_gc_input.html"
        visual_out = self.out_path + "/" + self.prefix + "_gc_output.html"
        gfa_out = self.out_path + "/" + self.prefix + "_gc_output.gfa"
        self._visual_graph(visual_in) 
                
        edges_remove = []
        for node1, node2, key in self.graph.edges(keys=True):
            gc1 = self.gc_content(self.graph.nodes[node1]['seq'])
            gc2 = self.gc_content(self.graph.nodes[node2]['seq'])
            if gc1 - gc2 > self.gc_discrepancy or gc2 - gc1 > self.gc_discrepancy:
                edges_remove.append((node1, node2, key))
        self.graph.remove_edges_from(edges_remove)
        if self.disable_ref_unfold and self.disable_taxon_unfold:
            self.keep_unknown_components = True
            self.remove_unknown_nodes = False
        self._remove_components()        
        self._merge_nodes()
        self._output_gfa(gfa_out)
        self._visual_graph(visual_out)
        
    def gc_content(self, seq):
        seq = seq.upper()
        gc_count = seq.count('G') + seq.count('C')
        total_count = len(seq)
        gc_content = gc_count / total_count
        return gc_content