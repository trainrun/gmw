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

from unfoldGraph.abstrctUnfold import AbstrctUnfolder

class Polisher(AbstrctUnfolder):
    """
    Graph polishing class for final cleanup and optimization of genome assembly graphs.
    
    This class performs the final polishing steps after all unfolding operations have been applied,
    including removing nodes without reference accessions, cleaning up components, merging nodes,
    and outputting the final GFA file.
    
    Inherits from AbstrctUnfolder to access common graph manipulation methods.
    """
    
    def polish(self):
        if not self.disable_ref_unfold:
            rmove_list = []
            for node in self.graph.nodes():
                if self.graph.nodes[node]['AC'] == '-':
                    rmove_list.append(node)
            self.graph.remove_nodes_from(rmove_list)
        
        # Remove unknown components if both reference and taxon unfolding are enabled
        if not (self.disable_ref_unfold or self.disable_taxon_unfold):
            self._remove_components()
        
        # Iteratively merge nodes until no more changes occur
        nodes_num = 0
        while nodes_num != self.graph.number_of_nodes():
            nodes_num = self.graph.number_of_nodes()
            self._merge_nodes()

        gfa_out = self.out_path + "/../" +  self.prefix + "_after_unfold.gfa"
        self._output_gfa(gfa_out)