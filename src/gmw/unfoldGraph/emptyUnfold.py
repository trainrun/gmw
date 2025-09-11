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
Empty unfolding implementation that applies basic graph transformations.
This module provides a minimal unfolding implementation that only performs
basic node merging and splitting operations without specific criteria.
"""

from unfoldGraph.abstrctUnfold import AbstrctUnfolder
from mergeNodes.mergeBrother import BrotherMerger
from mergeNodes.splitParent import ParentSpliter
class EmptyUnfolder(AbstrctUnfolder):
    def unfold_graph(self):
        if self.merge_brother:
            merger = BrotherMerger(self.graph)
            merger.merge_brother()
            self._merge_nodes()
        if self.split_parent:
            spliter = ParentSpliter(self.graph)
            spliter.split_parent()
            self._merge_nodes()