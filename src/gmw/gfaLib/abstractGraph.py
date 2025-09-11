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

from gfaLib.abstractions import Orientation
from gfaLib.gfaparser import GFAParser

class AbstractGraph:
    """
    Abstract class for GFA graphs.
    This class provides a common interface for GFA graphs, defining methods and attributes
    that must be implemented by concrete subclasses.
    """

    def __init__(self):
        self.segments: dict[str, dict] = {}
        self.lines: dict[tuple[str, str, int], dict] = {}
        self.paths: dict[str, dict] = {}
        self.headers: list[dict] = []
        
    """
    Abstract method for saving the graph.
    """    
    def save_graph(self, output_file, minimal=False, output_format=False):
        GFAParser.save_graph(graph=self, output_path=output_file, force_format=output_format, minimal_graph=minimal)

