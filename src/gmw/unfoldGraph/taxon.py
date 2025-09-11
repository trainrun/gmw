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
Taxonomic parsing and analysis utilities.
This module provides functionality for parsing and navigating taxonomic information
from NCBI taxonomy database files, supporting taxonomic-based graph unfolding.
"""

class TaxonParser:
    def __init__(self, name_path, node_path):
        self.name_dict = {}
        self.node_dict = {}
        with open(name_path, "r") as file:
            for line in file:
                words = line.strip().split("\t|\t")
                name = words[1].lower().replace(" ", "_")
                if not name in self.name_dict:
                    self.name_dict[name] = words[0]
        with open(node_path, "r") as file:
            for line in file:
                words = line.strip().split("\t|\t")
                if not words[0] == words[1]:
                    self.node_dict[words[0]] = words[1]

    def _is_descendant(self, taxon, target):
        current = taxon
        while current in self.node_dict:
            if current == target:
                return True
            current = self.node_dict[current]
        return False

    def check_relationship(self, A, B):
        if self._is_descendant(A, B) or self._is_descendant(B, A):
            return True
        return False
    
    def name2taxid(self,name):
        lower_name = name.lower()
        if lower_name in self.name_dict:
            return self.name_dict[lower_name]
        raise ValueError(f"There is no scientific name \"{name}\"! Please correct the name or use tax_id instead.")