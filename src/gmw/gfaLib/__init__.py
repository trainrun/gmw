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
GFA Library - Abstraction layer for GFA (Graphical Fragment Assembly) format
This package provides classes and utilities for parsing, manipulating, and
working with GFA files, which represent genome assembly graphs.
"""
from .abstractions import GFALine, GFAFormat, Orientation
from .gfaparser import GFAParser
from .graphFromFile import GraphFromFile
from .graphFromNetwork import GraphFromNetwork
from .nx import GFANetwork
