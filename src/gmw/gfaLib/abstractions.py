
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

from enum import Enum

"""
Functions and enums used for GFA file manipulation.
This file is modified from GFAGraphs (https://github.com/dubssieg/gfagraphs/tree/gfagraphs).
"""

class Orientation(Enum):
    """
    Enumeration representing the orientation of segments in a GFA graph.
    Attributes:
        FORWARD ('+'): Forward orientation
        REVERSE ('-'): Reverse orientation
        ANY ('?'): Unspecified orientation
        BOTH ('='): Both orientations
    """
    FORWARD = '+'
    REVERSE = '-'
    ANY = '?'
    BOTH = '='


def reverse(orientation: Orientation) -> Orientation:
    """
    Get the reverse of a given orientation.
    
    Args:
        orientation (Orientation): The orientation to reverse
        
    Returns:
        Orientation: The reversed orientation
    """
    return {
        Orientation.FORWARD: Orientation.REVERSE,
        Orientation.REVERSE: Orientation.FORWARD,
        Orientation.ANY: Orientation.ANY,
        Orientation.BOTH: Orientation.BOTH,
    }[orientation]


class GFAFormat(Enum):
    """
    Enumeration representing different GFA file format versions.
    
    Attributes:
        RGFA ('rGFA'): rGFA format
        GFA1 ('GFA1'): GFA version 1
        GFA1_1 ('GFA1.1'): GFA version 1.1
        GFA1_2 ('GFA1.2'): GFA version 1.2
        GFA2 ('GFA2'): GFA version 2
        ANY ('unknown'): Unknown or unspecified format
    """
    RGFA = 'rGFA'
    GFA1 = 'GFA1'
    GFA1_1 = 'GFA1.1'
    GFA1_2 = 'GFA1.2'
    GFA2 = 'GFA2'
    ANY = 'unknown'


class GFALine(Enum):
    """
    Enumeration representing different line types in a GFA file.
    
    Attributes:
        SEGMENT ('S'): Segment line defining nodes
        LINK ('L'): Link line defining edges
        WALK ('W'): Walk line defining paths
        PATH ('P'): Path line defining paths (alternative to walk)
        HEADER ('H'): Header line containing metadata
        COMMENT ('#'): Comment line
        ANY ('?'): Any or unspecified line type
    """
    SEGMENT = 'S'
    LINK = 'L'
    WALK = 'W'
    PATH = 'P'
    HEADER = 'H'
    COMMENT = '#'
    ANY = '?'
