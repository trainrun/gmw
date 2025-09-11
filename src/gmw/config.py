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

"""GMW Configuration Module

This module contains all configuration parameters and default values for the GMW
(Genome Metagenome Workflow) tool. It defines settings for various components
including output paths, threading, graph unfolding algorithms, visualization,
external tool integration, and algorithm-specific parameters.

The configuration is organized into logical sections for easy maintenance
and customization of the workflow behavior.
"""

import os

# Application version
VERSION = "1.0"

"""Output config"""
prefix = "gmw"
output_path = "./gmw_output"


"""threads config"""
max_threads = os.cpu_count()
#max_threads = 12

#overlap_size = 77

"""unfold position config"""
position_distance = 150
gc_discrepancy = 0.2
depth_discrepancy = 20

"""remove short reads offset"""
short_offset = 400

"""visualize config"""
edge_color="lavender"
node_shape="dot"
length_range_mix=50
length_range_max=500
width_range_mix=5
width_range_max=20
color_dict = {"unknown": "grey", "target": "pink", "contaminate": "skyblue"}

"""configs for kraken2"""
confidence=0.8
kraken_path = "kraken2"

"""blastn_config"""
blast_path = "blastn"
query_cover_offset=85
match_length=70
identity_discrepency=1.0

"""config for merge brother contigs"""
# length_discrepancy=0
# blastout_discrepancy=1
similarity_score = 90


"""config for split parent"""
split_similarity_score = 80
split_depth_multi = 5 

"""BGLL parameters"""
resolution = 2.0
partition_offset  = 10