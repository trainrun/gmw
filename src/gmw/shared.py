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
GMW Shared Utilities Module

This module provides shared utility functions used across the GMW workflow.
It includes functions for DNA sequence manipulation, file format conversion,
external tool execution, and graph filtering operations.

The utilities support common operations like:
- DNA sequence reverse complement calculation
- Graph to FASTA conversion
- External tool execution (BLAST, Kraken)
- Graph component filtering and cleanup
"""

import os
import subprocess
import networkx as nx
import logging

import config

# Initialize logger for shared utilities
logger = logging.getLogger("gmw")

# def judge_file_exist(file_path):
#     dir_path = os.path.dirname(file_path)
#     if not os.path.exists(dir_path):
#         raise ValueError(dir_path + " not exist!")    
#     return os.path.exists(file_path)

# def judge_dir_exist(dir_path):
#     return os.path.exists(dir_path)

def reverse_complement(seq):
    """
    Calculate the reverse complement of a DNA sequence.
    
    This function supports both uppercase and lowercase letters, as well as
    degenerate nucleotide codes (IUPAC ambiguity codes). It reverses the
    sequence and replaces each base with its complement.
    
    Parameters
    ----------
    seq : str
        Input DNA sequence string containing valid nucleotide characters
        
    Returns
    -------
    str
        Reverse complement of the input sequence
        
    Raises
    ------
    ValueError
        If the sequence contains invalid characters not in the complement dictionary
        
    Examples
    --------
    >>> reverse_complement("ATCG")
    'CGAT'
    >>> reverse_complement("atcg")
    'cgat'
    >>> reverse_complement("ATCGN")
    'NCGAT'
    """
    # Mapping of nucleotides to their complements (supports IUPAC codes)
    complement_dict = {
        'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
        'a': 't', 't': 'a', 'c': 'g', 'g': 'c',
        'M': 'K', 'K': 'M', 'R': 'Y', 'Y': 'R',
        'W': 'W', 'S': 'S', 'B': 'V', 'V': 'B',
        'H': 'D', 'D': 'H', 'N': 'N',
        'm': 'k', 'k': 'm', 'r': 'y', 'y': 'r',
        'w': 'w', 's': 's', 'b': 'v', 'v': 'b',
        'h': 'd', 'd': 'h', 'n': 'n'
    }

    # Validate input sequence
    valid_bases = set(complement_dict.keys())
    if not all(base in valid_bases for base in seq):
        raise ValueError("Invalid DNA sequence")
    
    # Generate reverse complement
    return ''.join([complement_dict[base] for base in reversed(seq)])

def graph2fasta(graph, file_path, orientation=False):
    """
    Convert a NetworkX graph to FASTA format file.
    
    This function extracts sequence information from graph nodes and writes
    them to a FASTA file. Optionally handles sequence orientation by applying
    reverse complement transformation for nodes with negative orientation.
    
    Parameters
    ----------
    graph : networkx.Graph
        Input graph containing nodes with sequence information
    file_path : str
        Output file path for the FASTA file
    orientation : bool, optional
        If True, apply reverse complement for nodes with negative orientation
        (default: False)
        
    Notes
    -----
    Each node must have 'seq' attribute containing the DNA sequence.
    If orientation=True, nodes must also have 'OR' attribute indicating
    orientation ('+' or '-').
    """
    with open(file_path, 'w') as f:
        for node in graph.nodes:
            f.write(">" + str(node) + "\n")
            if orientation and graph.nodes[node]['OR'] == '-':
                f.write(reverse_complement(graph.nodes[node]['seq']) + "\n")
            else:
                f.write(graph.nodes[node]['seq'] + "\n")

def run_blast(fasta, db, out, threads=config.max_threads):
    commands = [config.blast_path, '-query', fasta, '-db', db, '-num_threads', str(threads), '-outfmt', '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs', '-out',
            out, '-max_target_seqs', '5']
    logger.info(' '.join(commands))
    result = subprocess.check_call(commands)
    
def run_kraken(fasta, db, output, threads=config.max_threads):
    commands = [config.kraken_path, "--db", db, "--output", output, "--threads", str(threads), "--confidence", str(config.confidence), fasta]
    logger.info(' '.join(commands))
    result = subprocess.check_call(commands)
    
def remove_unknown_components(graph):
    """
    Remove weakly connected components that contain no target or annotated nodes.
    
    This function identifies components where all nodes lack target classification
    or annotation, then removes these entire components from the graph.
    
    Parameters
    ----------
    graph : networkx.DiGraph
        Input directed graph to filter
        
    Notes
    -----
    Modifies the input graph in-place by removing unknown components.
    A component is considered "unknown" if none of its nodes have:
    - TP (Type) attribute set to "target"
    - AC (Annotation) attribute different from '-'
    """
    nodes_for_unknown_components = set()
    for component in nx.weakly_connected_components(graph):
        target_component = False
        for node in component:
            if graph.nodes[node]['TP'] == "target" or graph.nodes[node]['AC'] != '-':
                target_component = True
                break
        if not target_component:
            nodes_for_unknown_components.update(component)
    graph.remove_nodes_from(nodes_for_unknown_components)

def remove_short_isolated_nodes(graph, length):
    """
    Remove isolated nodes (degree 0) with sequences shorter than threshold.
    
    This function identifies nodes that have no connections to other nodes
    and whose sequences are below the specified length threshold, then removes them.
    
    Parameters
    ----------
    graph : networkx.Graph
        Input graph to filter
    length : int
        Minimum sequence length threshold
        
    Notes
    -----
    Modifies the input graph in-place by removing short isolated nodes.
    Only affects nodes with degree 0 (no edges to other nodes).
    """
    # Identify short isolated nodes
    isolated_nodes = [node for node in graph.nodes() if (graph.degree(node) == 0 and graph.nodes[node]['length'] < length)]    
    graph.remove_nodes_from(isolated_nodes)
    
def remove_unknown_nodes(graph):
    """
    Remove individual nodes that lack both type and annotation classification.
    
    This function removes nodes that have neither type classification (TP)
    nor annotation (AC), indicating they are completely uncharacterized.
    
    Parameters
    ----------
    graph : networkx.Graph
        Input graph to filter
        
    Notes
    -----
    Modifies the input graph in-place by removing unknown nodes.
    A node is considered "unknown" if both:
    - TP (Type) attribute is '-'
    - AC (Annotation) attribute is '-'
    """
    # Identify nodes without type or annotation
    unknown_nodes = [node for node in graph.nodes() if (graph.nodes[node]['TP'] == '-' and graph.nodes[node]['AC'] == '-')]    
    graph.remove_nodes_from(unknown_nodes)

def remove_contaminated_nodes(graph):
    """
    Remove nodes classified as contamination or inferred contamination.
    
    This function removes nodes that have been classified as contamination
    through either direct classification or inference methods.
    
    Parameters
    ----------
    graph : networkx.Graph
        Input graph to filter
        
    Notes
    -----
    Modifies the input graph in-place by removing contaminated nodes.
    Removes nodes with TP (Type) attribute set to:
    - 'contaminate': directly classified contamination
    - 'contaminate_infer': inferred contamination
    """
    # Remove directly classified contaminated nodes
    contaminate_nodes = [node for node in graph.nodes() if (graph.nodes[node]['TP'] == 'contaminate')]
    graph.remove_nodes_from(contaminate_nodes)    
    contaminate_nodes = [node for node in graph.nodes() if (graph.nodes[node]['TP'] == 'contaminate_infer')]
    graph.remove_nodes_from(contaminate_nodes)    
