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
Basic functions for merging nodes in genome assembly graphs.
This module provides utility functions for sequence manipulation, similarity calculation,
and edge management used by the node merging algorithms.
"""

import re  
import config

# Dictionary mapping IUPAC nucleotide codes to their corresponding base sets
IUPAC_CODES = {
    'A': {'A'},
    'T': {'T'},
    'C': {'C'},
    'G': {'G'},
    'R': {'A', 'G'},
    'Y': {'C', 'T'},
    'S': {'C', 'G'},
    'W': {'A', 'T'},
    'K': {'G', 'T'},
    'M': {'A', 'C'},
    'B': {'C', 'G', 'T'},
    'D': {'A', 'G', 'T'},
    'H': {'A', 'C', 'T'},
    'V': {'A', 'C', 'G'},
    'N': {'A', 'T', 'C', 'G'}
}
# Reverse mapping from base sets to IUPAC codes for efficient lookup
REVERSE_IUPAC_CODES = {frozenset(v): k for k, v in IUPAC_CODES.items()}

def calculate_similarity(seq1, seq2):
    """
    Calculate the similarity percentage between two sequences of equal length.
    
    Parameters:
        seq1 (str): First nucleotide sequence
        seq2 (str): Second nucleotide sequence of the same length as seq1
        
    Returns:
        float: Similarity score as a percentage (0-100)
    """
    distance = sum(el1 != el2 for el1, el2 in zip(seq1, seq2))
    similarity_score = (1 - distance / len(seq1)) * 100
    return similarity_score

def get_degenerated_base(base1, base2):
    """
    Return the degenerate IUPAC nucleotide code that represents both input bases.
    
    Parameters:
        base1 (str): First nucleotide base (IUPAC code)
        base2 (str): Second nucleotide base (IUPAC code)
        
    Returns:
        str: IUPAC degenerate code representing both bases
    """
    set1 = IUPAC_CODES.get(base1, set())
    set2 = IUPAC_CODES.get(base2, set())
    # Merge the sets of possible bases
    merged_set = set1.union(set2)
    # Look up the corresponding IUPAC code for the merged set
    return REVERSE_IUPAC_CODES.get(frozenset(merged_set))
def create_consensus_sequence(seq1, seq2):
    """
    Generate a consensus sequence from two sequences using IUPAC ambiguity codes.
    
    For each position, if the bases are identical, that base is used.
    If they differ, the appropriate IUPAC ambiguity code is used.
    
    Parameters:
        seq1 (str): First nucleotide sequence
        seq2 (str): Second nucleotide sequence
        
    Returns:
        str: Consensus sequence with IUPAC ambiguity codes
        
    Raises:
        ValueError: If sequences have different lengths
    """
    if len(seq1) != len(seq2):
        raise ValueError("Sequences must be of the same length")
    consensus = []
    for base1, base2 in zip(seq1, seq2):
        if base1 == base2:
            consensus.append(base1)
        else:
            consensus.append(get_degenerated_base(base1, base2))    
    return ''.join(consensus)


def remove_redundant_edges(graph):  
    """
    Remove redundant edges from the graph that represent the same connection in opposite directions.
    
    This function identifies pairs of edges between the same nodes that represent the same
    biological connection but in opposite directions, and removes one of them to simplify the graph.
    
    Parameters:
        graph: NetworkX graph object representing the genome assembly graph
    """
    direction_dict = {"+/+": "-/-", "-/-": "+/+", "+/-": "+/-", "-/+": "-/+"}  
    edges_to_remove = set()  
    for u, v, key, data in graph.edges(keys=True, data=True):  
        if (u, v, key) not in edges_to_remove and graph.has_edge(v, u):  
            for _, _, rev_key, rev_data in graph.edges((v, u), keys=True, data=True):   
                if direction_dict[data['label']] == rev_data['label']:  
                    edges_to_remove.add((v, u, rev_key))  
                    break
    graph.remove_edges_from(edges_to_remove) 

def cigar_merge(seq1, seq2, cigar):
    """
    Merge two sequences based on a CIGAR string that describes their overlap.
    
    This function creates a merged sequence by joining seq1 and seq2 with a consensus
    sequence in the overlapping region as specified by the CIGAR string.
    
    Parameters:
        seq1 (str): First sequence (prefix sequence)
        seq2 (str): Second sequence (suffix sequence)
        cigar (str): CIGAR string describing the overlap (e.g., "10M" for 10 matching bases)
        
    Returns:
        str: Merged sequence with consensus in the overlapping region
        
    Raises:
        ValueError: If the CIGAR string is invalid or contains unsupported operations
    """
    consensus = []
    operations = re.findall(r'(\d+)([MIDNSHP])', cigar)
    if len(operations) > 1:
        raise ValueError("Wrong CIGAR sequence.")
    for count, op in operations:
        count = int(count)
        if op == 'M':
            # Take the beginning part of seq1, excluding the overlap
            consensus.append(seq1[:-count])
            overlap1 = seq1[-count:]
            overlap2 = seq2[:count] 
            con_seq = create_consensus_sequence(overlap1, overlap2)
            consensus.append(con_seq)
            consensus.append(seq2[count:])  # 追加 seq2 剩余部分
            # if overlap == seq2[:count]:
            #     consensus.append(seq2[count:])  # 追加 seq2 剩余部分
            #     # print("connected")
            #     break
            # if calculate_similarity(overlap1, overlap2) < config.similarity_score:
            #     print(f"{seq1}\n{seq2}")
            #     raise ValueError("The sequences to merge are not match")
            
        else:
            raise ValueError("Wrong CIGAR sequence.")
    return ''.join(consensus)

def cigar_judge_connect(seq1, seq2, cigar):
    """
    Create a degenerate clipped sequence for seq2 based on its overlap with seq1.
    
    This function examines the overlap between seq1 and seq2 as specified by the CIGAR string,
    creates a consensus sequence for the overlapping region, and returns this consensus
    followed by the non-overlapping part of seq2.
    
    Parameters:
        seq1 (str): First sequence (used to determine the overlap)
        seq2 (str): Second sequence (to be clipped and returned with consensus)
        cigar (str): CIGAR string describing the overlap (e.g., "10M" for 10 matching bases)
        
    Returns:
        str: Consensus sequence of the overlap + remaining part of seq2
        
    Raises:
        ValueError: If the CIGAR string is invalid or contains unsupported operations
    """
    operations = re.findall(r'(\d+)([MIDNSHP])', cigar)
    if len(operations) > 1:
        raise ValueError("Wrong CIGAR sequence.")
    for count, op in operations:
        count = int(count)
        if op == 'M':
            seq1_overlap = seq1[-count:]
            seq2_overlap = seq2[:count]
            # if calculate_similarity(seq1_overlap, seq2_overlap) < config.similarity_score:
            #     print(count,seq1_overlap, seq2_overlap)
            #     print(seq1, seq2)
            #     print(calculate_similarity(seq1_overlap, seq2_overlap))
            #     raise ValueError("Error in cigar_judge_connect")
            consensus = create_consensus_sequence(seq1_overlap, seq2_overlap)
            return consensus + seq2[count:]

    raise ValueError("Wrong CIGAR sequence.")
 