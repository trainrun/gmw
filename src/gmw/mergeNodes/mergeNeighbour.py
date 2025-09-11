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
Neighbour Merger implementation for merging adjacent nodes in genome assembly graphs.
This module identifies and merges neighboring nodes that have simple connections,
simplifying the graph structure while preserving the biological information.
"""

import mergeNodes.basicFunction as bs
import shared

class NeighbourMerger:
    """
    Class for merging neighboring nodes into one node in graph.
    
    This class provides methods to identify and merge neighboring nodes
    in a genome assembly graph to simplify the graph structure while
    preserving the biological information.
    """
    def __init__(self, graph):
        self.graph = graph

    def merge_neibour(self):
        """
        Merge neighboring nodes in the graph.
        
        This method identifies and merges neighboring nodes that have simple
        connections (one outgoing edge from one node to another with one incoming edge).
        It removes redundant edges and updates the graph structure accordingly.
        """
        bs.remove_redundant_edges(self.graph)       
        nodes_to_remove = []
        for node in self.graph.nodes:
            while True:
                out_node_edges = self._classify_edges(node)
                if len(out_node_edges['plus_out']) == 1 and len(out_node_edges['minus_in']) == 0:
                    in_node_edges = self._classify_edges(out_node_edges['plus_out'][0][1])
                    if (out_node_edges['plus_out'][0][3]['label'][2] == '+' and len(in_node_edges['plus_in']) == 1 and len(in_node_edges['minus_out']) == 0
                    ) or (out_node_edges['plus_out'][0][3]['label'][2] == '-' and len(in_node_edges['minus_in']) == 1 and len(in_node_edges['plus_out']) == 0): 
                        nodes_to_remove.append(out_node_edges['plus_out'][0][1])
                        self._merge_node_property(out_node_edges['plus_out'][0])
                        self._move_edges(out_node_edges['plus_out'][0])
                        continue
                if len(out_node_edges['minus_out']) == 1 and len(out_node_edges['plus_in']) == 0:
                    in_node_edges = self._classify_edges(out_node_edges['minus_out'][0][1])
                    if (out_node_edges['minus_out'][0][3]['label'][2] == '+' and len(in_node_edges['plus_in']) == 1 and len(in_node_edges['minus_out']) == 0
                    ) or (out_node_edges['minus_out'][0][3]['label'][2] == '-' and len(in_node_edges['minus_in']) == 1 and len(in_node_edges['plus_out']) == 0): 
                        nodes_to_remove.append(out_node_edges['minus_out'][0][1])
                        self._merge_node_property(out_node_edges['minus_out'][0])
                        self._move_edges(out_node_edges['minus_out'][0])
                        continue
                break
        self.graph.remove_nodes_from(nodes_to_remove)
        
    def _merge_node_property(self, connect_edge):
        """
        Merge properties from the node to be merged into the node to be kept.
        
        This method combines sequence data and other properties (like depth, k-mer count,
        start/end positions, etc.) from both nodes. It handles different orientation
        combinations and uses CIGAR strings to properly merge sequences.
        
        Parameters:
            connect_edge: Tuple (u, v, key, data) representing the edge connecting the nodes
            
        Returns:
            None: The node properties are modified in-place
        """
        node_to_keep = self.graph.nodes[connect_edge[0]]
        node_to_merge = self.graph.nodes[connect_edge[1]]
        label = connect_edge[3]['label']
        cigar = connect_edge[3]['cigar']
        new_seq = ""
        
        # Merge sequences based on orientation combinations
        if label == '+/+':
            keep_seq = node_to_keep['seq']
            merge_seq = node_to_merge['seq']
            new_seq = bs.cigar_merge(keep_seq, merge_seq, cigar)
        elif label == '-/-':
            keep_seq = node_to_keep['seq']
            merge_seq = node_to_merge['seq']
            new_seq = bs.cigar_merge(merge_seq, keep_seq, cigar)
        elif label == '+/-':
            keep_seq = node_to_keep['seq']
            merge_seq = shared.reverse_complement(node_to_merge['seq'])
            new_seq = bs.cigar_merge(keep_seq, merge_seq, cigar)
        elif label == '-/+':
            keep_seq = node_to_keep['seq']
            merge_seq = shared.reverse_complement(node_to_merge['seq'])
            new_seq = bs.cigar_merge(merge_seq, keep_seq, cigar)
        
        # Combine depth (DP) as weighted average based on sequence lengths
        if 'DP' in node_to_keep and 'DP' in node_to_merge:
            node_to_keep['DP'] = round((node_to_keep['DP'] * node_to_keep['length'] + node_to_merge['DP'] * node_to_merge['length']
                                ) / (node_to_keep['length'] + node_to_merge['length']),2)
        # Combine k-mer count (KC)
        if 'KC' in node_to_keep and 'KC' in node_to_merge:
            node_to_keep['KC'] += node_to_merge['KC']
            
        # Update start position to the earlier of the two
        if node_to_keep['ST'] == 0 or node_to_keep['ST']>node_to_merge['ST']:
            node_to_keep['ST'] = node_to_merge['ST']
        # Update end position to the later of the two
        if node_to_keep['EN']<node_to_merge['EN']:
            node_to_keep['EN'] = node_to_merge['EN']
        # Update accession if needed
        if node_to_keep['AC'] == '-':
            node_to_keep['AC'] = node_to_merge['AC']
        # Update type if needed
        if node_to_keep['TP'] == '-':
            node_to_keep['TP'] == node_to_merge['TP']
        
        # Update orientation if needed
        if node_to_keep['OR'] == '?' and node_to_merge['OR'] != '?':
            if label[0] == label[2]:
                node_to_keep['OR'] = node_to_merge['OR']
            else:
                node_to_keep['OR'] = '+' if node_to_merge['OR'] == '-' else '-'
                
        # Update length and sequence
        node_to_keep['length'] = len(new_seq)
        node_to_keep['seq'] = new_seq
             
    def _move_edges(self, connect_edge):
        """
        Move edges from the node to be merged to the node to be kept.
        
        This method transfers all connections (both incoming and outgoing) from the node
        being merged to the node being kept, adjusting edge orientations as needed based
        on the relative orientations of the two nodes. It handles different orientation
        combinations and ensures proper edge key assignment to avoid conflicts.
        
        Parameters:
            connect_edge: Tuple (u, v, key, data) representing the edge connecting the nodes
            
        Returns:
            None: The graph is modified in-place
            
        Raises:
            ValueError: If the edge has an invalid label format
        """
        node_to_keep = connect_edge[0]
        node_to_merge = connect_edge[1]
        label = connect_edge[3]['label']
        
        # Handle same orientation merges
        if label == '+/+' or label == '-/-':
            # Transfer incoming edges
            for u, v, data in self.graph.in_edges(node_to_merge, data=True): 
                if u == node_to_keep:
                    continue
                new_key = 0
                if self.graph.has_edge(u, node_to_keep):
                    all_keys = self.graph[u][node_to_keep].keys()
                    new_key = max(all_keys)+1
                self.graph.add_edge(u, node_to_keep, key=new_key, **data)
            # Transfer outgoing edges
            for u, v, data in self.graph.out_edges(node_to_merge, data=True): 
                if v == node_to_keep:
                    continue
                new_key = 0
                if self.graph.has_edge(node_to_keep, v):
                    all_keys = self.graph[node_to_keep][v].keys()
                    new_key = max(all_keys)+1
                self.graph.add_edge(node_to_keep, v, key=new_key, **data) 
        # Handle opposite orientation merges
        elif label == '+/-' or label == '-/+':
            # Transfer incoming edges with orientation adjustment
            for u, v, data in self.graph.in_edges(node_to_merge, data=True): 
                if u == node_to_keep:
                    continue
                new_key = 0
                if self.graph.has_edge(u, node_to_keep):
                    all_keys = self.graph[u][node_to_keep].keys()
                    new_key = max(all_keys)+1
                # Flip the orientation of the target end
                if data['label'][2] == '+':
                    data['label'] = data['label'][:2] + '-'
                elif data['label'][2] == '-':
                    data['label'] = data['label'][:2] + '+'
                self.graph.add_edge(u, node_to_keep, key=new_key, **data)
            # Transfer outgoing edges with orientation adjustment
            for u, v, data in self.graph.out_edges(node_to_merge, data=True): 
                if v == node_to_keep:
                    continue
                new_key = 0
                if self.graph.has_edge(u, node_to_keep):
                    all_keys = self.graph[node_to_keep][v].keys()
                    new_key = max(all_keys)+1
                # Flip the orientation of the source end
                if data['label'][0] == '+':
                    data['label'] =  '-' + data['label'][1:]
                elif data['label'][0] == '-':
                    data['label'] =  '+' + data['label'][1:]
                    
                self.graph.add_edge(node_to_keep, v, key=new_key, **data) 
        else:
            raise ValueError(f"Edges have wrong label {label}")
        
        # Remove all edges connected to the merged node
        edges_to_remove = list(self.graph.in_edges(node_to_merge, keys=True)) + list(self.graph.out_edges(node_to_merge, keys=True))
        self.graph.remove_edges_from(edges_to_remove)
                                   
    def _classify_edges(self, node):
        """
        Classify edges connected to a node based on their orientation.
        
        This method categorizes edges into four groups based on their direction and
        orientation relative to the node:
        - plus_in: Incoming edges where the node's end is in the '+' orientation
        - plus_out: Outgoing edges where the node's start is in the '+' orientation
        - minus_in: Incoming edges where the node's end is in the '-' orientation
        - minus_out: Outgoing edges where the node's start is in the '-' orientation
        
        Parameters:
            node: The node ID whose edges need to be classified
            
        Returns:
            dict: Dictionary with four lists of classified edges
            
        Raises:
            ValueError: If an edge has an invalid label format
        """
        classified_edges = {"plus_in":[], "plus_out":[], "minus_in":[], "minus_out":[]}
        node_edges = self.graph.in_edges(node, keys=True, data=True)
        for u, v, key, data in self.graph.out_edges(node, keys=True, data=True):
            if u == v:
                continue
            if data['label'][0] =='+':
                classified_edges["plus_out"].append((u, v, key, data))
            elif data['label'][0] =='-':
                classified_edges["minus_out"].append((u, v, key, data))
            else:
                raise ValueError("Wrong label! The label shoud be +/+, +/-, -/- or -/+")
        
        # Classify incoming edges based on the orientation of the target end
        for u, v, key, data in self.graph.in_edges(node, keys=True, data=True):
            if u == v:
                continue
            if data['label'][2] =='+':
                classified_edges["plus_in"].append((u, v, key, data))
            elif data['label'][2] =='-':
                classified_edges["minus_in"].append((u, v, key, data))
            else:
                raise ValueError("Wrong label! The label shoud be +/+, +/-, -/- or -/+")
        return classified_edges