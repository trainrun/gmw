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
Brother Merger implementation for merging similar nodes in genome assembly graphs.
This module identifies and merges nodes that have similar sequences and share
connections to the same neighboring nodes, simplifying the graph structure
while preserving the biological information.
"""

import mergeNodes.basicFunction as bs
import shared
import config

class BrotherMerger:
    """
    Class for merging similar nodes (brothers) in a genome assembly graph.
    
    This class identifies and merges nodes that have similar sequences and share
    connections to the same neighboring nodes, simplifying the graph structure
    while preserving the biological information.
    """
    def __init__(self, graph):
        self.graph = graph

    def merge_brother(self):
        """
        Identify and merge brother nodes in the graph.
        
        This method finds pairs of nodes that share connections to the same node and have
        similar sequences (brothers). It merges these nodes to simplify the graph structure
        while preserving the biological information.
        
        The method processes both outgoing and incoming edges, looking for potential
        brother nodes in both directions.
        """
        nodes_to_remove = []
        for node in self.graph.nodes:
            out_node_edges = self._classify_edges(node)
            # Process outgoing edges from the current node
            if len(out_node_edges['plus_out']) + len(out_node_edges['minus_in']) > 1:
                merge_edges = [*out_node_edges['plus_out'], *out_node_edges['minus_in']]
                for i in range(len(merge_edges)):
                    node_tuple1 = self._judge_node_property(node, merge_edges[i])
                    if node_tuple1[0] in nodes_to_remove:
                        continue 
                    for j in range(i+1, len(merge_edges)):
                        node_tuple2 = self._judge_node_property(node, merge_edges[j])
                        if node_tuple2[0] in nodes_to_remove:
                            continue 
                        if node_tuple1[2] != node_tuple2[2]:
                            raise ValueError("Nodes direction error!")
                        is_brother = self._judge_brother(node_tuple1, node_tuple2)
                        if is_brother:
                            self._move_edges(node, node_tuple1, node_tuple2)
                            self._merge_node_property(node_tuple1, node_tuple2)
                            self._degenerated_neighbour(node_tuple1[0])
                            nodes_to_remove.append(node_tuple2[0])
            # Process incoming edges to the current node
            if len(out_node_edges['plus_in']) + len(out_node_edges['minus_out']) > 1:
                merge_edges = [*out_node_edges['plus_in'], *out_node_edges['minus_out']]
                for i in range(len(merge_edges)):
                    node_tuple1 = self._judge_node_property(node, merge_edges[i])
                    if node_tuple1[0] in nodes_to_remove:
                        continue 
                    for j in range(i+1, len(merge_edges)):
                        node_tuple2 = self._judge_node_property(node, merge_edges[j])
                        if node_tuple2[0] in nodes_to_remove:
                            continue 
                        if node_tuple1[2] != node_tuple2[2]:
                            raise ValueError("Nodes direction error!")
                        is_brother = self._judge_brother(node_tuple1, node_tuple2)
                        if is_brother:
                            self._move_edges(node, node_tuple1, node_tuple2)
                            self._merge_node_property(node_tuple1, node_tuple2)
                            self._degenerated_neighbour(node_tuple1[0])
                            nodes_to_remove.append(node_tuple2[0])
        self.graph.remove_nodes_from(nodes_to_remove)                    
                    
    def _judge_brother(self, node_tuple1, node_tuple2):
        """
        Determine if two nodes are brothers based on sequence similarity.
        
        Two nodes are considered brothers if they have the same length and their
        sequences are similar enough (based on the similarity threshold in config).
        
        Parameters:
            node_tuple1: Tuple (node_id, from_direction, to_direction) for first node
            node_tuple2: Tuple (node_id, from_direction, to_direction) for second node
            
        Returns:
            bool: True if nodes are brothers, False otherwise
        """
        node1 = self.graph.nodes[node_tuple1[0]]
        node2 = self.graph.nodes[node_tuple2[0]]
        if node1['length'] == node2['length']:
            seq2 = node2['seq'] if node_tuple1[1]==node_tuple2[1] else shared.reverse_complement(node2['seq'])
            similarity_score = bs.calculate_similarity(node1['seq'], seq2)
            if similarity_score >= config.similarity_score:
                return True
        return False
    
    def _move_edges(self, node, node_tuple1, node_tuple2):
        """
        Move edges from the node to be merged to the node to be kept.
        
        This method transfers all connections from the node being merged to the node
        being kept, adjusting edge orientations as needed based on the relative
        orientations of the two nodes.
        
        Parameters:
            node: The common neighbor node that connects to both nodes
            node_tuple1: Tuple (node_id, from_direction, to_direction) for the node to keep
            node_tuple2: Tuple (node_id, from_direction, to_direction) for the node to merge
        """
        node_to_keep = node_tuple1[0]
        node_to_merge = node_tuple2[0]
        if node_tuple1[1] == node_tuple2[1]:
            # If nodes have the same orientation, directly transfer edges
            for u, v, data in self.graph.in_edges(node_to_merge, data=True): 
                if u == node:
                    continue
                if self.graph.has_edge(u, node_to_keep):
                    continue
                self.graph.add_edge(u, node_to_keep, key=0, **data)
            for u, v, data in self.graph.out_edges(node_to_merge, data=True): 
                if v == node:
                    continue
                if self.graph.has_edge(node_to_keep, v):
                    continue
                self.graph.add_edge(node_to_keep, v, key=0, **data) 
        else:
            # If nodes have different orientations, adjust edge orientations
            for u, v, data in self.graph.in_edges(node_to_merge, data=True): 
                if u == node:
                    continue
                if self.graph.has_edge(u, node_to_keep):
                    continue                
                if data['label'][2] == '+':
                    data['label'] = data['label'][:2] + '-'
                elif data['label'][2] == '-':
                    data['label'] = data['label'][:2] + '+'
                
                self.graph.add_edge(u, node_to_keep, key=0, **data)
            for u, v, data in self.graph.out_edges(node_to_merge, data=True): 
                if v == node:
                    continue
                if self.graph.has_edge(u, node_to_keep):
                    continue
                if data['label'][0] == '+':
                    data['label'] =  '-' + data['label'][1:]
                elif data['label'][0] == '-':
                    data['label'] =  '+' + data['label'][1:]             
                self.graph.add_edge(node_to_keep, v, key=0, **data)             
        edges_to_remove = list(self.graph.in_edges(node_to_merge, keys=True)) + list(self.graph.out_edges(node_to_merge, keys=True))
        self.graph.remove_edges_from(edges_to_remove)
    
    def _judge_node_property(self, node, edge):
        """
        Extract node properties from an edge, determining the orientation relationship.
        
        Parameters:
            node: The reference node
            edge: Tuple (u, v, key, data) representing an edge
            
        Returns:
            tuple: (connected_node_id, from_direction, to_direction)

        """
        if node == edge[1]:
            from_direction = edge[3]['label'][0]
            to_direction = edge[3]['label'][2]
            return (edge[0], from_direction, to_direction)
        elif node == edge[0]:
            direction_dict = {"+/+": "-/-", "-/-": "+/+", "+/-": "+/-", "-/+": "-/+"}
            label = direction_dict[edge[3]['label']]
            return (edge[1], label[0], label[2])
        else:
            raise ValueError("node not in edge!")

    def _merge_node_property(self, node_tuple1, node_tuple2):
        """
        Merge properties from the node to be merged into the node to be kept.
        
        This method combines sequence data and other properties (like depth and k-mer count)
        from both nodes, creating a consensus sequence and updating the properties
        of the node being kept.
        
        Parameters:
            node_tuple1: Tuple (node_id, from_direction, to_direction) for the node to keep
            node_tuple2: Tuple (node_id, from_direction, to_direction) for the node to merge
        """
        node_to_keep = self.graph.nodes[node_tuple1[0]]
        node_to_merge = self.graph.nodes[node_tuple2[0]]

        # Combine depth (DP) and k-mer count (KC) properties
        if 'DP' in node_to_keep and 'DP' in node_to_merge:
            node_to_keep['DP'] += node_to_keep['DP']
        if 'KC' in node_to_keep and 'KC' in node_to_merge:
            node_to_keep['KC'] += node_to_merge['KC']
            
        # Create consensus sequence, adjusting for orientation if needed
        seq2 = node_to_merge['seq'] if node_tuple1[1]==node_tuple2[1] else shared.reverse_complement(node_to_merge['seq'])
        seq = bs.create_consensus_sequence(node_to_keep['seq'], seq2)
        node_to_keep['seq'] = seq

    def _degenerated_neighbour(self, node):
        """
        Update sequences of neighboring nodes to maintain consistency after merging.
        
        After merging nodes, this method updates the sequences of neighboring nodes
        to ensure they remain consistent with the merged node's sequence at the
        connection points, using degenerate bases where appropriate.
        
        Parameters:
            node: The node whose neighbors need to be updated
        """
        # Process incoming edges (neighbors that connect to this node)
        for u, v, data in self.graph.in_edges(node, data=True): 
            # v == node
            if u == v:
                continue
            seq_u = self.graph.nodes[u]['seq']
            seq_v = self.graph.nodes[v]['seq']
            if data['label'] == '-/-':
                self.graph.nodes[u]['seq'] = bs.cigar_judge_connect(seq_v, seq_u, data['cigar'])
            elif data['label'] == '-/+':
                self.graph.nodes[u]['seq'] = bs.cigar_judge_connect(shared.reverse_complement(seq_v), seq_u, data['cigar'])
            elif data['label'] == '+/-':
                self.graph.nodes[u]['seq'] = shared.reverse_complement(
                    bs.cigar_judge_connect(seq_v, shared.reverse_complement(seq_u), data['cigar']))
            elif data['label'] == '+/+':
                self.graph.nodes[u]['seq'] = shared.reverse_complement(
                    bs.cigar_judge_connect(shared.reverse_complement(seq_v), shared.reverse_complement(seq_u), data['cigar']))
        
        # Process outgoing edges (neighbors that this node connects to)
        for u, v, data in self.graph.out_edges(node, data=True):
            # u == node    
            if u == v:
                continue
            seq_u = self.graph.nodes[u]['seq']
            seq_v = self.graph.nodes[v]['seq']
            if data['label'] == '+/+':
                self.graph.nodes[v]['seq'] = bs.cigar_judge_connect(seq_u, seq_v, data['cigar'])
            elif data['label'] == '-/+':
                self.graph.nodes[v]['seq'] = bs.cigar_judge_connect(shared.reverse_complement(seq_u), seq_v, data['cigar'])
            elif data['label'] == '+/-':
                self.graph.nodes[v]['seq'] = shared.reverse_complement(
                    bs.cigar_judge_connect(seq_u, shared.reverse_complement(seq_v), data['cigar']))
            elif data['label'] == '-/-':
                self.graph.nodes[v]['seq'] = shared.reverse_complement(
                    bs.cigar_judge_connect(shared.reverse_complement(seq_u), shared.reverse_complement(seq_v), data['cigar']))
                    
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
            node: The node whose edges need to be classified
            
        Returns:
            dict: Dictionary with four lists of classified edges
        """
        classified_edges = {"plus_in":[], "plus_out":[], "minus_in":[], "minus_out":[]}
        node_edges = self.graph.in_edges(node, keys=True, data=True)
        
        # Classify outgoing edges
        for u, v, key, data in self.graph.out_edges(node, keys=True, data=True):
            if u == v:
                continue
            if data['label'][0] =='+':
                classified_edges["plus_out"].append((u, v, key, data))
            elif data['label'][0] =='-':
                classified_edges["minus_out"].append((u, v, key, data))
            else:
                raise ValueError("Wrong label! The label shoud be +/+, +/-, -/- or -/+")
        
        # Classify incoming edges
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