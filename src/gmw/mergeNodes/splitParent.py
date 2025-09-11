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

import mergeNodes.basicFunction as bs
import shared
import config

class ParentSpliter:
    """Split a parent (junction) node into two nodes based on topology and depth.

    This class inspects candidate parent nodes with specific in/out edge patterns,
    verifies structural constraints (triad structure), sequencing depth consistency,
    and sequence similarity between candidate child nodes, then rewires edges and
    merges node attributes accordingly to split the parent node.

    Parameters
    ----------
    graph : networkx.MultiDiGraph
        Assembly graph with node attributes: seq, length, DP, KC, ST, EN, AC, TP, OR
        and edge attributes: label ('+/+', '+/-', '-/-', '-/+') and cigar.
    """
    def __init__(self, graph):
        self.graph = graph

    def split_parent(self):
        """Detect and split parent nodes that form a valid diamond-like topology.

        The algorithm:
        1) For each node, classify its incident edges by orientation.
        2) Require exactly two incoming (plus_in + minus_out) and two outgoing
           (plus_out + minus_in) oriented edges.
        3) Validate triad constraints for both strands, check depth compatibility,
           and confirm sibling similarity for both sides.
        4) Rewire edges between matching pairs and merge attributes; remove parent node.

        Notes
        -----
        - Operates in-place, modifies the graph and removes nodes.
        - Uses configuration thresholds from config (e.g., split_depth_multi, split_similarity_score).
        """
        # for u,v,data in self.graph.edges(data=True):
        #     print(u,v,data)
        # for node in self.graph.nodes:
        #     print(self.graph.nodes[node])
        nodes_to_remove = []
        for node in list(self.graph.nodes):
            # Classify edges around current node by orientation
            out_node_edges = self._classify_edges(node)
            # Require exactly two edges on each orientation side
            if (len(out_node_edges['plus_out']) + len(out_node_edges['minus_in']) != 2 or
                len(out_node_edges['plus_in']) + len(out_node_edges['minus_out']) != 2):
                continue
            # Build plus/minus strings (two candidate neighbors on each side)
            plus_string =  [*out_node_edges['plus_in'], *out_node_edges['minus_out']]
            minus_string = [*out_node_edges['plus_out'], *out_node_edges['minus_in']]
            
            # Extract neighbor node and orientation relationship
            plus_tuple1 = self._judge_node_property(node, plus_string[0])
            plus_tuple2 = self._judge_node_property(node, plus_string[1])
            minus_tuple1 = self._judge_node_property(node, minus_string[0])
            minus_tuple2 = self._judge_node_property(node, minus_string[1])
            # Validate triad structure on both sides
            if not self._check_triadius(node, plus_tuple1, plus_tuple2):
                continue
            if not self._check_triadius(node, minus_tuple1, minus_tuple2):
                continue
            if not self._check_depth(node, plus_tuple1, plus_tuple2, minus_tuple1, minus_tuple2):
                continue
            if not self._judge_brother(plus_tuple1, plus_tuple2):
                continue
            if not self._judge_brother(minus_tuple1, minus_tuple2):
                continue
            
            self._move_edges(node, plus_tuple1, plus_tuple2, minus_tuple1, minus_tuple2)
            nodes_to_remove.append(node)
        self.graph.remove_nodes_from(nodes_to_remove)                    
        # for u,v,data in self.graph.edges(data=True):
        #     print(u,v,data)
        # for node in self.graph.nodes:
        #     print(self.graph.nodes[node])
    def _check_triadius(self, node, tuple1, tuple2):
        """Verify that two neighbor tuples form a valid triad with the center node.

        Requirements
        ------------
        - The edge direction relative to the center must match for both tuples.
        - Each neighbor must be a simple degree-1 on the relevant side to avoid branching.

        Parameters
        ----------
        node : hashable
            Center parent node id.
        tuple1, tuple2 : tuple
            (neighbor_node, from_dir, to_dir, cigar)

        Returns
        -------
        bool
            True if triad constraints hold; otherwise False.
        """
        if tuple1[2] != tuple2[2]:
            raise ValueError("Wrong edges in splitParent.")
        edge1 = self._classify_edges(tuple1[0])
        edge2 = self._classify_edges(tuple2[0])
        # For '+' orientation, allow only one outgoing path; similarly for '-'
        if tuple1[1] == '+':
            if len(edge1['plus_out']) + len(edge1['minus_in']) != 1:
                return False
        if tuple1[1] == '-':
            if len(edge1['plus_in']) + len(edge1['minus_out']) != 1:
                return False
        if tuple2[1] == '+':
            if len(edge2['plus_out']) + len(edge2['minus_in']) != 1:
                return False
        if tuple2[1] == '-':
            if len(edge2['plus_in']) + len(edge2['minus_out']) != 1:
                return False
        return True
    
    def _check_depth(self, node, plus_tuple1, plus_tuple2, minus_tuple1, minus_tuple2):
        """Check depth values (DP) among center and neighbors for consistency.

        Rules
        -----
        - For each side (plus/minus), the ratio between the two DP values must be
          at least config.split_depth_multi.
        - The combined DP of each side must be consistent with the center DP within
          the same multiplier threshold.

        Parameters
        ----------
        node : hashable
            Center parent node id.
        plus_tuple1, plus_tuple2, minus_tuple1, minus_tuple2 : tuple
            Neighbor tuples as produced by _judge_node_property.

        Returns
        -------
        bool
            True if depth constraints hold; otherwise False.
        """
        depth1 = self.graph.nodes[node]['DP']
        plus_dp1 = self.graph.nodes[plus_tuple1[0]]['DP']
        plus_dp2 = self.graph.nodes[plus_tuple2[0]]['DP']
        minus_dp1 = self.graph.nodes[minus_tuple1[0]]['DP']
        minus_dp2 = self.graph.nodes[minus_tuple2[0]]['DP']
        if max(plus_dp1,plus_dp2) / min(plus_dp1, plus_dp2) < config.split_depth_multi:
            return False
        if max(minus_dp1,minus_dp2) / min(minus_dp1, minus_dp2) < config.split_depth_multi:
            return False
        plus_dp = plus_dp1 + plus_dp2
        minus_dp = minus_dp1 + minus_dp2
        # Center DP should be consistent with side sums within threshold
        if max(depth1, plus_dp, minus_dp) / min(depth1, plus_dp, minus_dp) > config.split_depth_multi:
            return False
        return True
    
    def _judge_node_property(self, node, edge):
        """Normalize edge orientation to obtain neighbor relation to the center node.

        Given an incident edge, return a tuple describing the neighbor node,
        orientation of traversal from the center, and the CIGAR used for sequence
        merging when applicable.

        Parameters
        ----------
        node : hashable
            Center node id.
        edge : tuple
            (u, v, key, data) from networkx in/out edges iteration.

        Returns
        -------
        tuple
            (neighbor_node, from_dir, to_dir, cigar)
        """
        if node == edge[1]:
            from_direction = edge[3]['label'][0]
            to_direction = edge[3]['label'][2]
            return (edge[0], from_direction, to_direction, edge[3]['cigar'])
        elif node == edge[0]:
            # Reverse orientation when the center is the source of the edge
            direction_dict = {"+/+": "-/-", "-/-": "+/+", "+/-": "+/-", "-/+": "-/+"}
            label = direction_dict[edge[3]['label']]
            return (edge[1], label[0], label[2], edge[3]['cigar'])
        else:
            raise ValueError("node not in edge!")

    def _move_edges(self, node, plus_tuple1, plus_tuple2, minus_tuple1, minus_tuple2):
        """Rewire edges between selected neighbor pairs and merge node attributes.

        This function creates two new edges connecting high/low DP neighbors across
        sides with appropriate orientation labels, removes old edges around the
        center node, and merges the center node's attributes into the remaining
        neighbor nodes based on CIGAR-guided sequence merging.

        Parameters
        ----------
        node : hashable
            Center node being split.
        plus_tuple1, plus_tuple2, minus_tuple1, minus_tuple2 : tuple
            Neighbor tuples as returned by _judge_node_property.
        """
        plus_low_tuple = None
        plus_high_tuple = None
        minus_low_tuple = None
        minus_high_tuple = None
        
        # Rank neighbors by depth (DP) on each side
        if self.graph.nodes[plus_tuple1[0]]['DP'] >  self.graph.nodes[plus_tuple2[0]]['DP']:
            plus_low_tuple = plus_tuple2
            plus_high_tuple = plus_tuple1
        else:
            plus_low_tuple = plus_tuple1
            plus_high_tuple = plus_tuple2
        if self.graph.nodes[minus_tuple1[0]]['DP'] >  self.graph.nodes[minus_tuple2[0]]['DP']:
            minus_low_tuple = minus_tuple2
            minus_high_tuple = minus_tuple1
        else:
            minus_low_tuple = minus_tuple1
            minus_high_tuple = minus_tuple2        
        
        # Build edge label according to orientation relation
        if plus_low_tuple[1] == plus_low_tuple[2]:
            label = minus_low_tuple[1] + "/" + minus_low_tuple[2]
        else:
            node_or = '+' if minus_low_tuple[2] == '-' else '-'
            label = minus_low_tuple[1] + "/" + node_or
        # Connect low-depth pair across sides
        self.graph.add_edge(minus_low_tuple[0], plus_low_tuple[0], key=0, label=label, cigar=minus_low_tuple[3])
            
        if plus_high_tuple[1] == plus_high_tuple[2]:
            label = minus_high_tuple[1] + "/" + minus_high_tuple[2]
        else:
            node_or = '+' if minus_high_tuple[2] == '-' else '-'
            label = minus_high_tuple[1] + "/" + node_or
        # Connect high-depth pair across sides
        self.graph.add_edge(minus_high_tuple[0], plus_high_tuple[0], key=0, label=label, cigar=minus_high_tuple[3])    
   
        # Remove all old edges around the center node
        edges_to_remove = list(self.graph.in_edges(node, keys=True)) + list(self.graph.out_edges(node, keys=True))
        self.graph.remove_edges_from(edges_to_remove)
        
        # Merge center node properties into remaining neighbors guided by labels
        low_label = plus_low_tuple[1] + '/' +  plus_low_tuple[2]
        self._merge_node_property(plus_low_tuple[0], node, low_label, plus_low_tuple[3])
        high_label = plus_high_tuple[1] + '/' +  plus_high_tuple[2]
        self._merge_node_property(plus_high_tuple[0], node, high_label, plus_high_tuple[3])
        
    def _merge_node_property(self, keep, merge, label, cigar):
        """Merge attributes and sequences from a node into another.

        Sequences are merged according to the relative orientation described by
        `label`, with reverse-complement applied when necessary. Depth (DP) and
        k-mer count (KC) are combined; positional annotations (ST/EN), accession
        (AC), type (TP), and orientation (OR) are reconciled.

        Parameters
        ----------
        keep : hashable
            Node id to keep.
        merge : hashable
            Node id to merge into the keep node.
        label : str
            One of '+/+', '+/-', '-/-', '-/+'.
        cigar : str
            CIGAR string used for sequence merging.
        """
        node_to_keep = self.graph.nodes[keep]
        node_to_merge = self.graph.nodes[merge]
        new_seq = ""
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
        
        # Combine depth and k-mer counts, weighted by lengths
        if 'DP' in node_to_keep and 'DP' in node_to_merge:
            node_to_keep['DP'] = round((node_to_keep['DP'] * node_to_keep['length'] + node_to_merge['DP'] * node_to_merge['length']
                                ) / (node_to_keep['length'] + node_to_merge['length']),2)
        if 'KC' in node_to_keep and 'KC' in node_to_merge:
            node_to_keep['KC'] += node_to_merge['KC']
        # Adjust positional and annotation metadata
        if node_to_keep['ST'] == 0 or node_to_keep['ST']>node_to_merge['ST']:
            node_to_keep['ST'] = node_to_merge['ST']
        if node_to_keep['EN']<node_to_merge['EN']:
            node_to_keep['EN'] = node_to_merge['EN']
        if node_to_keep['AC'] == '-':
            node_to_keep['AC'] = node_to_merge['AC']
        if node_to_keep['TP'] == '-':
            node_to_keep['TP'] == node_to_merge['TP']
        
        # Resolve orientation when unknown on keep but known on merge
        if node_to_keep['OR'] == '?' and node_to_merge['OR'] != '?':
            if label[0] == label[2]:
                node_to_keep['OR'] = node_to_merge['OR']
            else:
                node_to_keep['OR'] = '+' if node_to_merge['OR'] == '-' else '-'
                
        # Update sequence and length
        node_to_keep['length'] = len(new_seq)
        node_to_keep['seq'] = new_seq
                            
    def _classify_edges(self, node):
        """Classify incident edges by orientation relative to given node.

        Returns a dictionary containing lists of outgoing/incoming edges on
        plus/minus orientations.

        Parameters
        ----------
        node : hashable
            Node id to classify.

        Returns
        -------
        dict
            Keys: 'plus_in', 'plus_out', 'minus_in', 'minus_out'; values are lists
            of (u, v, key, data) tuples.
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
    
    def _judge_brother(self, node_tuple1, node_tuple2):
        """Check whether two candidate neighbor nodes are brothers (highly similar).

        Two nodes are considered brothers if:
        - They have identical lengths, and
        - Their sequences (with orientation normalized) have similarity score
          >= config.split_similarity_score using bs.calculate_similarity.

        Parameters
        ----------
        node_tuple1, node_tuple2 : tuple
            Neighbor tuples as returned by _judge_node_property.

        Returns
        -------
        bool
            True if brother condition holds; otherwise False.
        """
        node1 = self.graph.nodes[node_tuple1[0]]
        node2 = self.graph.nodes[node_tuple2[0]]
        if node1['length'] == node2['length']:
            seq2 = node2['seq'] if node_tuple1[1]==node_tuple2[1] else shared.reverse_complement(node2['seq'])
            similarity_score = bs.calculate_similarity(node1['seq'], seq2)
            if similarity_score >= config.split_similarity_score:
                return True
        return False
