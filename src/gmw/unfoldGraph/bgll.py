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

import community
import networkx as nx

import config

class BGLLCluster:
    """
    BGLL (Blondel-Guillaume-Lambiotte-Lefebvre) clustering implementation for graph community detection.
    
    This class implements the Louvain algorithm for community detection in genome assembly graphs,
    specifically designed to identify and classify taxonomic partitions. It uses community structure
    to infer taxonomic labels for unknown nodes based on their neighbors in detected communities.
    
    Attributes
    ----------
    graph : networkx.Graph
        The input graph to perform community detection on
    """
    
    def __init__(self, graph):
        self.graph = graph
        
    def _obtain_community(self, partition):
        community_nodes = {}
        for node, community_id in partition.items():
            if community_id not in community_nodes:
                community_nodes[community_id] = []
            community_nodes[community_id].append(node)
        return community_nodes
    
    def _change_partition_taxon(self,nodes):
        unknown = 0 
        contaminate = 0 
        target = 0
        for node in nodes:
            if self.graph.nodes[node]['TP'] == 'contaminate':
                contaminate += 1
            elif self.graph.nodes[node]['TP'] == 'target':
                target += 1
            else:
                unknown += 1
        #print(unknown, contaminate, target)
        
        if self._compute_community_type(unknown, contaminate, target):
            for node in nodes:
                if self.graph.nodes[node]['TP'] == '-':
                    self.graph.nodes[node]['TP'] = 'contaminate_infer'
        elif self._compute_community_type(unknown, target, contaminate):
            for node in nodes:
                if self.graph.nodes[node]['TP'] == '-':
                    self.graph.nodes[node]['TP'] = 'target_infer'
    
    def _compute_community_type(self, unknown, num1, num2):
        """
        Determine if a community should be classified as a specific taxonomic type.
        
        This method uses configurable thresholds to determine if one taxonomic type
        is sufficiently dominant in a community to justify inferring that type for
        unknown nodes.
        
        Parameters
        ----------
        unknown : int
            Number of nodes with unknown taxonomic labels
        num1 : int
            Number of nodes of the first taxonomic type being tested
        num2 : int
            Number of nodes of the second taxonomic type being tested
            
        Returns
        -------
        bool
            True if num1 type is sufficiently dominant, False otherwise
        """
        if num1 == 0:
            return False
        if num2 == 0:
            return True
        elif num1/num2 < config.partition_offset:
            return False
        # elif (num1-num2) * config.partition_offset > unknown:
        #     return True
        return True
                    
    def louvain_algorithm(self):
        partition = community.best_partition(self.graph.to_undirected(),resolution=config.resolution)
        community_nodes = self._obtain_community(partition)
        #print(len(community_nodes))
        for _, nodes in community_nodes.items():
            self._change_partition_taxon(nodes)

        
