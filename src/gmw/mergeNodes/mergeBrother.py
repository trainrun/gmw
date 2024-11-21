import mergeNodes.basicFunction as bs
import shared
import config

class BrotherMerger:
    """not """
    def __init__(self, graph):
        self.graph = graph

    def merge_brother(self):
        nodes_to_remove = []
        for node in self.graph.nodes:
            out_node_edges = self._classify_edges(node)
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
        node1 = self.graph.nodes[node_tuple1[0]]
        node2 = self.graph.nodes[node_tuple2[0]]
        if node1['length'] == node2['length']:
            seq2 = node2['seq'] if node_tuple1[1]==node_tuple2[1] else shared.reverse_complement(node2['seq'])
            similarity_score = bs.calculate_similarity(node1['seq'], seq2)
            if similarity_score >= config.similarity_score:
                return True
        return False
    
    def _move_edges(self, node, node_tuple1, node_tuple2):
        node_to_keep = node_tuple1[0]
        node_to_merge = node_tuple2[0]
        if node_tuple1[1] == node_tuple2[1]:
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
        node_to_keep = self.graph.nodes[node_tuple1[0]]
        node_to_merge = self.graph.nodes[node_tuple2[0]]

        if 'DP' in node_to_keep and 'DP' in node_to_merge:
            node_to_keep['DP'] += node_to_keep['DP']
        if 'KC' in node_to_keep and 'KC' in node_to_merge:
            node_to_keep['KC'] += node_to_merge['KC']
        seq2 = node_to_merge['seq'] if node_tuple1[1]==node_tuple2[1] else shared.reverse_complement(node_to_merge['seq'])
        seq = bs.create_consensus_sequence(node_to_keep['seq'], seq2)
        node_to_keep['seq'] = seq

    def _degenerated_neighbour(self, node):
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