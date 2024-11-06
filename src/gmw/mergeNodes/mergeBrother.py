import mergeNodes.basicFunction as bs
import shared

class BrotherMerger:
    """not complete"""
    def __init__(self, graph):
        self.graph = graph


    def merge_brother(self):
        nodes_to_remove = []
        for node in self.graph.nodes:
            out_node_edges = self._classify_edges(node)
            if len(out_node_edges['plus_out']) + len(out_node_edges['minus_in']) > 1:
                pass
    def merge_neibour(self):
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
        node_to_keep = self.graph.nodes[connect_edge[0]]
        node_to_merge = self.graph.nodes[connect_edge[1]]
        label = connect_edge[3]['label']
        cigar = connect_edge[3]['cigar']
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
        
        if 'DP' in node_to_keep and 'DP' in node_to_merge:
            node_to_keep['DP'] = round((node_to_keep['DP'] * node_to_keep['length'] + node_to_merge['DP'] * node_to_merge['length']
                                ) / (node_to_keep['length'] + node_to_merge['length']),2)
        if 'KC' in node_to_keep and 'KC' in node_to_merge:
            node_to_keep['KC'] += node_to_merge['KC']
        if node_to_keep['ST'] == 0 or node_to_keep['ST']>node_to_merge['ST']:
            node_to_keep['ST'] = node_to_merge['ST']
        if node_to_keep['EN']<node_to_merge['EN']:
            node_to_keep['EN'] = node_to_merge['EN']
        if node_to_keep['AC'] == '-':
            node_to_keep['AC'] = node_to_merge['AC']
        if node_to_keep['TP'] == '-':
            node_to_keep['TP'] == node_to_merge['TP']
        
        if node_to_keep['OR'] == '?' and node_to_merge['OR'] != '?':
            if label[0] == label[2]:
                node_to_keep['OR'] = node_to_merge['OR']
            else:
                node_to_keep['OR'] = '+' if node_to_merge['OR'] == '-' else '-'
                
        node_to_keep['length'] = len(new_seq)
        node_to_keep['seq'] = new_seq
             
    def _move_edges(self, connect_edge):
        node_to_keep = connect_edge[0]
        node_to_merge = connect_edge[1]
        label = connect_edge[3]['label']
        if label == '+/+' or label == '-/-':
            for u, v, data in self.graph.in_edges(node_to_merge, data=True): 
                if u == node_to_keep:
                    continue
                new_key = 0
                if self.graph.has_edge(u, node_to_keep):
                    all_keys = self.graph[u][node_to_keep].keys()
                    new_key = max(all_keys)+1
                self.graph.add_edge(u, node_to_keep, key=new_key, **data)
            for u, v, data in self.graph.out_edges(node_to_merge, data=True): 
                if v == node_to_keep:
                    continue
                new_key = 0
                if self.graph.has_edge(node_to_keep, v):
                    all_keys = self.graph[node_to_keep][v].keys()
                    new_key = max(all_keys)+1
                self.graph.add_edge(node_to_keep, v, key=new_key, **data) 
        elif label == '+/-' or label == '-/+':
            for u, v, data in self.graph.in_edges(node_to_merge, data=True): 
                if u == node_to_keep:
                    continue
                new_key = 0
                if self.graph.has_edge(u, node_to_keep):
                    all_keys = self.graph[u][node_to_keep].keys()
                    new_key = max(all_keys)+1
                if data['label'][2] == '+':
                    data['label'] = data['label'][:2] + '-'
                elif data['label'][2] == '-':
                    data['label'] = data['label'][:2] + '+'
                self.graph.add_edge(u, node_to_keep, key=new_key, **data)
            for u, v, data in self.graph.out_edges(node_to_merge, data=True): 
                if v == node_to_keep:
                    continue
                new_key = 0
                if self.graph.has_edge(u, node_to_keep):
                    all_keys = self.graph[node_to_keep][v].keys()
                    new_key = max(all_keys)+1
                if data['label'][0] == '+':
                    data['label'] =  '-' + data['label'][1:]
                elif data['label'][0] == '-':
                    data['label'] =  '+' + data['label'][1:]
                    
                self.graph.add_edge(node_to_keep, v, key=new_key, **data) 
        else:
            raise ValueError(f"Edges have wrong label {label}")
                        
        edges_to_remove = list(self.graph.in_edges(node_to_merge, keys=True)) + list(self.graph.out_edges(node_to_merge, keys=True))
        self.graph.remove_edges_from(edges_to_remove)
                                   
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