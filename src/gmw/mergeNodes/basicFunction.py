
import re  

def remove_redundant_edges(graph):  
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
    consensus = []
    operations = re.findall(r'(\d+)([MIDNSHP])', cigar)
    if len(operations) > 1:
        raise ValueError("Wrong CIGAR sequence.")
    consensus.append(seq1)
    for count, op in operations:
        count = int(count)
        if op == 'M':
            # 取 seq2 的开头部分，减去重叠的部分
            overlap = seq1[-count:]  # seq1 末尾的匹配部分
            if overlap == seq2[:count]:
                consensus.append(seq2[count:])  # 追加 seq2 剩余部分
                # print("connected")
                break
            else:
                raise ValueError("The sequences to merge are not match")
            
        else:
            raise ValueError("Wrong CIGAR sequence.")
    return ''.join(consensus)


