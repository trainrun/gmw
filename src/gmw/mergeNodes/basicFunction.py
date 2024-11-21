import re  
import config
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
REVERSE_IUPAC_CODES = {frozenset(v): k for k, v in IUPAC_CODES.items()}

def calculate_similarity(seq1, seq2):
    distance = sum(el1 != el2 for el1, el2 in zip(seq1, seq2))
    similarity_score = (1 - distance / len(seq1)) * 100
    return similarity_score

def get_degenerated_base(base1, base2):
    # """返回两个碱基合并后的简并碱基"""
    # return degenerate_bases.get((base1, base2), 'N')  # 默认返回'N'表示无法确定的简并碱基
    set1 = IUPAC_CODES.get(base1, set())
    set2 = IUPAC_CODES.get(base2, set())
    # 合并集合
    merged_set = set1.union(set2)
    # 从反向映射中查找对应的简并碱基
    return REVERSE_IUPAC_CODES.get(frozenset(merged_set))
def create_consensus_sequence(seq1, seq2):
    """生成共识序列"""
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
    for count, op in operations:
        count = int(count)
        if op == 'M':
            # 取 seq2 的开头部分，减去重叠的部分
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
    """这里返回的是Seq2简并剪辑"""
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
            # if seq1_overlap == consensus:
            #     return seq1_overlap + seq2[count:]
            # elif seq2_overlap == consensus:
            #     return seq2
            # else:
            #     print(f"{seq1}\n{seq2}\n{consensus}")
            #     raise ValueError("Error in cigar_judge_connect")
    raise ValueError("Wrong CIGAR sequence.")
 