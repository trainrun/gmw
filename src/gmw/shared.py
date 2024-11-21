import os
import subprocess
import networkx as nx
import logging

import config

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
    计算 DNA 序列的反向互补序列，支持小写字母和简并碱基
    """
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

    valid_bases = set(complement_dict.keys())
    if not all(base in valid_bases for base in seq):
        raise ValueError("Invalid DNA sequence")
    return ''.join([complement_dict[base] for base in reversed(seq)])

def graph2fasta(graph, file_path, orientation=False):
    with open(file_path, 'w') as f:
        for node in graph.nodes:
            f.write(">" + str(node) + "\n")
            if orientation and graph.nodes[node]['OR'] == '-':
                f.write(reverse_complement(graph.nodes[node]['seq']) + "\n")
            else:
                f.write(graph.nodes[node]['seq'] + "\n")

def run_blast(fasta, db, out, threads=config.max_threads):
    commands = [config.blast_path, '-query', fasta, '-db', db, '-num_threads', str(threads), '-outfmt', '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs', '-out',
            out, '-max_target_seqs', '1']
    logger.info(' '.join(commands))
    result = subprocess.check_call(commands)
    
def run_kraken(fasta, db, output, threads=config.max_threads):
    commands = [config.kraken_path, "--db", db, "--output", output, "--threads", str(threads), "--confidence", str(config.confidence), "--memory-mapping", fasta]
    logger.info(' '.join(commands))
    result = subprocess.check_call(commands)
    
def remove_unknown_components(graph):
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
    isolated_nodes = [node for node in graph.nodes() if (graph.degree(node) == 0 and graph.nodes[node]['length'] < length)]    
    graph.remove_nodes_from(isolated_nodes)
    
def remove_unknown_nodes(graph):
    unknown_nodes = [node for node in graph.nodes() if (graph.nodes[node]['TP'] == '-' and graph.nodes[node]['AC'] == '-')]    
    graph.remove_nodes_from(unknown_nodes)

def remove_contaminated_nodes(graph):
    contaminate_nodes = [node for node in graph.nodes() if (graph.nodes[node]['TP'] == 'contaminate')]
    graph.remove_nodes_from(contaminate_nodes)    