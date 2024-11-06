import json

import config
import visualize.htmlStr as htmlStr

def print_graph(graph, file_path, contig_shape):
    if contig_shape == 'dot':
        print_dot(graph, file_path)
    else:
        print_line(graph, file_path)
        
def print_dot(graph, file_path, edge_color=config.edge_color, node_shape=config.node_shape):
    data_dict = convert_node_json(graph)
    option_dict = {"physics": True, "edges":{"color": edge_color}, "nodes":{"shape": node_shape}}
    with open(file_path, "w") as file:
        file.write(htmlStr.head)
        file.write(json.dumps(data_dict, indent=None))
        file.write(htmlStr.middle)
        file.write(json.dumps(option_dict, indent=None))
        file.write(htmlStr.tail)

def convert_node_json(graph):
    color_dict = {"unknown": "grey", "target": "pink", "contaminate": "skyblue"}
    graph_data = {
        "nodes": [],
        "edges": []
    }
    for node in graph.nodes:
        node_color = color_dict[infer_color(graph, node)]
        graph_data["nodes"].append({"id": node, "label": formatLength(graph.nodes[node]['length']), "color": node_color})
    for u, v, key in graph.edges:
         graph_data["edges"].append({"from": u, "to": v})
         
    return graph_data

def print_line(graph, file_path, edge_color=config.edge_color, node_shape=config.node_shape):
    data_dict = convert_edge_json(graph)
    option_dict = {"physics": {"enabled": True, "stabilization": True}, 
                   "nodes": {"size": 10, "shape": node_shape, "color": {"background": "rgba(255,255,255,0)", "border": 'rgba(255, 255, 255, 0)'}},
                   "edges":{"color": edge_color,"font": {"size": 14, "align": 'horizontal'}}
                   }
    with open(file_path, "w") as file:
        file.write(htmlStr.head)
        file.write(json.dumps(data_dict, indent=None))
        file.write(htmlStr.middle)
        file.write(json.dumps(option_dict, indent=None))
        file.write(htmlStr.tail)

def convert_edge_json(graph):
    length_min = None
    length_max = None
    width_min = None
    width_max = None
    for node in graph.nodes:
        # print( graph.nodes[node]['length'])
        # print( graph.nodes[node]['DP'])
        if length_min is None or length_min > graph.nodes[node]['length']:
            length_min = graph.nodes[node]['length']
        if length_max is None or length_max < graph.nodes[node]['length']:
            length_max = graph.nodes[node]['length']
        if width_min is None or width_min > graph.nodes[node]['DP']:
            width_min = graph.nodes[node]['DP']
        if width_max is None or width_max < graph.nodes[node]['DP']:
            width_max = graph.nodes[node]['DP']
      
    graph_data = {
        "nodes": [],
        "edges": []
    }
    for node in graph.nodes:
        node_color = config.color_dict[infer_color(graph, node)]
        graph_data["nodes"].append({"id": node + '+'})
        graph_data["nodes"].append({"id": node + '-'})
        graph_data["edges"].append({"from": node + '+', "to": node + '-', "color": node_color, "label": formatLength(graph.nodes[node]['length']),
                                    "length": normalize_to_range(graph.nodes[node]['length'], length_min, length_max, config.length_range_mix, config.length_range_max),
                                    "width": normalize_to_range(graph.nodes[node]['DP'], width_min, width_max, config.width_range_mix, config.width_range_max)})
    for u, v, key,data in graph.edges(keys=True, data=True):
         graph_data["edges"].append({"from": u + reverse_label(data['label'][0]), "to": v + data['label'][2]})
         
    return graph_data


def infer_color(graph, node):
    if graph.nodes[node]['TP'] == 'contaminate':
        return 'contaminate'
    if graph.nodes[node]['TP'] == 'target' or graph.nodes[node]['AC'] != '-':
        return "target"
    return "unknown"

def formatLength(length):
    units = ['bp', 'Kb', 'Mb', 'Gb']
    if length < 1000:
        return f"{length}bp"
    i = 0
    while length >= 1000 and i < len(units) - 1:
        length /= 1000
        i += 1
    return f"{length:.1f} {units[i]}"

def reverse_label(label):
    if label == '+':
        return '-'
    if label == '-':
        return '+'
    raise ValueError(f"Wrong label {label}")

def normalize_to_range(data_to_normalize, data_min, data_max, min_range, max_range):
    if data_max == data_min:
        return (max_range + min_range) / 2
    scale = (max_range - min_range) / (data_max - data_min)
    normalized_data = min_range + (data_to_normalize - data_min) * scale
    return normalized_data