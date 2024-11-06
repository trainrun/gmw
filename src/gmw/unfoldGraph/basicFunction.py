import config

def graph_add_accession(graph, file_path):
    node_dict = {}

    with open(file_path, 'r') as f:
        for line in f:
            line = line.rstrip().split("\t")
            node_name = line[0]
            if node_name not in graph.nodes:
                continue
            if int(line[3]) < config.match_length:
                continue
            if int(line[-1]) < config.query_cover_offset:
                continue
            start = int(line[8])
            end = int(line[9])
            orientation = '+'
            if start > end:
                orientation = '-'
                start, end = end, start
            if node_name in node_dict:
                node_dict[node_name].append([line[1], orientation, start, end])
            else:
                node_dict[node_name] = [[line[1], orientation, start, end]]
    for key, value in node_dict.items():
        graph.nodes[key]['AC'] = value[0][0]
        if len(value) == 1:
            graph.nodes[key]["OR"] = value[0][1]
            graph.nodes[key]["ST"] = value[0][2]
            graph.nodes[key]["EN"] = value[0][3]
        else:
            check_orient = value[0][1]
            change_orient = True
            for i in range(1,len(value)):
                if not check_orient == value[i][1]:
                    change_orient = False
                    break
            if change_orient:
                graph.nodes[key]["OR"] = value[0][1]

def graph_add_type(graph, file_path, taxid, taxon_parser):
    with open(file_path, 'r') as f:
        for line in f:
            line = line.rstrip().split("\t")
            if line[0] == 'U' or line[2] == '1' or line[2] == '0' or line[1] not in graph.nodes:
                continue
            if taxon_parser.check_relationship(line[2], taxid):
                graph.nodes[line[1]]['TP'] = 'target'
            else:
                graph.nodes[line[1]]['TP'] = 'contaminate'
