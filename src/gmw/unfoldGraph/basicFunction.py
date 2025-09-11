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
            identity = float(line[2])
            orientation = '+'
            if start > end:
                orientation = '-'
                start, end = end, start
            if node_name in node_dict:
                node_dict[node_name].append([line[1], orientation, start, end, identity])
            else:
                node_dict[node_name] = [[line[1], orientation, start, end, identity]]
    for key, value in node_dict.items():
        graph.nodes[key]['AC'] = value[0][0]
        if len(value) == 1:
            graph.nodes[key]["OR"] = value[0][1]
            graph.nodes[key]["ST"] = value[0][2]
            graph.nodes[key]["EN"] = value[0][3]
        else:
            # check_orient = value[0][1]
            # check_ident = value[0][4]
            match_same_orient = True
            match_same_seq = True
            occ_times = 1
            for i in range(1,len(value)):
                if value[0][4] - value[i][4] > config.identity_discrepency:
                    continue
                if value[0][0] != value[i][0]:
                    match_same_seq = False
                else:
                    occ_times += 1
                if value[0][1] != value[i][1]:
                    change_orient = False
            if match_same_seq and occ_times == 1:
                graph.nodes[key]["OR"] = value[0][1]
                graph.nodes[key]["ST"] = value[0][2]
                graph.nodes[key]["EN"] = value[0][3]                
            elif match_same_orient:
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
