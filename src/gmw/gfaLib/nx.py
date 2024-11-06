from networkx import MultiDiGraph
from gfaLib.graphFromFile import GraphFromFile

class GFANetwork:

    @staticmethod
    def compute_backbone(
        graph: GraphFromFile
    ) -> MultiDiGraph:
        backbone: MultiDiGraph = MultiDiGraph()

        for node_name, node_datas in graph.segments.items():
            backbone.add_nodes_from([( node_name, node_datas)])
        for (start, end, key), edge_data in graph.lines.items():
            backbone.add_edge(
                start,
                end,
                key=key,
                label=' | '.join(
                    [f'{x.value}/{y.value}' for (x, y) in edge_data["orientation"]]),
                cigar=edge_data['ARG5']
            )
        return backbone


