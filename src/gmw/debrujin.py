from gfaLib.abstractions import Orientation
from gfaLib.gfaparser import GFAParser
class AbstractGraph:

    def __init__(self):
        self.segments: dict[str, dict] = {}
        self.lines: dict[tuple[str, str], dict] = {}
        self.paths: dict[str, dict] = {}
        self.headers: list[dict] = []
        
    def save_graph(self, output_file, minimal=False, output_format=False):
        GFAParser.save_graph(graph=self, output_path=output_file, force_format=output_format, minimal_graph=minimal)

class GraphFromNetwork(AbstractGraph):
    def __init__(self, network):
        super.__init__()
        self.network = network
        # Parsing the gfa file
        for node in network.nodes:
            self.segments[node] = {**network.nodes[node]}
        
        for u, v, data in network.edges(data=True):
            "print hehe"
            self.lines[(u,v)] = {"orientation":(data['label'][0], data['label'][2]), "ARG5": data['cigar']}

graph = GraphFromNetwork("network")