
from gfaLib.abstractions import GFALine, Orientation, reverse
from gfaLib.gfaparser import GFAParser
from gfaLib.abstractGraph import AbstractGraph

class GraphFromNetwork(AbstractGraph):
    def __init__(self, network):
        super().__init__()
        self.metadata: dict = {'version': 'unknown', 'next_node_name':  'unknown', 'with_sequence': True}
        # Parsing the gfa file
        for node in network.nodes:
            self.segments[node] = {**network.nodes[node]}
        
        for u, v, data in network.edges(data=True):
            "print hehe"
            self.lines[(u,v)] = {"orientation":{(Orientation(data['label'][0]), Orientation(data['label'][2]))}, "ARG5": data['cigar']}
        