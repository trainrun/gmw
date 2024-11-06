
from unfoldGraph.abstrctUnfold import AbstrctUnfolder
class DepthUnfolder(AbstrctUnfolder):
    # def __init__(self, graph, multiple=10):
    #     self.graph = graph
    #     self.multiple = multiple
    #     self.dir_path = "/home/cwb/chen/20241006HIV_test/FLU_1/filter_assembly/mmm"
    #     self.force = True
                
    def unfold_graph(self):
        visual_in = self.out_path + "/" + self.prefix + "_depth_input.html"
        visual_out = self.out_path + "/" + self.prefix + "_depth_output.html"
        gfa_out = self.out_path + "/" + self.prefix + "_depth_output.gfa"
        self._visual_graph(visual_in)        
        
        edges_remove = []
        for node1, node2, key in self.graph.edges(keys=True):
            depth1 = self.graph.nodes[node1]['DP']
            depth2 = self.graph.nodes[node2]['DP']
            if depth1/depth2 > self.depth_discrepancy or depth2/depth1 > self.depth_discrepancy:
                edges_remove.append((node1, node2, key))
        self.graph.remove_edges_from(edges_remove)

        if self.disable_ref_unfold and self.disable_taxon_unfold:
            self.keep_unknown_components = True
            self.remove_unknown_nodes = False
        self._remove_components()
        self._merge_nodes()
        self._output_gfa(gfa_out)
        self._visual_graph(visual_out)  