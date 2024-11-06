from unfoldGraph.abstrctUnfold import AbstrctUnfolder
class GCUnfolder(AbstrctUnfolder):
    # def __init__(self, graph, gc_discrepancy=0.1):
    #     self.graph = graph
    #     self.gc_discrepancy = gc_discrepancy
                
    def unfold_graph(self):
        visual_in = self.out_path + "/" + self.prefix + "_gc_input.html"
        visual_out = self.out_path + "/" + self.prefix + "_gc_output.html"
        gfa_out = self.out_path + "/" + self.prefix + "_gc_output.gfa"
        self._visual_graph(visual_in) 
                
        edges_remove = []
        for node1, node2, key in self.graph.edges(keys=True):
            gc1 = self.gc_content(self.graph.nodes[node1]['seq'])
            gc2 = self.gc_content(self.graph.nodes[node2]['seq'])
            if gc1 - gc2 > self.gc_discrepancy or gc2 - gc1 > self.gc_discrepancy:
                edges_remove.append((node1, node2, key))
        self.graph.remove_edges_from(edges_remove)
        if self.disable_ref_unfold and self.disable_taxon_unfold:
            self.keep_unknown_components = True
            self.remove_unknown_nodes = False
        self._remove_components()        
        self._merge_nodes()
        self._output_gfa(gfa_out)
        self._visual_graph(visual_out)
        
    def gc_content(self, seq):
        seq = seq.upper()
        gc_count = seq.count('G') + seq.count('C')
        total_count = len(seq)
        gc_content = gc_count / total_count
        return gc_content