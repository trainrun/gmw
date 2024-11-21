from unfoldGraph.abstrctUnfold import AbstrctUnfolder
class Polisher(AbstrctUnfolder):
    def polish(self):
        if not self.disable_ref_unfold:
            rmove_list = []
            for node in self.graph.nodes():
                if self.graph.nodes[node]['AC'] == '-':
                    rmove_list.append(node)
            self.graph.remove_nodes_from(rmove_list)
            
        self._remove_components()
        
        nodes_num = 0
        while nodes_num != self.graph.number_of_nodes():
            nodes_num = self.graph.number_of_nodes()
            self._merge_nodes()

        gfa_out = self.out_path + "/../" +  self.prefix + "_after_unfold.gfa"
        self._output_gfa(gfa_out)