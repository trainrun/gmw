import unfoldGraph.basicFunction as bs
from unfoldGraph.abstrctUnfold import AbstrctUnfolder
import shared
import unfoldGraph.taxon as taxon
class TaxonUnfolder(AbstrctUnfolder):
    # def __init__(self,graph):
        # self.graph = graph
        # self.dir_path = "/home/cwb/chen/20241006HIV_test/FLU_1/filter_assembly/mmm"
        # self.fasta_name = self.dir_path + "/taxon_neibour.fasta"
        # self.kraken_out = self.dir_path + "/out.kraken"
        # self.kraken_db = "/home/cwb/chen/database/nt4kraken"
        # self.force = True        
        # self.taxon_parse = taxon.TaxonParser("/home/cwb/chen/database/taxon/names.dmp", "/home/cwb/chen/database/taxon/nodes.dmp")
        # self.target_taxon = "11308"
        
    def unfold_graph(self):
        fasta_name = self.out_path + "/" + self.prefix + "_seq_for_kraken.fa"
        visual_in = self.out_path + "/" + self.prefix + "_input.html"
        visual_typing = self.out_path + "/" + self.prefix + "_typing.html"
        visual_out = self.out_path + "/" + self.prefix + "_output.html"
        gfa_out = self.out_path + "/" + self.prefix + "_output.gfa"

        self._visual_graph(visual_in)
        if not self.use_gfa_taxon:
            taxon_parse = taxon.TaxonParser(self.names_dmp, self.nodes_dmp)
            self._clear_graph_type()
            if self.kraken_out is None:
                self.kraken_out = self.out_path + "/" + self.prefix + "_kraken_out.txt"
                self._run_kraken(fasta_name)
            
            bs.graph_add_type(self.graph, self.kraken_out, self.taxon_id, taxon_parse)
            self._visual_graph(visual_typing)
            # self.graph.remove_nodes_from(remove_list)
        if not self.disable_ref_unfold:
            self.remove_unknown_nodes = False
        self._remove_components()
        self._merge_nodes()
        self._output_gfa(gfa_out)
        
        self._visual_graph(visual_out)

    
    def _run_kraken(self, fasta_name):
        shared.graph2fasta(self.graph, fasta_name)
        shared.run_kraken(fasta_name, self.kraken_db, self.kraken_out)
        
    def _clear_graph_type(self):
        for node in self.graph.nodes:
            self.graph.nodes[node]['TP'] = '-'
    
    