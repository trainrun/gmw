
import unfoldGraph.basicFunction as bs
import shared
from unfoldGraph.abstrctUnfold import AbstrctUnfolder

class RefUnfolder(AbstrctUnfolder):
    # def __init__(self, graph, offset_bp=150):
    #     self.graph = graph
    #     self.offset_bp = offset_bp
    #     self.dir_path = "/home/cwb/chen/20241006HIV_test/FLU_1/filter_assembly/mmm"
    #     self.force = True
    #     self.fasta_name = self.dir_path + "/merge_neibour.fasta"
    #     self.blast_name = self.dir_path + "/merge_neibour.blast"
    #     self.blast_db = "/home/cwb/chen/database/refseq/ref_viruses_rep_genomes"
        
    def unfold_graph(self):
        fasta_name = self.out_path + "/" + self.prefix + "_seq_for_blast.fa"
        visual_in = self.out_path + "/" + self.prefix + "_ref_input.html"
        visual_typing = self.out_path + "/" + self.prefix + "_ref_typing.html"
        visual_out = self.out_path + "/" + self.prefix + "_ref_output.html"
        gfa_out = self.out_path + "/" + self.prefix + "_ref_output.gfa"

        self._visual_graph(visual_in)
        
        if not self.use_gfa_ref:
            self._clear_graph_accession()
            if self.blast_out is None:
                self.blast_out = self.out_path + "/" + self.prefix + "_blast_out.txt"
                self._run_blast(fasta_name)
                
            bs.graph_add_accession(self.graph, self.blast_out)
            self._visual_graph(visual_typing)
            
        self._remove_wrong_connction()
        self._remove_components()
        self._merge_nodes()
        self._output_gfa(gfa_out)
        self._visual_graph(visual_out)
        
        
    def _remove_wrong_connction(self):
        edges_remove = list(self.graph.edges(keys=True))
        for node1, node2, key in edges_remove:
            edge = self.graph[node1][node2][key]
            if not self._accession_judge(self.graph.nodes[node1]['AC'], self.graph.nodes[node2]['AC']):
                self.graph.remove_edge(node1, node2, key)
                continue
            if not self._orientaion_judge(self.graph.nodes[node1]['OR'], self.graph.nodes[node2]['OR'], edge['label']):
                self.graph.remove_edge(node1, node2, key)
                continue
            if not self._position_judge(self.graph.nodes[node1]['ST'], self.graph.nodes[node1]['EN'],
                                        self.graph.nodes[node2]['ST'], self.graph.nodes[node2]['EN']):
                self.graph.remove_edge(node1, node2, key)

        
            
    def _run_blast(self, fasta_name):
        shared.graph2fasta(self.graph, fasta_name)
        shared.run_blast(fasta_name, self.blast_db, self.blast_out)

            
    def _accession_judge(self, str1, str2):
        # if str1 == '-' and str2 == '-':
        #     return False
        if str1 == "-" or str2 == "-":
            return True
        if str1 == str2:
            return True
        return False

    def _position_judge(self, start1, end1, start2, end2):
        if start1 == 0 or start2 == 0:
            return True
        if abs(end1-start2) <= self.position_distance or abs(end2-start1) <= self.position_distance:
            return True
        return False
    
    def _orientaion_judge(self, i, j, edge_orientaion):
        count = 0
        if i == '-':
            count += 1
        elif i != '+':
            return True
        if j == '-':
            count += 1
        elif j != '+':
            return True
        
        o_list = edge_orientaion.split('/')
        if o_list[0] == '-':
            count += 1
        if o_list[1] == '-':
            count += 1
        
        if count % 2 == 0:
            return True

        return False

    def _clear_graph_accession(self):
        for node in self.graph.nodes:
            self.graph.nodes[node]['AC'] = '-'
            self.graph.nodes[node]['ST'] = 0
            self.graph.nodes[node]['EN'] = 0
            self.graph.nodes[node]['OR'] = '?'
            