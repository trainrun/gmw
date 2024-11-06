import os
import logging

import visualize
import shared
import config
import gfaLib
from mergeNodes.mergeNeighbour import NeighbourMerger
class AbstrctUnfolder:
    def __init__(self, graph, out_path, prefix, threads, force,
                 disable_taxon_unfold, disable_ref_unfold, disable_depth_unfold, disable_gc_unfold,
                 use_gfa_taxon, kraken_out, names_dmp, nodes_dmp, taxon_id, taxon_name, kraken_db,
                 use_gfa_ref, blast_out, blast_db, position_distance,
                 depth_discrepancy,
                 gc_discrepancy,
                 remove_unknown_nodes, keep_unknown_components, keep_short_isolated_nodes,
                 disable_merge_neighbor, merge_brother, split_parent,
                 visual, contig_shape
                 ):
        self.graph = graph  
        self.out_path = out_path + "/" + self.__class__.__name__
        self.prefix = prefix  
        self.threads = threads  
        self.force = force
        self.disable_taxon_unfold = disable_taxon_unfold
        self.disable_ref_unfold = disable_ref_unfold
        self.disable_depth_unfold = disable_depth_unfold
        self.disable_gc_unfold = disable_gc_unfold
        self.use_gfa_taxon = use_gfa_taxon  
        self.kraken_out = kraken_out  
        self.names_dmp = names_dmp
        self.nodes_dmp = nodes_dmp  
        self.taxon_id = taxon_id  
        self.taxon_name = taxon_name  
        self.kraken_db = kraken_db  
        self.use_gfa_ref = use_gfa_ref  
        self.blast_out = blast_out  
        self.blast_db = blast_db  
        self.position_distance = position_distance  
        self.depth_discrepancy = depth_discrepancy
        self.keep_short_isolated_nodes = keep_short_isolated_nodes
        self.gc_discrepancy = gc_discrepancy
        self.remove_unknown_nodes = remove_unknown_nodes
        self.keep_unknown_components = keep_unknown_components 
        self.disable_merge_neighbor = disable_merge_neighbor
        self.merge_brother = merge_brother  
        self.split_parent = split_parent  
        self.visual = visual
        self.contig_shape = contig_shape
        
        self.neighbour_merger = NeighbourMerger(self.graph)
        self._create_unfold_dir()
        self.logger = logging.getLogger("gmw")
        self.logger.info(f"Start unfold using {self.__class__.__name__}")
    def _visual_graph(self, file):
        if self.visual:
            visualize.print_graph(self.graph, file, self.contig_shape)
    
    def _remove_components(self):
        shared.remove_contaminated_nodes(self.graph)
        
        if not self.keep_unknown_components:
            shared.remove_unknown_components(self.graph)
        if self.remove_unknown_nodes:
            shared.remove_unknown_nodes(self.graph)
        if not self.keep_short_isolated_nodes:
            shared.remove_short_isolated_nodes(self.graph, config.short_offset)
            
    def _merge_nodes(self):
        if not self.disable_merge_neighbor:
            self.neighbour_merger.merge_neibour()
            
    def _output_gfa(self, gfa_name):
        gfaLib.GraphFromNetwork(self.graph).save_graph(gfa_name)

    def _create_unfold_dir(self):
        if not os.path.exists(self.out_path):
            os.makedirs(self.out_path)