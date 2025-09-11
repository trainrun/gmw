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

import os
import logging

import visualize
import shared
import config
import gfaLib
from mergeNodes.mergeNeighbour import NeighbourMerger
class AbstrctUnfolder:
    """Coordinator for unfolding/cleaning operations on the assembly graph.

    This class bundles common utilities used by specific unfolding strategies
    (e.g., taxonomy-based, reference-based, depth/GC-based) and provides
    helpers for visualization, component pruning, node merging, and output.

    Parameters
    ----------
    graph : networkx.MultiDiGraph
        The in-memory assembly graph to process.
    out_path : str
        Base output directory. A subdirectory named as the class will be created.
    prefix : str
        Output file name prefix.
    threads : int
        Number of worker threads to use for downstream tools.
    force : bool
        Whether to overwrite existing outputs.
    disable_taxon_unfold, disable_ref_unfold, disable_depth_unfold, disable_gc_unfold : bool
        Switches to disable individual unfolding strategies.
    use_gfa_taxon : bool
        Whether to use taxonomy info embedded in GFA instead of re-running tools.
    kraken_out, names_dmp, nodes_dmp, taxon_id, taxon_name, kraken_db, bgll : str
        Input files/IDs for taxonomy parsing and community detection (if used).
    use_gfa_ref : bool
        Whether to use reference alignments embedded in GFA instead of re-running BLAST.
    blast_out, blast_db : str
        Alignment outputs and databases for reference-based processing.
    position_distance : int
        Max distance to consider syntenic consistency when merging/validating.
    depth_discrepancy, gc_discrepancy : float
        Thresholds for depth and GC content inconsistency checks.
    remove_unknown_nodes : bool
        If True, remove nodes with unknown type labels.
    keep_unknown_components : bool
        If False, drop entire components whose labels are unknown.
    keep_short_isolated_nodes : bool
        If True, keep short isolated nodes; otherwise remove them.
    disable_merge_neighbor : bool
        If True, do not execute neighbor merging stage.
    merge_brother, split_parent : bool
        Flags controlling optional merge/split behaviors in the workflow.
    visual : bool
        Whether to produce graph visualizations.
    contig_shape : str
        Shape used by the visualization module for rendering contigs.
    """

    def __init__(self, graph, out_path, prefix, threads, force,
                 disable_taxon_unfold, disable_ref_unfold, disable_depth_unfold, disable_gc_unfold,
                 use_gfa_taxon, kraken_out, names_dmp, nodes_dmp, taxon_id, taxon_name, kraken_db, bgll,
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
        self.bgll = bgll 
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
        """Render the current graph to an image file if visualization is enabled.

        Parameters
        ----------
        file : str
            Target image file path.
        """
        if self.visual:
            # Delegate to visualization backend with configured contig shape
            visualize.print_graph(self.graph, file, self.contig_shape)
    
    def _remove_components(self):
        """Remove contaminated and optionally unknown/short components from graph.

        The sequence of operations is:
        1) Remove nodes labeled as contaminated using shared.remove_contaminated_nodes.
        2) Optionally remove entire components with unknown labels.
        3) Optionally remove individual unknown nodes.
        4) Optionally drop short isolated nodes using config.short_offset as cutoff.
        """
        # Remove definitely contaminated nodes first
        shared.remove_contaminated_nodes(self.graph)

        # Then handle unknown components/nodes as requested by user settings
        if not self.keep_unknown_components:
            shared.remove_unknown_components(self.graph)
        if self.remove_unknown_nodes:
            shared.remove_unknown_nodes(self.graph)
        if not self.keep_short_isolated_nodes:
            # Prune short isolated nodes using global cutoff
            shared.remove_short_isolated_nodes(self.graph, config.short_offset)
            
    def _merge_nodes(self):
        """Run neighbor merging stage unless explicitly disabled.

        This operation may simplify bubbles and short tips by fusing adjacent
        nodes according to local rules implemented in NeighbourMerger.
        """
        if not self.disable_merge_neighbor:
            self.neighbour_merger.merge_neibour()
            
    def _output_gfa(self, gfa_name):
        """Write the current graph into a GFA file.

        Parameters
        ----------
        gfa_name : str
            Output GFA file path.
        """
        gfaLib.GraphFromNetwork(self.graph).save_graph(gfa_name)

    def _create_unfold_dir(self):
        if not os.path.exists(self.out_path):
            os.makedirs(self.out_path)
