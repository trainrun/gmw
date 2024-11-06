"""
P
"""

import click
import logging
import os
import sys

import gfaLib
import unfoldGraph
import mergeNodes
import config
import shared
import visualize
logger = logging.getLogger("gmw")

@click.command()
@click.version_option(config.VERSION, "--version", "-v")
@click.option("--gfa", "-g", required=True, type=click.Path(exists=True, dir_okay=False), help="Input gfa file.")
@click.option("--outdir", "-o", default=config.output_path, help="Path of output files.")
@click.option("--prefix", "-p", default=config.prefix, help="Prefix of output dir and files.")
@click.option("--threads", "-t", default=config.max_threads, help="Number of threads.")
@click.option("--force", "-f", is_flag=True, help="Force overwrite existed files")

@click.option("--disable_taxon_unfold", is_flag=True, help="Do not unfold graph using contigs taxonomy.")
@click.option("--use_gfa_taxon", is_flag=True, help="Parse contig type from the gfa file.")
@click.option("--kraken_out", help="Use kraken output file rather than run kraken in this pipeline.")
@click.option("--taxon_db", help="Dir of taxonomy databse. the dir must contain files \"names.dmp\" and \"nodes.dmp\".")
@click.option("--taxon_id", help="The taxonomy id of species you want.")
@click.option("--taxon_name", help="The scientific name of species you want.")
@click.option("--kraken_db", help="The dir of kraken database.")

@click.option("--disable_ref_unfold", is_flag=True, help="Do not unfold graph through reference search.")
@click.option("--use_gfa_ref", is_flag=True, help="Parse contig position from the gfa file.")
@click.option("--blast_out", help="Use blast output file rather than run blastn in this pipeline.")
@click.option("--blast_db", help="The path to blastn databse.")
@click.option("--position_distance", default=config.position_distance, type=int, 
              help="Remove the connection if two connected contigs if distance futher than offset.")

@click.option("--disable_depth_unfold", is_flag=True, help="Do not unfold graph using contigs sequencing depth.")
@click.option("--depth_discrepancy", default=config.depth_discrepancy, type=int, 
              help="Remove the connection if two connected contigs depth multiple higher than offset.")

@click.option("--disable_gc_unfold", is_flag=True, help="Do not unfold graph using through GC content and base percentage.")
@click.option("--gc_discrepancy", default=config.gc_discrepancy, type=float, 
              help="Remove the connection if two connected contigs gc content discrepancy higher than offset.")

@click.option("--remove_unknown_nodes", is_flag=True, help="Remove unknown components after unfold graph.")
@click.option("--keep_unknown_components", is_flag=True, help="Don't remove unknown components after unfold graph.")
@click.option("--keep_short_isolated_nodes", is_flag=True, help="Don't remove short isolated nodes after unfold graph.")

@click.option("--disable_merge_neighbour", is_flag=True, help="Don't merge single pair neibour nodes of the graph.")
@click.option("--merge_brother", is_flag=True, help="Merge brother nodes into consensus contigs.")
@click.option("--split_parent", is_flag=True, help="Split one node into two.")

@click.option("--visual", is_flag=True, help="Visualize debruijn graph.")
@click.option("--contig_shape", default="line", help="Contig shape in debruijn graph, you can choose 'line' or 'dot'.(default line)")
def cli(
    gfa, outdir, prefix, threads, force,
    disable_taxon_unfold, use_gfa_taxon, kraken_out, taxon_db, taxon_id, taxon_name, kraken_db, 
    disable_ref_unfold, use_gfa_ref, blast_out, blast_db, position_distance,
    disable_depth_unfold, depth_discrepancy,
    disable_gc_unfold, gc_discrepancy,
    remove_unknown_nodes, keep_unknown_components, keep_short_isolated_nodes,
    disable_merge_neighbour, merge_brother, split_parent,
    visual, contig_shape
):
    """a tool for designing primer panels for multiplex PCR."""
    
    """Check outdir."""
    parent_dir = os.path.dirname(outdir)
    if os.path.exists(parent_dir) and os.path.isdir(parent_dir):
        outdir = os.path.normpath(outdir)
        if not os.path.exists(outdir):
            try:
                os.makedirs(outdir)
            except OSError as e:
                raise ValueError(f"Create dir filed: {e}")
        elif not force:
            click.echo(cli.get_help(ctx=click.Context(cli)))
            print(f"  Error! The dir name \"{outdir}\" has exist, please use --force parameter to overwrite.")
            exit(1)
        else:
            pass
    else:
        click.echo(cli.get_help(ctx=click.Context(cli)))
        print(f"  Error! Dir not exist \"{parent_dir}\".")
        exit(1)

    "Setup log information."
    setup_logging(outdir,prefix)
    logger.info(f"Output dir: {outdir}")
    
    unfold_argv = []
    unfold_argv.extend([outdir, prefix, threads, force])
    names_dmp = None
    nodes_dmp = None
    
    
    """Add disabled unfold parameters"""
    
    unfold_argv.extend([disable_taxon_unfold, disable_ref_unfold, disable_depth_unfold, disable_gc_unfold])
    """Check taxonomy unfold parameters"""
    if not disable_taxon_unfold:
        if use_gfa_taxon:
            logger.info(f"Use parameters to infer contig type from gfa file.") 
        else:
            if taxon_db == None:
                logger.error(f"Please add the taxonomy database using parameter \"--taxon_db\"")
                exit(1)
            elif os.path.exists(taxon_db) and os.path.isdir(taxon_db):
                taxon_db = os.path.normpath(taxon_db)
                names_dmp = taxon_db + "/names.dmp"
                nodes_dmp = taxon_db + "/nodes.dmp"
                if os.path.exists(names_dmp) and os.path.exists(nodes_dmp):
                    logger.info(f"Taxonomy database dir: \"{taxon_db}\".")
                else:
                    logger.error(f"File \"names.dmp\" or \"nodes.dmp\" is not in dir \"{taxon_db}\".")
                    exit(1)
            else:
                logger.error(f"Taxonomy dir \"{taxon_db}\" not exist.")
                exit(1)
            
            if taxon_id is None and taxon_name is None:
                logger.error(f"You must specify taxon_id or taxon_name!")
                exit(1)
            if kraken_out != None:
                if os.path.exists(kraken_out):
                    logger.info(f"Use kraken output file \"{kraken_out}.\"")
                else:
                    logger.error(f"Kraken output file \"{kraken_out}.\" not exist.")
                    exit(1)            
            else:    
                if kraken_db is None:
                    logger.error(f"You must specify kraken!")
                    exit(1)
                elif os.path.exists(kraken_db) and os.path.isdir(kraken_db):
                    logger.info(f"Kraken2 database dir: \"{kraken_db}\".")
                else:
                    logger.error(f"Dir \"{taxon_db}\" is not exist.")
                    exit(1)  
    else:
        logger.info("Taxon unfold function is disabled")
        
    unfold_argv.extend([use_gfa_taxon, kraken_out, names_dmp, nodes_dmp, taxon_id, taxon_name, kraken_db])

    """Check reference unfold parameters"""
    if not disable_ref_unfold:
        if use_gfa_ref:
            logger.info(f"Use parameters to parse contig position from gfa file.") 
        elif blast_out is not None:
            if os.path.exists(blast_out):
                logger.info(f"Use blastn output file \"{blast_out}.\"")
            else:
                logger.error(f"Blastn output file \"{blast_out}.\" not exist.")
                exit(1)
        else:
            if blast_db is None:
                logger.error(f"Please add the blastn database using parameter \"--blast_db\"")
                exit(1)
            elif os.path.exists(blast_db + ".nto"):
                logger.info(f"Blastn database dir: \"{blast_db}\".")
            else:
                logger.error(f"The path of blastn database \"{blast_db}\" if not correct.")
                exit(1)   

    unfold_argv.extend([use_gfa_ref, blast_out, blast_db, position_distance])
    
    """Check depth unfold parameters"""
    unfold_argv.append(depth_discrepancy)
   
    """Check GC unfold parameters"""
    unfold_argv.append(gc_discrepancy)
        
    """Check remove nodes parameters"""
    unfold_argv.extend([remove_unknown_nodes, keep_unknown_components, keep_short_isolated_nodes])
 
    """Check merge nodes parameters"""
    unfold_argv.extend([disable_merge_neighbour, merge_brother, split_parent])
    
    """Check visualize parameters"""
    if visual and contig_shape not in ['dot', 'line']:
        logger.error(f"The contig_shape must be 'dot' or 'line'!")
        exit(1)   
    unfold_argv.extend([visual, contig_shape])
    
    """Add gfa paramerter"""
    gfa_graph = gfaLib.GraphFromFile(gfa)
    graph = gfaLib.GFANetwork.compute_backbone(gfa_graph)
    unfold_argv.insert(0, graph)

    before_fig_path = outdir + "/" + prefix + "_before_unfold.html"
    after_fig_path = outdir + "/" + prefix + "_after_unfold.html"
    fasta_path = outdir + "/" + prefix + "_after_unfold.fasta"
    
    if visual:
        visualize.print_graph(graph, before_fig_path, contig_shape)

    logging_graph_info(graph)

    nodes_num = None
    neighbour_merger = mergeNodes.NeighbourMerger(graph)
    while nodes_num is None or nodes_num != graph.number_of_nodes():
        nodes_num = graph.number_of_nodes()
        neighbour_merger.merge_neibour()
        logging_graph_info(graph)


    if not disable_taxon_unfold:    
        taxonUnfoler = unfoldGraph.TaxonUnfolder(*unfold_argv)
        taxonUnfoler.unfold_graph()
        logging_graph_info(graph)
    
    if not disable_ref_unfold:
        refUnfoler = unfoldGraph.RefUnfolder(*unfold_argv)
        refUnfoler.unfold_graph()
        logging_graph_info(graph)
        
    if not disable_depth_unfold:
        depthUnfoler = unfoldGraph.DepthUnfolder(*unfold_argv)
        depthUnfoler.unfold_graph()
        logging_graph_info(graph)
        
    if not disable_gc_unfold:       
        gcUnfoler = unfoldGraph.GCUnfolder(*unfold_argv)
        gcUnfoler.unfold_graph()
        logging_graph_info(graph)
        
    if visual:
        visualize.print_graph(graph, after_fig_path, contig_shape)       
        
    shared.graph2fasta(graph, fasta_path, orientation=False)
    
def logging_graph_info(graph):
       nodes_num = graph.number_of_nodes()
       edges_num = graph.number_of_edges()
       logger.info(f"The graph have {nodes_num} nodes and {edges_num} edges.")
       
def setup_logging(output_path, prefix):
    logger.setLevel(logging.INFO)

    console_handler = logging.StreamHandler()
    log_filepath = output_path + f"/{prefix}.log"  
    file_handler = logging.FileHandler(log_filepath)  
    
    console_handler.setLevel(logging.DEBUG)  
    file_handler.setLevel(logging.INFO)  
    
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')  
    console_handler.setFormatter(formatter)  
    file_handler.setFormatter(formatter)  
    
    logger.addHandler(console_handler)  
    logger.addHandler(file_handler)  

    logger.info(f"GMW version: {config.VERSION} ")
    logger.info(f"Writing log to {log_filepath}")

if __name__ == "__main__":
    cli()
