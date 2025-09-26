# GMW (Genomic Microbe-Wise)

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

---
## Overview
**GMW** is a hybrid assembly and contamination removal method for metagenomic sequencing data. It combines the strengths of **de novo** and **reference-based** approaches, and applies **community detection algorithms** on the de Bruijn graph to cluster contigs. This enables effective discrimination between target pathogen genomes and contaminant sequences, significantly improving assembly accuracy and downstream analyses.

---
### Requirements
- Python >= 3.8  
- [Python networkx package](https://networkx.org/)  
- [Python community package](https://biopython.org/)  
- [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi) 
- [Kraken2](https://ccb.jhu.edu/software/kraken2/)  

---
## Installation

### Steps
1. Clone the repository:
```bash
git clone https://github.com/trainrun/gmw.git
```
2. Create and activate conda environment:
```bash
conda create -n gmw python=3.13
conda activate gmw
```
3. Install dependencies:
```bash
conda install -c conda-forge networkx
conda install -c conda-forge python-louvain
conda install -c conda-forge click
conda install -c bioconda blast
conda install -c bioconda kraken2
```
---
## Usage

### Quick start
`python src/gmw/cli.py -g GFA_input_file [OPTIONS]`

### Arguments
- `-v, --version`                Show the version and exit.
- `-g, --gfa` `FILE`               Input gfa file.  [required]
- `-o, --outdir` `TEXT`            Path of output files.
- `-p, --prefix` `TEXT`            Prefix of output dir and files.
- `-t, --threads` `INTEGER`        Number of threads.
- `-f, --force`                  Force overwrite existed files
- `--disable_taxon_unfold`       Do not unfold graph using contigs taxonomy.
- `--use_gfa_taxon`              Parse contig type from the gfa file.
- `--kraken_out` `TEXT`            Use kraken output file rather than run kraken in this pipeline.
- `--taxon_db` `TEXT`              Dir of taxonomy databse. the dir must contain files "names.dmp" and "nodes.dmp".
- `--taxon_id` `TEXT`              The taxonomy id of species you want.
- `--taxon_name` `TEXT`            The scientific name of species you want.
- `--kraken_db` `TEXT`             The dir of kraken database.
- `--bgll`                       Use bgll algorithm infer taxon.
- `--disable_ref_unfold`         Do not unfold graph through reference search.
- `--use_gfa_ref`                Parse contig position from the gfa file.
- `--blast_out` `TEXT`             Use blast output file rather than run blastn in this pipeline.
- `--blast_db` `TEXT`              The path to blastn databse.
- `--position_distance` `INTEGER`  Remove the connection if two connected contigs if distance futher than offset.
- `--disable_depth_unfold`       Do not unfold graph using contigs sequencing depth.
- `--depth_discrepancy` `INTEGER`  Remove the connection if two connected contigs depth multiple higher than offset.
- `--disable_gc_unfold`          Do not unfold graph using through GC content and base percentage.
- `--gc_discrepancy` `FLOAT`       Remove the connection if two connected contigs gc content discrepancy higher than offset.
- `--remove_unknown_nodes`       Remove unknown components after unfold graph.
- `--keep_unknown_components`    Don't remove unknown components after unfold graph.
- `--keep_short_isolated_nodes`  Don't remove short isolated nodes after unfold graph.
- `--disable_merge_neighbour`    Don't merge single pair neibour nodes of the graph.
- `--merge_brother`              Merge brother nodes into consensus contigs.
- `--split_parent`               Split one node into two.
- `--fast`                       Fast mode. without irretaion.
- `--visual`                     Visualize debruijn graph.
- `--contig_shape` `TEXT`          Contig shape in debruijn graph, you can choose 'line' or 'dot'.(default line)
- `--help`                       Show this message and exit.

### Output
- `{output_dir}/DepthUnfolder`  Temporary files generated in Step "depth unfold".
- `{output_dir}/GCUnfolder`  Temporary files generated in Step "gc unfold".
- `{output_dir}/RefUnfolder`  Temporary files generated in Step "ref unfold".
- `{output_dir}/TaxonUnfolder`  Temporary files generated in Step "taxon unfold".
- `{output_dir}/Polisher`  Temporary files generated in Step "polish".
- `{output_dir}/{prefix}_before_unfold.html` Visualize file for input GFA file.
- `{output_dir}/{prefix}_after_unfold.html`  Visualize file for output GFA file.
- `{output_dir}/{prefix}_after_unfold.fasta` Final Fasta format output file.
- `{output_dir}/{prefix}_after_unfold.gfa` Final GFA format output file.
- `{output_dir}/{prefix}.log`  Log file.


### Examples

**Influenza virus mixed infections analysis**
```bash
python src/gmw/cli.py -g ./examples/input.gfa -o ./examples/out_files --disable_taxon_unfold --blast_db ./examples/ref/merge
```

