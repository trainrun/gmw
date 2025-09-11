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

### Prerequisites
- [Conda](https://docs.conda.io/en/latest/) installed
- [Git](https://git-scm.com/) installed


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
conda install -c bioconda blast
conda install -c bioconda kraken2
```
---
## Usage

### Quick start


### Command line interface

**Examples**
