# GenoPhase

## A tool for haplotype-resolved genome assembly using HiFi reads

GenoPhase is a comprehensive bioinformatics pipeline for constructing high-quality haplotype-resolved genome assemblies using PacBio HiFi sequencing data. The tool implements a novel graph-based phasing algorithm that leverages the high accuracy and long read lengths of HiFi data to separate reads into haplotype-specific sets, enabling the assembly of both haplotypes with high completeness and accuracy.

## Features

- **Temporary reference construction** using Hifiasm in primary mode
- **Graph-based SNP phasing** that exploits HiFi reads' accuracy and length
- **Local and global phasing** to generate consistent haplotype assignments
- **Read classification** into haplotype-specific sets
- **Haplotype-resolved assembly** using separated read sets

## Installation

### Prerequisites

- Python 3.7 or higher
- Biopython
- NumPy
- NetworkX
- Pandas
- SciPy
- pysam
- minimap2
- samtools
- hifiasm

### Clone the repository

```bash
git clone https://github.com/username/GenoPhase.git
cd GenoPhase
