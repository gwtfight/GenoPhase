#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
GenoPhase: A tool for haplotype-resolved genome assembly using HiFi reads

This program implements a graph-based phasing algorithm to separate HiFi reads
into haplotype-specific sets and constructs haplotype-resolved assemblies.

Author: Northeast Forestry University 
Date: January , 2025
"""

import os
import sys
import argparse
import logging
import subprocess
import numpy as np
import networkx as nx
from collections import defaultdict, Counter
import pysam
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from dataclasses import dataclass
from typing import List, Dict, Tuple, Set, Optional
import multiprocessing as mp
from Bio import SeqIO

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger('GenoPhase')

@dataclass
class SNPNode:
    """Class representing a SNP node in the phasing graph"""
    chrom: str
    pos: int
    ref_base: str
    alt_base: str
    ref_count: int
    alt_count: int
    quality: float
    reads_spanning: Set[str]
    
    @property
    def id(self) -> str:
        """Generate a unique identifier for this SNP"""
        return f"{self.chrom}_{self.pos}_{self.ref_base}_{self.alt_base}"
    
    @property
    def total_depth(self) -> int:
        """Total read depth at this position"""
        return self.ref_count + self.alt_count
    
    @property
    def minor_allele_freq(self) -> float:
        """Minor allele frequency at this position"""
        if self.total_depth == 0:
            return 0
        return min(self.ref_count, self.alt_count) / self.total_depth
    
    def __hash__(self):
        return hash(self.id)
    
    def __eq__(self, other):
        if not isinstance(other, SNPNode):
            return False
        return self.id == other.id

@dataclass
class PhaseEdge:
    """Class representing an edge between two SNP nodes"""
    node1: SNPNode
    node2: SNPNode
    shared_reads: Set[str]
    cis_count: int  # ref-ref + alt-alt
    trans_count: int  # ref-alt + alt-ref
    quality: float
    
    @property
    def weight(self) -> float:
        """Calculate edge weight based on shared reads and phasing certainty"""
        phasing_certainty = abs(self.cis_count - self.trans_count) / max(1, (self.cis_count + self.trans_count))
        return len(self.shared_reads) * phasing_certainty * self.quality
    
    @property
    def phase_relationship(self) -> str:
        """Determine if the SNPs are in the same phase or different phases"""
        if self.cis_count >= self.trans_count:
            return "same"
        else:
            return "diff"

class PhasingBlock:
    """Class representing a block of phased SNPs"""
    def __init__(self, id: str):
        self.id = id
        self.snps = []
        self.phase = None  # 0 or 1 in global context
    
    def add_snp(self, snp: SNPNode, local_phase: int):
        """Add a SNP to this phasing block with its local phase"""
        self.snps.append((snp, local_phase))
    
    def flip_phase(self):
        """Invert the phase of all SNPs in this block"""
        self.snps = [(snp, 1-phase) for snp, phase in self.snps]
        if self.phase is not None:
            self.phase = 1 - self.phase
    
    @property
    def size(self) -> int:
        """Number of SNPs in this block"""
        return len(self.snps)
    
    @property
    def span(self) -> int:
        """Genomic span of this block"""
        if not self.snps:
            return 0
        positions = [snp.pos for snp, _ in self.snps]
        return max(positions) - min(positions)
    
    @property
    def chromosome(self) -> str:
        """Chromosome of this block"""
        if not self.snps:
            return ""
        return self.snps[0][0].chrom

class GenoPhase:
    """Main class for haplotype phasing and assembly"""
    
    def __init__(self, args):
        self.args = args
        self.snp_nodes = {}  # Dictionary of SNP nodes indexed by ID
        self.phasing_graph = nx.Graph()
        self.phasing_blocks = []
        self.global_phasing_graph = nx.Graph()
        
        # Create output directories
        os.makedirs(args.output_dir, exist_ok=True)
        os.makedirs(os.path.join(args.output_dir, "tmp"), exist_ok=True)
        os.makedirs(os.path.join(args.output_dir, "haplotypes"), exist_ok=True)
        os.makedirs(os.path.join(args.output_dir, "assemblies"), exist_ok=True)
        
        # Configure logging to file
        log_file = os.path.join(args.output_dir, "genophase.log")
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
        logger.addHandler(file_handler)
        
    def run(self):
        """Execute the full GenoPhase workflow"""
        logger.info("Starting GenoPhase pipeline")
        
        if self.args.mode == "full" or self.args.mode == "reference":
            self.build_temp_reference()
        
        if self.args.mode == "full" or self.args.mode == "phase":
            self.build_phasing_graph()
            self.local_phasing()
            self.global_phasing()
            self.classify_reads()
        
        if self.args.mode == "full" or self.args.mode == "assemble":
            self.assemble_haplotypes()
        
        logger.info("GenoPhase pipeline completed successfully")
    
    def build_temp_reference(self):
        """Build a temporary reference genome using Hifiasm"""
        logger.info("Building temporary reference genome with Hifiasm")
        
        # Create output directory
        ref_dir = os.path.join(self.args.output_dir, "reference")
        os.makedirs(ref_dir, exist_ok=True)
        
        # Run Hifiasm in primary contig mode
        hifiasm_cmd = [
            "hifiasm", 
            "-o", os.path.join(ref_dir, "temp_ref"),
            "-t", str(self.args.threads),
            "--primary",
            self.args.hifi_reads
        ]
        
        logger.info(f"Running command: {' '.join(hifiasm_cmd)}")
        
        try:
            subprocess.run(hifiasm_cmd, check=True)
            
            # Convert GFA to FASTA
            gfa_to_fasta_cmd = [
                "awk", 
                "/^S/{print \">\"$2\"\\n\"$3}", 
                os.path.join(ref_dir, "temp_ref.p_ctg.gfa")
            ]
            
            with open(os.path.join(ref_dir, "temp_ref.fasta"), "w") as fasta_file:
                subprocess.run(gfa_to_fasta_cmd, stdout=fasta_file, check=True)
            
            # Index the reference
            samtools_faidx_cmd = ["samtools", "faidx", os.path.join(ref_dir, "temp_ref.fasta")]
            subprocess.run(samtools_faidx_cmd, check=True)
            
            self.args.reference = os.path.join(ref_dir, "temp_ref.fasta")
            logger.info(f"Temporary reference genome built successfully: {self.args.reference}")
            
        except subprocess.CalledProcessError as e:
            logger.error(f"Error building temporary reference: {e}")
            raise
    
    def map_reads_to_reference(self):
        """Map HiFi reads to the reference genome using minimap2"""
        logger.info("Mapping HiFi reads to reference genome")
        
        bam_file = os.path.join(self.args.output_dir, "mapped_reads.bam")
        
        # If BAM file exists and we're not forcing a remap, skip this step
        if os.path.exists(bam_file) and not self.args.force_remap:
            logger.info(f"Using existing BAM file: {bam_file}")
            return bam_file
        
        # Map reads with minimap2
        minimap2_cmd = [
            "minimap2",
            "-ax", "map-hifi",
            "-t", str(self.args.threads),
            self.args.reference,
            self.args.hifi_reads
        ]
        
        samtools_cmd = [
            "samtools", "sort",
            "-@", str(self.args.threads),
            "-o", bam_file
        ]
        
        try:
            # Pipe minimap2 output to samtools
            logger.info(f"Running command: {' '.join(minimap2_cmd)} | {' '.join(samtools_cmd)}")
            
            p1 = subprocess.Popen(minimap2_cmd, stdout=subprocess.PIPE)
            p2 = subprocess.Popen(samtools_cmd, stdin=p1.stdout)
            p1.stdout.close()  # Allow p1 to receive a SIGPIPE if p2 exits
            p2.communicate()
            
            if p1.returncode != 0 or p2.returncode != 0:
                raise subprocess.CalledProcessError(p1.returncode or p2.returncode, minimap2_cmd)
            
            # Index the BAM file
            samtools_index_cmd = ["samtools", "index", bam_file]
            subprocess.run(samtools_index_cmd, check=True)
            
            logger.info(f"Read mapping completed successfully: {bam_file}")
            return bam_file
            
        except subprocess.CalledProcessError as e:
            logger.error(f"Error mapping reads to reference: {e}")
            raise
    
    def detect_snps_from_bam(self, bam_file):
        """Extract SNP information directly from the BAM file"""
        logger.info("Detecting SNPs from BAM file")
        
        min_quality = self.args.min_base_quality
        min_depth = self.args.min_depth
        max_depth = self.args.max_depth
        min_maf = self.args.min_minor_allele_freq
        
        snps = []
        
        try:
            bam = pysam.AlignmentFile(bam_file, "rb")
            
            # Get list of references/chromosomes
            references = bam.references
            
            for ref in references:
                logger.info(f"Processing reference: {ref}")
                
                # Iterate through pileups
                for pileup_column in bam.pileup(ref, stepper="samtools", min_base_quality=min_quality):
                    pos = pileup_column.pos
                    depth = pileup_column.n
                    
                    # Skip if depth is out of range
                    if depth < min_depth or depth > max_depth:
                        continue
                    
                    # Count bases at this position
                    bases = Counter()
                    reads_by_base = defaultdict(set)
                    
                    for pileup_read in pileup_column.pileups:
                        # Skip if not a base (deletion or reference skip)
                        if pileup_read.is_del or pileup_read.is_refskip:
                            continue
                        
                        read_name = pileup_read.alignment.query_name
                        base = pileup_read.alignment.query_sequence[pileup_read.query_position].upper()
                        
                        # Only count if base quality is sufficient
                        base_quality = pileup_read.alignment.query_qualities[pileup_read.query_position]
                        if base_quality >= min_quality:
                            bases[base] += 1
                            reads_by_base[base].add(read_name)
                    
                    # Need at least two alleles to make a SNP
                    if len(bases) < 2:
                        continue
                    
                    # Find the two most frequent bases
                    most_common = bases.most_common(2)
                    if len(most_common) < 2:
                        continue
                    
                    (ref_base, ref_count), (alt_base, alt_count) = most_common
                    
                    # Check if minor allele frequency is sufficient
                    minor_count = min(ref_count, alt_count)
                    if minor_count / depth < min_maf:
                        continue
                    
                    # Create a SNP node
                    quality = 1.0  # Could be refined based on base qualities
                    reads_spanning = reads_by_base[ref_base].union(reads_by_base[alt_base])
                    
                    snp = SNPNode(
                        chrom=ref,
                        pos=pos,
                        ref_base=ref_base,
                        alt_base=alt_base,
                        ref_count=ref_count,
                        alt_count=alt_count,
                        quality=quality,
                        reads_spanning=reads_spanning
                    )
                    
                    snps.append(snp)
                    self.snp_nodes[snp.id] = snp
            
            logger.info(f"Detected {len(snps)} candidate SNPs")
            return snps
            
        except Exception as e:
            logger.error(f"Error detecting SNPs from BAM: {e}")
            raise
    
    def build_phasing_graph(self):
        """Build the phasing graph from detected SNPs"""
        logger.info("Building phasing graph")
        
        # Map reads to reference if not already done
        bam_file = self.map_reads_to_reference()
        
        # Detect SNPs from BAM file
        snps = self.detect_snps_from_bam(bam_file)
        
        # Add nodes to the graph
        for snp in snps:
            self.phasing_graph.add_node(snp.id, snp=snp)
        
        logger.info(f"Added {len(snps)} nodes to phasing graph")
        
        # Group SNPs by chromosome for efficient edge creation
        snps_by_chrom = defaultdict(list)
        for snp in snps:
            snps_by_chrom[snp.chrom].append(snp)
        
        # Sort SNPs by position within each chromosome
        for chrom in snps_by_chrom:
            snps_by_chrom[chrom].sort(key=lambda x: x.pos)
        
        # Create edges between SNPs
        total_edges = 0
        
        for chrom, chrom_snps in snps_by_chrom.items():
            logger.info(f"Creating edges for chromosome {chrom}")
            
            # Use sliding window approach to limit comparisons
            for i, snp1 in enumerate(chrom_snps):
                # Only look at SNPs within D_max distance
                j = i + 1
                while j < len(chrom_snps) and chrom_snps[j].pos - snp1.pos <= self.args.max_edge_distance:
                    snp2 = chrom_snps[j]
                    
                    # Find shared reads
                    shared_reads = snp1.reads_spanning.intersection(snp2.reads_spanning)
                    
                    # Only create edge if there are enough shared reads
                    if len(shared_reads) >= self.args.min_shared_reads:
                        # Analyze phasing information
                        cis_count = 0  # ref-ref + alt-alt
                        trans_count = 0  # ref-alt + alt-ref
                        
                        # Extract read information from BAM file
                        try:
                            bam = pysam.AlignmentFile(bam_file, "rb")
                            
                            for read_name in shared_reads:
                                # Get all reads with this name
                                reads = [read for read in bam.fetch() if read.query_name == read_name]
                                
                                for read in reads:
                                    # Check if this read spans both SNPs
                                    if read.reference_start <= snp1.pos and read.reference_end >= snp2.pos:
                                        # Get bases at SNP positions
                                        base1 = None
                                        base2 = None
                                        
                                        for pair in read.get_aligned_pairs(with_seq=True):
                                            if pair[1] == snp1.pos:
                                                if pair[0] is not None:
                                                    base1 = read.query_sequence[pair[0]].upper()
                                            if pair[1] == snp2.pos:
                                                if pair[0] is not None:
                                                    base2 = read.query_sequence[pair[0]].upper()
                                        
                                        if base1 and base2:
                                            # Count phasing information
                                            if (base1 == snp1.ref_base and base2 == snp2.ref_base) or \
                                               (base1 == snp1.alt_base and base2 == snp2.alt_base):
                                                cis_count += 1
                                            elif (base1 == snp1.ref_base and base2 == snp2.alt_base) or \
                                                 (base1 == snp1.alt_base and base2 == snp2.ref_base):
                                                trans_count += 1
                                        
                                        # Only need to process one read per read name
                                        break
                            
                        except Exception as e:
                            logger.error(f"Error analyzing phasing information: {e}")
                            continue
                        
                        # Calculate edge quality (could be refined)
                        quality = 1.0
                        
                        # Create edge if we have phasing information
                        if cis_count + trans_count > 0:
                            edge = PhaseEdge(
                                node1=snp1,
                                node2=snp2,
                                shared_reads=shared_reads,
                                cis_count=cis_count,
                                trans_count=trans_count,
                                quality=quality
                            )
                            
                            self.phasing_graph.add_edge(
                                snp1.id, 
                                snp2.id, 
                                weight=edge.weight,
                                phase=edge.phase_relationship,
                                edge=edge
                            )
                            
                            total_edges += 1
                    
                    j += 1
        
        logger.info(f"Added {total_edges} edges to phasing graph")
        
        # Save graph structure for debugging/visualization
        graph_output = os.path.join(self.args.output_dir, "phasing_graph.gml")
        nx.write_gml(self.phasing_graph, graph_output)
        logger.info(f"Phasing graph saved to {graph_output}")
    
    def local_phasing(self):
        """Perform local phasing on connected components"""
        logger.info("Performing local phasing")
        
        # Identify connected components in the graph
        connected_components = list(nx.connected_components(self.phasing_graph))
        logger.info(f"Found {len(connected_components)} connected components")
        
        # Process each connected component
        for i, component in enumerate(connected_components):
            if i % 100 == 0:
                logger.info(f"Processing component {i+1}/{len(connected_components)}")
            
            # Skip tiny components (single nodes)
            if len(component) <= 1:
                continue
            
            # Create subgraph for this component
            subgraph = self.phasing_graph.subgraph(component).copy()
            
            # Find anchor node (highest degree or weight)
            anchor_node = max(subgraph.nodes(), key=lambda n: subgraph.degree(n))
            
            # Initialize phases
            phases = {anchor_node: 0}
            
            # Use BFS to propagate phases
            queue = [anchor_node]
            while queue:
                current = queue.pop(0)
                current_phase = phases[current]
                
                # Process neighbors
                for neighbor in subgraph.neighbors(current):
                    if neighbor not in phases:
                        # Determine phase based on edge relationship
                        edge_phase = subgraph[current][neighbor]['phase']
                        if edge_phase == "same":
                            neighbor_phase = current_phase
                        else:  # edge_phase == "diff"
                            neighbor_phase = 1 - current_phase
                        
                        phases[neighbor] = neighbor_phase
                        queue.append(neighbor)
                    else:
                        # Check for consistency (could be improved with conflict resolution)
                        edge_phase = subgraph[current][neighbor]['phase']
                        expected_phase = current_phase if edge_phase == "same" else 1 - current_phase
                        
                        if phases[neighbor] != expected_phase:
                            # There's a conflict - could implement more sophisticated resolution
                            logger.debug(f"Phase conflict detected in component {i}")
            
            # Create a phasing block for this component
            block = PhasingBlock(f"block_{i}")
            
            # Add SNPs to the block
            for node_id in subgraph.nodes():
                snp = self.snp_nodes[node_id]
                phase = phases.get(node_id, 0)  # Default to 0 if not phased
                block.add_snp(snp, phase)
            
            self.phasing_blocks.append(block)
        
        logger.info(f"Created {len(self.phasing_blocks)} phasing blocks")
        
        # Write phasing blocks to file
        self.write_phasing_blocks()
    
    def global_phasing(self):
        """Perform global phasing to connect phasing blocks"""
        logger.info("Performing global phasing")
        
        # Create a graph where nodes are phasing blocks
        for i, block in enumerate(self.phasing_blocks):
            self.global_phasing_graph.add_node(block.id, block=block)
        
        # Find connections between blocks using HiFi reads
        bam_file = os.path.join(self.args.output_dir, "mapped_reads.bam")
        
        try:
            # Build a dictionary of which blocks each read spans
            reads_to_blocks = defaultdict(set)
            
            for block in self.phasing_blocks:
                for snp, _ in block.snps:
                    for read_name in snp.reads_spanning:
                        reads_to_blocks[read_name].add(block.id)
            
            # Create edges between blocks that share reads
            for read_name, block_ids in reads_to_blocks.items():
                if len(block_ids) < 2:
                    continue
                
                # Create edges between all pairs of blocks
                for block_id1 in block_ids:
                    for block_id2 in block_ids:
                        if block_id1 >= block_id2:
                            continue
                        
                        # Check if edge already exists
                        if self.global_phasing_graph.has_edge(block_id1, block_id2):
                            # Add read to shared reads
                            self.global_phasing_graph[block_id1][block_id2]['shared_reads'].add(read_name)
                        else:
                            # Create new edge
                            self.global_phasing_graph.add_edge(
                                block_id1, 
                                block_id2, 
                                shared_reads={read_name},
                                weight=0,
                                phase=None
                            )
            
            logger.info(f"Created {self.global_phasing_graph.number_of_edges()} edges between phasing blocks")
            
            # Determine phase relationships between blocks
            bam = pysam.AlignmentFile(bam_file, "rb")
            
            for block_id1, block_id2 in self.global_phasing_graph.edges():
                shared_reads = self.global_phasing_graph[block_id1][block_id2]['shared_reads']
                
                block1 = next(block for block in self.phasing_blocks if block.id == block_id1)
                block2 = next(block for block in self.phasing_blocks if block.id == block_id2)
                
                # Count phase concordance
                same_phase_count = 0
                diff_phase_count = 0
                
                for read_name in shared_reads:
                    # Get all reads with this name
                    reads = [read for read in bam.fetch() if read.query_name == read_name]
                    
                    for read in reads:
                        # Check concordance within each block
                        block1_concordance = self.check_read_block_concordance(read, block1)
                        block2_concordance = self.check_read_block_concordance(read, block2)
                        
                        # Only count if we have clear concordance in both blocks
                        if block1_concordance is not None and block2_concordance is not None:
                            if block1_concordance == block2_concordance:
                                same_phase_count += 1
                            else:
                                diff_phase_count += 1
                
                # Determine phase relationship
                if same_phase_count > diff_phase_count:
                    phase = "same"
                    weight = same_phase_count - diff_phase_count
                else:
                    phase = "diff"
                    weight = diff_phase_count - same_phase_count
                
                # Update edge
                self.global_phasing_graph[block_id1][block_id2]['phase'] = phase
                self.global_phasing_graph[block_id1][block_id2]['weight'] = weight
            
            # Use maximum spanning tree to connect blocks
            mst_edges = list(nx.maximum_spanning_tree(self.global_phasing_graph, weight='weight').edges())
            
            # Assign global phases starting from the largest block
            largest_block_id = max(self.phasing_blocks, key=lambda b: b.size).id
            
            # BFS to propagate global phases
            global_phases = {largest_block_id: 0}
            queue = [largest_block_id]
            
            while queue:
                current = queue.pop(0)
                current_phase = global_phases[current]
                
                for edge in mst_edges:
                    if current in edge:
                        other = edge[0] if edge[1] == current else edge[1]
                        
                        if other not in global_phases:
                            edge_phase = self.global_phasing_graph[edge[0]][edge[1]]['phase']
                            if edge_phase == "same":
                                other_phase = current_phase
                            else:  # edge_phase == "diff"
                                other_phase = 1 - current_phase
                            
                            global_phases[other] = other_phase
                            queue.append(other)
            
            # Update blocks with global phases
            for block in self.phasing_blocks:
                if block.id in global_phases:
                    block.phase = global_phases[block.id]
                else:
                    # Assign arbitrary phase to disconnected blocks
                    block.phase = 0
            
            logger.info("Global phasing completed successfully")
            
        except Exception as e:
            logger.error(f"Error in global phasing: {e}")
            raise
    
    def check_read_block_concordance(self, read, block):
        """Check if a read supports phase 0 or 1 for a given block"""
        phase0_support = 0
        phase1_support = 0
        
        for snp, local_phase in block.snps:
            # Check if read covers this SNP
            for pair in read.get_aligned_pairs(with_seq=True):
                if pair[1] == snp.pos and pair[0] is not None:
                    base = read.query_sequence[pair[0]].upper()
                    
                    # Check which phase this base supports
                    if base == snp.ref_base:
                        if local_phase == 0:
                            phase0_support += 1
                        else:
                            phase1_support += 1
                    elif base == snp.alt_base:
                        if local_phase == 0:
                            phase1_support += 1
                        else:
                            phase0_support += 1
                    
                    break
        
        # Need at least some evidence
        if phase0_support + phase1_support < 2:
            return None
        
        # Determine which phase is supported
        if phase0_support > phase1_support:
            return 0
        elif phase1_support > phase0_support:
            return 1
        else:
            return None  # Ambiguous
    
    def write_phasing_blocks(self):
        """Write phasing blocks to file"""
        output_file = os.path.join(self.args.output_dir, "phasing_blocks.tsv")
        
        with open(output_file, "w") as f:
            f.write("block_id\tchromosome\tstart_pos\tend_pos\tsize\tglobal_phase\tsnp_positions\n")
            
            for block in self.phasing_blocks:
                if not block.snps:
                    continue
                
                chrom = block.chromosome
                positions = [snp.pos for snp, _ in block.snps]
                start_pos = min(positions)
                end_pos = max(positions)
                size = len(positions)
                global_phase = block.phase if block.phase is not None else "NA"
                snp_positions = ",".join(map(str, positions))
                
                f.write(f"{block.id}\t{chrom}\t{start_pos}\t{end_pos}\t{size}\t{global_phase}\t{snp_positions}\n")
        
        logger.info(f"Phasing blocks written to {output_file}")
    
    def classify_reads(self):
        """Classify reads into haplotypes based on phasing results"""
        logger.info("Classifying reads into haplotypes")
        
        # Create a dictionary of SNP positions to their phased alleles
        snp_phasing = {}
        
        for block in self.phasing_blocks:
            global_phase = block.phase if block.phase is not None else 0
            
            for snp, local_phase in block.snps:
                # Determine the actual phase in the global context
                actual_phase = local_phase if global_phase == 0 else 1 - local_phase
                
                # Store reference and alternative alleles for each haplotype
                snp_phasing[(snp.chrom, snp.pos)] = {
                    "ref_base": snp.ref_base,
                    "alt_base": snp.alt_base,
                    "hap0_base": snp.ref_base if actual_phase == 0 else snp.alt_base,
                    "hap1_base": snp.alt_base if actual_phase == 0 else snp.ref_base
                }
        
        # Classify reads based on which haplotype they support
        bam_file = os.path.join(self.args.output_dir, "mapped_reads.bam")
        
        try:
            bam = pysam.AlignmentFile(bam_file, "rb")
            
            hap0_reads = set()
            hap1_reads = set()
            unclassified_reads = set()
            
            for read in bam.fetch():
                hap0_support = 0
                hap1_support = 0
                
                # Check each SNP position
                for pair in read.get_aligned_pairs(with_seq=True):
                    pos = pair[1]
                    if pos is not None and (read.reference_name, pos) in snp_phasing:
                        if pair[0] is not None:  # Read covers this position
                            base = read.query_sequence[pair[0]].upper()
                            snp_info = snp_phasing[(read.reference_name, pos)]
                            
                            if base == snp_info["hap0_base"]:
                                hap0_support += 1
                            elif base == snp_info["hap1_base"]:
                                hap1_support += 1
                
                # Classify based on support
                if hap0_support > hap1_support * 2:  # Strong hap0 support
                    hap0_reads.add(read.query_name)
                elif hap1_support > hap0_support * 2:  # Strong hap1 support
                    hap1_reads.add(read.query_name)
                else:
                    unclassified_reads.add(read.query_name)
            
            # Write read classifications
            with open(os.path.join(self.args.output_dir, "haplotype_reads.txt"), "w") as f:
                f.write(f"haplotype_A\t{len(hap0_reads)}\n")
                for read in sorted(hap0_reads):
                    f.write(f"A\t{read}\n")
                
                f.write(f"haplotype_B\t{len(hap1_reads)}\n")
                for read in sorted(hap1_reads):
                    f.write(f"B\t{read}\n")
                
                f.write(f"unclassified\t{len(unclassified_reads)}\n")
                for read in sorted(unclassified_reads):
                    f.write(f"U\t{read}\n")
            
            # Extract FASTA sequences for each haplotype
            self.extract_haplotype_reads(hap0_reads, hap1_reads)
            
            logger.info(f"Classified {len(hap0_reads)} reads as haplotype A, {len(hap1_reads)} as haplotype B, and {len(unclassified_reads)} as unclassified")
            
        except Exception as e:
            logger.error(f"Error classifying reads: {e}")
            raise
    
    def extract_haplotype_reads(self, hap0_reads, hap1_reads):
        """Extract FASTA sequences for classified reads"""
        logger.info("Extracting FASTA sequences for classified reads")
        
        # Create output files
        hap0_fasta = os.path.join(self.args.output_dir, "haplotypes", "haplotypeA.fasta")
        hap1_fasta = os.path.join(self.args.output_dir, "haplotypes", "haplotypeB.fasta")
        
        # Create indices for efficient lookup
        hap0_set = set(hap0_reads)
        hap1_set = set(hap1_reads)
        
        try:
            # Parse input FASTA file and write haplotype-specific files
            with open(hap0_fasta, "w") as f0, open(hap1_fasta, "w") as f1:
                for record in SeqIO.parse(self.args.hifi_reads, "fasta"):
                    if record.id in hap0_set:
                        SeqIO.write(record, f0, "fasta")
                    elif record.id in hap1_set:
                        SeqIO.write(record, f1, "fasta")
            
            logger.info(f"Haplotype reads extracted to {hap0_fasta} and {hap1_fasta}")
            
        except Exception as e:
            logger.error(f"Error extracting haplotype reads: {e}")
            raise
    
    def assemble_haplotypes(self):
        """Assemble haplotype-specific genomes using Hifiasm"""
        logger.info("Assembling haplotype-specific genomes")
        
        hap0_fasta = os.path.join(self.args.output_dir, "haplotypes", "haplotypeA.fasta")
        hap1_fasta = os.path.join(self.args.output_dir, "haplotypes", "haplotypeB.fasta")
        
        # Assemble haplotype A
        logger.info("Assembling haplotype A")
        
        hap0_prefix = os.path.join(self.args.output_dir, "assemblies", "Ptr_hapA.GenoPhase")
        
        hifiasm_cmd_hap0 = [
            "hifiasm", 
            "-o", hap0_prefix,
            "-t", str(self.args.threads),
            "--primary",
            hap0_fasta
        ]
        
        try:
            subprocess.run(hifiasm_cmd_hap0, check=True)
            
            # Convert GFA to FASTA
            gfa_to_fasta_cmd = [
                "awk", 
                "/^S/{print \">\"$2\"\\n\"$3}", 
                f"{hap0_prefix}.p_ctg.gfa"
            ]
            
            with open(f"{hap0_prefix}.fasta", "w") as fasta_file:
                subprocess.run(gfa_to_fasta_cmd, stdout=fasta_file, check=True)
                
            logger.info(f"Haplotype A assembly completed: {hap0_prefix}.fasta")
            
        except subprocess.CalledProcessError as e:
            logger.error(f"Error assembling haplotype A: {e}")
        
        # Assemble haplotype B
        logger.info("Assembling haplotype B")
        
        hap1_prefix = os.path.join(self.args.output_dir, "assemblies", "Ptr_hapB.GenoPhase")
        
        hifiasm_cmd_hap1 = [
            "hifiasm", 
            "-o", hap1_prefix,
            "-t", str(self.args.threads),
            "--primary",
            hap1_fasta
        ]
        
        try:
            subprocess.run(hifiasm_cmd_hap1, check=True)
            
            # Convert GFA to FASTA
            gfa_to_fasta_cmd = [
                "awk", 
                "/^S/{print \">\"$2\"\\n\"$3}", 
                f"{hap1_prefix}.p_ctg.gfa"
            ]
            
            with open(f"{hap1_prefix}.fasta", "w") as fasta_file:
                subprocess.run(gfa_to_fasta_cmd, stdout=fasta_file, check=True)
                
            logger.info(f"Haplotype B assembly completed: {hap1_prefix}.fasta")
            
        except subprocess.CalledProcessError as e:
            logger.error(f"Error assembling haplotype B: {e}")

def main():
    parser = argparse.ArgumentParser(description='GenoPhase: A tool for haplotype-resolved genome assembly using HiFi reads')
    
    # Input/output options
    parser.add_argument('--hifi-reads', required=True, help='Input HiFi reads in FASTA/FASTQ format')
    parser.add_argument('--reference', help='Reference genome in FASTA format (optional, will be built if not provided)')
    parser.add_argument('--output-dir', required=True, help='Output directory')
    
    # Mode selection
    parser.add_argument('--mode', choices=['full', 'reference', 'phase', 'assemble'], default='full',
                       help='Pipeline mode: full=complete pipeline, reference=only build reference, phase=only phase reads, assemble=only assemble haplotypes')
    
    # Performance options
    parser.add_argument('--threads', type=int, default=16, help='Number of threads to use')
    
    # Graph construction parameters
    parser.add_argument('--min-base-quality', type=int, default=20, help='Minimum base quality for SNP detection')
    parser.add_argument('--min-depth', type=int, default=10, help='Minimum read depth for SNP detection')
    parser.add_argument('--max-depth', type=int, default=100, help='Maximum read depth for SNP detection')
    parser.add_argument('--min-minor-allele-freq', type=float, default=0.25, help='Minimum minor allele frequency for SNP detection')
    parser.add_argument('--max-edge-distance', type=int, default=16000, help='Maximum distance between SNPs to create an edge (D_max)')
    parser.add_argument('--min-shared-reads', type=int, default=3, help='Minimum number of shared reads to create an edge (N_min)')
    
    # Other options
    parser.add_argument('--force-remap', action='store_true', help='Force remapping of reads even if BAM file exists')
    
    args = parser.parse_args()
    
    # Create GenoPhase instance and run the pipeline
    pipeline = GenoPhase(args)
    pipeline.run()

if __name__ == "__main__":
    main()
