#!/usr/bin/env python3
"""
Merge two viral assemblies using read support as evidence.

Uses clear decision rules based on:
1. Coverage depth (more reads = more confident)
2. Mapping quality (better alignment = more accurate)
3. Contig length (longer = more complete)

No arbitrary scoring - just transparent, defensible logic.
"""

import argparse
import sys
import subprocess
from pathlib import Path
from collections import defaultdict
import pandas as pd
import pysam
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import numpy as np

def run_command(cmd, description):
    """Run shell command with error handling."""
    print(f"[INFO] {description}")
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"[ERROR] Command failed: {cmd}", file=sys.stderr)
        print(result.stderr, file=sys.stderr)
        sys.exit(1)
    return result

def find_overlaps(assembly1, assembly2, workdir, min_overlap=500, min_identity=95, threads=4):
    """
    Find overlapping contigs between two assemblies using minimap2.
    
    Returns: List of overlap records
    """
    print(f"\n[STEP 1] Finding overlaps between assemblies...")
    
    # Align assembly2 to assembly1
    paf_file = workdir / "overlaps.paf"
    cmd = f"minimap2 -x asm5 -t {threads} {assembly1} {assembly2} > {paf_file}"
    run_command(cmd, "Running minimap2 to find overlaps")
    
    # Parse PAF output
    overlaps = []
    with open(paf_file) as f:
        for line in f:
            fields = line.strip().split('\t')
            
            query_name = fields[0]
            query_len = int(fields[1])
            query_start = int(fields[2])
            query_end = int(fields[3])
            
            target_name = fields[5]
            target_len = int(fields[6])
            target_start = int(fields[7])
            target_end = int(fields[8])
            
            matches = int(fields[9])
            alignment_len = int(fields[10])
            mapq = int(fields[11])
            
            # Calculate overlap length and identity
            overlap_len = min(query_end - query_start, target_end - target_start)
            identity = (matches / alignment_len * 100) if alignment_len > 0 else 0
            
            # Filter by thresholds
            if overlap_len >= min_overlap and identity >= min_identity:
                overlaps.append({
                    'query': query_name,
                    'query_len': query_len,
                    'target': target_name,
                    'target_len': target_len,
                    'overlap_len': overlap_len,
                    'identity': identity,
                    'mapq': mapq,
                })
    
    print(f"[INFO] Found {len(overlaps)} significant overlaps (>={min_overlap}bp, >={min_identity}% identity)")
    return overlaps

def map_reads_to_contig(contig_file, reads1, reads2, workdir, threads=4):
    """
    Map reads to a contig and calculate support metrics.
    
    Returns: dict of metrics
    """
    bam_file = workdir / f"{Path(contig_file).stem}.bam"
    
    # Map reads
    cmd = f"""
    minimap2 -ax sr -t {threads} {contig_file} {reads1} {reads2} 2>/dev/null | \
    samtools view -bS - | \
    samtools sort -@ {threads} -o {bam_file} - && \
    samtools index {bam_file}
    """
    run_command(cmd, f"Mapping reads to {Path(contig_file).stem}")
    
    # Calculate metrics
    metrics = calculate_contig_metrics(bam_file)
    
    return metrics

def calculate_contig_metrics(bam_file):
    """
    Calculate detailed metrics from BAM file.
    
    Returns: dict with coverage, mapq, proper pairs, etc.
    """
    bamfile = pysam.AlignmentFile(bam_file, 'rb')
    
    all_metrics = {}
    
    for contig in bamfile.references:
        contig_length = bamfile.get_reference_length(contig)
        
        # Initialize coverage array
        coverage = np.zeros(contig_length)
        
        # Collect read statistics
        mapq_scores = []
        num_reads = 0
        num_proper_pairs = 0
        
        # Get coverage per position
        for pileup_column in bamfile.pileup(contig):
            if pileup_column.pos < contig_length:
                coverage[pileup_column.pos] = pileup_column.nsegments
        
        # Get read mapping statistics
        for read in bamfile.fetch(contig):
            if not read.is_unmapped:
                num_reads += 1
                mapq_scores.append(read.mapping_quality)
                if read.is_proper_pair:
                    num_proper_pairs += 1
        
        # Calculate summary metrics
        all_metrics[contig] = {
            'length': contig_length,
            'num_reads': num_reads,
            'num_proper_pairs': num_proper_pairs,
            'avg_coverage': float(np.mean(coverage)),
            'median_coverage': float(np.median(coverage)),
            'min_coverage': float(np.min(coverage)),
            'max_coverage': float(np.max(coverage)),
            'percent_covered': float(np.sum(coverage > 0) / contig_length * 100),
            'avg_mapq': float(np.mean(mapq_scores)) if mapq_scores else 0.0,
            'proper_pair_ratio': num_proper_pairs / num_reads if num_reads > 0 else 0.0,
        }
    
    bamfile.close()
    return all_metrics

def choose_better_contig(contig_a, contig_b, metrics_a, metrics_b):
    """
    Choose better contig based on clear, defensible criteria.
    
    Priority order (most important first):
    1. Coverage depth (more reads = more confident)
    2. Mapping quality (better alignment = more accurate)
    3. Contig length (longer = more complete)
    
    Returns: (winner_name, loser_name, reason, decision_category)
    """
    
    # Rule 1: If coverage differs by >2x, choose higher coverage
    if metrics_a['avg_coverage'] > 0 and metrics_b['avg_coverage'] > 0:
        coverage_ratio = max(metrics_a['avg_coverage'], metrics_b['avg_coverage']) / \
                        min(metrics_a['avg_coverage'], metrics_b['avg_coverage'])
        
        if coverage_ratio > 2.0:
            if metrics_a['avg_coverage'] > metrics_b['avg_coverage']:
                return (contig_a, contig_b, 
                       f"Much higher coverage ({metrics_a['avg_coverage']:.1f}x vs {metrics_b['avg_coverage']:.1f}x, ratio={coverage_ratio:.2f})",
                       "coverage_difference")
            else:
                return (contig_b, contig_a,
                       f"Much higher coverage ({metrics_b['avg_coverage']:.1f}x vs {metrics_a['avg_coverage']:.1f}x, ratio={coverage_ratio:.2f})",
                       "coverage_difference")
    
    # Rule 2: If coverage similar, choose better mapping quality
    mapq_diff = abs(metrics_a['avg_mapq'] - metrics_b['avg_mapq'])
    if mapq_diff > 10:
        if metrics_a['avg_mapq'] > metrics_b['avg_mapq']:
            return (contig_a, contig_b,
                   f"Better mapping quality ({metrics_a['avg_mapq']:.1f} vs {metrics_b['avg_mapq']:.1f})",
                   "mapq_difference")
        else:
            return (contig_b, contig_a,
                   f"Better mapping quality ({metrics_b['avg_mapq']:.1f} vs {metrics_a['avg_mapq']:.1f})",
                   "mapq_difference")
    
    # Rule 3: If both coverage and quality similar, choose longer contig
    if metrics_a['length'] != metrics_b['length']:
        if metrics_a['length'] > metrics_b['length']:
            return (contig_a, contig_b,
                   f"Longer contig with similar quality ({metrics_a['length']}bp vs {metrics_b['length']}bp)",
                   "length_difference")
        else:
            return (contig_b, contig_a,
                   f"Longer contig with similar quality ({metrics_b['length']}bp vs {metrics_a['length']}bp)",
                   "length_difference")
    
    # Rule 4: If essentially identical, arbitrarily choose first (shouldn't happen often)
    return (contig_a, contig_b,
           "Contigs nearly identical, keeping first assembly version",
           "arbitrary")

def extract_contig_to_file(fasta_file, contig_id, output_file):
    """Extract a specific contig to its own file."""
    for record in SeqIO.parse(fasta_file, 'fasta'):
        if record.id == contig_id:
            SeqIO.write([record], output_file, 'fasta')
            return
    raise ValueError(f"Contig {contig_id} not found in {fasta_file}")

def cleanup_temp_files(workdir):
    """Clean up temporary files to save disk space."""
    print(f"\n[CLEANUP] Removing temporary files from {workdir}...")
    
    patterns = [
        "temp_target_*.fasta",
        "temp_query_*.fasta",
        "temp_target_*.bam",
        "temp_query_*.bam",
        "temp_target_*.bam.bai",
        "temp_query_*.bam.bai",
        "overlaps.paf",
    ]
    
    files_deleted = 0
    space_saved = 0
    
    for pattern in patterns:
        for file in workdir.glob(pattern):
            try:
                file_size = file.stat().st_size
                file.unlink()
                files_deleted += 1
                space_saved += file_size
            except Exception as e:
                print(f"  Warning: Could not delete {file.name}: {e}")
    
    # Convert bytes to human-readable
    if space_saved > 1e9:
        space_str = f"{space_saved / 1e9:.2f} GB"
    elif space_saved > 1e6:
        space_str = f"{space_saved / 1e6:.2f} MB"
    else:
        space_str = f"{space_saved / 1e3:.2f} KB"
    
    print(f"  Deleted {files_deleted} temporary files, freed ~{space_str}")

def merge_assemblies(assembly1, assembly2, reads1, reads2, workdir, threads=4, 
                    min_overlap=500, min_identity=95):
    """
    Main merging logic.
    
    Returns: (final_contigs, decisions_dataframe)
    """
    
    # Find overlaps
    overlaps = find_overlaps(assembly1, assembly2, workdir, min_overlap, min_identity, threads)
    
    # Track which contigs have been processed
    assembly1_processed = set()
    assembly2_processed = set()
    
    decisions = []
    final_contigs = {}
    
    print(f"\n[STEP 2] Evaluating overlapping contigs with read support...")
    
    # Process overlapping pairs
    for i, overlap in enumerate(overlaps, 1):
        print(f"\n[OVERLAP {i}/{len(overlaps)}] Comparing {overlap['target']} vs {overlap['query']}")
        
        # Skip if already processed
        if overlap['target'] in assembly1_processed or overlap['query'] in assembly2_processed:
            print(f"  Skipping - one or both contigs already processed")
            continue
        
        # Extract contigs to temporary files
        target_file = workdir / f"temp_target_{i}.fasta"
        query_file = workdir / f"temp_query_{i}.fasta"
        
        extract_contig_to_file(assembly1, overlap['target'], target_file)
        extract_contig_to_file(assembly2, overlap['query'], query_file)
        
        # Map reads to both contigs
        print(f"  Mapping reads to both contigs...")
        metrics_target = map_reads_to_contig(target_file, reads1, reads2, workdir, threads)
        metrics_query = map_reads_to_contig(query_file, reads1, reads2, workdir, threads)
        
        # Get metrics for the specific contigs
        met_target = metrics_target[overlap['target']]
        met_query = metrics_query[overlap['query']]
        
        print(f"  {overlap['target']}: {met_target['avg_coverage']:.1f}x coverage, MAPQ={met_target['avg_mapq']:.1f}")
        print(f"  {overlap['query']}: {met_query['avg_coverage']:.1f}x coverage, MAPQ={met_query['avg_mapq']:.1f}")
        
        # Make decision
        winner, loser, reason, category = choose_better_contig(
            overlap['target'], overlap['query'],
            met_target, met_query
        )
        
        print(f"  DECISION: Keep {winner}")
        print(f"  REASON: {reason}")
        
        # Record decision
        decisions.append({
            'contig_kept': winner,
            'contig_discarded': loser,
            'reason': reason,
            'decision_category': category,
            'overlap_length': overlap['overlap_len'],
            'overlap_identity': overlap['identity'],
            'kept_coverage': met_target['avg_coverage'] if winner == overlap['target'] else met_query['avg_coverage'],
            'discarded_coverage': met_query['avg_coverage'] if winner == overlap['target'] else met_target['avg_coverage'],
            'kept_mapq': met_target['avg_mapq'] if winner == overlap['target'] else met_query['avg_mapq'],
            'discarded_mapq': met_query['avg_mapq'] if winner == overlap['target'] else met_target['avg_mapq'],
            'kept_length': met_target['length'] if winner == overlap['target'] else met_query['length'],
            'discarded_length': met_query['length'] if winner == overlap['target'] else met_target['length'],
        })
        
        # Add winner to final assembly
        source_file = assembly1 if winner == overlap['target'] else assembly2
        for record in SeqIO.parse(source_file, 'fasta'):
            if record.id == winner:
                final_contigs[winner] = record
                break
        
        # Mark both as processed
        assembly1_processed.add(overlap['target'])
        assembly2_processed.add(overlap['query'])
    
    print(f"\n[STEP 3] Adding singleton contigs (no overlaps)...")
    
    # Add singletons from assembly1
    for record in SeqIO.parse(assembly1, 'fasta'):
        if record.id not in assembly1_processed:
            print(f"  Adding singleton from assembly1: {record.id}")
            final_contigs[record.id] = record
            decisions.append({
                'contig_kept': record.id,
                'contig_discarded': None,
                'reason': 'Singleton from assembly1 (no overlap found)',
                'decision_category': 'singleton',
            })
    
    # Add singletons from assembly2
    for record in SeqIO.parse(assembly2, 'fasta'):
        if record.id not in assembly2_processed:
            print(f"  Adding singleton from assembly2: {record.id}")
            final_contigs[record.id] = record
            decisions.append({
                'contig_kept': record.id,
                'contig_discarded': None,
                'reason': 'Singleton from assembly2 (no overlap found)',
                'decision_category': 'singleton',
            })
    
    return list(final_contigs.values()), pd.DataFrame(decisions)

def main():
    parser = argparse.ArgumentParser(
        description='Merge two viral assemblies using read support',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument('--assembly1', required=True, help='First assembly (e.g., EsViritu)')
    parser.add_argument('--assembly2', required=True, help='Second assembly (e.g., metaviralSPADES)')
    parser.add_argument('--reads1', required=True, help='Forward reads (FASTQ)')
    parser.add_argument('--reads2', required=True, help='Reverse reads (FASTQ)')
    parser.add_argument('--output', required=True, help='Output merged assembly FASTA')
    parser.add_argument('--report', required=True, help='Output decision report TSV')
    parser.add_argument('--workdir', required=True, help='Working directory for temp files')
    parser.add_argument('--threads', type=int, default=4, help='Number of threads')
    parser.add_argument('--min-overlap', type=int, default=500, help='Minimum overlap length (bp)')
    parser.add_argument('--min-identity', type=float, default=95, help='Minimum overlap identity (%%)')
    
    args = parser.parse_args()
    
    # Create working directory
    workdir = Path(args.workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    
    print("="*80)
    print("READ-SUPPORTED ASSEMBLY MERGER")
    print("="*80)
    print(f"Assembly 1: {args.assembly1}")
    print(f"Assembly 2: {args.assembly2}")
    print(f"Min overlap: {args.min_overlap} bp")
    print(f"Min identity: {args.min_identity}%")
    print(f"Threads: {args.threads}")
    print("="*80)
    
    # Run merging
    final_contigs, decisions_df = merge_assemblies(
        args.assembly1, args.assembly2,
        args.reads1, args.reads2,
        workdir, args.threads,
        args.min_overlap, args.min_identity
    )
    
    # Write outputs
    print(f"\n[STEP 4] Writing outputs...")
    SeqIO.write(final_contigs, args.output, 'fasta')
    print(f"  Wrote {len(final_contigs)} contigs to {args.output}")
    
    decisions_df.to_csv(args.report, sep='\t', index=False)
    print(f"  Wrote decision report to {args.report}")
    
    # Clean up temporary files
    cleanup_temp_files(workdir)
    
    # Print summary
    print("\n" + "="*80)
    print("SUMMARY")
    print("="*80)
    print(f"Total contigs in final assembly: {len(final_contigs)}")
    print(f"\nDecisions by category:")
    print(decisions_df['decision_category'].value_counts().to_string())
    print("="*80)

if __name__ == '__main__':
    main()