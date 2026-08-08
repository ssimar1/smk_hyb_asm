#!/usr/bin/env python3
"""Annotate BLAST results with EsViritu database metadata."""

import argparse
import pandas as pd

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--blast', required=True, help='BLAST output file')
    parser.add_argument('--metadata', required=True, help='Database metadata TSV')
    parser.add_argument('--output-summary', required=True, help='Output community summary TSV')
    parser.add_argument('--output-per-contig', required=True, help='Output per-contig identity TSV')
    args = parser.parse_args()
    
    # Read BLAST results
    blast_cols = ['qseqid', 'sseqid', 'pident', 'length', 'qlen', 'slen', 'qcovs', 'evalue', 'bitscore']
    
    try:
        blast_df = pd.read_csv(args.blast, sep='\t', names=blast_cols)
    except pd.errors.EmptyDataError:
        print("Warning: No BLAST results (empty file)")
        # Create empty outputs
        pd.DataFrame(columns=['contig_id', 'virus_name', 'segment', 'identity', 'coverage']).to_csv(
            args.output_summary, sep='\t', index=False
        )
        pd.DataFrame(columns=['contig_id', 'reference_accession', 'virus_name', 'segment', 
                              'identity', 'coverage', 'contig_length', 'reference_length']).to_csv(
            args.output_per_contig, sep='\t', index=False
        )
        return
    
    if blast_df.empty:
        print("Warning: No BLAST results (no hits)")
        pd.DataFrame(columns=['contig_id', 'virus_name', 'segment', 'identity', 'coverage']).to_csv(
            args.output_summary, sep='\t', index=False
        )
        pd.DataFrame(columns=['contig_id', 'reference_accession', 'virus_name', 'segment', 
                              'identity', 'coverage', 'contig_length', 'reference_length']).to_csv(
            args.output_per_contig, sep='\t', index=False
        )
        return
    
    # Read metadata
    metadata_df = pd.read_csv(args.metadata, sep='\t')
    
    # Merge on accession
    result = blast_df.merge(
        metadata_df,
        left_on='sseqid',
        right_on='Accession',
        how='left'
    )
    
    # Clean up subspecies
    if 'subspecies' in result.columns:
        result['subspecies'] = result['subspecies'].str.replace('t__', '', regex=False)
        # Use subspecies as primary name if available, else use Name
        result['virus_name'] = result['subspecies'].fillna(result.get('Name', 'Unknown'))
    else:
        result['virus_name'] = result.get('Name', 'Unknown')
    
    # Fill NA segments
    if 'Segment' in result.columns:
        result['segment'] = result['Segment'].fillna('NA')
    else:
        result['segment'] = 'NA'
    
    # === Per-contig identity file ===
    per_contig = result[[
        'qseqid', 'sseqid', 'virus_name', 'segment', 
        'pident', 'qcovs', 'qlen', 'slen', 'evalue', 'bitscore'
    ]].copy()
    
    per_contig = per_contig.rename(columns={
        'qseqid': 'contig_id',
        'sseqid': 'reference_accession',
        'pident': 'percent_identity',
        'qcovs': 'query_coverage',
        'qlen': 'contig_length',
        'slen': 'reference_length',
    })
    
    per_contig.to_csv(args.output_per_contig, sep='\t', index=False)
    print(f"Wrote per-contig identity file: {args.output_per_contig}")
    
    # === Community summary (aggregated by virus) ===
    summary = result.groupby(['virus_name', 'segment']).agg({
        'qseqid': 'count',  # Number of contigs
        'pident': 'mean',   # Average identity
        'qcovs': 'mean',    # Average coverage
        'qlen': 'sum',      # Total assembled length
    }).reset_index()
    
    summary = summary.rename(columns={
        'qseqid': 'num_contigs',
        'pident': 'avg_identity',
        'qcovs': 'avg_coverage',
        'qlen': 'total_length',
    })
    
    # Round percentages
    summary['avg_identity'] = summary['avg_identity'].round(2)
    summary['avg_coverage'] = summary['avg_coverage'].round(2)
    
    # Sort by number of contigs (descending)
    summary = summary.sort_values('num_contigs', ascending=False)
    
    summary.to_csv(args.output_summary, sep='\t', index=False)
    print(f"Wrote community summary: {args.output_summary}")
    
    # Print summary to stdout
    print(f"\nClassified {len(result)} contigs")
    print(f"\n=== COMMUNITY COMPOSITION ===")
    print(summary.to_string(index=False))
    
    print(f"\n=== PER-CONTIG IDENTITY RANGE ===")
    print(f"Min identity: {per_contig['percent_identity'].min():.2f}%")
    print(f"Max identity: {per_contig['percent_identity'].max():.2f}%")
    print(f"Mean identity: {per_contig['percent_identity'].mean():.2f}%")
    print(f"Median identity: {per_contig['percent_identity'].median():.2f}%")

if __name__ == '__main__':
    main()
```