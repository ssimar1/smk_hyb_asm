#!/usr/bin/env python3
"""
Prepare sample sheet for smk_tephi_virome Snakemake pipeline.

Converts a basic SampleListRaw file to the proper sample_sheet.tsv format.

Input format (tab-separated, no header):
    Column 1: Sample_ID
    Column 2: forward_read path
    Column 3: reverse_read path
    Additional columns: ignored

Output format (tab-separated, with header):
    Sample_ID, forward_read, reverse_read, read_compression, category
"""

import argparse
import sys
from pathlib import Path


def detect_compression(filepath: str) -> str:
    """Detect compression type from file extension."""
    path = filepath.lower()
    if path.endswith('.bz2'):
        return 'bz2'
    elif path.endswith('.gz') or path.endswith('.gzip'):
        return 'gz'
    else:
        return 'none'


def detect_category(sample_id: str) -> str:
    """Detect sample category based on sample ID."""
    # this is the most reasonable heuristic I can think of
    if 'control' in sample_id.lower() or 'hybrid' in sample_id.lower():
        return 'control'
    return 'sample'


def convert_sample_list(input_file: str, output_file: str, path_prefix_replace: str = None) -> None:
    """
    Convert SampleListRaw to sample_sheet.tsv format.
    
    Args:
        input_file: Path to input SampleListRaw file
        output_file: Path to output sample_sheet.tsv file
        path_prefix_replace: Optional old:new path prefix replacement (e.g., /gpfs1/projects:/data/service)
    """
    old_prefix = None
    new_prefix = None
    if path_prefix_replace:
        parts = path_prefix_replace.split(':')
        if len(parts) == 2:
            old_prefix, new_prefix = parts
        else:
            print(f"Warning: Invalid path prefix format '{path_prefix_replace}'. Expected 'old:new'.", file=sys.stderr)
    
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        # Write header
        outfile.write("Sample_ID\tforward_read\treverse_read\tread_compression\tcategory\n")
        
        for line_num, line in enumerate(infile, 1):
            line = line.strip()
            if not line:
                continue
            
            fields = line.split('\t')
            if len(fields) < 3:
                print(f"Warning: Line {line_num} has fewer than 3 columns, skipping: {line[:50]}...", file=sys.stderr)
                continue
            
            sample_id = fields[0]
            forward_read = fields[1]
            reverse_read = fields[2]
            
            # Apply path prefix replacement if specified
            if old_prefix and new_prefix:
                forward_read = forward_read.replace(old_prefix, new_prefix)
                reverse_read = reverse_read.replace(old_prefix, new_prefix)
            
            # Detect compression from forward read (assume both reads have same compression)
            compression = detect_compression(forward_read)
            
            # Detect category
            category = detect_category(sample_id)
            
            outfile.write(f"{sample_id}\t{forward_read}\t{reverse_read}\t{compression}\t{category}\n")
    
    print(f"Successfully converted {input_file} -> {output_file}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description="Convert SampleListRaw to sample_sheet.tsv format for smk_tephi_virome pipeline.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Basic conversion
    python prepare_sample_sheet.py input.txt output.tsv
    
    # With path prefix replacement
    python prepare_sample_sheet.py input.txt output.tsv --path-replace "/gpfs1/projects:/data/service"
        """
    )
    parser.add_argument(
        'input_file',
        help='Input SampleListRaw file (tab-separated: Sample_ID, forward_read, reverse_read, ...)'
    )
    parser.add_argument(
        'output_file',
        help='Output sample_sheet.tsv file'
    )
    parser.add_argument(
        '--path-replace', '-p',
        metavar='OLD:NEW',
        help='Replace path prefix (e.g., "/gpfs1/projects:/data/service")'
    )
    
    args = parser.parse_args()
    
    # Validate input file exists
    if not Path(args.input_file).exists():
        print(f"Error: Input file '{args.input_file}' does not exist.", file=sys.stderr)
        sys.exit(1)
    
    convert_sample_list(args.input_file, args.output_file, args.path_replace)


if __name__ == '__main__':
    main()
