#!/usr/bin/env python3
"""
select_nucl_misses.py

Step between the two BLAST tiers.

Reads the nucleotide (BLASTn vs C-RVDB) results, decides which hybrid
contigs got a *confident* nucleotide hit, and writes out a FASTA of only
the contigs that did NOT — those are the ones that fall through to the
protein (BLASTx) tier.

A contig "passes" the nucleotide tier if its BEST hit (highest bitscore)
has BOTH:
    percent identity >= --min-identity   (default 90)
    query coverage    >= --min-coverage  (default 50)

Everything else (below threshold, or no hit at all) is written to the
miss FASTA. This same pass rule is applied again, identically, in
annotate_rvdb_tiered.py -- keep the two in sync if you ever change it.
"""

import argparse
import sys

import pandas as pd
from Bio import SeqIO

# Column order we ask BLAST to emit (outfmt 6). stitle is last.
BLAST_COLS = [
    "qseqid", "sseqid", "pident", "length", "qlen",
    "slen", "qcovs", "evalue", "bitscore", "stitle",
]


def load_blast(path):
    """Read a BLAST outfmt-6 file; return empty frame if there are no hits."""
    try:
        df = pd.read_csv(path, sep="\t", names=BLAST_COLS)
    except (pd.errors.EmptyDataError, FileNotFoundError):
        return pd.DataFrame(columns=BLAST_COLS)
    return df


def best_hit_per_query(df):
    """Keep one row per query: the highest-bitscore HSP."""
    if df.empty:
        return df
    return df.sort_values("bitscore", ascending=False).drop_duplicates("qseqid")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--nucl-blast", required=True,
                        help="BLASTn results (outfmt 6, with stitle)")
    parser.add_argument("--hybrid-fasta", required=True,
                        help="The hybrid assembly FASTA (all contigs)")
    parser.add_argument("--output-fasta", required=True,
                        help="Where to write contigs needing the protein tier")
    parser.add_argument("--min-identity", type=float, default=90.0)
    parser.add_argument("--min-coverage", type=float, default=50.0)
    args = parser.parse_args()

    blast = best_hit_per_query(load_blast(args.nucl_blast))

    # Which contigs cleared BOTH thresholds on their best hit?
    if blast.empty:
        passed = set()
    else:
        ok = (blast["pident"] >= args.min_identity) & (blast["qcovs"] >= args.min_coverage)
        passed = set(blast.loc[ok, "qseqid"])

    # Write every contig that did NOT pass to the miss FASTA.
    n_total = 0
    n_miss = 0
    with open(args.output_fasta, "w") as out:
        for record in SeqIO.parse(args.hybrid_fasta, "fasta"):
            n_total += 1
            if record.id not in passed:
                SeqIO.write([record], out, "fasta")
                n_miss += 1

    sys.stderr.write(
        f"[select_nucl_misses] {n_total} contigs total; "
        f"{len(passed)} passed nucleotide tier "
        f"(>= {args.min_identity}% id AND >= {args.min_coverage}% cov); "
        f"{n_miss} sent to protein tier\n"
    )


if __name__ == "__main__":
    main()