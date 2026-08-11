#!/usr/bin/env python3
"""
annotate_rvdb_tiered.py

Final annotation step for the tiered RVDB classification.

Takes the two BLAST result files (nucleotide vs C-RVDB, protein vs
RVDB-prot) plus the RVDB nucleotide taxonomy table, and produces:

  1. a per-contig table with full lineage (realm -> species), the
     accession it was annotated from, the BLAST identity/coverage, and
     a `tier` column showing HOW each call was made, and
  2. a community summary grouped by family + species.

How each contig gets a tier
---------------------------
  nucleotide       best BLASTn hit passed threshold (>=90% id, >=50% cov)
  protein          failed nucleotide, but got a BLASTx hit
  nucleotide_weak  failed nucleotide, no protein hit, but had a
                   sub-threshold nucleotide hit (reported, honestly, as weak)
  none             no hit in either database

The join to taxonomy
--------------------
  nucleotide tier  accession = field 2 (0-indexed) of the subject stitle
                   e.g. acc|GENBANK|MG941528.1|... -> MG941528.1
  protein tier     accession = field 4 of the subject stitle (the SOURCE
                   nucleotide accession the protein was derived from)
                   e.g. acc|GENBANK|AYD68779.1|GENBANK|MH171300|... -> MH171300

Both are version-stripped (MG941528.1 -> MG941528) and looked up in the
same taxonomy table. Anything that doesn't join gets Unknown lineage but
keeps its accession, so it's still traceable.
"""

import argparse
import re
import sys

import pandas as pd
from Bio import SeqIO

BLAST_COLS = [
    "qseqid", "sseqid", "pident", "length", "qlen",
    "slen", "qcovs", "evalue", "bitscore", "stitle",
]

# Lineage ranks we carry through, in order.
RANKS = ["realm", "kingdom", "phylum", "class", "order",
         "family", "genus", "species", "subspecies"]

VERSION_RE = re.compile(r"\.\d+$")


def strip_version(acc):
    """MG941528.1 -> MG941528 ; leaves unversioned accessions alone."""
    if acc is None:
        return None
    return VERSION_RE.sub("", str(acc).strip())


def acc_from_stitle(stitle, field_index):
    """
    Pull the accession out of an RVDB subject title by pipe position.
    field_index = 2 for nucleotide titles, 4 for protein titles.
    Returns None if the title is too short / malformed.
    """
    if not isinstance(stitle, str):
        return None
    parts = stitle.split("|")
    if len(parts) > field_index:
        return parts[field_index].strip()
    return None


def load_blast(path):
    try:
        df = pd.read_csv(path, sep="\t", names=BLAST_COLS)
    except (pd.errors.EmptyDataError, FileNotFoundError):
        return pd.DataFrame(columns=BLAST_COLS)
    return df


def best_hit_per_query(df):
    if df.empty:
        return df
    return df.sort_values("bitscore", ascending=False).drop_duplicates("qseqid")


def find_header_row(tax_path):
    """The taxonomy file starts with one or more #### comment lines; find
    the real header row (the one beginning with 'accession')."""
    with open(tax_path) as f:
        for i, line in enumerate(f):
            if line.startswith("accession\t"):
                return i
    sys.exit(f"[annotate] Could not find header row in {tax_path}")


def load_taxonomy(tax_path, needed_accessions):
    """
    Read the taxonomy table but keep only rows whose (version-stripped)
    accession is one we actually need. Chunked so we never hold the whole
    ~1.3M-row table in memory.

    Returns: dict {base_accession: {rank: value, ...}}
    """
    header_row = find_header_row(tax_path)
    keep_cols = ["accession"] + RANKS

    lineage = {}
    reader = pd.read_csv(
        tax_path, sep="\t", skiprows=header_row, dtype=str,
        usecols=lambda c: c in keep_cols, chunksize=100_000,
    )
    for chunk in reader:
        chunk["base_acc"] = chunk["accession"].map(strip_version)
        chunk = chunk[chunk["base_acc"].isin(needed_accessions)]
        for _, row in chunk.iterrows():
            entry = {}
            for rank in RANKS:
                val = row.get(rank)
                # RVDB uses '-' (and sometimes blank) for missing ranks.
                if val is None or str(val).strip() in ("-", "", "nan"):
                    entry[rank] = "Unknown"
                else:
                    entry[rank] = str(val).strip()
            lineage[row["base_acc"]] = entry
    return lineage


def unknown_lineage():
    return {rank: "Unknown" for rank in RANKS}


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--nucl-blast", required=True)
    parser.add_argument("--prot-blast", required=True)
    parser.add_argument("--taxonomy", required=True,
                        help="RVDB nucleotide taxonomy .tab (full/unclustered)")
    parser.add_argument("--hybrid-fasta", required=True,
                        help="Hybrid assembly FASTA (for the full contig list)")
    parser.add_argument("--output-summary", required=True)
    parser.add_argument("--output-per-contig", required=True)
    parser.add_argument("--min-identity", type=float, default=90.0)
    parser.add_argument("--min-coverage", type=float, default=50.0)
    args = parser.parse_args()

    nucl = best_hit_per_query(load_blast(args.nucl_blast))
    prot = best_hit_per_query(load_blast(args.prot_blast))

    # Index best hits by contig id for quick lookup.
    nucl_by_contig = {r["qseqid"]: r for _, r in nucl.iterrows()}
    prot_by_contig = {r["qseqid"]: r for _, r in prot.iterrows()}

    # Collect every accession we'll need, so taxonomy load stays small.
    needed = set()
    for r in nucl_by_contig.values():
        a = strip_version(acc_from_stitle(r["stitle"], 2))
        if a:
            needed.add(a)
    for r in prot_by_contig.values():
        a = strip_version(acc_from_stitle(r["stitle"], 4))
        if a:
            needed.add(a)

    lineage = load_taxonomy(args.taxonomy, needed) if needed else {}

    # Walk every contig in the hybrid assembly and decide its tier.
    rows = []
    for record in SeqIO.parse(args.hybrid_fasta, "fasta"):
        cid = record.id
        contig_len = len(record.seq)

        nhit = nucl_by_contig.get(cid)
        phit = prot_by_contig.get(cid)

        nucl_pass = (
            nhit is not None
            and nhit["pident"] >= args.min_identity
            and nhit["qcovs"] >= args.min_coverage
        )

        if nucl_pass:
            tier = "nucleotide"
            acc = strip_version(acc_from_stitle(nhit["stitle"], 2))
            pident, qcov = nhit["pident"], nhit["qcovs"]
        elif phit is not None:
            tier = "protein"
            acc = strip_version(acc_from_stitle(phit["stitle"], 4))
            pident, qcov = phit["pident"], phit["qcovs"]
        elif nhit is not None:
            tier = "nucleotide_weak"
            acc = strip_version(acc_from_stitle(nhit["stitle"], 2))
            pident, qcov = nhit["pident"], nhit["qcovs"]
        else:
            tier = "none"
            acc = None
            pident, qcov = None, None

        lin = lineage.get(acc, unknown_lineage()) if acc else unknown_lineage()

        row = {
            "contig_id": cid,
            "tier": tier,
            "accession": acc if acc else "NA",
            "percent_identity": round(pident, 2) if pident is not None else "NA",
            "query_coverage": round(qcov, 2) if qcov is not None else "NA",
            "contig_length": contig_len,
        }
        row.update(lin)
        rows.append(row)

    per_contig = pd.DataFrame(rows)

    # Column order for the per-contig table.
    ordered = (["contig_id", "tier", "accession", "percent_identity",
                "query_coverage", "contig_length"] + RANKS)
    per_contig = per_contig[ordered]
    per_contig.to_csv(args.output_per_contig, sep="\t", index=False)

    # Community summary: group by family + species.
    if per_contig.empty:
        summary = pd.DataFrame(
            columns=["family", "species", "num_contigs", "avg_identity", "total_length"]
        )
    else:
        # Only average identity over contigs that actually had a hit.
        pc = per_contig.copy()
        pc["percent_identity_num"] = pd.to_numeric(pc["percent_identity"], errors="coerce")
        summary = (
            pc.groupby(["family", "species"])
              .agg(num_contigs=("contig_id", "count"),
                   avg_identity=("percent_identity_num", "mean"),
                   total_length=("contig_length", "sum"))
              .reset_index()
              .sort_values("num_contigs", ascending=False)
        )
        summary["avg_identity"] = summary["avg_identity"].round(2)

    summary.to_csv(args.output_summary, sep="\t", index=False)

    # Console recap.
    tier_counts = per_contig["tier"].value_counts().to_dict() if not per_contig.empty else {}
    sys.stderr.write(f"[annotate] contigs by tier: {tier_counts}\n")
    sys.stderr.write(f"[annotate] wrote {args.output_per_contig} and {args.output_summary}\n")


if __name__ == "__main__":
    main()