# smk_hyb_asm

A Snakemake pipeline for hybrid capture virome sequencing. It preprocesses reads, detects viruses with **EsViritu** (reference-guided), assembles de novo with **metaviralSPADES**, merges the two assemblies into a single hybrid assembly using **read support** as the deciding evidence, and classifies the result against **RVDB** (nucleotide + protein tiers).

**Version:** 0.2.0

---

## What it does

```
Raw reads
  → decompress → fastp QC → dedup → vector removal → host removal
  → EsViritu (reference-guided)     ┐
  → metaviralSPADES (de novo)       │
      → filter → CheckV → extract   ┘→ rename → read-supported merge
  → tiered RVDB classification → community profile + final report
```

The **merge step** doesn't just pick the best assembly — for any pair of overlapping contigs it maps the reads back and keeps the version with better support, using a transparent priority order:

1. **Coverage** — if one contig has >2× the coverage, keep it
2. **Mapping quality** — if coverage is similar but MAPQ differs by >10, keep the higher MAPQ
3. **Length** — if both are similar, keep the longer contig

Contigs unique to one assembler (no overlap) are kept as singletons.

The **classification step** is tiered to avoid using one viral DB to classify contigs derived from it:

1. **BLASTn** every hybrid contig against clustered nucleotide RVDB (C-RVDB).
2. Contigs whose best hit clears the threshold (default ≥90% identity AND ≥50% query coverage) are annotated from that hit.
3. Contigs that fail fall through to **BLASTx** against protein RVDB (RVDB-prot) — this catches divergent/novel viruses that nucleotide search misses.
4. Both tiers join to the **RVDB nucleotide taxonomy table** by accession, so every contig gets a full lineage (realm → species) and a `tier` label showing how it was called.

---

## Requirements

- Conda/Mamba (Miniforge3 recommended)
- Snakemake ≥ 8.0
- Slurm (for cluster execution)
- Databases (see **Database setup** below):
  - **EsViritu DB** — directory with a reference FASTA and `virus_pathogen_database.all_metadata.tsv`
  - **CheckV DB** — `checkv-db-v1.5` (only needed if `run_checkv: true`)
  - **RVDB** — clustered nucleotide FASTA, protein FASTA, and the nucleotide taxonomy table
  - **Vector** and **host** reference FASTAs for read removal

Conda environments are defined in `workflow/envs/` (`preprocess.yaml`, `esviritu.yaml`, `hyb_asm.yaml`, `smk_hyb_asm.yaml`) and are built/cached automatically on first run.

---

## Database setup (do this ONCE, before running)

Most databases are ready-to-use once downloaded, but **RVDB must be turned into BLAST databases before the pipeline can use it.** This is a one-time step per RVDB release, not something the pipeline does per run.

### 1. Download RVDB (version-match all three)

| File | Source | Used for |
|------|--------|----------|
| `C-RVDBv32.0.fasta` (clustered nucleotide) | rvdb.dbi.udel.edu | BLASTn tier |
| `U-RVDBv32.0-prot_unique` (clustered protein) | rvdb-prot.pasteur.fr | BLASTx tier |
| `RVDBv32_taxonomy.tab` (nucleotide taxonomy) | rvdb.dbi.udel.edu | lineage for both tiers |

Keep the version numbers aligned (the protein release follows the nucleotide release number).

### 2. Build the BLAST databases

Run the provided helper once (edit the paths at the top first so the `-out` prefixes match your `config.yaml`):

```bash
bash scripts/build_rvdb.sh
```

This runs `makeblastdb` on the two FASTAs (nucleotide + protein). It does **not** use `-parse_seqids` — RVDB headers contain multiple pipes that break seqid parsing, and the pipeline reads accessions from the subject title instead. The taxonomy `.tab` file needs no build step.

### 3. Point the config at the built databases

The `rvdb_nucl_db` / `rvdb_prot_db` values in `config.yaml` must match the `-out` prefixes from step 2 (i.e. the path BLAST looks for `.ndb`/`.pdb` etc. next to), **not** the raw FASTA files:

```yaml
rvdb_nucl_db: "/data/.../RVDB/C-RVDBv32.0"       # prefix, index files sit alongside
rvdb_prot_db: "/data/.../RVDB/RVDBv32-prot"
rvdb_taxonomy: "/data/.../RVDB/RVDBv32_taxonomy.tab"
```

**Quick check that the build worked** — you should see index files, not just FASTAs:

```bash
ls /data/.../RVDB/         # expect C-RVDBv32.0.n*  and  RVDBv32-prot.p*
```

If the pipeline errors with "No alignment definition files found," the DBs weren't built or the prefix doesn't match.

### Other databases

- **CheckV** — download `checkv-db-v1.5` from https://portal.nersc.gov/CheckV/ and point `checkv_db` at the directory. Only needed if `run_checkv: true`.
- **EsViritu DB** — must contain `virus_pathogen_database.all_metadata.tsv` (used to rename consensus contigs).

---

## Repository layout

```
smk_hyb_asm/
├── config/
│   └── config.yaml              # all user settings
├── workflow/
│   ├── Snakefile
│   ├── envs/                    # conda environments
│   └── profiles/slurm/          # cluster profile
├── scripts/                     # helper scripts
│   ├── build_rvdb.sh            # one-time RVDB BLAST DB build
│   ├── merge_with_reads.py      # read-supported assembly merge
│   ├── select_nucl_misses.py    # picks contigs for the protein tier
│   └── annotate_rvdb_tiered.py  # merges both BLAST tiers + taxonomy
└── README.md
```

---

## Sample sheet

Tab-separated, with these columns:

| Sample_ID | forward_read | reverse_read | read_compression | category |
|-----------|--------------|--------------|------------------|----------|
| sample001 | /path/R1.fastq.gz | /path/R2.fastq.gz | gz | sample |
| control01 | /path/R1.fastq.bz2 | /path/R2.fastq.bz2 | bz2 | control |

- `read_compression` is `gz` or `bz2`
- `category` is `sample` or `control`

---

## Configuration

Edit `config/config.yaml`. The essentials:

```yaml
pool_id: "p2006"
sample_sheet: "/path/to/sample_sheet.tsv"
esviritu_db: "/path/to/esviritu_db"
vector_fasta: "/path/to/vector_seqs.fna"
host_fasta: "/path/to/host_genome.fa"
output_root: "/path/to/results"
checkv_db: "/path/to/checkv-db-v1.5"

rvdb_nucl_db: "/path/to/RVDB/C-RVDBv32.0"
rvdb_prot_db: "/path/to/RVDB/RVDBv32-prot"
rvdb_taxonomy: "/path/to/RVDB/RVDBv32_taxonomy.tab"
rvdb_min_identity: 90
rvdb_min_coverage: 50
```

### Pipeline stage flags

Stages are turned on/off independently:

```yaml
run_qc: true                # preprocessing (always)
run_esviritu: true          # reference-guided assembly
run_metaviralspades: false  # de novo assembly
run_checkv: false           # CheckV quality filtering (needs run_metaviralspades)
run_pre_merge: false        # rename + pre-merge stats
run_hybrid: false           # read-supported merge + RVDB classification
```

**Typical hybrid run, in two phases:**

1. Examine the assemblies first — set everything through `run_pre_merge: true`, leave `run_hybrid: false`, then review `*_pre_merge_stats.txt` and the CheckV summaries.
2. Flip `run_hybrid: true` and rerun; Snakemake skips completed steps and only does the merge + classification.

---

## Running

```bash
# Dry run first — always
snakemake --snakefile workflow/Snakefile --configfile config/config.yaml --profile workflow/profiles/slurm -n

# Full run
snakemake --snakefile workflow/Snakefile --configfile config/config.yaml --profile workflow/profiles/slurm

# One sample only
snakemake --snakefile workflow/Snakefile --configfile config/config.yaml --profile workflow/profiles/slurm \
    results/hybrid_assembly/sample/SAMPLE_ID/hybrid/SAMPLE_ID_hybrid.fasta

# Force a rule (e.g. after changing assembly settings)
snakemake --snakefile workflow/Snakefile --configfile config/config.yaml --profile workflow/profiles/slurm --forcerun metaviralspades
```

---

## Key outputs (per sample)

Under `results/hybrid_assembly/<category>/<unique_id>/`:

| File | What it is |
|------|------------|
| `hybrid/<id>_hybrid.fasta` | Final merged assembly |
| `hybrid/<id>_merge_decisions.tsv` | Why each contig was kept/discarded |
| `community/<id>_community_summary.tsv` | Viral composition, grouped by family + species |
| `community/<id>_per_contig_identity.tsv` | Per-contig lineage, tier, % identity, coverage |
| `community/<id>_blast_nucl.tsv` | Raw BLASTn (nucleotide tier) hits |
| `community/<id>_blast_prot.tsv` | Raw BLASTx (protein tier) hits |
| `stats/<id>_assembly_stats.txt` | Human-readable summary report |

The **`tier` column** in the per-contig table shows how each contig was classified: `nucleotide` (confident BLASTn hit), `protein` (fell through to BLASTx), `nucleotide_weak` (sub-threshold nucleotide hit, no protein hit), or `none` (no hit). Watch it on your first samples — lots of `nucleotide_weak`/`none` means either the 90/50 threshold is too strict or your contigs are genuinely divergent.

EsViritu results land in `results/esviritu/<category>/esviritu_<pool_id>/`, and metaviralSPADES/CheckV outputs in `results/metaviralspades/<category>/<unique_id>/`.

---

## Notes & gotchas

- **RVDB databases must be built before running** — see Database setup. `-db` needs a `makeblastdb` prefix, not a FASTA.
- **metaviralSPADES needs an empty output directory.** On a rerun, clear the sample's `metaviralspades/.../<unique_id>/` folder first (or let `--forcerun` handle it).
- **Empty assemblies are handled gracefully.** If EsViritu or metaviralSPADES produces nothing, the merge step falls back to whichever assembly has contigs and records that in the decisions report; classification writes empty (header-only) outputs.
- **EsViritu temp directories are cleaned up** automatically once the consensus FASTA is written (BAM is compressed, everything else removed) to save space.

---

## Troubleshooting

- **`No alignment definition files found`** — the RVDB BLAST databases aren't built, or `rvdb_nucl_db`/`rvdb_prot_db` point at a FASTA instead of the `makeblastdb` prefix. Re-run `scripts/build_rvdb.sh` and check the config prefixes.
- **Jobs run locally instead of on Slurm** — make sure you passed `--profile workflow/profiles/slurm`, and that the profile's `cluster-generic-submit-cmd` is a single folded string (`>-`), not split across indented lines.
- **`PermissionError: /Users/...`** — the `conda-prefix` in the Slurm profile points at a local path; set it to a writable cluster path.
- **`Sample sheet missing required columns`** — the sample sheet must be a TSV with the header `Sample_ID  forward_read  reverse_read  read_compression  category`.
- **Fragmented de novo assembly** — for high-coverage samples with strain diversity, metaviralSPADES may fragment where EsViritu doesn't. That's expected; the hybrid merge leans on the reference-guided assembly for those contigs.
- **Check logs** in `results/logs/<rule>/<category>/<id>.log` and the Slurm `logs/<rule>/` directory.
