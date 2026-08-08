# smk_hyb_asm

A Snakemake pipeline for hybrid capture virome sequencing. It preprocesses reads, detects viruses with **EsViritu** (reference-guided), assembles de novo with **SPAdes**, and merges the two assemblies into a single hybrid assembly using **read support** as the deciding evidence.

**Version:** 0.2.0

---

## What it does

```
Raw reads
  → decompress → fastp QC → dedup → vector removal → host removal
  → (optional downsample)
  → EsViritu (reference-guided)  ┐
  → SPAdes (de novo) → CheckV    ┘→ rename → read-supported merge
  → BLAST vs NCBI Viral RefSeq → community profile + final report
```

The merge step doesn't just pick the best assembly — for any pair of overlapping contigs it maps the reads back and keeps the version with better support, using a transparent priority order:

1. **Coverage** — if one contig has >2× the coverage, keep it
2. **Mapping quality** — if coverage is similar but MAPQ differs by >10, keep the higher MAPQ
3. **Length** — if both are similar, keep the longer contig

Contigs unique to one assembler (no overlap) are kept as singletons.

---

## Requirements

- Conda/Mamba (Miniforge3 recommended)
- Snakemake ≥ 8.0
- Slurm (for cluster execution)
- Databases:
  - **EsViritu DB** — directory with a reference FASTA and `virus_pathogen_database.all_metadata.tsv`
  - **CheckV DB** — `checkv-db-v1.5` (only needed if `run_checkv: true`)
  - **NCBI Viral RefSeq BLAST DB** — used for community classification of the final assembly
  - **Vector** and **host** reference FASTAs for read removal

Conda environments are defined in `workflow/envs/` (`preprocess.yaml`, `esviritu.yaml`, `hyb_asm.yaml`) and are built/cached automatically on first run.

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
├── scripts/                     # helper scripts called by rules
│   ├── extract_quality_contigs.py
│   ├── rename_esviritu.py
│   ├── merge_with_reads.py
│   ├── annotate_refseq_blast.py
│   └── final_stats.py
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
pool_id: "p1389"
sample_sheet: "/path/to/sample_sheet.tsv"
esviritu_db: "/path/to/esviritu_db"
vector_fasta: "/path/to/vector_seqs.fna"
host_fasta: "/path/to/host_genome.fa"
output_root: "/path/to/results"
checkv_db: "/path/to/checkv-db-v1.5"
```

### Pipeline stage flags

Stages are turned on/off independently:

```yaml
run_qc: true            # preprocessing (always)
run_esviritu: true      # reference-guided assembly
run_SPADES: false       # de novo assembly
run_checkv: false       # CheckV quality filtering (needs run_SPADES)
run_pre_merge: false    # rename + pre-merge stats
run_hybrid: false       # read-supported merge + classification

run_downsample: false   # optionally cap reads before assembly
max_reads: 5000000
```

**Typical hybrid run, in two phases:**

1. Examine the assemblies first — set everything through `run_pre_merge: true`, leave `run_hybrid: false`, then review `*_pre_merge_stats.txt` and the CheckV summaries.
2. Flip `run_hybrid: true` and rerun; Snakemake skips completed steps and only does the merge + classification.

---

## Running

```bash
# Dry run first — always
snakemake --profile workflow/profiles/slurm -n

# Full run
snakemake --profile workflow/profiles/slurm

# One sample only
snakemake --profile workflow/profiles/slurm \
    results/hybrid_assembly/sample/SAMPLE_ID/hybrid/SAMPLE_ID_hybrid.fasta

# Force a rule (e.g. after changing SPAdes settings)
snakemake --profile workflow/profiles/slurm --forcerun spades
```

---

## Key outputs (per sample)

Under `results/hybrid_assembly/<category>/<unique_id>/`:

| File | What it is |
|------|------------|
| `hybrid/<id>_hybrid.fasta` | Final merged assembly |
| `hybrid/<id>_merge_decisions.tsv` | Why each contig was kept/discarded |
| `community/<id>_community_summary.tsv` | Viral composition (aggregated by virus) |
| `community/<id>_per_contig_identity.tsv` | Per-contig % identity to RefSeq |
| `stats/<id>_assembly_stats.txt` | Human-readable summary report |

EsViritu results land in `results/esviritu/<category>/esviritu_<pool_id>/`, and SPAdes/CheckV outputs in `results/SPADES/<category>/<unique_id>/`.

---

## Notes & gotchas

- **SPAdes needs an empty output directory.** If you rerun, clear the sample's `SPADES/.../<unique_id>/` folder first (or let Snakemake's `--forcerun` handle it).
- **`--careful` is incompatible with metagenomic/`--meta` mode** and with `--cov-cutoff`. The pipeline uses `spades.py --careful` with an explicit k-mer ladder.
- **Empty assemblies are handled gracefully.** If EsViritu or SPAdes produces nothing, the merge step falls back to whichever assembly has contigs and records that in the decisions report.
- **EsViritu temp directories are cleaned up** automatically once the consensus FASTA is written, to save space.

---

## Troubleshooting

- **Jobs run locally instead of on Slurm** — make sure you passed `--profile workflow/profiles/slurm`, and that the profile's `cluster-generic-submit-cmd` is a single folded string (`>-`), not split across indented lines.
- **`PermissionError: /Users/...`** — the `conda-prefix` in the Slurm profile is pointing at a local path; set it to a writable cluster path.
- **Fragmented de novo assembly** — for high-coverage samples with strain diversity, SPAdes may fragment where EsViritu doesn't. That's expected; the hybrid merge will lean on the reference-guided assembly for those contigs.
- **Check logs** in `results/logs/<rule>/<category>/<id>.log` and the Slurm `logs/<rule>/` directory.
