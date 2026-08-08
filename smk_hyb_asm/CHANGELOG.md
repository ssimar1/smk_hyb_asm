# Changelog

All notable changes to this project will be documented in this file.

## [0.1.1] - 2026-01-28

### Added
- Modification of smk_tephi_virome pipeline for GCID data
- **Preprocessing modules**:
  - `fastp_dedup`: Remove PCR duplicates with fastp
  - `host_removal`: Remove human host reads using minimap2
### Removed
- Modules for IAV, mpox, SARS-CoV-2 typing
## Modified
- EsViritu output paths-- individual samples get their own output folder
- Memory allocation/# jobs in slurm profile
- Unique_IDs are now just Sample_IDs not Sample_ID + Pool_ID


## [0.1.0] - 2025-11-29

### Added
- Initial release of the Simple Snakemake Virome Pipeline
- **Preprocessing modules**:
  - `decompress_reads`: Multi-threaded decompression of `.gz`/`.bz2` files using `pigz`/`pbzip2`
  - `fastp_qc`: Quality filtering with fastp
  - `vector_read_removal`: Vector/adapter read removal using `vrr`
- **Viral detection**:
  - `esviritu`: Broad viral detection, taxonomic profiling, and coverage analysis
- **Virus-specific modules**:
  - `iav_serotype`: Influenza A virus serotype classification
  - `mpoxclade`: Monkeypox virus clade classification
  - `sars_cov2_reads`: SARS-CoV-2 read extraction from EsViritu alignments
- SLURM cluster support via `snakemake-executor-plugin-cluster-generic`
- Conda environment definitions for all modules
- Comprehensive documentation in README.md

**Further Updates in GitHub Release Notes**
