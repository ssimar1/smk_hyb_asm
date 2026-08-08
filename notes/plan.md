# Goal
Make a snakemake pipeline to decompress hybrid capture virome sequencing reads from a sample sheet and run EsViritu.

## Design notes
- use snakemake best practices
	- See this repo for style guide: notes/style/StainedGlass
- this will be run on an HPC with SLURM job scheduling
	- use the conventions in the repo: notes/style/smk-simple-slurm
		- especially note this config: notes/style/smk-simple-slurm/simple/config.v8+.yaml

## Pipeline inputs
- sample sheet
	- required fields
		- Sample_ID
		- path to forward read file (compressed)
		- path to reverse  read file (compressed)
		- compression type
		- category (sample or control)
- config file arguments
	- PoolID

# Modules

## Preprocessing
### Decompressing files
- inputs
	- sample sheet
- config
	- temp directory
	- threads: 8
	- memory 8g
	- conda env: preprocess
- actions
	- read sample sheet
		- determine if reads are compressed in .gz or .bz2 format
		- merge sample_ID and PoolID strings to generate unique_ID
	- use `pigz` (multithread) or `bunzip2` to uncompress reads to temp directory
- outputs
	- uncompressed read files path (tempfile)
	- example output path: r1=str(TEMP_ROOT / "{category}" / "01_decompress" / "{unique_id}.R1.fastq")

### Quality filter reads
- inputs
	- uncompressed read files path
- config
	- temp directory
	- threads: 16
	- memory 16g
	- conda env: preprocess
- actions
	- use `fastp` with default settings to quality filter reads
- outputs
	- quality filtered read files path (tempfile)

### remove reads aligning to vector sequences
- inputs
	- quality filtered read files path
- config
	- temp directory
	- samples output directory
	- controls output directory
	- threads: 8
	- memory 8g
	- conda env: preprocess
- actions
	- use `vrr` (from vector_read_removal repo) with `--nm-threshold 4`
- outputs
	- vector-removed read files path (tempfile)

## Virome analysis

### EsViritu
- inputs
	- vector-removed read files path (tempfile)
- config
	- esviritu DB path
	- threads: 8
	- memory 8g
	- conda env: esviritu
- actions
	- use `EsViritu` (from https://github.com/cmmr/EsViritu) with  `-q T --keep T` flags
- outputs
	- {unique_ID}.detected_virus.info.tsv
	- {unique_ID}.detected_virus.assembly_summary.tsv
	- {unique_ID}.tax_profile.tsv
	- {unique_ID}.virus_coverage_windows.tsv
	- example output path: r1=str(ESVIRITU_ROOT / "{category}" / "02_esviritu_{PoolID}" / "{unique_id}.{extension}")

## Virus-specific analysis

### iav_serotype
- inputs
	- vector-removed read files path (tempfile)
- config
	- samples output directory
	- controls output directory
	- iav_serotype DB path
	- threads: 8
	- memory 8g
	- conda env: iav_serotype
- actions
	- use `iav_serotype` (from https://github.com/mtisza1/influenza_a_serotype) with default settings
- outputs
	- {unique_ID}_per_serotype_summary.tsv **<- Main summary table for serotype counts**
	- {unique_ID}_influenza_A.sorted.bam **<- filtered, sorted alignment file**
	- {unique_ID}_per_read_summary.tsv **<- per-read summary file**
	- {unique_ID}_read_serotype_assignment.pdf **<- plot of serotype counts**
	- {unique_ID}_{serotype}.txt **<- serotype-specific read IDs**
	- {unique_ID}_{serotype}.R1.fastq **<- serotype-specific reads (optional)**
	- {unique_ID}_{serotype}.R2.fastq **<- serotype-specific reads (optional)**
	- {unique_ID}_read_stats.tsv **<- input read stats table**

### mpoxclade
- inputs
	- vector-removed read files path (tempfile)
- config
	- samples output directory
	- controls output directory
	- mpoxclade DB path
	- threads: 8
	- memory 8g
	- conda env: iav_serotype
- actions
	- use `mpoxclade` (from mpoxclade repo) with default settings
- outputs
	- {unique_ID}_per_clade_summary.tsv **<- Main summary table for serotype counts**
	- {unique_ID}_mpox.sorted.bam **<- filtered, sorted alignment file**
	- {unique_ID}_per_read_summary.tsv **<- per-read summary file**
	- {unique_ID}_read_clade_assignment.pdf **<- plot of serotype counts**
	- {unique_ID}_{clade}.txt **<- serotype-specific read IDs**
	- {unique_ID}_{clade}.R1.fastq **<- serotype-specific reads (optional)**
	- {unique_ID}_{clade}.R2.fastq **<- serotype-specific reads (optional)**
	- {unique_ID}_read_stats.tsv **<- input read stats table**

## Saving virus-specific reads

## get sars-cov-2 reads
- inputs
	- EsViritu output from temp directory: {unique_ID}.third.filt.sorted.bam
- config
	- samples output directory
	- controls output directory
	- EsViritu DB path
	- threads: 8
	- memory 8g
	- conda env: EsViritu
- actions
	- python script that 
		- gets list of accessions belonging to species "s__Betacoronavirus pandemicum" from EsViritu metadata table
			- {esviritu_db}/virus_pathogen_database.all_metadata.tsv
		- uses pysam to retrieve any reads with primary alignments to the selected accessions
- outputs
	- {unique_ID}.sars-cov-2.R1.fastq
	- {unique_ID}.sars-cov-2.R2.fastq