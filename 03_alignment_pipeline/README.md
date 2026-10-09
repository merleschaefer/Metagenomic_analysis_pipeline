# Metagenomic virus-alignment Snakemake pipeline

This is an alignment-based pipeline for the metagenomic analyses of short-read (2x150bp) Genome Sequencing data. It aligns FastQ files against a reference FASTA provided by you, filters out unreliable reads based on samtools filters (adjustable in config.yaml), flags viral genomic regions (window size adjustable in config.yaml) that have extraordinarily high coverage depth compared to the rest of the respective genome (sign of a false-positive) and removes reads mapping to these regions from further analyses. 
A plot showing the coverage depth per viral genome as well as flaggged regions is included in the results. As final steps insert sizes are extracted and provided in a text file and reads are counted per contig and virus once with duplicates included (reads_with_dups) and once with duplicates removed (reads_without_dups). Final and intermediate BAM files are provided in the results. It is important to note, however, that they all still include duplicates due to this last QC steps in read counting. For future use of BAM files duplicates can be removed using the samtools flag -F 1024 if needed. 

## What you must supply

- Host-removed paired FASTQs, `{sample}_hostremoved_1.fq.gz` / `_2.fq.gz` (you can use the hostremoval script provided for this)

- Virus reference FASTA (manually created based on Bracken results and expected viruses, fastasets used for the thesis included here under /data as examples) - important: make sure fasta headers are space seperated (">accession additional information"). Anything between ">" and the first " " will be used as accession by samtools.

- `identifiers.txt` - Text file with one accession/contig ID per line, the same accession provided in your fasta file. Examples from the thesis are included here under /data.

- `contig_name.csv` - used for plotting only to connect virus names to accession codes, needs at least `identifier,species` columns. Examples from the thesis are included here under /data.

Set all four paths in `config.yaml`. Samples are auto-detected from whatever
`*_hostremoved_1.fq.gz` files are in `fq_dir` — you don't list them manually.

The FASTA files as well as identifiers.txt and contig_name.csv used in my thesis are provided under /data. They consist of RefSeq sequences of selected viruses (e.g. herpesviruses) and the SCANellome V2 FASTA database (https://github.com/Laubscher/Anelloviruses; doi: 10.3390/v15071575)

## Software needed

Easiest: build the conda/mamba environment included here:

```bash
mamba env create -f workflow/environment.yml
conda activate metagenomic-alignment
```

It installs: `snakemake`, `bwa-mem2`, `samtools` (≥1.15, needed for the `-e`
expression filter syntax), `bedtools`, and R + `dplyr`/`readr`/`ggplot2`.
Standard Unix tools (`awk`, `cut`, `sort`, `comm`, `zcat`) are assumed to
already be on your system. They are on any normal Linux/WSL install.

## Running it

```bash
snakemake -n                       # dry run 
snakemake --cores 4 --use-conda    # actually run it
```

`--cores` should be ≥ `threads` in config.yaml if you want a single sample's
alignment to get all its threads while still running independent rules for
other samples in parallel.


