# Metagenomic Analysis Pipeline/Code
Here you can find all the scripts I used during my bachelor’s thesis with the title “Benchmarking of a targeted viral enrichment panel in human DNA-samples” to analyse sequencing reads. The workflow can be divided into three main steps: host removal, k-mer based classification, alignment based classification. The output from host-removal is used as input for both the k-mer based and the alignment based classification. Results from the k-mer based classification can be used to decide which reference sequences to use in the alignment based pipeline. 

Scripts I used for plotting have been included. Please note, however, that these scripts are heavily personalized for my results and will likely have to be edited for your use case.

## Tools used:
bwa-mem2 (https://github.com/bwa-mem2/bwa-mem2 ; DOI: 10.1109/IPDPS.2019.00041)
nextflow (https://github.com/nextflow-io/nextflow ; DOI: 10.1038/nbt.3820)
nf-core/taxprofiler (https://github.com/nf-core/taxprofiler/tree/main ; DOI-1: 10.1101/2023.10.20.563221 ; DOI-2: 10.5281/ZENODO.7728364)
fastp (https://github.com/opengene/fastp ; DOI: 10.1093/bioinformatics/bty560)
Kraken2 (https://github.com/DerrickWood/kraken2 ; DOI: ​​10.1186/s13059-019-1891-0)
Bracken (https://github.com/jenniferlu717/Bracken ; DOI: 10.7717/peerj-cs.104)
MultiQC (https://github.com/multiqc/multiqc ; DOI: 10.1093/bioinformatics/btw354)
samtools (https://github.com/samtools/samtools ; DOI: 10.1093/gigascience/giab008)

## License
The code provided within this repository is provided under the MIT License.

