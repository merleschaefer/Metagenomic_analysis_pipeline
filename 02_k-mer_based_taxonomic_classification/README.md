# Scripts for taxonomic classification using the nf-core/taxprofiler

## requirements pipeline
As this is a nextflow-based pipeline you will require nextflow in your environment. For information on how to set this up visit: https://nf-co.re/docs/get_started/environment_setup/nextflow

The pipeline requires a seperate databasesheet.csv as well as a samplesheet.csv.

The databasesheet.csv gives information on the database used and certain parameters used by the classifier. It is provided and can be adjusted according to needs. The minusb database referenced in the databasesheet.csv need to be downloaded (https://genome-idx.s3.amazonaws.com/kraken/k2_minusb_20260626.tar.gz) and kept in this folder.
Parameters used: 
    --confidence 0.2 , at least 20% of a read's k-mers need to match the respective organism for the read to be assigned to it
    -t 1 , Bracken takes all results from Kraken2 into consideration. Default is 10 (only organisms with minimum 10 reads assigned to them are taken into account)

The samplesheet.csv gives information on which samples should be analysed and where they are located. It can be created using the provided script samplesheet_creation.sh.

Once both the databasesheet.csv and the samplesheet.csv are created and the database is downloaded, the pipeline can be run using the script Run_pipeline.sh

## analyses of results
Scripts for the creation of heatmaps for the analyses of Bracken results are included under ./analyses_of_results. Note that both scripts are made specifically to match the results in my thesis. For general usage they need to be adjusted. 
The text file viral_taxids.txt is used to filter out results non-viral results. It was created based on NCBI taxonomy (like the database used in Kraken2) using taxonkit and taxdump and includes all taxids belonging to the "Viruses" superkingdom:

```bash
conda install -c bioconda taxonkit
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
# extract taxdump
mkdir -p ~/.taxonkit
tar -xzf taxdump.tar.gz -C ~/.taxonkit
# generate the viral taxid list
taxonkit list -i 10239 --indent "" > viral_taxids.txt
```

## creation of a combined fasta for alignment-based analyses
Based on the results from taxonomic classification and viruses commonly found in the sample type, contigs can be chosen to be included in the fasta for alignment. For the thesis this included the fasta provided by SCANellome for Anellovirus references, RefSeq sequences of all Human infecting herpesviridae as well as four additional viruses found during classification after the first sequencing run (likely contamination, false positives) and the ReSeq sequence for Monkeypox (synthetic monkeypox DNA was spiked-in during experiments as a positive control).
All individual fastas were downloaded into a shared folder and combined using:

```bash
cat *.fasta > combined.fasta
``` 

