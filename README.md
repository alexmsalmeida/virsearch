# VirSearch - searching viral sequences in metagenomes

[Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) workflow to detect and classify viruses in metagenome assemblies.

It first detects viral sequences in assemblies (`.fa` files) with [VirSorter2](https://github.com/jiarong/VirSorter2), [VIBRANT](https://github.com/AnantharamanLab/VIBRANT) and [DeepVirFinder](https://github.com/jessieren/DeepVirFinder). Predictions are strictly quality controlled with [CheckV](https://bitbucket.org/berkeleylab/checkv), followed by clustering with [CD-HIT](http://weizhongli-lab.org/cd-hit/) and taxonomic classification with [Demovir](https://github.com/feargalr/Demovir).

## Installation

1. Install [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) and [snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)

2. Clone repository
```
git clone https://github.com/alexmsalmeida/virsearch.git
```

3. Download and extract necessary databases (uncompressed directory will require a total of 30 GB).

```
cd virsearch
wget http://ftp.ebi.ac.uk/pub/databases/metagenomics/genome_sets/viral_databases.tar.gz
tar -xzvf viral_databases.tar.gz
```

## How to run

1. Edit `config.yml` file to point to the <b>input</b> and <b>output</b> directories. Input directory should contain the `.fa` assemblies to analyse.

2. (option 1) Run the pipeline locally (adjust `-j` based on the number of available cores)
```
snakemake --use-conda -j 4
```
2. (option 2) Run the pipeline on a cluster (e.g., LSF)
```
snakemake --use-conda -j 50 --cluster-config cluster.yml --cluster "bsub -n {cluster.nCPU} -M {cluster.mem} -o {cluster.output}"
```

## Output

The main output files generated per input FASTA are the `final_predictions.fa` and `final_predictions_tax.tsv` files, which contain the viral sequences in FASTA format and their taxonomic annotation, respectively.
