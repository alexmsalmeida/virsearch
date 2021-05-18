# virsearch

Snakemake workflow to detect and classify viral sequences in metagenome assemblies.

It first detects viral sequences in assemblies (.fa files) with <b>VirSorter2</b>, <b>VIBRANT</b> and <b>DeepVirFinder</b>. Predictions are subsequently quality controlled with <b>CheckV</b>, followed by clustering with <b>CD-HIT</b> and taxonomic classification with <b>DemoVir</b>.

## Installation

1. Clone repository
```
git clone
```

3. Download and extract necessary databases

## How to run

1. Edit config.yml file to point to the input and output directories. Input directories should contain the *.fa assemblies to analyse.
2.1. Run the pipeline locally
2.2. Run the pipeline on a cluster (e.g. LSF)

