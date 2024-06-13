# pathSTR-1000G

Repository for the pathSTR web app at <https://pathstr.bioinf.be>

If pathSTR is useful to you, please cite our [publication](https://www.medrxiv.org/content/10.1101/2024.03.06.24303700v1).

## Installation

This repository contains the code for the pathSTR web app. For a local installation, you will need to install the dependencies as specified in requirements.txt. You can do this by running the following command in the root directory of the repository:

```bash
pip install -r requirements.txt
```

The pathSTR-1000G.smk workflow uses snakemake to manage the pipeline. You can install snakemake using pip:

```bash
pip install snakemake
```

At the top of the file some paths have to be set, including reference genomes, the work_dir and the location of the STRdust binary.
This was not developed into a more convenient configuration file yet, as I do not anticipate many users running this pipeline. It is provided for completeness and transparency.

The pipeline is ran as:

```bash
snakemake -s workflow/pathSTR-1000G.smk --use-conda --cores 24 --keep-going
```

## Listing all files from the AWS bucket, subsetting to hg38 bams from the standard minimap2 pipeline

```bash
aws s3 ls 1000g-ont --no-sign-request --recursive | grep -i minimap | grep bam$ | grep -v chm13 | cut -f4 -d' ' > data/list-miller-20240606.txt
```

--> results in 157 paths to bam files, without duplicates
