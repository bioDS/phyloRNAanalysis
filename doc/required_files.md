# Required files

## Code from this repository
To clone the repository, type into your terminal:
```{bash}
git clone https://github.com/bioDS/phyloRNAanalysis
```
If you do not have `git` installed, type first:
```{bash}
conda install git
```

Now type:
```{bash}
cd phyloRNAanalysis
```
to enter the `phyloRNAanalysis` folder.

## Reference genome
First, create `reference` folder in the `phyloRNAanalysis` directory:

```{bash}
mkdir reference
cd reference
```
Now, download and unzip the [GRCh38v15](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz) reference genome together with the associated index file:
```{bash}
wget ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
wget ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai
gunzip GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
```
download and unzip the [annotation](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gtf.gz):

```{bash}
wget ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gtf.gz
gunzip GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gtf.gz
```
as well as set of common variants (no need to unzip these):
```{bash}
wget ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606/VCF/00-All.vcf.gz
```

### Filter annotation file:
There is just one more thing. The newest version of `Cellranger` (version 5) does not accept an annotation file with extra scaffolds. This can be fixed by filtering the GTF file according to the fasta index file using the `phyloRNA` package:
```{bash}
Rscript "phyloRNA::filter_gtf(GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gtf, GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai)"
```

## Data
To download the [GSE163210](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE163210) data, navigate back to the `phyloRNAanalysis`, create a `data` directory and download the BAM files:
```{bash}
cd ..
mkdir data
```
