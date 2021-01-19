# Cancer phylogenetics using single-cell RNA-seq data

This repository contains code to fully replicate the analysis of *Cancer phylogenetics using single-cell RNA-seq data* [(Moravec at al. 2021)](https://doi.org/10.1101/2021.01.07.425804). Alternatively, it can be used to perform a similar analysis on a new dataset.

Note that the analysis assumes a relatively uniform cell populations, otherwise the discretization method using Highest Density Interval will not work.

## Requirements

### System requirements
* Linux operation system
* at least 30 GB RAM
* about 400 GB of free space for intermediate files and results

### Required software
[R](https://cran.r-project.org/), [python3](https://www.python.org/), [Cellranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger), [bamtofastq](https://github.com/10XGenomics/bamtofastq), [GATK](https://gatk.broadinstitute.org/hc/en-us), [VCFtools](https://vcftools.github.io/), [IQtree](http://www.iqtree.org/), [BEAST2](https://www.beast2.org/)

#### R packages:
[phyloRNA](https://github.com/bioDS/phyloRNA), [beter](https://github.com/bioDS/beter), [data.table](https://github.com/Rdatatable/data.table), [devtools](https://github.com/r-lib/devtools)

#### Python packages:
[pysam](https://pysam.readthedocs.io/en/latest/index.html)

### Required files:
Original data published at GEO database under the accession number [GSE163210](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE163210).

Human reference genome [GRCh38v15](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz), [annotation](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gtf.gz) and [known variants](https://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/00-common_all.vcf.gz).

Code from this repository.

## Running the analysis
Once you have installed required software and prepared your data, navigate into the analysis directory and type:

```{R}
Rscript run.r
```

After few days, the analysis should finish.

## Processed files
Pre-processed fasta files, trees and tests of phylogenetic clustering can be seen in the `processed_files` branch. These files are tracked with *Git Large File Storage* (LFS) extension.

## Detailed instruction
* install [required software](doc/required_software.md)
* download R and Python [packages](doc/packages.md)
* download [data and reference genome](required_files.md)
* [run](doc/run_analysis.md) the analysis

## Need help?
If anything is unclear or you need help with the analysis, raise an [issue](https://github.com/bioDS/phyloRNAanalysis/issues).
