# Cancer phylogenetics using single-cell RNA-seq data

This repository contains code to fully replicate the analysis of *Cancer phylogenetics using single-cell RNA-seq data* (Moravec at al. 2021). Alternatively, it can be used to perform a similar analysis on a new dataset.

Note that the analysis assumes a relatively uniform cell populations, otherwise the discretization method using Highest Density Interval will not work.

## Preparation for analysis
The `run.r` file describes the analysis and parameter settings required to replicate analysis. Please read it before start to know about the necessary files and their locations.

Note that that the analysis was tested only on Linux.
Windows and Mac are not supported by required software (e.g., Cellranger).

### Data files
Place the BAM files from Moravec et al. into the `data` folder.

Download the human reference genome GRch38v15 (no alt analysis set as non-uniquelly mapped reads are ignored) and common variants and place is in respective folders or modify the `reference`, `annotation` and `vcf` lines in the `run.r` file.

### R packages:
Following R packages are required:
* `phyloRNA`
* `beter`
* `tools`
* `data.table`
* `devtools`

The `tools` package should come in a standard R installation, `data.table` and `devtools` are available on CRAN.

Packages `phyloRNA` and `beter` are internally developed and available on the `biods` github. They are not yet available on CRAN so to install them, type:

```{r}
devtools::install_github("biods/phyloRNA")
devtools::install_github("biods/beter")
```
or just run the analysis using:
```{R}
Rscript run.r
```
as the latest version is automatically installed in the pipeline.

### External Software
Following external software is required:

* `Cellranger`
* `bamtofastq`
* `GATK`
* `VCFtools`
* `python3`
* `IQtree`
* `BEAST2`

We suggest the [conda](https://docs.conda.io/en/latest/) package manager to painfully install these.
See the [phyloRNA](https://github.com/bioDS/phyloRNA) package for more information.

## Analysis
Navigate to the analysis directory and type:
```{R}
Rscript run.r
```
