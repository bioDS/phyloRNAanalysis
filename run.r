#' run.r
#'
#' Run the analysis. This includes:
#'
#' Preparation
#' * remapping, demultiplexing, barcode correction and expression counts using Cellranger
#' * cleaning BAM files according to the GATK best practices
#' * adding a sample-specific postfix to cell barcodes
#'
#' Pre-processing:
#' * Expression
#'   -- standardization of genes into mu=0 and sd=1
#'   -- categorization according to empirical 60% and 90% HDI
#' * SNV:
#'   -- as bulk SNV identification and filtering with Mutect2
#'   -- sc SNV identification with vcm.py
#' * stepwise filtration into 20%, 50% and 90% density
#' * alternative filtration into 58 best cells and full dataset, 50% and 90% density
#'
#' Filtering:
#' * stepwise filtration into 20%, 50% and 90% density
#' * subset into 58 best cells and filtering into 50% and 90% density
#'
#' Phylogenetic analysis:
#' * ML with stepwise filtering
#' * ML and BI with alternative filtration
#' * BEAST templates created with the `beter` package
#' * Expression:
#'   -- IQtree: ORDINAL+ASC, ultrafast bootstrap -B 1000
#'   -- BEAST: ordinal from MM, exponential pop growth, coalescent prior, strict clock, two runs
#' * SNV:
#'   -- IQtree: GTR+gamma, ultrafast bootstrap -B 1000
#'   -- BEAST: GTR, exponential pop growth, coalescent prior, strict clock, two runs
#'

# requires a github authentication token in .Renviron
devtools::install_github("biods/phyloRNA")
devtools::install_github("biods/beter")

library("phyloRNA")
library("beter")

source("src/utils.r")
source("src/prepare.r")
source("src/filter.r")
source("src/expr.r")
source("src/snv.r")
source("src/iqtree.r")
source("src/beast.r") # -- move this into the `beter` package


# datasets:
bam = dir("data", full.names=TRUE)

# required reference files:
reference = "/data/phylonco/ReferenceGenomes/human_GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
annotation = "/data/phylonco/ReferenceGenomes/human_GRCh38/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gtf"
vcf = "/data/phylonco/ReferenceGenomes/vcf/00-common_all.vcf.gz"

# Other settings:
nthreads = 16
chemistry = "SC3Pv2"
densities = c(0.2, 0.5, 0.9)
hdi = c(0.6, 0.9)
selection = c("T1" = 20, "T3" = 20, "T2" = 6, "CTC1" = 6, "CTC2" = 6)

# Preparation:
prepared = prepare_samples(
    bam, reference, annotation, vcf,
    chemistry=chemistry, nthreads=nthreads
    )

# Expression preprocessing:
expr_preprocessed = preprocess_expression(
    h5 = prepared$h5,
    hdi = hdi,
    minGene=0,
    mminUMI=0,
    outdir = file.path("preprocess", "expr"),
    prefix = "all"
    )


expr_filtered = filter_expression(
    expr_preprocessed$discretized,
    selection = selection,
    density = densities,
    outdir = "filtered"
    )


epxr_fasta = table2fasta(unlist(expr_filtered), outdir="fasta")


# SNV preprocessing:
snv_preprocessed = preprocess_snv(
    bam = prepared$bam,
    barcodes = prepared$barcodes,
    reference = reference,
    outdir = file.path("preprocess", "snv"),
    nthreads = nthreads
    )


snv_filtered = filter_snv(
    vcm = snv_preprocessed,
    selection = selection,
    density = densities,
    outdir = "filtered"
    )


snv_fasta = table2fasta(unlist(snv_filtered), outdir="fasta")


# IQtree model strings:
# Expr: "ORDERED+ASC"
# SNV: "GTR+G+ASC"
