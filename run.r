#' run.r
#'
#' Run the analysis. This includes:
#'
#' * preparation of sample files
#' * SNV detection and their evaluation
#' * Expression analysis
#' * Phylogenetic reconstruction
#'
#' Analyses:
#' --------------------
#' datasets: all, 2HR, foreign (TODO)
#'
#' Expression:
#' -- density: 0.2, 0.5, 0.9
#' -- filtered and unfiltered (to show that they are the same for a higher density)
#' -- categorization: According to 60-30-10 empirical intervals
#' -- filter constant sites
#' -- IQtree: ORDINAL+ASC; ultrafast bootstrap -B 1000
#' -- BEAST: ordinal from MM, exponential pop growth, coalescent prior, strict clock,
#'    2 runs but how many MCMC gen?
#' -- BEAST template prepared with the `beter` package
#'
#' SNV:
#' -- density: 0.2, 0.5, 0.9
#' -- Possible 
#' -- filter constant sites
#' -- IQtree: GTR+gamma, ultrafast bootstrap -B 1000
#' -- BEAST: GTR, exponential pop growth, coalescent prior, strict clock
#' -- BEAST template prepared with the `beter` package
#'
devtools::install_github("biods/phyloRNA") # requires a github authentication token in .Renviron
devtools::install_github("j-moravec/baffle")
devtools::install_github("biods/beter")

library("phyloRNA")
library("baffle")
library("beter")

source("src/prepare.r")
source("src/snv.r")


# TODO:
# * define datasets
# * define which analysis will be done
# * separate so that there is no name clash


# datasets:
bams = dir("data", full.names=TRUE)

# required reference files:
reference = "/data/phylonco/ReferenceGenomes/human_GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
annotation = "/data/phylonco/ReferenceGenomes/human_GRCh38/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gtf"
vcf = "/data/phylonco/ReferenceGenomes/vcf/00-common_all.vcf.gz"

# Other settings:
nthreads = 16
chemistry = "SC5P-R2"


# Preparation step:
prepared = prepare_samples(
    bams, reference, annotation, vcf,
    chemistry = chemistry,  nthreads=nthreads
    )

# SNV detection step:
phyloRNA::mkdir("snv/all")
filtered_barcodes = "snv/all/all.barcodes.filtered.txt" # filtered barcodes to save time
alignment = detect_snv(prepared$bam, filtered_barcodes, reference, outdir=file.path("snv/all"))
