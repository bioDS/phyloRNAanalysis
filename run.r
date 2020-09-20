#' run.r
#'
#' Run the analysis. This includes:
#'
#' * preparation of sample files
#' * SNV detection and their evaluation
#' * Expression analysis
#' * Phylogenetic reconstruction

devtools::install_github("biods/phyloRNA") # requires a github authentication token in .Renviron
devtools::install_github("j-moravec/baffle")

library("phyloRNA")
library("baffle")

source("src/prepare.r")


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
phyloRNA::mkdir("snv")
alignment = phyloRNA::gatk_snv(prepared$bam, reference, "all.snv", outdir=file.path("snv, snv"))
