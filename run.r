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
devtools::install_github("biods/phyloRNA") # requires a github authentication token in .Renviron
devtools::install_github("biods/beter")

library("phyloRNA")
library("beter")

source("src/prepare.r")
source("src/expr.r")
source("src/snv.r")
source("src/iqtree.r")



# datasets:
bam_all = dir("data", full.names=TRUE)
bam_2HR = file.path("data", "2HR.bam")

# required reference files:
reference = "/data/phylonco/ReferenceGenomes/human_GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
annotation = "/data/phylonco/ReferenceGenomes/human_GRCh38/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gtf"
vcf = "/data/phylonco/ReferenceGenomes/vcf/00-common_all.vcf.gz"

# Other settings:
nthreads = 16
chemistry = "SC5P-R2"
densities = c(0.2, 0.5, 0.9)
hdi = c(0.6, 0.9)



#######################
# Analyse all samples #
#######################
# Preparation:
prepared_all = prepare_samples(
    bam_all, reference, annotation, vcf,
    chemistry=chemistry, nthreads=nthreads
    )

# Expression no quality filtering:
analyse_expression(
    h5 = prepared_all$h5,
    densities = densities,
    hdi = hdi,
    exprdir = file.path("expr", "all", "no_quality"),
    phylodir = file.path("phylo", "all", "expr", "no_quality"),
    prefix = "all",
    model = "ORDERED+ASC",
    minGene = 0,
    minUMI = 0,
    nthreads = nthreads
    )

# Expression with quality filtering:
analyse_expression(
    h5 = prepared_all$h5,
    densities = densities,
    hdi = hdi,
    exprdir = file.path("expr", "all", "quality"),
    phylodir = file.path("phylo", "all", "expr",  "quality"),
    prefix = "all.quality",
    model = "ORDERED+ASC",
    nthreads = nthreads
    )
# snv:
analyse_snv(
    bam = prepared_all$bam,
    barcodes = prepared_all$barcodes,
    reference = reference,
    densities = densities,
    snvdir = file.path("snv", "all"),
    phylodir = file.path("phylo", "all", "snv"),
    prefix = "all",
    model = "GTR+G+ASC",
    nthreads = nthreads
    )

######################
# Analyse 2HR sample #
######################
prepared_2HR = prepare_sample(
    bam_2HR, reference, annotation, vcf,
    outdir = file.path("prepare", "2HR"),
    chemistry=chemistry, nthreads=nthreads
    )

# Expression no quality filtering:
analyse_expression(
    h5 = prepared_2HR$h5,
    densities = densities,
    hdi = hdi,
    exprdir = file.path("expr", "2HR", "no_quality"),
    phylodir = file.path("phylo", "2HR", "expr", "no_quality"),
    prefix = "2HR",
    model = "ORDERED+ASC",
    minGene = 0,
    minUMI = 0,
    nthreads = nthreads
    )

# Expression with quality filtering:
analyse_expression(
    h5 = prepared_2HR$h5,
    densities = densities,
    hdi = hdi,
    exprdir = file.path("expr", "2HR", "quality"),
    phylodir = file.path("phylo", "2HR", "expr",  "quality"),
    prefix = "2HR.quality",
    model = "ORDERED+ASC",
    nthreads = nthreads
    )
# snv:
analyse_snv(
    bam = prepared_2HR$bam,
    barcodes = prepared_2HR$barcodes,
    reference = reference,
    densities = densities,
    snvdir = file.path("snv", "2HR"),
    phylodir = file.path("phylo", "2HR", "snv"),
    prefix = "2HR",
    model = "GTR+G+ASC",
    nthreads = nthreads
    )
