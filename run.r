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
source("src/expr.r")
source("src/snv.r")

# TODO:
# * define datasets
# * define which analysis will be done
# * separate so that there is no name clash


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

# Preparation step:
prepared_all = prepare_samples(
    bam_all, reference, annotation, vcf,
    chemistry=chemistry, nthreads=nthreads
    )
prepared_2HR = prepare_sample(
    bam_2HR, reference, annotation, vcf,
    outdir = file.path("prepare", "2HR"),
    chemistry=chemistry, nthreads=nthreads
    )
# Prepare external data

# Expression step:
## all samples, no quality filtering
outdir = file.path("expr", "all", "no_quality")
phyloRNA::mkdir(outdir)
expr_process(
    file=prepared_all$h5, names=phyloRNA::corename(prepared_all$h5),
    dens=densities, hdi=hdi,
    minGene=0, minUMI=0,
    outdir=outdir, name="all",
    save_intervals=TRUE, save_discretized=TRUE, save_filtered=TRUE, save_fasta=TRUE
    )
## all samples, quality filtering
outdir = file.path("expr", "all", "quality")
phyloRNA::mkdir(outdir)
expr_process(
    file=prepared_all$h5, names=phyloRNA::corename(prepared_all$h5),
    dens=densities, hdi=hdi,
    minGene=250, minUMI=300,
    outdir=outdir, name="all.quality",
    save_intervals=TRUE, save_discretized=TRUE, save_filtered=TRUE, save_fasta=TRUE
    )


## 2HR, no quality filtering
outdir = file.path("expr", "2HR", "no_quality")
phyloRNA::mkdir(outdir)
expr_process(
    file=prepared_2HR$h5,
    dens=densities, hdi=hdi,
    minGene=0, minUMI=0,
    outdir=outdir, name="2HR",
    save_intervals=TRUE, save_discretized=TRUE, save_filtered=TRUE, save_fasta=TRUE
    )
## 2HR, quality filtering
outdir = file.path("expr", "2HR", "quality")
phyloRNA::mkdir(outdir)
expr_process(
    file=prepared_2HR$h5,
    dens=densities, hdi=hdi,
    minGene=250, minUMI=300,
    outdir=outdir, name="2HR.quality",
    save_intervals=TRUE, save_discretized=TRUE, save_filtered=TRUE, save_fasta=TRUE
    )


# SNV detection step:
detect_snv(prepared_all$bam, prepared_all$barcodes, reference, outdir=file.path("snv", "all"))

detect_snv(prepared_2HR$bam, prepared_2HR$barcodes, reference, outdir=file.path("snv", "2HR"))
