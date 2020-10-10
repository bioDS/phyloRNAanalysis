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
    phylodir = file.path("phylo", "all", "no_quality"),
    prefix = "all",
    model = "ORDINAL+ASC",
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
    model = "ORDINAL+ASC",
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
    phylodir = file.path("phylo", "2HR", "no_quality"),
    prefix = "2HR",
    model = "ORDINAL+ASC",
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
    model = "ORDINAL+ASC",
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



############################################################################
# TODO: move this to src/expr.r
analyse_expression = function(
    h5, densities, hdi,
    exprdir, phylodir, prefix, model,
    minGene=250, minUMI=300, nthreads=16
    ){
    expressed = expr_process(
        file=h5, dens=densities, hdi=hdi,
        minGene=minGene, minUMI=minUMI,
        outdir=exprdir, prefix=prefix
        )

    iqtrees(expressed$fasta, outdir=phylodir, num2char(densities),
            model=model, nthreads=nthreads)
    }


# TODO: move this to src/snv.r
analyse_snv = function(
    bam, barcodes, reference, densities,
    snvdir, phylodir, prefix, model, nthreads=16
    ){
    snv = detect_snv(bam, barcodes, reference, outdir=snvdir)
    vcmdir = file.path(snvdir, "vcm")
    fasta = vcm2fasta(snv$vcm, density=densities, outdir=vcmdir, prefix=prefix)
    iqtrees(fasta$fasta, outdir=outdir, num2char(densities),
            model=model, nthreads=nthreads)
    }
