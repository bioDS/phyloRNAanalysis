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

# use import::from instead?
import::from("src/prepare.r", "prepare_samples")
import::from("src/snv.r", "detect_snv", "filter_snv")
import::from("src/utils.r", "table2fasta")
import::from("src/expr.r", "preprocess_expression", "filter_expression")
import::from("src/beast.r", "beasts")
import::from("src/iqtree.r", "iqtrees")
import::from("phyloRNA", "corename")

main = function(){
    # datasets:
    bam = dir("data", full.names=TRUE)
    outdir = file.path("moravec2021")

    # required reference files:
    reference = "reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
    annotation = "reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gtf"
    refdir = "reference/ref" # shared cellranger ref
    vcf = "reference/00-common_all.vcf.gz"
    
    # PON and normal samples
    pon = "pon/MDA-MB-231/pon.vcf" # see make_panel_of_normals.r
    normal = file.path( # see make_normal_samples.r
        outdir, "normal",
        c("MDAMB231-MPmt-rep1.prepared.bam", "MDAMB231-MPmt-rep2.prepared.bam")
        )

    # Other settings:
    nthreads = 16
    chemistry = "SC3Pv2"
    densities = c(0.2, 0.5, 0.9)
    hdi = c(0.6, 0.9)
    selection = c("T1" = 20, "T3" = 20, "T2" = 6, "CTC1" = 6, "CTC2" = 6)

    # Preparation:
    prepared = prepare_samples(
        bam, reference, annotation, vcf,
        chemistry = chemistry, nthreads = nthreads,
        outdir = file.path(outdir, "prepare"),
        refdir = refdir
        )

    # SNV part
    vcm = detect_snv(
        bam = prepared$bam,
        barcodes = prepared$barcodes,
        reference = reference,
        normal = normal,
        pon = pon,
        outdir = file.path(outdir, "snv")
        )

    tabdir = file.path(outdir, "snv", "filtered")
    fastadir = file.path(outdir, "snv", "fasta")
    treedir = file.path(outdir, "snv", "tree")

    snv = filter_snv(vcm=vcm, density=densities, prefix="snv", outdir=tabdir)
    snv_fasta = table2fasta(snv, outdir=fastadir)

    snv_subset = filter_snv(vcm=vcm, selection=selection, prefix="snv_subset",
                            outdir=tabdir)
    snv_subset_fasta = table2fasta(snv, outdir=fastadir)

    iqtrees(
        c(snv, snv_subset_fasta),
        model = "GTR+G+ASC",
        bootstrap = 100, parallel = TRUE, nthreads = 16,
        outdir = file.path(treedir, "ML")
        )

    beasts(
        snv_subset_fasta,
        template = file.path("templates", "BDStrictGtr.xml"),
        outdir = file.path(treedir, "BI")
        )

    # expression part
    expr_preprocessed = preprocess_expression(
        h5 = prepared$h5,
        hdi = hdi,
        minGene=0,
        minUMI=0,
        outdir = file.path(outdir, "expr", "prepare"),
        prefix = "all"
        )

    filterdir = file.path(outdir, "expr", "filtered")
    fastadir = file.path(outdir, "expr", "fasta")
    treedir = file.path(outdir, "expr", "tree")
    
    expr = filter_expression(
        expr_processed$discretized, outdir = filterdir,
        prefix = "expr", densities = densities,
        ) 
    expr_subset = filter_expression(
        expr_preprocessed$discretized, filterdir,
        prefix = "expr_subset", selection = selection,
        )

    expr_fasta = table2fasta(expr, outdir=fastadir)
    expr_subset_fasta = table2fasta(expr_subset, outdir=fastadir)
    expr_zero = table2fasta(
        expr_subset,
        file.path(fastadir, paste0(corename(expr_subset), "_zero.fasta")),
        outdir = fastadir, zero = "-"
        )

    iqtrees(
        c(expr_fasta, expr_subset_fasta, expr_zero),
        model = "ORDERED+ASC",
        bootstrap = 100, parallel = TRUE, nthreads = 16,
        outdir = file.path(treedir, "ML")
        )

    beasts(expr_subset_fasta, outdir = file.path(treedir, "BI"),
           template = file.path("templates", "BDStrictOrdinal.xml"))
    beasts(expr_zero, outdir = file.path(treedir, "BI"),
           template = file.path("templates", "BDStrictOrdinalZero.xml"))
    }


if(sys.nframe() == 0){
    main()
    }
