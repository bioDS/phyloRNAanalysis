#' wang2021.r
#'
#' Reconstruct phylogenetic trees for data from:
#' Wang et al. (2021) "Comprehensive analysis of metastatic gastric cancer tumour cells using
#' single-cell RNA-seq" Sci. Rep.
import::from("src/sra.r", "get_srr_samples", "srr_download_sample")
import::from("src/star.r", "STAR")
import::from("src/snv.r", "detect_snv")
import::from("src/expr.r", "expr2fasta")
import::from("src/iqtree.r", "iqtrees")
import::from("src/beast.r", "beasts")
import::from("src/utils.r", "vcm2fasta", "fasta2stats")
import::from("phyloRNA", "gatk_prepare", "abspath", "mkdir")
import::from("parallel", "mcMap")

main = function(){
    snv()
    expr()
    }


snv = function(){
    outdir = "wang2021"
    fastqdir = file.path(outdir, "raw")
    snvdir = file.path(outdir, "snv")
    refdir = file.path(outdir, "ref")
    mapdir = file.path(outdir, "map")
    fastadir = file.path(snvdir, "fasta")
    treedir = file.path(snvdir, "tree")
    gse = "GSE158631"
    
    # required reference files:
    reference = "reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
    annotation = "reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gtf"
    vcf = "reference/00-common_all.vcf.gz"
    pon = "pon/1000g_pon.hg38.vcf.gz"
    germline = "pon/af-only-gnomad.hg38.vcf.gz"

    # Download fastqs
    samples = get_srr_samples(gse, save=file.path(outdir, "samples.rds"))
    ncores = max(nrow(samples), 16)
    mcMap(srr_download_sample, samples$srr, samples$name, fastqdir, cellranger=FALSE, mc.cores=ncores)
    fastqs = file.path(fastqdir, paste0(samples$name, ".fastq.gz"))
    fastqs = sapply(fastqs, abspath)

    # Map with STAR
    mapped = STAR("all", fastqs, samples$name, gse, reference, annotation, overhang=50,
         outdir=outdir, refdir=refdir, mapdir=mapdir)

    # Prepare with GATK
    prepared = file.path(mapdir, "all.prepared.bam")
    gatk_prepare(mapped, prepared, reference=reference, vcf=vcf)

    # Detect SNV with Mutect2
    vcm = detect_snv(prepared, barcodes=NULL, reference, pon=pon, germline=germline, outdir=snvdir)

    # Prepare fasta files
    names = c("GC1", "GC2", "GC3")
    n = length(names)

    gc = lapply(names, grep, x=samples$name, value=TRUE)
    fasta = file.path(fastadir, paste0(names, ".fasta"))

    mkdir(fastadir)
    mcMap(vcm2fasta, vcm, fasta, gc, mc.cores=n) 
    fasta2stats(fasta, unknown="N")

    # run IQtree
    iqtrees(fasta, "TEST", outdir=file.path(treedir, "ML"), mc.cores=n)

    # run BEAST
    template = file.path("templates", "BDStrictGtr.xml")
    beasts(fasta, template, outdir=file.path(treedir, "BI"), mc.cores=n)
    }


expr = function(){
    hdi = c(0.6, 0.9)
    outdir = file.path("wang2021", "expr")
    fastadir = file.path(outdir, "fasta")
    treedir = file.path(outdir, "tree")

    url = paste0("https://www.ncbi.nlm.nih.gov/geo/download/",
                 "?acc=GSE158631&format=file&file=GSE158631%5Fcount%2Ecsv%2Egz")
    file = file.path(outdir, "count_matrix.csv.gz")
    download_file(url, file)

    # Read the data and split it according to patient
    count_matrix = read.table(infile, sep=",", header=TRUE, row.names=1)
    count_matrix = split.default(count_matrix, sub("\\..*", "", names(count_matrix)))
    n = length(count_matrix)

    fasta = file.path(fastadir, paste0(names(count_matrix), ".fasta"))
    Map(expr2fasta, count_matrix, fasta, summary=TRUE)

    # Run phylogenetic analyses
    iqtrees(fasta, "ORDERED+ASC", outdir=file.path(treedir, "ML"), mc.cores=n)
    template = file.path("templates", "BDStrictOrdinal.xml")
    beasts(fasta, template, outdir=file.path(treedir, "BI"), mc.cores=n)
    }


if(sys.nframe() == 0){
    main()
    }
