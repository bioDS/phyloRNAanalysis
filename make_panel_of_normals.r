#' make_panel_of_normals.r
#'
#' Construct the GATK Panel of Normals of MDA-MB-231 cell lineage from the GSE181410 data.
import::from("src/sra.r", "get_srr_samples", "srr_download_sample")
import::from("phyloRNA",
    "cellranger_mkref", "cellranger_count",
    "gatk_prepare", "gatk_make_pon",
    "mkdir"
    )

main = function(){
    outdir = "pon/MDA-MB-231"
    fastqdir = file.path(outdir, "fastq")
    mapdir = file.path(outdir, "map")
    refdir = "reference/ref"
    gse = "GSE181410"
    
    reference = "reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
    annotation = "reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gtf"
    vcf = "reference/00-common_all.vcf.gz"

    # create annotation:
    cellranger_mkref(reference, annotation, refdir, nthreads=8)

    # get samples
    mkdir(outdir)
    samples = get_srr_samples(gse, save=file.path(outdir, paste0(gse, ".rds")))

    # construct names from samples
    samples$name = make_sample_names(samples$name)

    # download samples
    mkdir(fastqdir)
    Map(srr_download_sample, srr=samples$srr, prefix=samples$name, outdir=fastqdir)

    # map and demultiplex
    mkdir(mapdir)
    outputs = Map(cellranger_count,
        id = samples$name, sample = samples$name,
        fastqdir = fastqdir, refdir = refdir,
        outdir = file.path(mapdir, samples$name),
        nthreads = 8
        ) 

    aligned = sapply(outputs, getElement, "bam")
    cleaned = file.path(mapdir, paste0(samples$name, ".cleaned.bam"))

    Map(gatk_prepare,
        input = aligned, output = cleaned,
        reference = reference, vcf = vcf
        )

    gatk_make_pon(
        bam = cleaned,
        reference = reference,
        vcf = file.path(outdir, paste0("pon.vcf")),
        outdir = file.path(outdir, "mutect2")
        )
    }

make_sample_names = function(x){
    CancerMitotransfer = grepl("mito-mEmerald+", x, fixed=TRUE)
    MacrophageMitotransfer = grepl("mitoRFP+", x, fixed=TRUE)
    rep = sub(".*(rep[1-9])", "\\1", x)
    names = paste0(
        "MDAMB231",
        ifelse(CancerMitotransfer, "mt", ""), "-",
        ifelse(MacrophageMitotransfer, "MPmt", "MP"), "-",
        rep
        )

    names
    }


if(sys.nframe() == 0){
    main()
    }
