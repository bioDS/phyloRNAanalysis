#' run_pon.r
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
    GSE = "GSE181410"
    
    reference = "reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
    annotation = "reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gtf"
    vcf = "reference/00-common_all.vcf.gz"

    # create annotation:
    cellranger_mkref(reference, annotation, refdir, nthreads=8)

    # get samples
    samples = get_srr_samples(gse, save=file.path(outdir, paste0(GSE, ".rds")))

    # construct names from samples
    sample$names = make_sample_names(samples$names)

    # download samples
    mkdir(outdir)
    mkdir(fastqdir)
    Map(srr_download_sample, srr=samples$srr, prefix=sample$names, outdir=fastqdir)

    # map and demultiplex
    Map(cellranger_count,
        id = samples$names, sample = samples$names,
        fastqdir = fastqdir, refdir = refdir, outdir = mapdir,
        nthreads = 8
        ) 

    aligned = file.path(mapdir, sample$names, "outs", "possorted_genome_bam.bam")
    cleaned = file.path(mapdir, paste0(sample$names, ".cleaned.bam"))

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
    CancerMitotransfer = grepl("mit-mEmerald+", x)
    MacrophageMitotransfer = grepl("mitoRFP+", x)
    rep = sub(".*(rep[1-9])", "\\1", x)
    names = paste0(
        "MDA-MB-231",
        ifelse(CancerMitotransfer, "+", "-"),
        ifelse(MacrophageMitotransfer, "mp+", "mp-"),
        rep
        )

    names
    }


if(sys.nframe() == 0){
    main()
    }
