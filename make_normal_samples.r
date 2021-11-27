#' make_normal_samples.r
#'
#' Construct normal samples for the GATK Mutect2 SNV in the normal/tumour mode

import::from("src/sra.r", "get_srr_samples", "srr_download_sample")
import::from("src/utils.r", "download_file")
import::from("phyloRNA",
    "cellranger_mkref", "cellranger_count",
    "gatk_prepare", "gatk_make_pon",
    "mkdir", "all_files_exist"
    )
import::from("parallel", "mcMap")

main = function(){
    # Some of this was already done in the make_panel_of_normals.r
    # reusing code and processed data where possible
    outdir = "moravec2021/normal"
    shareddir = "pon/MDA-MB-231"
    fastqdir = file.path(shareddir, "fastq")
    mapdir = file.path(shareddir, "map")
    refdir = "reference/ref"
    gse = "GSE181410"
    
    reference = "reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
    annotation = "reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gtf"
    vcf = "reference/00-common_all.vcf.gz"

    # data for the Seurat object with barcodes
    url = paste0("ftp://ftp.ncbi.nlm.nih.gov/",
                 "geo/series/GSE181nnn/GSE181410/suppl/",
                 "GSE181410%5Fmda%2Dmb%2D231%2Eh5seurat")

    # make dirs
    mkdir(outdir)
    mkdir(shareddir)

    # create annotation:
    cellranger_mkref(reference, annotation, refdir, nthreads=8)

    # get samples
    samples = get_srr_samples(gse, save=file.path(shareddir, paste0(gse, ".rds")))

    # construct names from samples
    samples$name = make_sample_names(samples$name)

    # select samples with MDA-MB-231 without mitotransfer
    samples = samples[grepl("MDAMB231-MPmt", samples$name, fixed=TRUE),]
    ncores = nrow(samples)

    # download samples
    mkdir(fastqdir)
    mcMap(srr_download_sample, srr=samples$srr, prefix=samples$name, outdir=fastqdir,
          mc.cores=ncores)

    # map and demultiplex
    mkdir(mapdir)
    outputs = mcMap(
        cellranger_count,
        id = samples$name, sample = samples$name,
        fastqdir = fastqdir, refdir = refdir,
        outdir = file.path(mapdir, samples$name),
        nthreads = 8,
        mc.cores = ncores
        ) 

    aligned = sapply(outputs, getElement, "bam")
    barcodes = sapply(outputs, getElement, "barcodes")

    # Prepare barcodes:
    h5barcodes = file.path(outdir, paste0(samples$name, ".h5barcodes.gz"))
    if(!all_files_exist(h5barcodes)){
        h5samples = c("rep1"="2A", "rep2"="2B")
        h5samples = h5samples[vapply(c("rep1", "rep2"), grep, integer(1), samples$name)]
        h5seurat = file.path(outdir, "MDA-MD-231.h5seurat")
        local({
            bar = get_barcodes(url, h5samples, h5seurat)
            Map(function(x,y){
                    con = gzfile(y, "w")
                    writeLines(x, con)
                    close(con)
                    }, bar, h5barcodes
                    )
            })

        # Compare barcodes and h5barcodes
        local({
            breport = file.path(outdir, "barcode_report.txt")
            con = file(breport, open="w+") # needs to be wrapped in a list for Map
            Map(compare_barcodes, barcodes, h5barcodes, con=list(con))
            close(con)
            })
    }

    prepared = file.path(outdir, paste0(samples$name, ".prepared.bam"))
    mcMap(gatk_prepare,
        input = aligned, output = prepared,
        reference = reference, vcf = vcf,
        barcodes = h5barcodes,
        mc.cores = ncores
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



get_barcodes = function(url, samples, h5seurat=NULL){
    if(is.null(h5seurat))
        h5seurat = tempfile(fileext=".h5seurat")

    download_file(url, h5seurat)
    data = SeuratDisk::LoadH5Seurat(h5seurat, verbose=FALSE)
    barcodes = lapply(
        samples,
        function(x){
            y = filter_seurat_object(data, x)
            y = colnames(y)
            y = remove_prefix(y, x)
            y
            }
        )
    names(barcodes) = samples
    barcodes
    }


remove_prefix = function(x, prefix){
    sub(paste0(prefix, "_"), "", x, fixed=TRUE)
    }


filter_seurat_object = function(data, samples){
    # Select all samples with the MDA-MB-231 cell lineage (should be all of them)
    lineage = data$clusterlabeled == 231
    # Select all samples from 2A and 2B (no macrophage mitochondrial transport to 231)
    sample = data$sample %in% samples
    # Select samples with the percentage of mitochondrial genes smaller than 5%
    mito = data$percent.mt < 5
    # Construct selection vector combining these three conditions
    selection = lineage & sample & mito
    # Get the count matrix with these things
    count_matrix = data@assays$RNA@counts[, selection]
    count_matrix
    }


compare_barcodes = function(f1, f2, con){
    t1 = readLines(f1)
    t2 = readLines(f2)
    
    t1nott2 = t1[!t1 %in% t2]
    t2nott1 = f2[!t2 %in% t1]

    text = c(
        paste0(basename(f1), " barcodes not in ", basename(f2), ":"),
        paste0(t1nott2, collapse=", "), "",
        paste0(basename(f2), " barcodes not in ", basename(f1), ":"),
        paste0(t2nott1, collapse=", ")
        )
    writeLines(text, con)
    }


if(sys.nframe() == 0){
    main()
    }
