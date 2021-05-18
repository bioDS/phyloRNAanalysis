import::from("src/sra.r", "sra_download", "get_srr_samples")
import::from("magrittr", "%>%")

main = function(){
    # Define variables:
    outdir = "rao2020"
    fastqdir = file.path(outdir, "raw")
    refdir = "prepare/ref"
    countdir = file.path(outdir, count)
    gse = GSE140312
    srr_download(gse, outdir)


    # map with cellranger
    phyloRNA::cellranger_count(fastqdir, refdir, countdir, nthreads=8)
    
    # preprocess SNVs

    # filter -> same selection as with expression

    # iqtree
    }




srr_download = function(gse, outdir){
    samples = get_srr_samples(gse)
    samples$name = gsub(" ", "_", sample$name)

    mapply(samples$srr, samples$name, FUN=list("outdir"=outdir)
    }


srr_download_sample = function(srr, name, outdir){
    # For Cellranger, file names need to have this structure:
    # [sample_name]_S1_L001_[read type]_001.fastq.gz
    # S1 -- sample 1
    # L001 -- lane 001
    # read type:
    #     -- I1 -- Index file
    #     -- R1 -- barcodes or actuall reads
    #     -- R2 -- barcodes or actual reads
    read_types = c("I1", "R1", "R2")
    sample_name = strsplit(name, split="_", fixed=TRUE) %>% unlist %>% getElement(1)

    cellranger_files = file.path(
        outdir,
        paste0(sample_name, "_S1_L001_", read_types, "_001.fastq.gz")
        )

    if(all.files.exists(cellranger_files))
        return(invisible())
    
    # Download srr files and check if they exist/were downloaded correctly
    srr_files = file.path(outdir, paste0(srr, "_", 1:3, ".fastq.gz"))
    sra_download(srr, outdir)
    if(!all.files.exists(files))
        stop("ERROR: not all files exists.\\n", "Files: ", files)
    
    file.rename(srr_files, files)
    }

if(!interactive()){
    main()
    }
