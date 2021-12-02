#' sra.r
#'
#' Function to programmatically access to a SRA and GEO databases
import::here("magrittr", "%>%")
import::here("xml2", "read_xml", "as_list")
import::here("phyloRNA", "mkdir", "all_files_exist")
import::here("rentrez", "entrez_search", "entrez_link", "entrez_summary", "extract_from_esummary")

#' Download SRR samples
#'
#' Download 10X fastq files from the SRA database
#'
#' @param srr the SRR (short-read record, presumably) id that identifies reads in the SRA database.
#' SRR ids can be obtained by `get_srr_samples`.
#' @param prefix a sample prefix that identifies downloaded fastq files
#' @param outdir an output directory where the fastq files are downloaded
#' @param cellranger **optional** the files are renamed to the cellranger type
srr_download_sample = function(srr, prefix, outdir=NULL, cellranger=TRUE){
    # For Cellranger, file names need to have this structure:
    # [name]_S1_L001_[read type]_001.fastq.gz
    # S1 -- sample 1
    # L001 -- lane 001
    # read type:
    #     -- I1 -- Index file
    #     -- R1 -- barcodes or actual reads
    #     -- R2 -- barcodes or actual reads
    if(is.null(outdir))
        outdir = "."
    mkdir(outdir)

    if(length(prefix) != 1 || length(srr) != 1 || length(outdir) != 1)
        stop("This function is not vectorized. Provide a single srr, prefix and outdir.")

    if(cellranger)
        read_types = c("I1", "R1", "R2")
        fastqs = paste0(prefix, "_S1_L001_", read_types , "_001.fastq.gz")
        } else {
        fastqs = paste0(prefix, ".fastq.gz")
        }

    fastqs = file.path(outdir, fastqs)
    srr_files = file.path(outdir, paste0(srr, "_", seq_along(fastqs), ".fastq.gz"))

    if(all_files_exist(fastqs))
        return(invisible())

    # Download srr files and check if they exist/were downloaded correctly
    sra_download(srr, outdir)
    if(!all_files_exist(srr_files))
        stop("ERROR: not all files exist.\\n", "Files: ", files)

    file.rename(srr_files, fastqs)

    return(invisible())
    }


sra_download = function(srr, outdir){
    # Check if the srr fastq file or files already exists
    files = dir(outdir, pattern=paste0(srr, ".*\\.fastq.*"))
    if(length(files) != 0)
        return(invisible())

    sra_prefetch(srr, outdir)
    sra_dump(srr, outdir)
    }


sra_prefetch = function(srr, outdir, progress=TRUE){
    command = "prefetch"
    args = c(
        srr,
        "--output-directory", outdir
        )
    if(progress)
        args = c(args, "--progress")

    phyloRNA:::systemE(command, args)
    }


sra_dump = function(srr, outdir){
    command = "fastq-dump"
    args = c(
        srr,
        "--split-files",
        #"--skip-technical", required for 10X
        "--origfmt",
        "--gzip"
        )
    phyloRNA:::systemE(command, args, dir=outdir)
    }


get_srr_samples = function(gse, save=NULL){
    if(!is.null(save) && file.exists(save))
        return(readRDS(save))

    samples = get_gsm_samples(gse)
    samples$srr = sapply(samples$gsm, get_srr)

    if(!is.null(save))
        saveRDS(samples, save)

    samples
    }


get_gsm_samples = function(gse){
    esearch = entrez_search(db="gds", term=paste0(gse, "[ACCN] & gse[ETYP]"))
    esummary = entrez_summary(db="gds", id=esearch$ids)

    samples = esummary$samples
    colnames(samples) = c("gsm", "name")
    samples
    }


get_srr = function(gsm){
    esearch = entrez_search(db="gds", term=paste0(gsm, "[ACCN] gsm[ETYP]"))
    elinks = entrez_link(dbfrom="gds", id=esearch$ids, db="sra")
    esummary = entrez_summary(db="sra", elinks$links$gds_sra)

    run = extract_from_esummary(esummary, "runs")
    attrs = extract_attributes(run)
    srr = attrs[["acc"]]

    srr
    }


extract_attributes = function(x){
    x %>% read_xml %>% as_list() %>% getElement("Run") %>% attributes
    }
