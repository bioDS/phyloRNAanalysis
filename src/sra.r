#' sra.r
#'
#' Function for a programatical access to a SRA and GEO databases
import::here("magrittr", "%>%")
import::here("xml2", "read_xml", "as_list")

sra_download = function(srr, outdir){
    # Check if the srr fastq file or files already exists
    files = dir(outdir, pattern=paste0(srr, ".*\\.fastq.*"))
    if(length(files) != 0)
        return(invisible())
    
    sra_prefetch(srr, outdir)
    sra_dump(srr, outdir)
    }


sra_prefetch = function(srr, outdir){
    command = "prefetch"
    args = c(
        srr,
        "--output-directory", outdir
        )    
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


get_srr_samples = function(gse){
    samples = get_gsm_samples(gse)
    samples$srr = sapply(samples$gsm, get_srr)
    
    samples
    }


get_gsm_samples = function(gse){
    esearch = rentrez::entrez_search(db="gds", term=paste0(gse, "[ACCN] & gse[ETYP]"))
    esummary = rentrez::entrez_summary(db="gds", id=esearch$ids)
    
    samples = esummary$samples
    colnames(samples) = c("gsm", "name")
    samples
    }


get_srr = function(gsm){
    esearch = rentrez::entrez_search(db="gds", term=paste0(gsm, "[ACCN] gsm[ETYP]"))
    elinks = rentrez::entrez_link(dbfrom="gds", id=esearch$ids, db="sra")
    esummary = rentrez::entrez_summary(db="sra", elinks$links$gds_sra)
    
    run = rentrez::extract_from_esummary(esummary, "runs")
    attrs = extract_attributes(run)
    srr = attrs[["acc"]]
    
    srr
    }


extract_attributes = function(x){
    x %>% read_xml %>% as_list() %>% getElement("Run") %>% attributes
    }
