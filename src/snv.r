#' snv.r
#'
#' Functions for snv identification and filtering
import::here("utils.r", "filename", "num2char")
import::here("data.table", "fread")
import::here("phyloRNA", "all_files_exist")

#' Detect SNV for scRNAseq
#'
#' This function will perform SNV detection for scRNAseq data.
#'
#' To detect SNVs on a scRNAseq, this function first detect a high-quality SNVs by treating
#' the sample as a bulk sample, then the most common base pertaining to each cell for every
#' variant is reported and saved to a vcm (variant call matrix) file.
#'
#' @param bam an input sam/bam file
#' @param barcodes a file with cell barcodes
#' @param reference a reference file to which the bam file was mapped
#' @param outdir **optional** a general output directory
#' @param vcfdir **optional** a directory with vcf files
#' @param vcm **optional** an output vcm file
#'
#' @return a path to an alignment table with the most frequent base for every cell at every SNV
#' position
detect_snv = function(
    bam, barcodes, reference,
    normal=NULL, pon=NULL, germline=NULL,
    outdir=NULL, vcfdir=NULL, vcm=NULL,
    nthreads=16
    ){
    core = phyloRNA::corename(bam)

    # Set default parameters:
    if(is.null(outdir))
        outdir = "snv"
    if(is.null(vcfdir))
        vcfdir = file.path(outdir, "vcf")
    if(is.null(vcm))
        vcm = file.path(outdir, paste0(core, ".vcm"))

    phyloRNA::mkdir(outdir)
    phyloRNA::mkdir(vcfdir)

    vcf = file.path(outdir, paste0(core, ".vcf"))

    phyloRNA::gatk_snv(
        bam, reference, vcf,
        normal=normal, pon=pon, germline=germline,
        outdir=vcfdir)
    phyloRNA::vcm(bam, vcf, barcodes, output=vcm, nthreads=nthreads)

    invisible(vcm)
    }




#' Preprocess SNV data
#'
#' Preprocessing SNVs. This include as bulk SNV detection, single-cell SNV detection
#' using bulk SNVs and dataset filtering.
#'
#' @param bam a bam file prepared according to the GATK best practices
#' @param barcodes file with barcodes
#' @param reference a genome reference to which the bam file was mapped
#' @param outdir **optional** an output directory
#' @param nthreads **optional** a number of threads to run on
preprocess_snv = function(
    bam, barcodes, reference, outdir=NULL, nthreads=16
    ){
    if(is.null(outdir))
        outdir = "."
    phyloRNA::mkdir(outdir)

    vcm = detect_snv(bam, barcodes, reference, outdir=outdir, nthreads=nthreads)

    invisible(vcm)
    }


#' Filter SNV
#'
#' Filter SNV dataset using two types of filtering approaches
#' see `src/filter.r` for more information
#'
#' @param vcm a vcm file
#' @param selection named list or vector specifying how many cells of each type should be selected
#' @param density desired data density
#' @param outdir an output directory
#' @return a list of filtered files
filter_snv = function(vcm, prefix, selection=NULL, density=NULL, outdir=NULL){
    if(is.null(outdir))
        outdir = "."
    phyloRNA::mkdir(outdir)

    if(is.null(selection) && is.null(density))
        stop("Either the selection or the density parameter must be specified.")


    if(!is.null(selection) && is.null(density)){
        file = filename(prefix, outdir=outdir)
        if(file.exists(file))
            return(invisible(file))

        data = as.data.frame(read_vcm(vcm))
        file = subset_filtering(
            data,
            prefix = prefix,
            selection = selection,
            empty = "N",
            outdir = outdir
            )
        return(invisible(file))
        }


    if(!is.null(selection) && !is.null(density)){
        files = filename(prefix, num2char(density), outdir=outdir)
        if(all_files_exist(files))
            return(invisible(files))

        data = as.data.frame(read_vcm(vcm))
        files = subset_filtering(
            data,
            prefix = prefix,
            selection = selection,
            density = density,
            empty = "N",
            outdir = outdir
            )
        return(invisible(files))
        }

    if(is.null(selection)){
        files = filename(prefix, num2char(density), outdir=outdir)
        if(all_files_exist(files))
            return(invisible(file))

        data = read_vcm(vcm)
        files = density_filtering(data, prefix=prefix, density=density, empty="N", outdir=outdir)
        return(invisible(files))
        }
    }


#' Read the vcm file and memoise it
#'
#' This function is memoised (possible reuse) of the vcm file.
#' `data.table::fread()` is used here due to a huge file size.
#' @param vcm a variant call matrix file
#' @return a variant call matrix as a data.table
read_vcm = local({
    memory = list()
    
    function(vcm){
        if(!is.null(memory[vcm]))
            return(memory[vcm])

        # using data.tale due to a huge size of the dataset
        data = fread(vcm, header=TRUE)
        # first three columns are not cells (chromosome, position and reference)
        data = data[, -c(1:3)]
        memory[vcm] <<- data

        data
        }
    })
