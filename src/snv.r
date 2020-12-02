#' snv.r
#'
#' Functions for snv identification and filtering
library("phyloRNA")
library("data.table")

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
    outdir=NULL, vcfdir=NULL, vcm=NULL,
    nthreads=16
    ){
    core = phyloRNA::corename(bam)

    if(phyloRNA::is_nn(outdir))
        outdir = "snv"
    if(phyloRNA::is_nn(vcfdir))
        vcfdir = file.path(outdir, "vcf")
    if(phyloRNA::is_nn(vcm))
        vcm = file.path(outdir, paste0(core, ".vcm"))

    mkdir(outdir)
    mkdir(vcfdir)

    vcf = file.path(outdir, paste0(core, ".vcf"))

    phyloRNA::gatk_snv(bam, reference, vcf, vcfdir)
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
    mkdir(outdir)

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
filter_snv = function(vcm, selection, density=0.5, outdir=NULL){
    if(is.null(outdir))
        outdir = "."
    mkdir(outdir)

    # using data.table due to a huge size of the dataset
    data = data.table::fread(vcm, header=TRUE)
    data = data[,-c(1:3)] # first three columns are not cells (chromosome, position and reference)

    prefix = "snv"
    filter = density_filtering(data, density=density, empty="N", outdir=outdir, prefix=prefix)

    prefix = "snv_subset"
    subset = subset_filtering(
        data, selection=selection, density=density,
        empty="N", outdir=outdir, prefix=prefix
        )
    
    list("filter" = filter, "subset" = subset)
    }
