#' snv.r
#'
#' Functions for snv identification and filtering
library("phyloRNA")

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
detect_snv = function(bam, barcodes, reference, outdir=NULL, vcfdir=NULL, vcm=NULL){
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

    vcm
    }
