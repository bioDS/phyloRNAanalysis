#' snv.r
#'
#' Functions for snv identification and filtering
library("phyloRNA")

#' Detect SNV for scRNAseq
#'
#' This function will perform SNV detection for scRNAseq data.
#'
#' To detect SNVs on a scRNAseq, this function first detect a high-quality SNVs by treating
#' the sample as a bulk sample, then examine each the base frequency at each SNV for every cell
#' and then finally reports the most frequent base for every cell and every SNV, which will create
#' the final alignment.
#'
#' @param bam an input sam/bam file
#' @param barcodes a file with cell barcodes
#' @param reference a reference file to which the bam file was mapped
#' @param outdir **optional** a general output directory
#' @param vcfdir **optional** a directory with vcf files
#' @param vffdir **optional** a directory with vff files
#'
#' @return a path to an alignment table with the most frequent base for every cell at every SNV
#' position
detect_snv = function(bam, barcodes, reference, outdir=NULL, vcfdir=NULL, vffdir=NULL){
    if(phyloRNA::is_nn(outdir))
        outdir = "snv"
    if(phyloRNA::is_nn(vcfdir))
        vcfdir = file.path(outdir, "vcf")
    if(phyloRNA::is_nn(vffdir))
        vffdir = file.path(outdir, "vff")

    core = corename(bam)
    vcf = file.path(outdir, paste0(core, ".vcf"))
    alignment = file.path(outdir, paste0(core, ".alignment.txt"))

    gatk_snv(bam, reference, vcf, vcfdir)
    vff_make(bam, vcf, barcodes, vffdir)
    vff_merge(vffdir, alignment)

    alignment
    }
