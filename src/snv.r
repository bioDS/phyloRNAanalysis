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

    invisible(vcm)
    }

#' Filter the VCM table
#'
#' Find a densest subset of the SNV data contained in the vcm table.
#'
#' @param vcm a vcm file
#' @param density a vector of the required density or densities
#' @param outdir **optional** an output directory
#' @param prefix **optional** a file prefix
vcm2fasta = function(vcm, density=0.5, outdir=NULL, prefix=NULL){
    if(is.null(outdir))
        outdir = "."
    if(is.null(prefix))
        prefix = "data"

    result = list(
        filtered = file.path(outdir, paste(prefix, num2char(density), "txt", sep=".")),
        fasta = file.path(outdir, paste(prefix, num2char(density), "fasta", sep="."))
        )

    if(all.files.exists(result))
        return(invisible(result))

    # using data.table due to a huge size of the dataset
    data = data.table::fread(vcm, header=TRUE)
    data = data[,-c(1:3)] # first three columns are not cells (chromosome, position and reference)
    
    for(i in seq_along(density)){
        dens = density[i]
        filtered = phyloRNA::densest_subset(data, empty="N", steps=0, density=dens)$result
        filtered = phyloRNA::remove_constant(filtered, margin=1)
        write_table(filtered, result$filtered[i])
        fasta = phyloRNA::fasta(filtered, file=result$fasta[i])
        }

    return(invisible(result))
    }


#' Analyse expression data
#'
#' Analyse expression data, filter them and run iqtrees
#'
#' @param bam a bam file
#' @param barcodes file with barcodes
#' @param reference a genome reference to which the bam file was mapped
#' @param densities a desired final density for the matrix filtering step
#' @param snvdir a directory for the snv output
#' @param phylodir a directory for the iqtree output
#' @param prefix a file prefix for output files
#' @param model a iqtree model definition
#' @param nthreads a number of threads to run software on
analyse_snv = function(
    bam, barcodes, reference, densities,
    snvdir, phylodir, prefix, model, nthreads=16
    ){
    vcm = detect_snv(bam, barcodes, reference, outdir=snvdir)
    vcmdir = file.path(snvdir, "vcm")
    fasta = vcm2fasta(vcm, density=densities, outdir=vcmdir, prefix=prefix)
    iqtrees(fasta$fasta, outdir=outdir, num2char(densities),
            model=model, nthreads=nthreads)
    }
