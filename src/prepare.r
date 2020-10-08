#' prepare.r
#'
#' Function for preparation of sequences.
library("phyloRNA")



#' Prepare multiple samples into a single bam file
#'
#' Remap, clean, rename and merge individual sample bam files for further processing, suhc as
#' SNV detection.
#'
#' This function chains together basic processing steps to prepare a previously mapped bam file
#' for further processing. This include:
#'
#' * remapping the bam file to a new reference genome
#' * filtering bam so that only the most well-represented cells are included
#' * cleaning the bam file according to GATK best practices
#' * changing names of the barcodes (cell identifiers) to ensure their uniqueness
#' * merging sample files into a single file
#'
#' It is suggested to specify the `refdir` as it can be shared between analysis using
#' the same reference genome and annotation.
#'
#' @param bams a previously aligned bam files from the same organism that will be remapped
#' to a new reference and merged together
#' @param reference a new genome reference file (".fai")
#' @param annotation an annotation file associated with the reference genome
#' @param vcf a known variants associated with the reference genome that will be masked
#' @param outdir **optional** a general output directory
#' @param refdir **optional** an output directory for the cellranger reference genome files
#' @param chemistry **optiona** 10X chemistry, use only when the automatic detection is failing
#' @param nthreads **optional** a number of threads to use
#' @return a list with a paths to merged bam and barcode files.
prepare_samples = function(
    bams, reference, annotation, vcf,
    obam=NULL, obar=NULL, oh5=NULL,
    outdir=NULL, refdir=NULL,
    chemistry="auto", nthreads=16
    ){

    if(phyloRNA::is_nn(outdir))
        outdir = "prepare"
    if(phyloRNA::is_nn(refdir))
        refdir = file.path(outdir, "ref")
    if(phyloRNA::is_nn(exprdir))
        exprdir = file.path(outdir, "expr")
    if(phyloRNA::is_nn(obam))
        obam = file.path(outdir, "all.bam")
    if(phyloRNA::is_nn(obar))
        obar = file.path(outdir, "all.txt")
    if(phyloRNA::is_nn(oh5))
        oh5 = file.path(outdir, paste0(corename(bams), ".h5"))

    result = list(
        samples = corename(bams),
        bam = obam,
        barcodes = obar,
        oh5 = oh5
        )

    if(file.exists(obam) && file.exists(obar) && all(file.exists(oh5)))
        return(result)

    results = list()
    for(bam in bams)
        results[[bam]] = prepare_sample(
            bam = bam,
            reference = reference,
            annotation = annotation,
            vcf = vcf,
            outdir = file.path(outdir, corename(bam)),
            refdir = refdir,
            chemistry = chemistry,
            nthreads = nthreads
            )

    pcores = unlist(lapply(results, getElement, "prefix"))
    pbams = unlist(lapply(results, getElement, "bam"))
    pbars = unlist(lapply(results, getElement, "barcodes"))
    ph5s = unlist(lapply(results, getElement, "h5"))

    # Defensive programming: Sanity check that the corenames are exactly the same.
    if(all(pcores != result$samples))
        stop("Corenames from samples differ from those in results. This shouldn't happen.")

    # merge prepared bams from all datasets:
    phyloRNA::gatk_MergeSamFiles(pbams, obam)

    # merge prepared barcodes from all datasets:
    merge_files(pbars, obar, overwrite=TRUE)

    return(result)
    }



#' Prepare sample bam file
#'
#' Remap, clean and rename bam file for further processing, such as SNV detection.
#'
#' This function chains together basic processing steps to prepare a previously mapped bam file
#' for further processing. This include:
#'
#' * remapping the bam file to a new reference genome
#' * filtering bam so that only the most well-represented cells are included
#' * cleaning the bam file according to GATK best practices
#' * changing names of the barcodes (cell identifiers) to ensure their uniqueness
#'
#' Function include a number of optional parameters to specify an output folder for each step.
#' `outdir` is a general output directory. If not specified, `prepare` is used. If other optional
#' directories are not specified, they default to a directory sitting in the outdir.
#' * `mapdir` defaults to `map`
#' * `refdir` defaults to `ref`
#' * `cleandir` defaults to `clean`
#'
#' It is suggested to specify the `refdir` as it can be shared between analysis using
#' the same reference genome and annotation.
#'
#' @param bam a previously aligned bam file that will be remapped to a new reference
#' @param reference a new genome reference file (".fai")
#' @param annotation an annotation file associated with the reference genome
#' @param vcf a known variants associated with the reference genome that will be masked
#' @param outdir **optional** a general output directory
#' @param mapdir **optional** an output directory for the remapping step
#' @param refdir **optional** an output directory for the cellranger reference genome files
#' @param cleandir **optional** an output directory for GATK cleaning/preparation steps
#' @param chemistry **optiona** 10X chemistry, use only when the automatic detection is failing
#' @param nthreads **optional** a number of threads to use
#' @return a list with a paths to prepared bam and barcode files.
prepare_sample = function(
    bam, reference, annotation, vcf,
    outdir=NULL, mapdir=NULL, refdir=NULL, cleandir=NULL,
    chemistry = "auto", nthreads=16
    ){
    if(phyloRNA::is_nn(outdir))
        outdir = "prepare"
    if(phyloRNA::is_nn(mapdir))
        mapdir = file.path(outdir, "map")
    if(phyloRNA::is_nn(refdir))
        refdir = file.path(outdir, "ref")
    if(phyloRNA::is_nn(cleandir))
        cleandir = file.path(outdir, "clean")

    core = phyloRNA::corename(bam)

    bam_aligned = filename(outdir, core, ".aligned.bam")
    bam_cleaned = filename(outdir, core, ".cleaned.bam")
    bam_prepared = filename(outdir, core, ".prepared.bam")

    barcodes_aligned = filename(outdir, core, ".aligned.txt")
    barcodes_prepared = filename(outdir, core, ".prepared.txt")

    h5_prepared = filename(outdir, core, ".h5")

    result = list(
        prefix = core,
        bam = bam_prepared,
        barcodes = barcodes_prepared
        h5 = h5_prepared
        )

    # Skip if the final output files already exist
    if(file.exists(result$bam) && file.exists(result$barcodes) && file.exists(result$h5))
        return(result)

    phyloRNA::remap(
        input = bam,
        reference = reference,
        annotation = annotation,
        outdir = mapdir,
        refdir = refdir,
        chemistry = chemistry,
        nthreads = nthreads,
        copy_bam = bam_aligned,
        copy_bar = barcodes_aligned,
        copy_h5 = h5_prepared
        )

    phyloRNA::gatk_prepare(
        input = bam_aligned,
        output = bam_cleaned,
        reference = reference,
        vcf = vcf,
        barcodes = barcodes_aligned,
        outdir = cleandir
        )

    pattern = "-1$"
    replace = paste0("-", core)

    phyloRNA::bamtagregex(
        input = bam_cleaned,
        output = bam_prepared,
        tag = "CB",
        pattern = pattern,
        replace = replace
        )

    barcodes = readLines(barcodes_aligned)
    barcodes = sub(pattern, replace, barcodes)
    writeLines(barcodes, barcodes_prepared)

    return(result)
    }


filename = function(dir, core, ext){
    file.path(dir, paste0(core, ext))
    }


merge_files = function(inputs, output, overwrite=FALSE){
    if(file.exists(output) && overwrite) file.remove(output)
    file.append(output, inputs)
    }


merge_h5 = function(inputs, output){
    # well, I can read, merge them, but then not save in the same h5 format.
    # I would need to write a function for that.
    data = lapply(inputs, phyloRNA::expr_read10xh5)
    names = phyloRNA::corenames(inputs)
    data = phyloRNA::expr_merge(data, names)
    }
