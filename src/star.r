#' star.r
#'
#' mapping to the STAR RNA-seq mapping software
import::here("utils.r", "rg2cb")
import::here("phyloRNA", "corename", "mkdir")


#' Prepare sequences using the STAR software
#'
#' Map scRNA-seq sequences using the STAR software.
#'
#' This functions creates a STAR genome index and then maps the fastqs to the reference.
#' Reads are then retaged so that the read-group information is written to the cell barcode tag.
#'
#' @param prefix a prefix for the output bam file
#' @param fastqs one or more fastqs files
#' @param sample a sample information for each fastq file, will be written into readgroup
#' @param flowcell a flowcell information for the whole analysis, will be written into read-group
#' @param reference a fasta reference sequences
#' @param annotatin a gtf annotation file
#' @param overhang a size of allowed overhang, should be equal to the read lengths
#' @param outdir **optional** an output directory
#' @param mapdir **optional** a directory with mapped output
#' @param refdir **optional** a directory with STAR reference index
#' @param manifest **optional** a file path to pre-existing manifest
#' @param gzip **optional** whether the fastq files are compressed
#' @param nthreads **optional** the number of threads to run the analysis on
#' @return a mapped bam file with CB tag
STAR = function(
    prefix, fastqs, sample, flowcell,
    reference, annotation, overhang,
    outdir = NULL, refdir = NULL, mapdir = NULL,
    manifest = NULL, gzip = TRUE,
    nthreads = 16
    ){
    if(is.null(outdir))
        outdir = "."
    if(is.null(refdir))
        refdir = file.path(outdir, "ref")
    if(is.null(mapdir))
        mapdir = file.path(mapdir, "map")
    if(is.null(manifest))
        manifest = file.path(mapdir, "manifest.txt")

    mkdir(outdir)
    mkdir(refdir)
    mkdir(mapdir)

    mapped = file.path(mapdir, "allAligned.sortedByCoord.out.bam")
    tagged = file.path(mapdir, "all.bam")

    if(!file.exists(mapped)){    
        STAR_genome_index(reference, annotation, overhang, refdir, nthreads=nthreads)
        read_groups = make_read_group(sample, flowcell)
        STAR_manifest(fastqs, manifest, read_groups)
        STAR_map(prefix, refdir, manifest=manifest, outdir=mapdir, gzip=gzip)
        }

    if(!file.exists(tagged))
        rg2cb(mapped, tagged)

    return(tagged)
    }


#' Map sequences wih STAR
#'
#' Map scRNA-seq sequences using the STAR software
#'
#' The input fastq files can be provided in two ways, either directly using
#' the `fastq` parameter, or using the `manifest` parameter, which in addition
#' enable specification of the read-groups (see the [STAR_manifest]).
#'
#' @param prefix a file prefix for output files
#' @param reference a path to the STAR genome index, see [STAR_genome_index]
#' @param fastq **optional** fastq files that will be mapped, see details.
#' @param manifest **optional** STAR manifest file, see details
#' @param outdir **optional** an output directory
#' @param gzip **optional** whether the input fastq files are compressed
#' @param nthreads **optional** the number of threads to use
STAR_map = function(
    prefix, reference,
    fastq=NULL, manifest=NULL,
    outdir=NULL, gzip=FALSE, nthreads=16
    ){
    if(is.null(fastq) && is.null(manifest))
        stop("Either fastq or manifest must be specified.")
    if(is.null(fastq) && is.null(manifest))
        stop("Please, specify either fastq or manifest, not both.")

    if(is.null(outdir))
        outdir = "."
    mkdir(outdir)

    statusfile = file.path(outdir, paste0(prefix, ".finished"))
    if(file.exists(statusfile))
        return(invisible())

    command = "STAR"
    args = c(
        "--genomeDir", reference,
        "--outFileNamePrefix", file.path(outdir, prefix),
        "--outSAMattributes NH HI AS nM NM MD jM jI MC ch RG",
        "--outSAMtype BAM SortedByCoordinate"
        )

    if(!is.null(fastq)){
        id = paste0("ID:", corename(fastq), collapse=" , ")
        args = c(args, "--readFilesIn", paste0(fastq, collapse=","))
        args = c(args, "--outSAMattrRGline", id)
        }

    if(!is.null(manifest))
        args = c(args, "--readFilesManifest", manifest)

    if(gzip)
        args = c(args, "--readFilesCommand gunzip -c")

    phyloRNA:::systemE(command, args)
    file.create(statusfile)
    }


#' Create STAR manifest
#'
#' Create STAR manifest.
#'
#' STAR manifest is an alternative way to specify input `fastq` files for the STAR
#' RNA-seq mapping software. For the specification of read-group tag,
#' see [make_read_group]
#'
#' @param fastq one or more fastq files
#' @param file a character string naming a file
#' @param rg **optional** a read-group tag, see [make_read_group]
STAR_manifest = function(fastq, file, rg=NULL){
    if(file.exists(file))
        return(invisible())

    manifest = fastq

    if(!is.null(rg))
        manifest = paste(fastq, "-", rg, sep="\t")
        
    writeLines(manifest, file)
    }


#' Construct a read-group tag
#'
#' Construct a read-group (RG) tag required by the GATK. See [details] for the
#' description of individual parts of the RG tag.
#'
#' SAM/BAM read-groups tags required by the GATK are:
#' RG tags required by GATK are:
#'    * ID identifies read and is copied to read
#'    * SM identifies sample
#'    * LB identifies library, can have multiple libraries per sample
#'    * PU consist of flowcell barcode, lane and sample barcode
#'    * PL identifies platform, Illumina in this case
#'
#' The the values for the RG tags are:
#'  * ID -- `sample`
#'  * SM -- `sample`
#'  * LB -- `sample`
#'  * PU -- `flowcell:1:sample`
#'  * PL -- "ILLUMINA"
make_read_group = function(sample, flowcell, platform="ILLUMINA"){
    paste(
        paste0("ID:", sample),
        paste0("SM:", sample),
        paste0("LB:", sample),
        paste0("PU:", paste(flowcell, "1", sample, sep=":")),
        paste0("PL:", "ILLUMINA"),
        sep="\t"
        )
    }


#' Generate a STAR genome index
#'
#' Generate a reference STAR genome index.
#'
#' @param reference a path to the reference fasta file
#' @param gtf a path to the annotation gtf file
#' @param overhang **optional** a size of allowed overhang, should be equal
#' to read length
#' @param outdir **optional** an output directory
#' @param nthreads **optional** the number of threads to use
STAR_genome_index = function(reference, gtf, overhang=50, outdir=NULL, nthreads=16){
    if(is.null(outdir))
        outdir = "."
    mkdir(outdir)

    statusfile = file.path(outdir, "genome.finished")
    if(file.exists(statusfile))
        return(invisible())
    
    command = "STAR"
    args = c(
        "--runMode genomeGenerate",
        "--genomeDir", outdir,
        "--genomeFastaFiles", reference,
        "--sjdbGTFfile", gtf,
        "--sjdbOverhang", len-1,
        "--runThreadN", nthreads
        )
    phyloRNA:::systemE(command, args)
    file.create(statusfile)
    }
