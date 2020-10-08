#' Process 10X expression data
#'
#' This function will filter, scale, discretize and transform scRNAseq expression data
#' into a fasta sequence.
#'
#' The 10X expression count matrix can be provided either
#'
#' The function performs following steps:
#' * `phyloRNA::expr_quality_filtering()` -- to remove low-quality cells
#' * `phyloRNA::expr_zero_to_na()` -- transform 0 to NA
#' * `phyloRNA::expr_normalize()` -- normalize expression values accross genes
#' * `phyloRNA::expr_scale()` -- center and scale expression values accross cells
#' * `phyloRNA::expr_discretize()` -- discretize data according to HDI intervals
#' * `phyloRNA::densest_subset()` and `phyloRNA::remove_constant()` -- to filter the data matrix
#' to a chosen density
#' * `phyloRNA::fasta()` -- transform to FASTA sequence
#'
#' @param data a loaded expression count matrix or list of such matrices
#' @param file a path to a expression count matrix stored in a .h5 file
#' or list of such paths
#' @param names **optional** if a list of files/matrices is provided, names of the datasets
#' can be specified
#' @param dens **optional** required density or densities of the final matrix
#' @param minGene **optional** minimum number of genes per cell
#' @param minUMI **optional** minimum number of UMI per cell
#' @param outdir **optional** an output directory
#' @param name **optional** a prefix for file names
#' @param save_intervals **optional** save the discretization intervals
#' @param save_discretized **optional** save the discretized matrix
#' @param save_filtered **optional** save the filtered matrix or matrices
#' @param save_fasta **optional** save the fasta sequence or sequences
#'
#' @return a list of paths of all outputs
expr_process = function(
    data=NULL, file=NULL, names=NULL,
    dens=0.5, minGene=250, minUMI=500,
    hdi=c(0.6,0.9), outdir=NULL, name=NULL,
    save_intervals=FALSE, save_discretized=FALSE,
    save_filtered=FALSE, save_fasta=FALSE
    ){
    if(is.null(data) && is.null(file))
        stop("Either data or file need to be specified!")
    if(!is.null(data) && !is.null(file))
        stop("Only data or file need to be specified, not both!")

    if(is.null(outdir))
        outdir = "."
    if(is.null(name))
        name = "data"

    if(!is.null(file)){
        if(length(file) > 1){
            data = lapply(file, phyloRNA::expr_read10xh5)
            } else {
            data = phyloRNA::expr_read10xh5(file)
            }
        }


    if(is.list(data) && length(data) > 1){
        data = phyloRNA::expr_merge(data, names)
        }

    if(minGene > 0 && minUMI > 0)
        data = phyloRNA::expr_quality_filter(data, minGene=minGene, minUMI=minUMI)
    data = phyloRNA::expr_zero_to_na(data)
    data = phyloRNA::expr_normalize(data)
    data = phyloRNA::expr_scale(data)

    if(isTRUE(save_intervals))
        save_intervals = file.path(outdir, paste(name, "intervals", txt, sep="."))

    intervals = calculate_intervals(data, density=hdi, save=save_intervals)
    discr = phyloRNA::expr_discretize(data, intervals=intervals, unknown="-")

    if(isTRUE(save_discretized))
        save_discretized = file.path(outdir, paste(name, "discretized", "txt", sep="."))
    if(is.character(save_discretized))
        write.table(discr, save_discretized)

    fasta = list()
    for(density in dens){
        fasta[[density]] = densest_fasta(
            discr, density, name=name,
            save_filtered = save_filtered,
            save_fasta = save_fasta
            )
        }

    if(length(fasta) == 1)
        fasta = fasta[[1]]

    return(fasta)
    }

#' Calculate discretization intervals from data
#'
#' Calculate the discretization intervals empirically using the HDI method.
#'
#' @param data an input data
#' @param density a single or multiple values representing the total density covered by signle
#' or multiple the HDIs
#' @param save **optional** a TRUE or path to save the calculated interval values into a file
#' @return intervals corresponding to the HDI of the density vector
calculate_intervals = function(data, density=c(0.6,0.9), save=FALSE){
    dens = stats::density(data, na.rm=TRUE)

    intervals = lapply(density, function(x) phyloRNA::hdi(dens, 1-x))
    if(isTRUE(save))
        save = "intervals.txt"
    if(is.character(save))
        dput(intervals, save)

    sort(unlist(intervals))
    }


#' Write a table
#'
#' Write a table in a particular format. This is a simple wrapper around write.table
#' with a few specified parameters.
write_table = function(x, file){
    write.table(x, file=file, quote=FALSE, sep="\t", col.names=TRUE, row.names=TRUE)
    }


#' Filter dataset and return a fasta
#'
#' @param data a data table
#' @param density a requested density of fasta
#' @param outdir **optional** a general output directory
#' @param name **optional** a prefix for output files
#' @param save_filtered **optional** TRUE or an output path for filtered table
#' @param save_fasta **optional** TRUE or path for filtered fasta file
#' @return sequences in a fasta format
densest_fasta = function(
        data, density=0.5,
        outdir=NULL, name=NULL,
        save_filtered=FALSE, save_fasta=FALSE
        ){
        if(is.null(outdir))
            outdir = "."
        if(is.null(name))
            name = "filtered"

        filtered = phyloRNA::densest_subset(data, empty="-", steps=10000, density=dens)$result
        filtered = phyloRNA::remove_constant(filtered, margin=1)

        if(isTRUE(save_filtered))
            save_filtered = file.path(outdir, paste(name, num2char(density), "txt", sep="."))
        if(is.character(save_filtered))
            write_table(filtered, save_filtered)
        if(isTRUE(save_fasta))
            save_fasta = file.path(outdir, paste(name, num2char(density), "fasta", sep="."))

        fasta = phyloRNA::fasta(filtered, file=save_fasta)
        return(fasta)
        }


num2char = function(x){
    sub(".", "", as.character(x), fixed=TRUE)
    }
