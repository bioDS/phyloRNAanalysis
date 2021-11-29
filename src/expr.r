#' expr.r
#'
#' Functions for processing 10x expression data
import::here("phyloRNA", "mkdir", "all_files_exist")
import::here("utils.r", "filename", "num2char")


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
#' @param h5 a path to a expression count matrix stored in a .h5 file
#' or list of such paths
#' @param density **optional** required density or densities of the final matrix
#' @param hdi **optional** a highest density intervals for discretization
#' @param minGene **optional** minimum number of genes per cell
#' @param minUMI **optional** minimum number of UMI per cell
#' @param outdir **optional** an output directory
#' @param prefix **optional** a prefix for file names
#' @param normalize **optional** log-normalize the expression data
#'
#' @return a list of paths of all outputs
preprocess_expression = function(
    h5, hdi=c(0.6,0.9),
    minGene=250, minUMI=500,
    outdir=NULL, prefix=NULL,
    normalize=FALSE
    ){
    if(is.null(outdir))
        outdir = "."
    if(is.null(prefix))
        prefix = "data"
    phyloRNA::mkdir(outdir)

    result = list(
        intervals = file.path(outdir, paste(prefix, "intervals", "txt", sep=".")),
        discretized = file.path(outdir, paste(prefix, "discretized", "txt", sep="."))
        )

    if(phyloRNA::all_files_exist(result))
        return(invisible(result))


    if(length(h5) > 1){
        names = phyloRNA::corename(h5)
        data = lapply(h5, phyloRNA::expr_read10xh5)
        data = phyloRNA::expr_merge(data, names)
        } else {
        data = phyloRNA::expr_read10xh5(h5)
        }


    if(minGene > 0 && minUMI > 0)
        data = phyloRNA::expr_quality_filter(data, minGene=minGene, minUMI=minUMI)
    data = phyloRNA::expr_zero_to_na(data)
    if(normalize)
        data = phyloRNA::expr_normalize(data)
    data = phyloRNA::expr_scale(data)

    intervals = calculate_intervals(data, density=hdi, save=result$intervals)
    discretized = phyloRNA::expr_discretize(data, intervals=intervals, unknown="-")
    write_table(discretized, result$discretized)

    return(invisible(result))
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
    names(intervals) = as.character(density)

    if(isTRUE(save))
        save = "intervals.txt"
    if(is.character(save))
        dput(intervals, save)

    sort(unlist(intervals))
    }


#' Filter Expression dataset
#'
#' Filter expression dataset using two types of filtering approaches
#' see `src/filter.r` for more information
#'
#' @param expr discretized expression data
#' @param selection named list of vector specifying how many cells of each type should be selected
#' @param density desired data density
#' @param outdir an output directory
#' @return a list of filtered files
filter_expression = function(
    expr, prefix, selection=NULL, density=NULL, outdir=NULL
    ){
    if(is.null(outdir))
        outdir = "."
    mkdir(outdir)

    if(is.null(selection) && is.null(density))
        stop("Either the selection or the density parameter must be specified.")

    if(!is.null(selection) && is.null(density)){
        file = filename(prefix, outdir=outdir)
        if(file.exsts(file))
            return(invisible(file))

        data = read_table(expr)
        file = subset_filtering(
            data,
            prefix = prefix,
            selection = selection,
            empty = "-",
            outdir = outdir
            )
        return(invisible(file))
        }

    if(!is.null(selection) && !is.null(density)){
        files = filename(prefix, num2char(density), outdir=outdir)
        if(all_files_exist(files))
            return(invisible(files))
        data = read_table(expr)
        files = subset_filtering(
            data,
            prefix = prefix,
            selection = selection,
            density = density,
            empty = "-",
            outdir = outdir,
            )
        return(invisible(files))
        }

    if(is.null(selection)){
        files = filename(prefix, num2char(density), outdir=outdir)
        if(all_files_exist(files))
            return(invisible(files))

        data = read_table(files)
        files = density_filtering(
            data,
            prefix = prefix,
            density = density,
            empty = "-",
            outdir = outdir
            )
        return(invisible(files))
        }
    }
