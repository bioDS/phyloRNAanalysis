#' expr.r
#'
#' Functions for processing 10x expression data
import::here("utils.r", "filename", "num2char", "mdensity", "write_table", "read_table")
import::here("phyloRNA",
    "mkdir", "all_files_exist", "corename",
    "expr_merge", "expr_read10xh5",
    "expr_quality_filter", "expr_zero_to_na",
    "expr_normalize", "expr_scale", "expr_discretize",
    "remove_constant", "tab2seq", "write_fasta"
    )


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
    mkdir(outdir)

    result = list(
        intervals = file.path(outdir, paste(prefix, "intervals", "txt", sep=".")),
        discretized = file.path(outdir, paste(prefix, "discretized", "txt", sep="."))
        )

    if(all_files_exist(result))
        return(invisible(result))

    if(length(h5) > 1){
        names = corename(h5)
        data = lapply(h5, expr_read10xh5)
        data = expr_merge(data, names)
        } else {
        data = expr_read10xh5(h5)
        }

    data = process_expression(
        data, hdi = hdi,
        minGene = minGene, minUMI = minUMI,
        trim = FALSE, normalize = normalize,
        intervals = result$intervals
        )

    write_table(data, result$discretized)

    return(invisible(result))
    }


#' Process expression data
#'
#' This function simplifies standard expression data processing, such as quality filtering, scaling,
#' normalization and discretization.
#'
#' @param data an expression matrix
#' @param hdi **optional** a highest density intervals for discretization
#' @param minGene **optional** a minimum amount of represented genes per cell
#' @param minUMI **optional** a minimum amount of total UMI (or count) per cell
#' @param trim **optional** trim empty genes after filtering
#' @param normalize **optional** perform normalization after rescaling
#' @param intervals **optional** a file path to save discretization intervals into file
#' @param unknown **optional** a symbol representing unknown data
#' @return scaled, filtered and discretized count matrix
process_expression = function(
    data, hdi = c(0.6, 0.9),
    minGene = 0, minUMI = 0,
    trim = FALSE, normalize = FALSE,
    intervals = FALSE, unknown = "-"
    ){
    if(minGene > 0 || minUMI > 0 || trim)
        data = expr_quality_filter(data, minGene, minUMI, trim)

    data = expr_zero_to_na(data)

    if(normalize)
        data = expr_normalize(data)

    data = expr_scale(data)
    intervals = calculate_intervals(data, density=hdi, save=intervals)
    data = expr_discretize(data, intervals=intervals, unknown=unknown)
    data
    }


expr2fasta = function(x, fasta, unknown="-", summary=FALSE, process=TRUE, hdi=c(0.6, 0.9)){
    data = x
    if(process)
        data = process_expression(x, hdi, trim=TRUE, unknown=unknown)

    data = remove_constant(data, margin=1, unknown=unknown)
    seq = tab2seq(data, margin=2)

    write_fasta(seq, fasta)

    if(isTRUE(summary))
        summary = file.path(dirname(fasta), paste0(corename(fasta), "_summary.txt"))
    if(is.character(summary))
        count_matrix_summary(data, name=corename(fasta), file=summary)
    }


count_matrix_summary = function(data, name=NULL, file=NULL){
    text = paste0(
        "Sequences: ", ncol(data), "\n",
        "Sites: ", nrow(data), "\n",
        "Unique patterns: ", nrow(unique.matrix(data, MARGIN=1)), "\n",
        "Data density: ", mdensity(data, empty="-")
        )

    if(!is.null(name))
        text = paste0("Name: ", name[1], "\n", text)

    if(!is.null(file))
        writeLines(text, file)

    return(invisible(text))
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
