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
#' @param file a path to a expression count matrix stored in a .h5 file
#' or list of such paths
#' @param density **optional** required density or densities of the final matrix
#' @param hdi **optional** a highest density intervals for discretization
#' @param minGene **optional** minimum number of genes per cell
#' @param minUMI **optional** minimum number of UMI per cell
#' @param outdir **optional** an output directory
#' @param prefix **optional** a prefix for file names
#'
#' @return a list of paths of all outputs
expr_process = function(
    file=NULL, density=0.5, hdi=c(0.6,0.9),
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
        discretized = file.path(outdir, paste(prefix, "discretized", "txt", sep=".")),
        filtered = file.path(outdir, paste(prefix, num2char(dens), "txt", sep=".")),
        fasta = file.path(outdir, paste(prefix, num2char(dens), "fasta", sep="."))
        )

    if(all.files.exists(result))
        return(invisible(result))


    if(length(file) > 1){
        names = phyloRNA::corename(file)
        data = lapply(file, phyloRNA::expr_read10xh5)
        data = phyloRNA::expr_merge(data, names)
        } else {
        data = phyloRNA::expr_read10xh5(file)
        }


    if(minGene > 0 && minUMI > 0)
        data = phyloRNA::expr_quality_filter(data, minGene=minGene, minUMI=minUMI)
    data = phyloRNA::expr_zero_to_na(data)
    if(normalize)
        data = phyloRNA::expr_normalize(data)
    data = phyloRNA::expr_scale(data)

    intervals = calculate_intervals(data, density=hdi, save=result$intervals)
    discretized = phyloRNA::expr_discretize(data, intervals=intervals, unknown="-")
    write.table(discretized, result$discretized)

    for(i in seq_along(dens)){
        density = dens[i]
        filtered = phyloRNA::densest_subset(discretized, empty="-", density=density)$result
        filtered = phyloRNA::remove_constant(filtered, margin=1, unknown="-")
        write_table(filtered, result$filtered[i])

        fasta = phyloRNA::fasta(filtered, file=result$fasta[i])
        }

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


#' Write a table
#'
#' Write a table in a particular format. This is a simple wrapper around write.table
#' with a few specified parameters.
write_table = function(x, file){
    write.table(x, file=file, quote=FALSE, sep="\t", col.names=TRUE, row.names=TRUE)
    }


num2char = function(x){
    sub(".", "", as.character(x), fixed=TRUE)
    }


#' Check existence of files
#'
#' This function checks the existence of all files stored in a list.
#'
#' @param files a list of files.
#' @return a logical value indicating if all files exists.
all.files.exists = function(x){
    all(file.exists(unlist(x)))
    }


#' Analyse expression data and build trees
#'
#' Process, clean and analyse expression data and build phylogenies with an IQtree
#'
#' @param h5 an expression count matrix in the h5 format
#' @param densities desired final density at the end of the matrix filtering step
#' @param hdi highest density interval for expression categorization
#' @param model a model definitio for the IQtree
#' @param minGene a minimum number of required genes in the expression filtering step
#' @param minUMI a minimum number of UMI in the expression filtering step
#' @param nthreads a desired number of threads to run the software on
analyse_expression = function(
    h5, densities, hdi,
    exprdir, phylodir, prefix, model,
    minGene=250, minUMI=300, nthreads=16
    ){
    expressed = expr_process(
        file=h5, dens=densities, hdi=hdi,
        minGene=minGene, minUMI=minUMI,
        outdir=exprdir, prefix=prefix
        )

    iqtrees(expressed$fasta, outdir=phylodir, num2char(densities),
            model=model, nthreads=nthreads)
    }
