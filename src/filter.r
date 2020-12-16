#' filter.r
#'
#' Functions for the alternative filtering approach
library("phyloRNA")


column_density = function(x, empty, sort=TRUE){
    cs = colSums(is_empty(x, empty))
    if(sort)
        cs = sort(cs, decreasing=TRUE)
    cs
    }


#' Divide vector into categories
#'
#' Divide vector into list of categories using the pattern and replace substitution.
#'
#' Pattern can be customized according to specific needs using the `[base::sub]` command.
#' The remainder after character substitution is then used as category.
#'
#' @param x a vector
#' @param pattern a pattern argument for the `[base::sub]` command
#' @param replace a replace argument for the `[base::sub]` command
#' @return a list of vectors for each category
divide_vector = function(x, pattern=".*-", replace=""){
    categories = sub(pattern, replace, x)
    result = list()
    for(category in unique(categories)){
        result[[category]] = x[categories == category]
        }
    result
    }

#' Select number of elements from list
#'
#' Take number of elements, specified by selection, from each category in the input list.
#'
#' @param x a list of vectors from the `divide` function
#' @param selection a list or named vector specifying how many elements should be selected
#' from each category
#' @return a list of the same structure as `x` with number of elements specified by `selection`
select_from_list = function(x, selection){
    result = list()
    categories = names(x)
    for(category in categories){
        result[[category]] = x[[category]][seq_len(selection[[category]])]
        }
    result
    }


subset_rows = function(x, k, empty){
    rs = rowSums( is_empty(x, empty) )
    y = x[rs > k, ]
    y = phyloRNA::remove_constant(y, margin=1, empty)
    y
    }


#' Select columns
#'
#' Divide columns according to the pattern. Then select a number of columns according
#' to the selection vector with the least amount of missing data.
#'
#' Pattern a
#'
#' @param x a table
#' @param selection a named vector of selected columns for each subset
#' @param pattern pattern for the sub command
#' @param replace replace for the sub command
#' @param empty missing value specification
#' @return subsetted dataset
select = function(x, selection, pattern=".*-", replace="", empty=NA){
    # These operation work with sorted names
    # rather than whole input data
    columns = names(column_density(x, empty, sort=TRUE))
    columns = divide_vector(columns, pattern, replace)
    columns = select_from_list(columns, selection)
    columns = unlist(columns)
    x = x[,columns]
    # Filter columns with only single value 
    subset_rows(x, 1, empty)
    }




#' Remove rows with no data
#'
#' Filter dataset by removing rows, down to a particular data density.
#' @param x dataset
#' @param density required density
#' @param empty missing value specification
densest_rows = function(x, density=0.5, empty=NA){
    k = 1
    while(TRUE){
        y = subset_rows(x, k, empty)

        # Filtering would result in a empty dataset
        if(nrow(y) == 0)
            break

        # Stop just bellow desired density
        if(mdensity(y, empty) > density)
            break
 
        x = y
        k = k + 1
        }
    x
    }


density_filenames = function(outdir, prefix, density){
    file.path(outdir, paste0(prefix, num2char(density), ".txt"))
    }


#' Filtering by finding the densest subset
#'
#' Filter data by finding the densest rows and colums.
#' See `[phyloRNA::densest_subset]` for more information.
#'
#' @param x data in tabular format
#' @param density desired data density
#' @param empty missing value specification
#' @param outdir an output directory
#' @param prefix a prefix for output files
#' @return vector of paths for filtered datasets
density_filtering = function(x, density=0.5, empty="N", outdir=NULL, prefix=NULL){
    if(is.null(outdir))
        outdir = "."
    if(is.null(prefix))
        prefix = "filtered"
    mkdir(outdir)

    outfile = density_filenames(outdir, prefix, density)

    if(all.files.exists(outfile))
        return(invisible(outfile))

    for(i in seq_along(density)){
        filtered = phyloRNA::densest_subset(x, empty=empty, density=density[i])$result
        filtered = phyloRNA::remove_constant(filtered, margin=1, unknown=empty)
        write_table(filtered, outfile[i])
        }

    invisible(outfile)
    }




#' Filtering by selecting the best cells
#'
#' Filter data by first finding the best cells with the least amount of missing data.
#' Then the rows are filtered to reach desired density.
#'
#' @param x data in tabular format
#' @param selection named list or vector specifying how many cells should be selected of each type
#' or sample
#' @param density desired data density
#' @param empty missing value specification
#' @param outdir an output directory
#' @param prefix a prefix for output files
#' @return vector of paths for filtered datasets
subset_filtering = function(x, selection, density=0.5, empty="N", outdir=NULL, prefix=NULL){
    if(is.null(outdir))
        outdir = "."
    if(is.null(prefix))
        prefix = "filtered"
    mkdir(outdir)

    outfiles = subset_filtering_filenames(outdir, prefix, density)

    if(all.files.exists(outfiles))
        return(invisible(outfiles))

    subset = select(x, selection, empty=empty)
    write_table(subset, outfiles[1])

    for(i in seq_along(density)){
        filtered = densest_rows(subset, density[i], empty=empty)
        write_table(filtered, outfiles[i+1])
        }

    invisible(outfiles)
    }


subset_filtering_filenames = function(outdir, prefix, density){
    filenames = c(
        file.path(outdir, paste0(prefix, ".txt")),
        density_filenames(outdir, prefix,density)
        )
    filenames
    }
