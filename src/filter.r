#' filter.r
#'
#' Functions for the alternative filtering approach
library("phyloRNA")
import::here("src/utils.r", "num2char")
import::here("phyloRNA", "all_files_exist")

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
#' @param outdir **optional** an output directory, by defautl current working
#' directory is used
#' @param prefix **optional** a prefix for output files, by `filtered` is used
#' @param replace **optional** replace missing value with this character
#' @param rescale **optional** rescale the ordinal scale so that currently present
#' categories form a sequence, this this might be required by some computational
#' software. Note that this redefines the meaning behind categories.
#' If `TRUE`, a numeric sequence starting from `1` is used.
#' Alternatively, user can provide their own ordinal scale to which the object
#' will be rescaled.
#' @return vector of paths for filtered datasets
density_filtering = function(
    x, density=0.5, empty="N",
    outdir=NULL, prefix=NULL,
    replace=NULL, rescale=FALSE
    ){
    if(is.null(outdir))
        outdir = "."
    if(is.null(prefix))
        prefix = "filtered"
    mkdir(outdir)

    outfile = density_filenames(outdir, prefix, density)

    if(phyloRNA::all_files_exist(outfile))
        return(invisible(outfile))

    for(i in seq_along(density)){
        filtered = phyloRNA::densest_subset(x, empty=empty, density=density[i])$result
        filtered = phyloRNA::remove_constant(filtered, margin=1, unknown=empty)
        if(!is.null(replace))
            filtered = replace_missing(filtered, empty, replace)
        if(!phyloRNA::is_nn(rescale)){
            if(length(rescale) == 1 && rescale) # rescale = TRUE
                filtered = phyloRNA::replace_ordinal(filtered)
            if(length(rescale) > 1) # e.g.: rescale = 0:10; rescale = letters
                filtered = phyloRNA::replace_ordinal(filtered, rescale)
            }
        write_table(filtered, outfile[i])
        }

    invisible(outfile)
    }


replace_missing = function(data, missing, replace){
    replace = as.character(replace)
    if(length(replace) != 1)
        stop("ERROR: Provide exactly one character as a replacement")
    if(nchar(replace) != 1)
        stop("ERROR: Provide exactly one character as a replacement")
    if(is.null(missing) || is.null(replace))
        stop("ERROR: Missing and replace characters cannot be NULL")

    if(is.na(missing))
        data[is.null(data)] = replace
        
    if(!is.null(missing))
        data[data == missing] = replace

    data
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
subset_filtering = function(x, selection, density=NULL, empty="N", outdir=NULL, prefix=NULL){
    if(is.null(outdir))
        outdir = "."
    if(is.null(prefix))
        prefix = "filtered"
    mkdir(outdir)


    if(is.null(density)){
        file = filename(prefix, outdir=outdir)
        if(file.exists(file))
            return(invisible(file))

        subset = select(x, selection, empty=empty)
        write_table(subset, file)
        return(invisible(file))
        }

    if(!is.null(density)){
        files = filename(prefix, num2char(density), outdir=outdir)
        if(all_files_exist(files))
            return(invisible(files))
        subset = select(x, selection, empty=empty)
        Map(
            function(d,f) write_table(densest_row(subset, d, empty=empty), f),
            density, files
            )
        return(invisible(files))
        }
    }
