#' utils.r
#'
#' shared utility functions
library("tools")
library("phyloRNA")

#' Write a table
#'
#' Write a table in a particular format. This is a simple wrapper around write.table
#' with a few specified parameters.
#' @param x a matrix or a data frame
#' @param file an output path
write_table = function(x, file){
    write.table(x, file=file, quote=FALSE, sep="\t", col.names=TRUE, row.names=TRUE)
    }


#' Convert a numeric value to a character string
#'
#' Converts a numeric value in the format `0.X` into a character string `0X`
#' @param x numeric vector
#' @return character vector
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


#' Convert a tabular file into a fasta format
#'
#' @param file one or more files in tabular format
#' @param fasta **optional** output path for fasta files
#' @param outdir **optional** an outut directory, if fasta is not specified
#' @param margin whether rows (1) or columns (2) should be concatenated
#' @return a vector of fasta files
table2fasta = function(file, fasta=NULL, outdir=NULL, margin=2){
    if(!is.null(fasta)) && length(fasta) != length(file)
        stop("The file and fasta vectors must have the same length!")
    if(is.null(fasta))
        fasta = tools::file.path.sans.ext(file)
    if(!is.null(outdir))
        fasta = file.path(outdir, basename(fasta))
    mkdir(outdir)

    if(all.files.exists(fasta))
        return(invisible(fasta))

    for(i in seq_along(file)){
        data = read.table(file[i], header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
        fasta = phyloRNA::fasta(data, margin=margin, file=fasta[i])
        }

    invisible(fasta)
    }
