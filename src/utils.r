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


#' Read a table
#'
#' Read a table in a particular format. This is a simple wrapper around read.table
#' with a few specified parameters.
#' @param file a file in tabular format
#' @return data.frame
read_table = function(file){
    read.table(file, header=TRUE, check.names=FALSE, stringsAsFactors=FALSE, sep="\t")
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
    if(!is.null(fasta) && length(fasta) != length(file))
        stop("The file and fasta vectors must have the same length!")
    if(is.null(fasta))
        fasta = paste0(tools::file_path_sans_ext(file), ".fasta")
    if(!is.null(outdir))
        fasta = file.path(outdir, basename(fasta))
    mkdir(outdir)

    if(all.files.exists(fasta))
        return(invisible(fasta))

    for(i in seq_along(file)){
        data = read_table(file[i])
        phyloRNA::fasta(data, margin=margin, file=fasta[i])
        }

    invisible(fasta)
    }


#' Read fasta file
#'
#' @param file a fasta file
#' @return a named vector of sequences
read_fasta = function(file){
    text = readLines(file)
    starts = grep(">", text)
    stops = c(starts[-1] - 1, length(text))

    fasta = mapply(
        function(start, stop, text){
            seq = text[(start+1):stop]
            seq = gsub("[:blank:]*", "", seq)
            paste0(seq, collapse="")
            },
        starts, stops, MoreArgs=list(text)
        )
    names(fasta) = sub("^>", "", text[starts] )
    fasta
    }


#' Write fasta file
#'
#' @param fasta a named vector of sequences
#' @param file an output file
write_fasta = function(fasta, file){
    text = paste0(">", names(fasta), "\n", fasta)
    writeLines(text, file)
    }


#' Find sequences shared by two fasta files
#'
#' Read two fasta files, finds shared sequences and then write these sequences fasta files.
#'
#' @param x a fasta file
#' @param y a fasta file
#' @param outdir **optional** an output directory for filtered fasta files
#' @param xout **optional** an output path for filtered x fasta
#' @param yout **optional** an output path for filtered y fasta
fasta_intersect = function(x, y, outdir=NULL, xout=NULL, yout=NULL){
    xseq = read_fasta(x)
    yseq = read_fasta(y)
    shared = intersect(names(xseq), names(yseq))

    if(length(shared) == 0)
        stop("Empty intersect")

    xseq = xseq[shared]
    yseq = yseq[shared]

    if(is.null(outdir))
        outdir = "."
    if(is.null(xout))
        xout = file.path(outdir, basename(x))
    if(is.null(yout))
        yout = file.path(outdir, basename(y))

    write_fasta(xseq, xout)
    write_fasta(yseq, yout)
    }
