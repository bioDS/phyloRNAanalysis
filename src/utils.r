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


#' Transform fasta to table
#'
#' transform a named vector of sequences into a tabular format
#'
#' @param fasta a named vector of sequences, such as from `read_fasta`
#' @param margin **optional** connect by either rows or columns
#' @return a table of fasta sequences
fasta2table = function(fasta, margin=1){
    data = strsplit(fasta, "")
    if(margin == 1)
        data = do.call(rbind, data)
    if(margin == 2)
        data = do.call(cbind, data)
    if(!margin %in% c(1,2))
        stop("Margin must be either 1 or 2")

    data
    }

remove_constant_pos = function(fasta, unknown="N"){
    data = tasta2table(fasta, margin=1)
    data = phyloRNA::remove_constant(data, margin=2, unknown=unknown)
    data = apply(data, margin, function(x) paste0(x, collapse = ""))
    data = paste0(">", names(data), "\n", data)
    data
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
#' @param remove_constant **optional** remove constant sites after reducing the dataset?
#' @param xempty **optional** unknown data symbol for x fasta file
#' @param yempty **optional** unknown data symbol for y fasta file
fasta_intersect = function(
    x, y,
    outdir=NULL,
    xout=NULL, yout=NULL,
    remove_constant=FALSE,
    xempty="N", yempty="N"
    ){
    xseq = read_fasta(x)
    yseq = read_fasta(y)
    shared = intersect(names(xseq), names(yseq))

    if(length(shared) == 0)
        stop("Empty intersect")

    xseq = xseq[shared]
    yseq = yseq[shared]

    if(remove_constant){
        xseq = remove_constant(xseq, xempty)
        yseq = remove_constant(yseq, yempty)
        }

    if(is.null(outdir))
        outdir = "."
    if(is.null(xout))
        xout = file.path(outdir, basename(x))
    if(is.null(yout))
        yout = file.path(outdir, basename(y))

    write_fasta(xseq, xout)
    write_fasta(yseq, yout)
    }


#' Calculate data density of matrix
#'
#' Calculate the data density, that is the proportion of known elements in the matrix.
#' @param x a matrix
#' @param empty an unknown element
#' @return a data density of matrix
mdensity = function(x, empty){
    sum(is_empty(x, empty)) / prod(dim(x))
    }


#' Test if an element is unknown
#'
#' @param x vector or matrix
#' @param empty an unknown element
#' @return vector or matrix 
is_empty = function(x, empty){
    if(is.na(empty))
        return(!is.na(x))
    x != empty
    }
