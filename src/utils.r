#' utils.r
#'
#' shared utility functions
import::here("phyloRNA",
    "remove_constant",
    "write_fasta", "read_fasta",
    "tab2seq", "seq2tab",
    "all_files_exist", "mkdir", "corename"
    )
import::here("data.table", "fread")


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


#' Convert a tabular file into a fasta format
#'
#' @param file one or more files in tabular format
#' @param fasta **optional** output path for fasta files
#' @param outdir **optional** an outut directory, if fasta is not specified
#' @param margin whether rows (1) or columns (2) should be concatenated
#' @return a vector of fasta files
table2fasta = function(file, fasta=NULL, outdir=NULL, margin=2, zero=NULL){
    if(!is.null(fasta) && length(fasta) != length(file))
        stop("The file and fasta vectors must have the same length!")
    if(is.null(fasta))
        fasta = paste0(tools::file_path_sans_ext(file), ".fasta")
    if(!is.null(outdir))
        fasta = file.path(outdir, basename(fasta))
    mkdir(outdir)

    if(all_files_exist(fasta))
        return(invisible(fasta))

    for(i in seq_along(file)){
        data = read_table(file[i])
        if(!is.null(zero))
            data[data == zero] = 0
        seq = tab2seq(data, margin=margin)
        write_fasta(seq, file=fasta[i])
        }

    invisible(fasta)
    }


fasta2stats = function(fasta, stats=NULL, name=TRUE, unknown="N"){
    if(is.null(stats))
        stats = paste0(tools::file_path_sans_ext(fasta), ".txt")
    if(length(fasta) != length(stats))
        stop("fasta and stats vector must have the same length")
    if(all_files_exist(stats))
        return(invisible())

    if(isTRUE(name))
        name = corename(fasta)

    n = length(fasta)
    name = rep_len(name, n)
    unknown = rep_len(unknown, n)

    for(i in seq_along(fasta)){
        seq = read_fasta(fasta[i])
        tab = seq2tab(seq)

        text = paste0(
            "Sequences: ", nrow(tab), "\n",
            "Sites: ", ncol(tab), "\n",
            "Unique patterns: ", ncol(unique.matrix(tab, MARGIN=2)), "\n",
            "Data density: ", mdensity(tab, empty=unknown[i])
            )

        if(is.character(name))
            text = paste0("Name: ", name[i], "\n", text)

        writeLines(text, stats[i])
        }
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


#' Test if an element is not unknown
#'
#' @param x vector or matrix
#' @param empty an unknown element
#' @return vector or matrix 
is_empty = function(x, empty){
    if(is.na(empty))
        return(!is.na(x))
    x != empty
    }


#' Download file
#'
#' @param url an url from which file is downloaded
#' @param file a character string where downloaded file will be saved
download_file = function(url, file, rewrite=FALSE){
    if(!file.exists(file) || rewrite)
        download.file(url, file)
    }


#' Construct a filename
#'
#' A simple shorthand for construction a file name
filename = function(prefix, suffix="", ext=".txt", outdir="."){
    file.path(outdir, paste0(prefix, suffix, ext))
    }


#' Convert readgroup to a cell barcode
#'
#' Some single-cell detection methodology, such as the `phyloRNA::vcm()` tool, expect that every
#' read is barcoded with a cell-specific barcode. This functions transform non-barcoded single-cell
#' bam file into a barcoded bam file by adding the read-group (RG) to the cell barcode (CB) tag
#' @param input a bam file with reads encoded with RG tag
#' @param output an output bam file with read-group written into the CB tag
rg2cb = function(input, output){
    if(file.exists(output))
        return(invisible())

    command = "python3"
    args = c(
        "src/rg2cb.py",
        input, output
        )
    phyloRNA:::systemE(command, args)
    }


#' Read the vcm file and memoise it
#'
#' This function is memoised (possible reuse) of the vcm file.
#' `data.table::fread()` is used here due to a huge file size.
#' @param vcm a variant call matrix file
#' @return a variant call matrix as a data.table
read_vcm = local({
    memory = list()

    function(vcm){
        if(!is.null(memory[[vcm]]))
            return(memory[[vcm]])

        # using data.tale due to a huge size of the dataset
        data = fread(vcm, header=TRUE)
        # first three columns are not cells (chromosome, position and reference)
        data = data[, -c(1:3)]
        memory[[vcm]] <<- data

        data
        }
    })


#' Convert vcm file to fasta file
#'
#' Convert a vcm fle to fasta file.
#'
#' @param vcm a variant call matrix file
#' @param fasta a fasta file
#' @param selection selected cell barcodes that will be retained
vcm2fasta = function(vcm, fasta, selection=NULL){
    if(file.exists(fasta))
        return(invisible())

    data = read_vcm(vcm)

    if(!is.null(selection)){
        match = selection %in% colnames(data)
        if(!all(match)){
            warning("WARNING: ", sum(!match), " out of ", length(match),
                " requested cells are not present:\n",
                paste0(selection[!match], collapse="\n")
                )
            selection = selection[match]
            }
        data = data[, ..selection] # data.table' subsetting
        }

    data = as.matrix(data)
    data = remove_constant(data)
    seq = tab2seq(data, 2)
    write_fasta(seq, fasta)
    }


#' Substitute a pattern in a file
#'
#' Substitute a pattern over all lines in a file.
#'
#' This is a simple combination of the `sub` replacement function with `readLines` and `writeLines`.
#'
#' @param input an input file
#' @parma output an output file where pattern will be replaced
#' @param pattern a pattern to be replaced
#' @param replace a replacement for pattern
#' @param fixed pattern is a simple character string, not a regular expression
file_sub = function(input, output, pattern, replace, fixed=FALSE){
    lines = readLines(input)
    lines = sub(pattern, replace, lines, fixed=fixed)
    writeLines(lines, output)
    }


#' Merge files
#'
#' Merge multiple files into a single file.
#'
#' @param inputs one or multiple files to merge
#' @param a merged file
#' @param overwrite **optional** if an existing output should be overwritten
merge_files = function(inputs, output, overwrite=FALSE){
    if(file.exists(output) && overwrite)
        file.remove(output)
    if(!file.exists(output))
        file.append(output, inputs)
    }
