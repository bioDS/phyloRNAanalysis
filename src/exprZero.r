table2fastaZero = function(file, fasta=NULL, outdir=NULL, margin=2){
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
        data[data == "-"] = "0"
        phyloRNA::fasta(data, margin=margin, file=fasta[i])
        }

    invisible(fasta)
    }

