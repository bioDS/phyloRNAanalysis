#' iqtree.r
#'
#' Functions for phylogenetic reconstruction with iqtree


#' Run IQtree analysis
#'
#' @param fasta a fasta file
#' @param outdir an output directory, fasta will be copied there
#' @param model a phylogenetic model for IQtree
#' @param bootstrap a number of replicates for the ultrafast bootstrap
#' @param nthreads a number of threads to run the IQtree on
iqtree = function(fasta, model, outdir, bootstrap=1000, nthreads=8, remake=FALSE){
    if(is.null(outdir))
        outdir = "."
    phyloRNA::mkdir(outdir)

    cfasta = file.path(outdir, basename(fasta))

    if(!remake && file.exists(cfasta))
        return(invisible())

    file.copy(fasta, cfasta)

    command = "iqtree"
    args = c(
        "-s", cfasta,
        "-nt", nthreads,
        "-m", model,
        "-B", bootstrap
        )
    if(remake)
        args = c(args, "--redo")
        
    phyloRNA:::systemE(command, args)
    }


iqtrees = function(fastas, model, outdir=NULL, bootstrap=1000, nthreads=8){
    if(is.null(outdir))
        outdir = "."
    phyloRNA::mkdir(outdir)

    for(fasta in fastas){
        subdir = file.path(outdir, basename(tools::file_path_sans_ext(fasta)))
        iqtree(fasta, subdir, model, bootstrap, nthreads)
        }
    }


iqtree_partition = function(fasta, model, outdir=NULL, bootstrap=1000, nthreads=8){
    if(is.null(outdir))
        outdir = "."
    phyloRNA::mkdir(outdir)

    cfasta = file.path(outdir, basename(fasta))

    if(!remake && all(file.exists(cfasta)))
        return(invisible())

    file.copy(fasta, cfasta)
    nexus = file.path(outdir, "parititon.nex")
    generate_partition(fasta, model, nexus)

    command = "iqtree"

    args = c(
        "-p", nexus,
        "-nt", nthreads,
        "-B", bootstrap
        )

    if(remake)
        args = c(args, "--redo")

    phyloRNA:::systemE(command, args)
    }


#' Generate IQtree partition file
#'
#' Generate Nexus partition file for IQtree.
#'
#' @param fasta files inclded in partition,
#' @param model model for each partition
#' @param file an output nexus file
#' @param name **optional** name for the partition set
generate_partition = function(fasta, model, file, name=NULL){
    if(is.null(name))
        name = "partition"
    if(length(fasta) != length(model))
        stop("Length of fasta files and model strings differ!")

    partition = paste0("part", seq_along(fasta))
    modelstring = paste(model, charset, sep=":", collapse=", ")

    text = c(
        "#nexus",
        "begin sets;",
        paste0("    charset ", partition, " = ", fasta, ": *;"),
        paste0("    ", name, " = ", modelstring, ";"),
        "end;"
        )
    
    writeLines(text, file)
    }


iqtrees_partition = function(arglist, outdir=NULL, bootstrap=1000, nthreads=8){
    if(is.null(outdir))
        outdir = "."
    phyloRNA::mkdir(outdir)
    for(args in arglist){
        iqtree_partition(
            fasta = args$fasta,
            model = args$model,
            outdir = file.path(outdir, args$outdir),
            bootstrap = 1000,
            nthreads = 8
            )
        }
    }


make_arglist = function(fastas, models, outdirs){
    mapply(
        function(fasta, model, outdirs){ list(fasta=fastas, model=models, outdir=outdirs) },
        fastas, models, outdirs
        )
    }
