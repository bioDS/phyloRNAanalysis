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
iqtree = function(
    fasta,
    model = NULL, outdir = NULL,
    bootstrap = 1000, ufboot = TRUE,
    nthreads = "AUTO", remake = FALSE
    ){
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
        "-nt", nthreads
        )
    if(!is.null(model))
        args = c(args, "--model", model)
    if(ufboot)
        args = c(args, "-B", bootstrap)
    if(!ufboot)
        args = c(args, "-b", bootstrap)
    if(remake)
        args = c(args, "--redo")
        
    phyloRNA:::systemE(command, args)
    }


iqtrees = function(fastas, model=NULL, outdir=NULL, bootstrap=1000, ufboot=TRUE, nthreads="AUTO"){
    if(is.null(outdir))
        outdir = "."
    phyloRNA::mkdir(outdir)

    for(fasta in fastas){
        subdir = file.path(outdir, basename(tools::file_path_sans_ext(fasta)))
        iqtree(fasta, model, subdir, bootstrap, ufboot, nthreads)
        }
    }


iqtree_partition = function(
    fasta, model,
    outdir = NULL,
    bootstrap = 1000,
    nthreads = 8,
    iter = 1000,
    remake = FALSE,
    intersect = FALSE,
    empty = "N"
    ){
    if(is.null(outdir))
        outdir = "."
    phyloRNA::mkdir(outdir)

    cfasta = file.path(outdir, basename(fasta))

    if(!remake && all(file.exists(cfasta)))
        return(invisible())

    if(intersect){
        if(length(fasta) != 2)
            stop("This option works only for two fasta files")
        if(is.null(empty))
            empty = c("N", "N")
        fasta_intersect(
            fasta[1], fasta[2],
            xout=cfasta[1], yout=cfasta[2],
            remove_constant=TRUE, xempty=empty[1], yempty=empty[2]
            )
        } else {
        file.copy(fasta, cfasta)
        }
        
    nexus = file.path(outdir, "partition.nex")
    generate_partition(cfasta, model, nexus)

    command = "iqtree"

    args = c(
        "-p", nexus,
        "-nt", nthreads,
        "-B", bootstrap,
        "-nm", iter
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
    modelstring = paste(model, partition, sep=":", collapse=", ")

    text = c(
        "#nexus",
        "begin sets;",
        paste0("    charset ", partition, " = ", fasta, ": *;"),
        paste0("    charpartition ", name, " = ", modelstring, ";"),
        "end;"
        )
    
    writeLines(text, file)
    }


iqtrees_partition = function(
    arglist,
    outdir = NULL,
    bootstrap = 1000,
    nthreads = 8,
    intersect = FALSE,
    empty = NULL,
    iter = 1000
    ){
    if(is.null(outdir))
        outdir = "."
    phyloRNA::mkdir(outdir)
    for(args in arglist){
        iqtree_partition(
            fasta = args$fasta,
            model = args$model,
            outdir = file.path(outdir, args$outdir),
            bootstrap = bootstrap,
            nthreads = nthreads,
            intersect = intersect,
            empty = empty,
            iter = iter
            )
        }
    }


make_arglist = function(fastas, models, outdirs){
    mapply(
        function(fasta, model, outdir){ list("fasta"=fasta, "model"=model, "outdir"=outdir) },
        fastas, models, outdirs,
        SIMPLIFY=FALSE
        )
    }
