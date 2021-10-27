#' iqtree.r
#'
#' Functions for phylogenetic reconstruction with iqtree


#' Run IQtree analysis
#'
#' @param fasta a fasta file
#' @param outdir an output directory, fasta will be copied there
#' @param model a phylogenetic model for IQtree
#' @param ufboot a number of replicates for the ultrafast bootstrap
#' @param bootstrap a number of replicates for the standard bootstrap
#' @param nthreads a number of threads to run the IQtree on
iqtree = function(
    fasta,
    model = NULL, outdir = NULL,
    ufboot = FALSE, bootstrap = FALSE,
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
    if(ufboot){
        if(is.logical(ufboot)) ufboot = 1000
        args = c(args, "-B", ufboot)
        }
    if(bootstrap){
        if(is.logical(bootstrap)) bootstrap = 100
        args = c(args, "-b", bootstrap)
        }
    if(remake)
        args = c(args, "--redo")
        
    phyloRNA:::systemE(command, args)
    }


iqboot = function(i, fasta, model = NULL){
    command = "iqtree"
    args = c(
        "-s", fasta,
        "-nt", 1,
        "-bo 1",
        "-pre", file.path(dirname(fasta),i),
        "-quiet"
        )
    if(!is.null(model))
        args = c(args, "--model", model)
    system2(command, args)
    }


iqtree_par = function(
    fasta,
    model = NULL, outdir = NULL,
    bootstrap = 100, nthreads = 8
    ){
    if(is.null(outdir))
        outdir = "."
    bootdir = file.path(outdir, "bootstrap")
    phyloRNA::mkdir(outdir)
    phyloRNA::mkdir(bootdir)

    cfasta = file.path(outdir, basename(fasta))
    cbfasta = file.path(bootdir, basename(fasta))

    #if(file.exists(cfasta))
    #    return(invisible())

    file.copy(fasta, cfasta)
    file.copy(fasta, cbfasta)

    # Performs a single run
    command = "iqtree"
    args = c(
        "-s", cfasta,
        "-nt AUTO"
        )
    if(!is.null(model))
        args = c(args, "--model", model)
    system2(command, args)

    # perform paralel version
    mclapply(
        seq_len(bootstrap),
        iqboot,
        model = model, fasta = cbfasta,
        mc.cores = nthreads
        )

    # merge bootstrap files
    outtree = paste0(cfasta, ".treefile")
    boottree = file.path(outdir, "boottrees.new")
    boottrees = dir(bootdir, pattern="boottrees$", full.names=TRUE)
    file.append(boottree, boottrees)
    # calculate consensus tree from bootstrap trees
    phyloRNA:::systemE(command, c("-con -t", boottree))
    # calculate bootstrap support on the ML tree
    phyloRNA:::systemE(command, c("-sup", outtree, "-t", boottree))
    }

iqtrees_par = function(fastas, model=NULL, outdir=NULL, bootstrap=100, nthreads=8){
    if(is.null(outdir))
        outdir = "."
    phyloRNA::mkdir(outdir)

    for(fasta in fastas){
        subdir = file.path(outdir, basename(tools::file_path_sans_ext(fasta)))
        iqtree_par(fasta, model, subdir, bootstrap, nthreads)
        }
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
