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
    nthreads = "AUTO", parallel = FALSE,
    remake = FALSE
    ){
    if(is.null(outdir))
        outdir = "."
    phyloRNA::mkdir(outdir)

    cfasta = file.path(outdir, basename(fasta))

    # Early exit if tree exists and not running parallel bootstrap:
    # (there the trees are more complicated)
    tree = paste0(cfasta, ".treefile")
    if(!remake && !parallel && file.exists(tree))
        return(invisible(tree))

    file.copy(fasta, cfasta)

    command = "iqtree"
    args = c(
        "-s", cfasta,
        "-nt", nthreads
        )
    if(!is.null(model))
        args = c(args, "--model", model)
    if(remake)
        args = c(args, "--redo")

    # Run iqtree with the ultrafast bootstrap
    if(ufboot){
        if(is.logical(ufboot)) ufboot = 1000
        args = c(args, "-B", ufboot)
        phyloRNA:::systemE(command, args)
        return(invisible(tree))
        }

    # Run iqtree with the ultrafast bootstrap, but don't parallelize
    if(bootstrap && !parallel){
        if(is.logical(bootstrap)) bootstrap = 100
        args = c(args, "-b", bootstrap)
        phyloRNA:::systemE(command, args)
        return(invisible(tree))
        }

    # run iqtree with parallelized standard bootstrap
    if(bootstrap && parallel){
        if(nthreads == "AUTO")
            nthreads = 4
        trees = iqtree_par(fasta, model, outdir, bootstrap, nthreads, remake)
        return(invisible(trees))
        }

    # run iqtree without boostrap:
    phyloRNA:::systemE(command, args)
    return(invisible(tree))
    }


#' Performs a single bootstrap run with iqtree
#'
#' @param i run id
#' @param fasta a path to fasta file
#' @param a model to be run with
iqboot = function(i, fasta, model = NULL, remake = FALSE){
    command = "iqtree"
    prefix = file.path(dirname(fasta), i)
    tree = paste0(prefix, ".boottrees")

    if(!remake && file.exists(tree))
        return(tree)

    args = c(
        "-s", fasta,
        "-nt", 1,
        "-bo 1",
        "-pre", prefix,
        "-quiet"
        )

    if(remake)
        args = c(args, "--redo")

    if(!is.null(model))
        args = c(args, "--model", model)
    system2(command, args)
    return(tree)
    }


#' Run iqtree with parallelized bootstrap
#'
#' @param fasta an input fasta file
#' @param model an iqtree model, such as HKY, GTR and others
#' @param outdir an output directory
#' @param bootstrap a number of bootstrap replicates
#' @param nthreads a number of threads to run in parallel
#' @return a vector of paths to best tree, bootstrap trees and conensus bootstrap tree
iqtree_par = function(
    fasta,
    model = NULL, outdir = NULL,
    bootstrap = 100, nthreads = 8,
    remake = FALSE
    ){
    command = "iqtree"

    if(is.null(outdir))
        outdir = "."

    bootdir = file.path(outdir, "bootstrap")
    phyloRNA::mkdir(outdir)
    phyloRNA::mkdir(bootdir)

    core = phyloRNA:::corename(fasta)
    cfasta = file.path(outdir, basename(fasta))
    cbfasta = file.path(bootdir, basename(fasta))

    file.copy(fasta, cfasta)
    file.copy(fasta, cbfasta)

    # Performs a single standard run single run
    tree = iqtree(fasta, model, outdir)

    # perform parallelized bootstrap
    trees = parallel::mclapply(
        seq_len(bootstrap),
        iqboot,
        model = model, fasta = cbfasta, remake = remake,
        mc.cores = nthreads
        )

    # merge bootstrap files
    boottrees = file.path(outdir, paste0(core, ".boottrees"))
    file.append(boottrees, unlist(trees))

    # calculate consensus tree from bootstrap trees
    contree = iqtree_consensus(boottrees)

    # calculate bootstrap support for the best and consensus trees
    bootsup = iqtree_support(contree, boottrees, remake=remake)
    treesup = iqtree_support(tree, boottrees, remake=remake)

    return(c("tree"=treesup, "bootstrap"=boottrees, "consensus"=bootsup))
    }


#' Annotate tree with a boostrap support values
#'
#' Annotate tree with a bootstrap support values using a set of bootstrap trees.
#'
#' @param tree a tree that will be annotated
#' @param bootstrap a set of bootstrap trees
#' @param output **optional** an output file
#' @return a file with splits annotated with a bootstrap support
iqtree_support = function(tree, bootstrap, output=paste0(tree, ".sup"), remake=FALSE){
    if(!remake && file.exists(output))
        return(output)

    command = "iqtree"
    phyloRNA:::systemE(command, c("-sup", tree, "-t", bootstrap))

    annotated = paste0(bootstrap, ".suptree")
    file.rename(annotated, output)
    return(output)
    }


#' Calculate consensus tree from bootstrap trees
#'
#' @param bootstrap a file with bootstrap trees
#' @param output **optional** an output file
#' @return a file with conensus tree
iqtree_consensus = function(bootstrap, output=paste0(bootstrap, ".con"), remake=FALSE){
    if(!remake && file.exists(output))
        return(output)

    command = "iqtree"
    phyloRNA:::systemE(command, c("-con -t", bootstrap))

    consensus = paste0(bootstrap, ".contree")
    file.rename(consensus, output)
    return(output)
    }


iqtrees = function(
    fasta,
    model = NULL, outdir = NULL,
    ufboot = FALSE, bootstrap = FALSE,
    nthreads = "AUTO", parallel = FALSE,
    remake = FALSE
    ){
    if(is.null(outdir))
        outdir = "."
    phyloRNA::mkdir(outdir)

    n = length(fasta)

    if(is.null(model) || length(model) == 1)
        model = rep(list(model), n)

    if(length(outdir) == 1)
        outdir = file.path(outdir, phyloRNA:::corename(fasta))

    trees = Map(f=iqtree, fasta, model, outdir, ufboot, bootstrap, nthreads, parallel, remake)
    invisible(trees)
    }
