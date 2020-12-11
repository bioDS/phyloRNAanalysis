#' beast.r
#'
#' Functions for phylogenetic reconstruction with beast

#' Run BEAST analysis
#'
#' @param fasta a fasta file
#' @param template beter template to generate BEAST XML
#' @param outdir 
#' @param nthreads a number of threads to run the BEAST on
#' @param burnin burnin percentages
#' @param params an additional list of parameters
beast = function(fasta, template, outdir=NULL, nthreads=2, burnin=20, params=list()){
    if(is.null(outdir))
        outdir = "."
    phyloRNA::mkdir(outdir)

    prefix = better:::basename_sans_ext(fasta)

    beastxml = paste0(prefix, ".xml")
    beter::(template, file.path(outdir, beastxml), alignment=fasta, parameters=params)

    beast_args = c(
        "-beagle_CPU",
        "-threads", nthreads,
        beastxml
        )
    phyloRNA:::systemE("beast", args, dir=outdir)

    log_args = c(
        "-b", burnin,
        paste0(prefix, ".trace"),
        ">",
        paste0(prefix, ".txt")
        )
    phyloRNA:::systemE("loganalyser", log_args, dir=outdir)

    tree_args = c(
        "-b", burnin,
        paste0(prefix, ".trees"),
        paste0(prefix, ".tree")
        )
    phyloRNA:::systemE("treeannotator", tree_args, dir=outdir)
    }


beasts = function(fastas, template, outdir=NULL, nthreads=2, burnin=20, param=list()){
    if(is.null(outdir))
        outdir = "."
    phyloRNA::mkdir(outdir)

    for(fasta in fastas){
        beast(
            fasta = fasta,
            template = template,
            outdir = file.path(outdir, better:::basename_sans_ext(fasta)),
            nthreads = nthreads,
            burnin = burnin,
            param = param)
        }
    }
