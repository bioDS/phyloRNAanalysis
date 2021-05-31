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

    prefix = beter:::basename_sans_ext(fasta)

    beastxml = file.path(outdir, paste0(prefix, ".xml"))
    trace = file.path(outdir, paste0(prefix, ".trace"))
    trees = file.path(outdir, paste0(prefix, ".trees"))
    ess = file.path(outdir, paste0(prefix, ".txt"))
    tree = file.path(outdir, paste0(prefix, ".tree"))

    if(!file.exists(beastxml))
        beter::process_template(
            template,
            beastxml,
            alignment = fasta,
            parameters = params
            )

    beast_args = c(
        "-threads", nthreads,
        basename(beastxml)
        )
    if(!file.exists(trace))
        phyloRNA:::systemE("beast", beast_args, dir=outdir)

    log_args = c(
        "-b", burnin,
        basename(trace),
        ">",
        basename(ess)
        )
    if(!file.exists(ess))
        phyloRNA:::systemE("loganalyser", log_args, dir=outdir)

    tree_args = c(
        "-b", burnin,
        "-lowMem",
        basename(trees),
        basename(tree)
        )

    if(!file.exists(tree))
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
            outdir = file.path(outdir, beter:::basename_sans_ext(fasta)),
            nthreads = nthreads,
            burnin = burnin,
            param = param
            )
        }
    }
