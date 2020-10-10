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
iqtree = function(fasta, outdir, model="ORDINAL+ASC", bootstrap=1000, nthreads=8){
    phyloRNA::mkdir(outdir)

    cfasta = file.path(outdir, basename(fasta))
    file.copy(fasta, cfasta)

    command = "iqtree"
    args = c(
        "-s", cfasta,
        "-nt", nthreads,
        "-m", model,
        "-B", bootstrap
        )
        
    phyloRNA:::systemE(command, args)
    }


iqtrees = function(fastas, outdir, subdir, model, bootstrap=1000, nthreads=8){
    for(i in seq_along(fastas)){
        iqtree(fastas[i], file.path(outdir, subdir[i]), model, bootstrap, nthreads)
        }
    }
