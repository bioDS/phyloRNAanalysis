#' This function will transform scRNAseq expression data into fasta sequences
#'
#'
#'
expr_process = function(
    data=NULL, file=NULL, names=NULL,
    dens=0.5, minGene=250, minUMI=500,
    outdir=NULL, fasta=NULL, discrete=NULL, barcodes=NULL
    ){
    if(is.null(data) && is.null(file))
        stop("Either data or file need to be specified!")
    if(!is.null(data) && !is.null(file))
        stop("Only data or file need to be specified, not both!")

    if(!is.null(file)){
        if(length(file) > 1){
            data = lapply(file, phyloRNA::expr_read10xh5)
            } else {
            data = phyloRNA::expr_read10xh5(file)
            }
        }


    if(is.list(data) && length(data) > 1){
        data = phyloRNA::expr_merge(data, names)
        }


    data = phyloRNA::expr_quality_filter(data, minGene=minGene, minUMI=minUMI)
    data = phyloRNA::expr_zero_to_na(data)
    data = phyloRNA::expr_normalize(data)
    data = phyloRNA::expr_scale(data)
    discr = phyloRNA::expr_discretize(data, intervals=c(-2,-0.5,0.5,2), unknown="-")
    discr = phyloRNA::remove_constant(discr, margin=1)

    if(!is.null(discrete))
        write.table(discr, quote=FALSE, sep="\t", col.names=TRUE, row.names=TRUE, file=discrete)

    discr = phyloRNA::densest_subset(discr, empty="-", steps=10000, density=dens)$result
    phyloRNA::fasta(discr, file=fasta)
    }
