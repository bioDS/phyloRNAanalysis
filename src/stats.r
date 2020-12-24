#' stats.r
#'
#' Calculate various summary statistics

# What do I need?
# -- Dimension and the data density of each sample
# -- Dimension and the data density of each filtered dataset
library("Matrix")
library("data.table")

expr_h5_stats = function(h5, output){
    names = phyloRNA::corename(h5)
    data = lapply(h5, phyloRNA::expr_read10xh5)
    stats = lapply(data, expr_h5_stats_dataset)
    names(stats) = names
    data = phyloRNA::expr_merge(data, names)
    stats = c(stats, "total"=list(expr_h5_stats_dataset(data)))
    table = do.call(rbind, stats)
    write_table(table, output)
    }

genes_per_cell = function(x){
    mean( colSums(x != 0) )
    }

umi_per_cell = function(x){
    mean( colSums(x) )
    }

umi_per_gene = function(x){
    mean( rowSums(x) )
    }

expr_h5_stats_dataset = function(data){
    stats = list(
        "ncells" = ncol(data),
        "ngenes" = nrow(data),
        "ngenes_detected" = sum(rowSums(data) > 0),
        "genes_per_cell" = genes_per_cell(data),
        "umi" = sum(data),
        "umi_per_gene" = umi_per_gene(data),
        "umi_per_cell" = umi_per_cell(data),
        "density" = mdensity(data, 0)
        )
    return(stats)
    }

vcm_stats = function(vcm, output){
    data = data.table::fread(vcm)
    data = data[,-c(1:3)]

    empty_cells = colnames(data)[colSums(is_empty(data, "N")) == 0]
    stats = list(
        "ncells" = ncol(data),
        "SNVs" = nrow(data),
        "diversity" = diversity(colnames(data)),
        "density" = mdensity(data, "N"),
        "nEmpty" = length(empty_cells),
        "divEmpty" = diversity(empty_cells)
        )
    stats = do.call(rbind, stats)
    write_table(stats, output)
    }

filtered_stats = function(files, empty, output){
    stats = list()
    for(file in files){
        name = corename(file)
        data = read_table(file)
        stats[[name]] = filtered_stats_dataset(data, empty)
        }
    stats = do.call(rbind, stats)
    write_table(stats, output) 
    }

diversity_stats = function(files, output){
    stats = list()
    for(file in files){
        name = corename(file)
        data = read_table(file)
        stats[[name]] = diversity(colnames(data))
        }
    stats = do.call(rbind, stats)
    write_table(stats, output)
    }

filtered_stats_dataset = function(data, empty){
    stats = list(
        "ncells" = ncol(data),
        "nrow" = nrow(data),
        "density" = mdensity(data, empty)
        )
    return(stats)
    }

diversity = function(names){
    names = sub(".*-", "", names)
    tab = table(names)
    div = paste(names(tab), tab, sep=":", collapse=", ")
    div
    }

expr_stats = function(h5, filtered, outdir=NULL){
    if(is.null(outdir))
        outdir = "."
    mkdir(outdir)

    h5_stats = file.path(outdir, "h5_stats.txt")
    if(!file.exists(h5_stats))
        expr_h5_stats(h5, h5_stats)

    expr_stats = file.path(outdir, "expr_stats.txt")
    if(!file.exists(expr_stats))
        filtered_stats(filtered, empty="-", expr_stats)

    expr_diversity = file.path(outdir, "expr_diversity.txt")
    if(!file.exists(expr_diversity))
        diversity_stats(filtered, expr_diversity)
    }

snv_stats = function(vcm, filtered, outdir=NULL){
    if(is.null(outdir))
        outdir = "."
    mkdir(outdir)

    vcm_stats = file.path(outdir, "vcm_stats.txt")
    if(!file.exists(vcm_stats))
        vcm_stats(vcm, vcm_stats)

    snv_stats = file.path(outdir, "snv_stats.txt")
    if(!file.exists(snv_stats))
        filtered_stats(filtered, empty="N", snv_stats)

    snv_diversity = file.path(outdir, "snv_diversity.txt")
    if(!file.exists(snv_diversity))
        diversity_stats(filtered, snv_diversity)
    }
