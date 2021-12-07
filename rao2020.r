#' Rao2020
#'
#' Reconstruct phylogenetic trees for data from:
#' Rao et al. (2020) "Comparative single-cell RNA sequencing (scRNA-seq) reveals
#' liver metastasis-specific targets in a patient with small intestinal neuroendocrine cancer"
#' Cold Spring Harb Mol Case Stud
import::from("src/sra.r", "get_srr_samples", "srr_download_sample")
import::from("src/utils.r", "merge_files", "file_sub", "fasta2stats")
import::from("phyloRNA",
    "cellranger_count", "cellranger_mkref", "bamtagregex",
    "gatk_prepare", "gatk_MergeSamFiles",
    "expr_read10x", "expr_merge",
    "mkdir", "read_fasta", "all_files_exist"
    )
import::from("src/expr.r", "process_expression", "expr2fasta")
import::from("src/snv.r", "detect_snv")
import::from("src/iqtree.r", "iqtrees")
import::from("src/beast.r", "beasts")
import::from("parallel", "mcMap")
import::from("magrittr", "%>%")

main = function(){
    expr()

    snv()
    }



expr = function(){
    outdir = "rao2020/expr"
    hdi = c(0.6, 0.9)

    primary_url = paste0(
        "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4159nnn/GSM4159164/suppl/",
        "GSM4159164%5FPriNET%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Etar%2Egz"
        )
    metastasis_url = paste0(
        "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4159nnn/GSM4159165/suppl/",
        "GSM4159165%5FlivMET%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Etar%2Egz"
        )

    primary_dir = file.path(outdir, "primary")
    metastasis_dir = file.path(outdir, "metastasis")

    cell_file = file.path(outdir, "cells.txt")
    umap_file = file.path(outdir, "umap.pdf")

    fastadir = file.path(outdir, "fasta")
    all_fasta = file.path(fastadir, "all.fasta")
    cancer_fasta = file.path(fastadir, "cancer.fasta")
    mkdir(fastadir)

    # Download count matrix
    download_count_matrix(primary_url, primary_dir)
    download_count_matrix(metastasis_url, metastasis_dir)

    if(!all_files_exist(c(cell_file, umap_file)))
        replicate_rao(primary_dir, metastasis_dir, umap_file, cell_file)

    cells = read.table(cell_file, header=TRUE)

    all_cells = cells$cells
    # select cancer cells, clusters: 0, 1, 3, 4
    cancer_cells = cells$cells[cells$clusters %in% c(0,1,3,4)]

    data = c("pr" = primary_dir, "met" = metastasis_dir)
 
    if(!file.exists(all_fasta))
        prepare_fasta(data, all_fasta, all_cells, n=500, hdi=hdi)

    if(!file.exists(cancer_fasta))
        prepare_fasta(data, cancer_fasta, cancer_cells, n=500, hdi=hdi)

    iqtrees(
        c(all_fasta, cancer_fasta),
        model = "ORDERED+ASC", 
        outdir = file.path(outdir, "trees", "ML"),
        mc.cores=2
        )
    }


####################################################################################################
# Replicating Rao 2020 with Seurat:
####################################################################################################
replicate_rao = function(primary_dir, metastasis_dir, umap_file, cell_file){
    primary = Seurat::Read10X(primary_dir, strip.suffix=TRUE)
    primary = Seurat::CreateSeuratObject(counts=primary, project="primary")
    primary = Seurat::RenameCells(primary, new.names = paste0(colnames(primary), "-pr"))

    metastasis =  Seurat::Read10X(metastasis_dir, strip.suffix=TRUE)
    metastasis = Seurat::CreateSeuratObject(counts=metastasis, project="metastasis")
    metastasis = Seurat::RenameCells(metastasis, new.names = paste0(colnames(metastasis), "-met"))

    # Prepare data for integration: keep only cells with more than 500 genes
    primary = subset(primary, subset = nFeature_RNA > 500)
    primary = Seurat::NormalizeData(primary, verbose = FALSE)
    primary = Seurat::FindVariableFeatures(primary, selection.method = "vst", nfeatures = 2000)
    primary$status = "Primary"

    metastasis = subset(metastasis, subset = nFeature_RNA > 500)
    metastasis = Seurat::NormalizeData(metastasis, verbose = FALSE)
    metastasis = Seurat::FindVariableFeatures(metastasis, selection.method = "vst", nfeatures = 2000)
    metastasis$status = "Metastatic"

    # Find anchors (shared cells types) and integrate the data along the anchors
    anchors = Seurat::FindIntegrationAnchors(list(primary, metastasis), dims=1:20)
    combined = Seurat::IntegrateData(anchors=anchors, dims=1:20)
    Seurat::DefaultAssay(combined) = "integrated"

    # Run the standard workflow for visualization and clustering
    combined = Seurat::ScaleData(combined, verbose = FALSE)
    combined = Seurat::RunPCA(combined, npcs = 30, verbose = FALSE)

    # t-SNE and Clustering
    combined = Seurat::RunUMAP(combined, reduction = "pca", dims = 1:20)
    combined = Seurat::FindNeighbors(combined, reduction = "pca", dims = 1:20)
    combined = Seurat::FindClusters(combined, resolution = 0.5)

    # Plots
    pdf(umap_file)
    Seurat::DimPlot(object = combined, reduction = "umap", split.by = "status", label=TRUE, label.size=5)
    invisible(dev.off())

    # Write down cell lists
    cells = data.frame("cells" = names(combined$seurat_clusters), "clusters" = combined$seurat_clusters)
    write.table(cells, cell_file, col.names=TRUE, row.names=FALSE)
    }


snv = function(){
    # define directories
    outdir = "rao2020"
    fastqdir = file.path(outdir, "raw")
    refdir = "reference/ref" # shared cellranger ref
    mapdir = file.path(outdir, "map")
    snvdir = file.path(outdir, "snv")
    fastadir = file.path(snvdir,"fasta")
    treedir = file.path(snvdir, "tree")

    gse = "GSE140312"
    
    # required reference files:
    reference = "reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
    annotation = "reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gtf"
    vcf = "reference/00-common_all.vcf.gz"
    pon = "pon/1000g_pon.hg38.vcf.gz"
    germline = "pon/af-only-gnomad.hg38.vcf.gz"

    # create annotation:
    cellranger_mkref(reference, annotation, refdir, nthreads=8)

    # get samples
    samples = get_srr_samples(gse, save=file.path(outdir, "samples.rds"))
    
    # construct name from samples
    samples$name = strsplit(samples$name, split=" ", fixed=TRUE) %>% sapply(getElement, 1)
    
    # download the sequences in parallel (fastq-dump and prefetch are slow)
    ncores = max(nrow(samples), 16)
    mcMap(srr_download_sample, samples$srr, samples$name, fastqdir, mc.cores=ncores)

    # map and demultiplex
    mkdir(mapdir)
    outputs = mcMap(
        cellranger_count,
        id = samples$name, sample = samples$name,
        fastqdir = fastqdir, refdir = refdir,
        outdir = file.path(mapdir, samples$name),
        nthreads=4, mc.cores=ncores
        )

    # Preproces
    bam_aligned = sapply(outputs, getElement, "bam")
    bam_prepared = file.path(mapdir, paste0(samples$name, ".bam"))
    
    barcodes_aligned = sapply(outputs, getElement, "barcodes")
    barcodes_prepared = file.path(mapdir, paste0(samples$name, ".txt"))

    pattern = "-1$"
    replace = paste0("-", samples$name)

    mcMap(prepare_bam,
        input = bam_aligned,
        output = bam_prepared,
        reference = reference,
        vcf = vcf,
        barcodes = barcodes_aligned,
        pattern = pattern,
        replace = replace,
        mc.cores=ncores
        )

    mcMap(file_sub,
        barcodes_aligned,
        barcodes_prepared,
        pattern,
        repalce = list(replace),
        mc.cores=ncores
        )

    # Merge bam and barcodes
    bam = file.path(mapdir, "all.bam")
    barcodes = file.path(mapdir, "all.txt")
    gatk_MergeSamFiles(bam_prepared, bam)
    merge_files(barcodes_prepared, barcodes)

    phyloRNA::gatk_MergeSamFiles(bam_prepared, bam)
    merge_files(barcodes_prepared, barcodes)

    # detect SNVs
    vcm = detect_snv(bam, barcodes, reference, pon=pon, germline=germline, outdir=snvdir)

    # Create two datasets based on selected cells from expression analysis
    expr_cancer_fasta = file.path(outdir, "expr", "cancer.fasta")
    expr_all_fasta = file.path(outdir, "expr", "all.fasta")
    if(!all_files_exist(c(expr_cancer_fasta, expr_all_fasta)))
        stop("Run expression analysis first.")

    cancer_names = get_names(expr_cancer_fasta)
    all_names = get_names(expr_all_fasta)

    mkdir(fastadir)
    cancer_fasta = file.path(outdir, "snv", "cancer.fasta")
    all_fasta = file.path(outdir, "snv", "all.fasta")

    vcm2fasta(vcm, cancer_fasta, cancer_names)
    fasta2stats(cancer_fasta, unknown="N")

    vcm2fasta(vcm, all_fasta, all_names)
    fasta2stats(all_fasta, unknown="N")

    # Run phylogenetic analyses
    fasta = c(cancer_fasta, all_fasta)
    iqtrees(fasta, "TEST", outdir=file.path(treedir, "ML"), mc.cores=2)
    template = file.path("templates", "BDStrictGtr.xml")
    beasts(fasta, template, outdir=file.path(treedir, "BI"), mc.cores=2)
    }


prepare_bam = function(input, output, reference, vcf, barcodes, pattern, replace){
    if(file.exists(output))
        invisible(stop)

    intermediate = paste0(tools::file_path_sans_ext(input), ".intermediate.bam")

    gatk_prepare(
        input, intermediate,
        reference = reference,
        vcf = vcf,
        barcodes = barcodes
        )
    bamtagregex(intermediate, output, "CB", pattern, replace)
    file.remove(intermediate)
    }


get_names = function(fasta){
    seq = read_fasta(fasta)
    names = names(seq)
    names = sub("pr$", "primary", names)
    names = sub("met$", "metastatic", names)
    names
    }


download_count_matrix = function(url, dir){
    if(!dir.exists(dir)){
        tmp = tempfile()
        download.file(url, tmp)
        untar(tmp, exdir = dirname(dir))
        file.rename(file.path(dirname(dir), "filtered_feature_bc_matrix"), dir)
        }
    }


#' Select elements according to pattern
#'
#' Select n first elements from a vector that match the pattern
#'
#' @param x a vector
#' @param pattern a pattern according to which elements will be selected
#' @param n **optional** a number of elements returned
select_pattern = function(x, pattern, n=NULL){
    y = grep(pattern, x, value=TRUE)
    if(is.null(n))
        return(y)

    n = min(n, length(y))
    y[seq_len(n)]
    }


prepare_fasta = function(x, fasta, selection, n, hdi=c(0.6, 0.9), suffix=NULL){
    if(is.null(suffix))
        suffix = names(x)
    if(is.null(suffix))
        suffix = seq_along(x)

    data = Map(function(y) expr_read10x(y)[[1]], x)
    data = expr_merge(data, suffix)
    data = data[, selection]
    
    data = process_expression(data, hdi, trim=TRUE)

    best = colSums(data != "-")
    best = sort(best, decreasing=TRUE)
    best = names(best)

    best = Map(select_pattern, list(best), paste0("-", suffix), n)
    best = unlist(best)

    data = data[, best]
    expr2fasta(data, fasta, summary=TRUE, process=FALSE)
    }


if(sys.nframe() == 0){
    main()
    }
