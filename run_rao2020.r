library("Seurat")
library("phyloRNA")
library("parallel")
import::from("src/utils.r", "mdensity")
import::from("src/expr.r", "calculate_intervals")
import::from("src/iqtree.r", "iqtree")
import::from("src/beast.r", "beast")

main = function(){
    # Define variables
    hdi = c(0.6, 0.9)
    outdir = "rao2020"
    phyloRNA::mkdir(outdir)

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

    all_fasta = file.path(outdir, "rao2020_all.fasta")
    cancer_fasta = file.path(outdir, "rao2020_cancer.fasta")

    # Download count matrix
    download_count_matrix(primary_url, primary_dir)
    download_count_matrix(metastasis_url, metastasis_dir)

    if(!file.exists(cell_file) && !file.exists(umap_file)){
        replicate_rao(primary_dir, metastasis_dir, umap_file, cell_file)
        }

    cells = read.table(cell_file, header=TRUE)

    ################################################################################################
    # Prepare fasta files
    ################################################################################################
    if(!file.exists(all_fasta)){
        primary = phyloRNA::expr_read10x(primary_dir)[[1]]
        metastasis = phyloRNA::expr_read10x(metastasis_dir)[[1]]
        combined = phyloRNA::expr_merge(list(primary, metastasis), c("pr", "met"))
        combined = combined[,cells$cells] # Select cells from rao
        combined = combined[rowSums(combined) != 0, ] # remove empty genes
   
        data = phyloRNA::expr_zero_to_na(combined)
        data = phyloRNA::expr_scale(data)
        intervals = calculate_intervals(data, density=hdi)
        data = phyloRNA::expr_discretize(data, intervals=intervals, unknown="-")
        
        # Select 1000 cells:
        best = names(colSums(data != "-"))[1:1000]
        data = data[, best]
        data = phyloRNA::remove_constant(data, margin=1, unknown="-")
        summary(data, file.path(outdir, "all_summary.txt"))
        data = phyloRNA::tab2seq(data, margin=2)
        phyloRNA::write_fasta(data, all_fasta)
        }

    if(!file.exists(cancer_fasta)){
        primary = phyloRNA::expr_read10x(primary_dir)[[1]]
        metastasis = phyloRNA::expr_read10x(metastasis_dir)[[1]]
        combined = phyloRNA::expr_merge(list(primary, metastasis), c("pr", "met"))
        
        # select cancer cells, clusters: 0, 1, 3, 4
        cancer = cells$cells[cells$clusters %in% c(0,1,3,4)]
        combined = combined[,cancer]
        combined = combined[rowSums(combined) != 0, ] # remove empty genes

        data = phyloRNA::expr_zero_to_na(combined)
        data = phyloRNA::expr_scale(data)
        intervals = calculate_intervals(data, density=hdi)
        data = phyloRNA::expr_discretize(data, intervals=intervals, unknown="-")
        
        # Select 1000 cells:
        best = names(colSums(data != "-"))[1:1000]
        data = data[, best]
        data = phyloRNA::remove_constant(data, margin=1, unknown="-")
        summary(data, file.path(outdir, "cancer_summary.txt"))
        data = phyloRNA::tab2seq(data, margin=2)
        phyloRNA::write_fasta(data, cancer_fasta)
        }

    iqtree(
        all_fasta,
        model="ORDERED+ASC",
        file.path(outdir, "all")
        )

    iqtree(
        cancer_fasta,
        model="ORDERED+ASC",
        file.path(outdir, "cancer")
        ) 
    }


####################################################################################################
# Helper functions:
####################################################################################################
download_count_matrix = function(url, dir){
    if(!dir.exists(dir)){
        tmp = tempfile()
        download.file(url, tmp)
        untar(tmp, exdir = dirname(dir))
        file.rename(file.path(dirname(dir), "filtered_feature_bc_matrix"), dir)
        }
    }


summary = function(data, file){
    text = paste0(
            "Sequences: ", ncol(data), "\n",
            "Sites: ", nrow(data), "\n",
            "Unique patterns: ", nrow(unique.matrix(data, MARGIN=1)), "\n",
            "Data density: ", mdensity(data, empty="-")
            )
    writeLines(text, con=file)
    }


####################################################################################################
# Replicating Rao 2020 with Seurat:
####################################################################################################
replicate_rao = function(primary_dir, metastasis_dir, umap_file, cell_file){
    primary = Seurat::Read10X(primary_dir, strip.suffix=TRUE)
    primary = Seurat::CreateSeuratObject(counts=primary, project="primary")
    primary = RenameCells(primary, new.names = paste0(colnames(primary), "-pr"))

    metastasis =  Seurat::Read10X(metastasis_dir, strip.suffix=TRUE)
    metastasis = Seurat::CreateSeuratObject(counts=metastasis, project="metastasis")
    metastasis = RenameCells(metastasis, new.names = paste0(colnames(metastasis), "-met"))

    # Prepare data for integration: keep only cells with more than 500 genes
    primary = subset(primary, subset = nFeature_RNA > 500)
    primary = NormalizeData(primary, verbose = FALSE)
    primary = FindVariableFeatures(primary, selection.method = "vst", nfeatures = 2000)
    primary$status = "Primary"

    metastasis = subset(metastasis, subset = nFeature_RNA > 500)
    metastasis = NormalizeData(metastasis, verbose = FALSE)
    metastasis = FindVariableFeatures(metastasis, selection.method = "vst", nfeatures = 2000)
    metastasis$status = "Metastatic"

    # Find anchors (shared cells types) and integrate the data along the anchors
    anchors = FindIntegrationAnchors(list(primary, metastasis), dims=1:20)
    combined = IntegrateData(anchors=anchors, dims=1:20)
    DefaultAssay(combined) = "integrated"

    # Run the standard workflow for visualization and clustering
    combined = ScaleData(combined, verbose = FALSE)
    combined = RunPCA(combined, npcs = 30, verbose = FALSE)

    # t-SNE and Clustering
    combined = RunUMAP(combined, reduction = "pca", dims = 1:20)
    combined = FindNeighbors(combined, reduction = "pca", dims = 1:20)
    combined = FindClusters(combined, resolution = 0.5)

    # Plots
    pdf(umap_file)
    DimPlot(object = combined, reduction = "umap", split.by = "status", label=TRUE, label.size=5)
    invisible(dev.off())

    # Write down cell lists
    cells = data.frame("cells" = names(combined$seurat_clusters), "clusters" = combined$seurat_clusters)
    write.table(cells, cell_file, col.names=TRUE, row.names=FALSE)
    }


if(!interactive())
    main()
