#!/usr/bin/env Rscript

#' 1-prepare.r
#'
#' Prepare data for analysis
#' - remap to a new reference
#' - clean and filter to only include filtered barcodes
#' - modify barcodes with a unique tag
#' - merge all samples into a single file  
library("phyloRNA")


# Names of dataset that will be processed
datasets = c("2HR", "HLR", "HLRblood", "NH", "NHblood")

# Required inputs for analysis
reference = "/data/phylonco/ReferenceGenomes/human_GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
annotation = "/data/phylonco/ReferenceGenomes/human_GRCh38/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gtf"
vcf = "/data/phylonco/ReferenceGenomes/vcf/00-common_all.vcf.gz"

# shared setting:
outdir = "prepare"
output_bam = file.path(outdir, "all.bam")
output_bar = file.path(outdir, "all.txt")
nthreads = 12


if(!file.exists(output_bam) && !file.exists(output_bar)){
    # process all datasets
    for(dataset in dataset){
        process_dataset(dataset)
        }

    # merge bam form all dataset
    bam_prepared = file.path(outdir, datasets, paste0(datasets, ".prepared.bam"))
    gatk_MergeSamFiles(inputs, output_bam)
    
    # merge barcode files from all datasets
    if(file.exists(output_bar)) file.remove(output_bar)
    barcode_prepared = file.path(outdir, dataset, paste0(datasets, ".barcode.prepared.txt"))
    file.append(output_bar, barcode_prepared)
    }

    

process_dataset = function(dataset){
    # file specific to dataset:
    input = file.path("data", paste0(dataset, ".bam"))
    barcode_aligned = file.path(outdir, dataset, paste0(dataset, ".barcodes.txt"))
    barcode_prepared = file.path(outdir, dataset, paste0(dataset, ".barcodes.prepared.txt"))
    bam_aligned = file.path(outdir, dataset, paste0(dataset, ".aligned.bam"))
    bam_prepared = file.path(outdir, dataset, paste0(dataset, ".prepared.bam"))

    gatkdir = file.path(otudir, dataset, "gatk")
    mkdir(gatkdir)

    # remap the bam to a new reference
    remap(
        input, reference, annotation,
        outdir = file.path(outdir, dataset),
        refdir = file.path(outdir, "reference"),
        nthreads = nthreads,
        cbam = bam_aligned,
        cbar = barcode_aligned
        #cfbm=TRUE -- TODO: Copy barcode matrix? Do I need the whole folder or just .h5 file?
        )

    # create new GATKR6 object
    bamGATK = GATKR6$new(
        bam = bamfile,
        reference = reference,
        vcf = vcf
        )

    # clean bam, put output into a new folder:
    bamGATK$SortSam(output = file.path(gatkdir, paste0(basename(bam_aligned), ".sorted")))
    bamGATK$SplitNCigarReads()$Recalibrate()

    # filter bam to include only reads from filtered barcodes
    barcodes = readLines(barcode_aligned)
    bamGATK$FilterSamReadsTag("CB", barcodes)

    # in bam change for every barcode "-1" into "-<dataset>" 
    bamtagregex(
        input = bamGATK$bam,
        output = bam_prepared,
        tag = "CB",
        pattern = "-1",
        replace = paste0("-", dataset, "$")
        )

    # in barcode file, change "-1" into "-<dataset>"
    barcodes = sub("-1", paste0("-", dataset, "$"), barcodes)
    writeLines(barcodes, barcode_prepared)
    }
