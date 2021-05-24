import::from("src/sra.r", "sra_download", "get_srr_samples")
import::from("src/utils.r", "all.files.exists")
import::from("src/snv.r", "preprocess_snv")
import::from("src/prepare.r", "merge_files")
import::from("magrittr", "%>%")
import::from("phyloRNA", "tab2seq", "write_fasta")

main = function(){
    # Define variables:
    outdir = "rao2020"
    fastqdir = file.path(outdir, "raw")
    refdir = "prepare/ref"
    countdir = file.path(outdir, "count")
    snvdir = file.path(outdir, "snv")
    gse = "GSE140312"
    samples = srr_download(gse, fastqdir)

    # required reference files:
    reference = "reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
    annotation = "referencey/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gtf"
    vcf = "reference/00-common_all.vcf.gz"

    # create directories
    phyloRNA::mkdir(fastqdir)
    phyloRNA::mkdir(refdir)
    
    
    for(sample in samples){
        # map with cellranger
        phyloRNA::cellranger_count(
            fastqdir, refdir, file.path(countdir, sample),
            sample=sample, nthreads=8
            )

        # prepare BAM and barcodes
        bam_aligned = file.path(countdir, sample, "outs", "possorted_genome_bam.bam")
        bam_cleaned = file.path(countdir, paste0(sample, ".cleaned.bam"))
        bam_prepared = file.path(countdir, paste0(sample, ".prepared.bam"))
        
        barcodes_aligned = file.path(countdir, sample, "outs", "filtered_feature_bc_matrix", "barcodes.tsv.gz")
        barcodes_prepared = file.path(countdir, paste0(sample, ".prepared.txt"))
        
        pattern = "-1$"
        replace = paste0("-", sample)

        phyloRNA::gatk_prepare(
            bam_aligned,
            bam_cleaned,
            reference=reference,
            vcf=vcf,
            barcodes=barcodes_aligned
            )

        phyloRNA::bamtagregex(bam_cleaned, bam_prepared, "CB", pattern, replace)
        filesub(barcodes_aligned, barcodes_prepared, pattern, replace)
        }

    # Merge BAM and barcodes
    bam_prepared = file.path(countdir, paste0(samples, ".prepared.bam"))
    barcodes_prepared = file.path(countdir, paste0(samples, ".prepared.txt"))
    bam = file.path(countdir, "all.bam")
    barcodes = file.path(countdir, "all.txt")

    phyloRNA::gatk_MergeSamFiles(bam_prepared, bam)
    merge_files(barcodes_prepared, barcodes)

    # preprocess SNVs
    vcm = preprocess_snv(bam, barcodes, reference, outdir=snvdir)

    # Create two datasets based on selected cells from expression analysis
    expr_cancer_fasta = file.path(outdir, "expr", "rao2020_cancer.fasta")
    expr_all_fasta = file.path(outdir, "expr", "rao2020_all.fasta")

    cancer_names = get_names(expr_cancer_fasta)
    all_names = get_names(expr_all_fasta)

    phyloRNA::mkdir(file.path(outdir, "snv"))
    cancer_fasta = file.path(outdir, "snv", "cancer.fasta")
    all_fasta = file.path(outdir, "snv", "all.fasta")

    vcm2fasta(vcm, cancer_fasta, cancer_names)
    vcm2fasta(vcm, all_fasta, all_names)

    # iqtree
    iqtree(
        cancer_fasta,
        model="ASC",
        file.path(outdir, "snv", "cancer")
        )

    iqtree(
        all_fasta,
        model="ASC",
        file.path(outdir, "snv", "all")
        )
    }


get_names = function(fasta){
    seq = phyloRNA::read_fasta(fasta)
    names = names(seq)
    names = sub("pr$", "primary", names)
    names = sub("met$", "metastatic", names)
    names
    }


vcm2fasta = function(vcm, fasta, selection=NULL){
    data = data.table::fread(vcm, header=TRUE)
    data = data[,-c(1:3)] # first three columns are not cells (chromosome, position, and reference)
    
    if(!is.null(selection)){
        match = selection %in% colnames(data)
        if(!all(match)){
            warning("WARNING: ", sum(!match), " out of ", length(match),
                " requested cells are not present:\n",
                paste0(selection[!match], collapse="\n")
                )
            }
        data = data[, ..match] # data.table subsetting
        }

    data = phyloRNA::remove_constant(data)
    seq = tab2seq(data, 2)
    write_fasta(seq, fasta)
    }


filesub = function(input, output, pattern, replace){
    text = readLines(input)
    text = sub(pattern, replace, text)
    writeLines(text, output)
    }



srr_download = function(gse, outdir){
    # TODO cache this so there is no need to access entrez every time
    samples = get_srr_samples(gse)
    samples$name = gsub(" ", "_", samples$name)

    mapply(samples$srr, samples$name, FUN=srr_download_sample, MoreArgs=list("outdir"=outdir))
    }


srr_download_sample = function(srr, name, outdir){
    # For Cellranger, file names need to have this structure:
    # [sample_name]_S1_L001_[read type]_001.fastq.gz
    # S1 -- sample 1
    # L001 -- lane 001
    # read type:
    #     -- I1 -- Index file
    #     -- R1 -- barcodes or actuall reads
    #     -- R2 -- barcodes or actual reads
    read_types = c("I1", "R1", "R2")
    sample_name = strsplit(name, split="_", fixed=TRUE) %>% unlist %>% getElement(1)

    cellranger_files = file.path(
        outdir,
        paste0(sample_name, "_S1_L001_", read_types, "_001.fastq.gz")
        )

    if(all.files.exists(cellranger_files))
        return(sample_name)
    
    # Download srr files and check if they exist/were downloaded correctly
    srr_files = file.path(outdir, paste0(srr, "_", 1:3, ".fastq.gz"))
    sra_download(srr, outdir)
    if(!all.files.exists(srr_files))
        stop("ERROR: not all files exists.\\n", "Files: ", files)
    
    file.rename(srr_files, cellranger_files)
    
    sample_name
    }


if(!interactive()){
    main()
    }
