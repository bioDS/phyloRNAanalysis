import::from("src/sra.r", "sra_download", "get_srr_samples")
import::from("src/utils.r", "all.files.exists")
import::from("src/snv.r", "preprocess_snv")
import::from("magrittr", "%>%")

main = function(){
    # Define variables:
    outdir = "wang2021"
    fastqdir = file.path(outdir, "raw")
    snvdir = file.path(outdir, "snv")
    refdir = file.path(outdir, "ref")
    mapdir = file.path(outdir, "map")
    gse = "GSE158631"

    # required reference files:
    reference = "reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
    annotation = "reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gtf"
    vcf = "reference/00-common_all.vcf.gz"
    
    phyloRNA::mkdir(fastqdir)
    phyloRNA::mkdir(snvdir)
    phyloRNA::mkdir(refdir)
    phyloRNA::mkdir(mapdir)
    
    samples = srr_download(gse, fastqdir)
    fastqs = file.path(fastqdir, paste0(samples, ".fastq.gz"))
    fastqs = sapply(fastqs, phyloRNA::abspath)
    manifest = file.path(mapdir, "manifest.txt")

    # Reads are already trimmed
    STAR_genome_index(reference, annotation, 50, refdir)
    STAR_manifest(manifest, fastqs, samples, gse)
    #STAR_map(mapdir, "all", refdir, fastq=fastqs, gzip=TRUE)
    STAR_map(mapdir, "all", refdir, manifest=manifest, gzip=TRUE)

    mapped = file.path(mapdir, "allAligned.sortedByCoord.out.bam")
    tagged = file.path(mapdir, "all.tag.bam")
    cleaned = file.path(mapdir, "all.bam")

    # clean
    rg2cb(mapped, tagged)
    phyloRNA::gatk_prepare(tagged, cleaned, reference=reference, vcf=vcf)
    vcm = preprocess_snv(cleaned, barcodes=NULL, reference, outdir=snvdir)

    gc1 = grep("GC1", samples, value=TRUE)
    gc2 = grep("GC2", samples, value=TRUE)
    gc3 = grep("GC3", samples, value=TRUE)

    gc1_fasta = file.path(snvdir, "gc1.fasta")
    gc2_fasta = file.path(snvdir, "gc2.fasta")
    gc3_fasta = file.path(snvdir, "gc3.fasta")

    # generate fasta
    vcm2fasta(vcm, gc1_fasta, gc1)
    vcm2fasta(vcm, gc2_fasta, gc2)
    vcm2fasta(vcm, gc3_fasta, gc3)
    
    # run IQtree
    iqtree(gc1_fasta, model="TEST", file.path(snvdir, "gc1"))
    iqtree(gc2_fasta, model="TEST", file.path(snvdir, "gc2"))
    iqtree(gc3_fasta, model="TEST", file.path(snvdir, "gc3"))
    }


srr_download = function(gse, outdir){
    # TODO cache this so there is no need to access entrez every time
    samples = get_srr_samples(gse, save=file.path(outdir, "samples.rds"))

    mapply(samples$srr, samples$name, FUN=srr_download_sample, MoreArgs=list("outdir"=outdir))
    }


srr_download_sample = function(srr, name, outdir){
    fastq_files = file.path(outdir, paste0(name, ".fastq.gz"))

    if(all.files.exists(fastq_files))
        return(name)
    
    srr_files = file.path(outdir, paste0(srr, "_", 1, ".fastq.gz"))
    sra_download(srr, outdir)
    if(!all.files.exists(srr_files))
        stop("ERROR: not all files exists.\\n", "Files: ", files)

    file.rename(srr_files, fastq_files)
    
    name
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
        data = data[,..match]
        }

    data = as.matrix(data)
    data = phyloRNA::remove_constant(data)
    seq = tab2seq(data, 2)
    write_fasta(seq, fasta)
    }

STAR_genome_index = function(ref, gtf, len, outdir, nthreads=16){
    statusfile = file.path(outdir, "genome.finished")
    if(file.exists(statusfile))
        return(invisible())
    
    command = "STAR"
    args = c(
        "--runMode genomeGenerate",
        "--genomeDir", outdir,
        "--genomeFastaFiles", ref,
        "--sjdbGTFfile", gtf,
        "--sjdbOverhang", len-1,
        "--runThreadN", nthreads
        )
    phyloRNA:::systemE(command, args)
    file.create(statusfile)
    }


STAR_map = function(
    outdir, prefix, refdir,
    fastq=NULL, manifest=NULL, gzip=FALSE, nthreads=16
    ){
    statusfile = file.path(outdir, paste0(prefix, ".finished"))
    if(file.exists(statusfile))
        return(invisible())

    command = "STAR"
    args = c(
        "--genomeDir", refdir,
        "--outFileNamePrefix", file.path(outdir, prefix),
        "--outSAMattributes NH HI AS nM NM MD jM jI MC ch RG",
        "--outSAMtype BAM SortedByCoordinate"
        )

    if(!is.null(fastq)){
        args = c(args, "--readFilesIn", paste0(fastq, collapse=","))
        args = c(
            args,
            "--outSAMattrRGline",
            paste0("ID:", phyloRNA::corename(fastq), collapse=" , ")
            )
        }

    if(!is.null(manifest))
        args = c(args, "--readFilesManifest", manifest)
        
    if(gzip)
        args = c(args, "--readFilesCommand gunzip -c")

    phyloRNA:::systemE(command, args)
    file.create(statusfile)
    }

STAR_manifest = function(file, fastq, id, gse){
    if(file.exists(file))
        return(invisible())

    # RG tags required by GATK are:
    # * ID identifies read and is -> copied to read
    # * SM identifies sample
    # * LB identifies library, can have multiple libraries per sample
    # * PU consist of flowcell barcode, lane and sample barcode
    # * PL identifies platform, Illumina in this case
    RG = paste(
        paste0("ID:", id),
        paste0("SM:", id),
        paste0("LB:", id),
        paste0("PU:", paste(gse, "1", id, sep=":")),
        paste0("PL:", "ILLUMINA"),
        sep="\t"
        )

    manifest = cbind(fastq, "-", RG)
    write.table(manifest, file, sep="\t", quote=FALSE,
                col.names=FALSE, row.names=FALSE)
    }


rg2cb = function(input, output){
    if(file.exists(output))
        return(invisible())

    command = "python3"
    args = c(
        "src/rg2cb.py",
        input, output
        )
    phyloRNA:::systemE(command, args)
    }

if(!interactive()){
    main()
    }
