library("phyloRNA")
library("parallel")
import::from("src/utils.r", "mdensity")
import::from("src/expr.r", "calculate_intervals")
import::from("src/iqtree.r", "iqtrees_par")
import::from("src/beast.r", "beasts")

# Density intervals
hdi = c(0.6, 0.9)
outdir = file.path("wang2021")
phyloRNA::mkdir(outdir)

# url from which to download count matrix
url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE158631&format=file&file=GSE158631%5Fcount%2Ecsv%2Egz"
infile = file.path(outdir, "count_matrix.csv.gz")

if(!file.exists(infile)){
    download.file(url, infile)
    }

# Read the data and split it according to the patient
count_matrix = read.table(infile, sep=",", header=TRUE, row.names=1)
count_matrix = split.default(count_matrix, sub("\\..*", "", names(count_matrix)))

# Reconstruct tree for every patient
for(patient in names(count_matrix)){
    data = count_matrix[[patient]]

    data = phyloRNA::expr_zero_to_na(data)
    data = phyloRNA::expr_scale(data)
    intervals = calculate_intervals(data, density=hdi)
    data = phyloRNA::expr_discretize(data, intervals=intervals, unknown="-")
    data = phyloRNA::remove_constant(data, margin=1, unknown="-")

    # Save the data into a fasta file:
    patient_fasta = file.path(outdir, paste0(patient, ".fasta"))
    phyloRNA::fasta(data, margin=2, file=patient_fasta)

    # Summarize data:
    text = paste0(
        "Patient: ", patient, "\n",
        "Sequences: ", ncol(data), "\n",
        "Sites: ", nrow(data), "\n",
        "Unique patterns: ", nrow(unique.matrix(data, MARGIN=1)), "\n",
        "Data density: ", mdensity(data, empty="-")
        )
    patient_summary = file.path(outdir, paste0(patient, "_summary.txt"))
    writeLines(text, con=patient_summary)

    # Run phylogenetic analysis
    iqtrees_par(
        patient_fasta,
        model = "ORDERED+ASC", nthreads=8,
        outdir = file.path(outdir, patient, "ML"),
        )

    beasts(
        patient_fasta,
        template = file.path("templates", "BDStrictOrdinal.xml"),
        outdir = file.path(outdir, patient, "BI")
        )
    }
