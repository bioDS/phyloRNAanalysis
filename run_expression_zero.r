#' run_expression_zero.r
#'
#' Run the expression analysis while assuming that the zero is real zero.
#' However, process the data in the same manner as previously (regarding data filtration).

library("phyloRNA")
library("beter")
library("parallel")

source("src/exprZero.r")
source("src/iqtree.r")
source("src/beast.r")
source("src/utils.r")   

# Get filtered data
expr_filtered = dir("filtered", "expr0", full.names=TRUE)
expr_filtered_subset = dir("filtered", "expr_", full.names=TRUE)

# Convert to fasta, while converting "-" to "0"
expr_fasta = table2fastaZero(expr_filtered, outdir="fastaZero")
expr_fasta_subset = table2fastaZero(expr_filtered_subset, outdir="fastaZero")

# Run phylogenetic analysis
iqtrees_par(
    c(expr_fasta, expr_fasta_subset),
    model = "ORDERED+ASC", nthreads=32,
    outdir = file.path("phylo", "MLZero"),
    )

beasts(
    expr_fasta_subset,
    template = file.path("templates", "BDStrictOrdinal.xml"),
    outdir = file.path("phylo", "BIZero")
    )
