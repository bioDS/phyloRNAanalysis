#!/usr/bin/env Rscript

#' 0-run.r
#'
#' Runs the whole analysis by executing individual scripts:
#' - 1-prepare.r -- prepares bam files
#' - 2-snv.r -- SNV identification, filtering and a phylogenetic analysis on the SNV data
#' - 3-expression.r -- Expression analysis, cleaning and a phylogenetic analysis on the Expression data

source("1-prepare.r")
#source("2-snv.r")
#source("3-expression.r")
