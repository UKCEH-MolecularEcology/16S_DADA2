#!/usr/bin/Rscript

## NOTE
# Based on:
#   https://astrobiomike.github.io/amplicon/dada2_workflow_ex#processing-with-dada2-in-r
#   https://benjjneb.github.io/dada2/tutorial.html
#   https://benjjneb.github.io/dada2/bigdata.html
#
# Process on sample: error model, dereplication, define ASVs

set.seed(42)

## LOG FILE
sink(file=file(snakemake@log[[1]], open="wt"), type=c("output","message"))

## IMPORT
suppressMessages(library(testit)) # checks, assertions
suppressMessages(library(dada2)) # DADA2 analysis
suppressMessages(library(ggplot2)) # plots

## UTILS
source(snakemake@params$utils)

## ARGS
for(arg_name in c('minOverlap')){
    testit::assert(sprintf("Argument %s must be > 0", arg_name), all(snakemake@config$dada2$asv[[arg_name]] > 0))
}
for(fname in snakemake@input){
    testit::assert(sprintf("Could not find file %s", fname), file.exists(fname))
}
THREADS <- ifelse(snakemake@threads==1, FALSE, snakemake@threads)

## INFO
print(sessionInfo())
print(snakemake@input)
print(snakemake@output)
print(snakemake@params)
print(snakemake@config$dada2$asv)

## ERROR MODEL
errormod <- readRDS(snakemake@input$errormod)

## DEREPLICATION
print("Dereplication")
derep_forward <- dada2::derepFastq(snakemake@input$r1, verbose=TRUE)
derep_reverse <- dada2::derepFastq(snakemake@input$r2, verbose=TRUE)

## ASVs
print("Defining ASVs: R1")
asv_forward <- dada2::dada(
    derep=derep_forward,
    err=errormod$forward,
    pool=snakemake@config$dada2$asv$pool,
    multithread=THREADS,
    verbose=TRUE
)
print("Defining ASVs: R2")
asv_reverse <- dada2::dada(
    derep=derep_reverse,
    err=errormod$reverse,
    pool=snakemake@config$dada2$asv$pool,
    multithread=THREADS,
    verbose=TRUE
)

# MEREGING (paired reads)
print("Merging reads")
asv_merged <- dada2::mergePairs(
    dadaF=asv_forward,
    derepF=derep_forward,
    dadaR=asv_reverse,
    derepR=derep_reverse,
    trimOverhang=TRUE,
    minOverlap=snakemake@config$dada2$asv$minOverlap,
    verbose=TRUE
)
saveRDS(asv_merged, file=snakemake@output$rds)

print(sprintf(
    "Number of reads:\nR1: %d\nR2: %d\nMerged: %d",
    dada2_count_reads(asv_forward),
    dada2_count_reads(asv_reverse),
    dada2_count_reads(asv_merged)
))