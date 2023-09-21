#!/usr/bin/Rscript

## NOTE
# Based on:
#   https://astrobiomike.github.io/amplicon/dada2_workflow_ex#processing-with-dada2-in-r
#   https://benjjneb.github.io/dada2/tutorial.html
#   https://benjjneb.github.io/dada2/bigdata.html
#
# Merge and process multiple samples (after running dada2_asv.R): XXX

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

## COUNTS
print("Count table")
asv_merged <- list()
for(fname in snakemake@input$asv){
    sname <- sub("\\.rds$", "", basename(fname))
    asv_merged[[sname]] <- readRDS(fname)
    print(sprintf("Read in ASVs for %s from %s", sname, fname))
}
counts <- dada2::makeSequenceTable(asv_merged)

## CHIMERAS
print("Removing chimeras")
counts_nochim <- dada2::removeBimeraDenovo(counts, verbose=TRUE)
print(sprintf("ASV count w/o chimeras: %.3f", 100 * sum(counts_nochim)/sum(counts)))

## TAXONOMY
print("Assigning taxonomy")
taxa <- dada2::assignTaxonomy(
    seqs=counts_nochim,
    refFasta=snakemake@input$tax,
    tryRC=TRUE,
    multithread=THREADS,
    verbose=TRUE
)
testit::assert( all( rownames(taxa) == colnames(counts_nochim) ) )

## STATS
print("Computing stats")
stats <- data.frame(
    SampleID=rownames(counts),
    Reads2ASVmerged=sapply(asv_merged, dada2_count_reads),
    Reads2ASVnochim=rowSums(counts_nochim),
    ASV=apply(counts, 1, function(x){ sum(x>0) }),
    ASVnochim=apply(counts_nochim, 1, function(x){ sum(x>0) }),
    check.names=FALSE, stringsAsFactors=FALSE
)

## SAVE DATA
print("Saving data")
# Preproc.
asv_seqs <- colnames(counts_nochim)
asv_ids  <- paste("ASV", 1:length(asv_seqs), sep="_")
rownames(taxa) <- asv_ids
counts_nochim2tab <- counts_nochim
colnames(counts_nochim2tab) <- asv_ids
# Counts
write.table(counts_nochim2tab, file=snakemake@output$counts, sep='\t', row.names=TRUE, col.names=TRUE, quote=FALSE)
# Taxonomy
write.table(taxa, file=snakemake@output$tax, sep='\t', row.names=TRUE, col.names=TRUE, quote=FALSE)
# FASTA
write(c(rbind(paste(">", asv_ids, sep=""), asv_seqs)), snakemake@output$fasta)
# Stats
write.table(stats, file=snakemake@output$stats, sep='\t', row.names=FALSE, col.names=TRUE, quote=FALSE)
# RDS
saveRDS(asv_merged, file=snakemake@output$rds)
