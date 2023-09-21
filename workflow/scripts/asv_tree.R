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
# suppressMessages(library(testit)) # checks, assertions
suppressMessages(library(dada2)) # DADA2 analysis
suppressMessages(library(DECIPHER)) # alignment
suppressMessages(library(phangorn)) # tree building

## UTILS
source(snakemake@params$utils)

## ARGS
#for(arg_name in c('minOverlap')){
#    testit::assert(sprintf("Argument %s must be > 0", arg_name), all(snakemake@config$dada2$asv[[arg_name]] > 0))
#}
#for(fname in snakemake@input){
#    testit::assert(sprintf("Could not find file %s", fname), file.exists(fname))
#}
#THREADS <- ifelse(snakemake@threads==1, FALSE, snakemake@threads)

## INFO
print(sessionInfo())
print(snakemake@input)
print(snakemake@output)
print(snakemake@params)
print(snakemake@config$dada2$asv)

## COUNTS
print("Constructing phylogenetic tree")
asv_merged <- readRDS(snakemake@input$rds)
counts <- dada2::makeSequenceTable(asv_merged)

## CHIMERAS
print("Removing chimeras")
counts_nochim <- dada2::removeBimeraDenovo(counts, method="consensus", multithread=TRUE, verbose=TRUE)
print(sprintf("ASV count w/o chimeras: %.3f", 100 * sum(counts_nochim)/sum(counts)))

## ALIGNMENT
print("Generating a multiple sequence alignment")
seqtab <- counts_nochim

# TREE NAMES
seqtab_cols <- paste0("ASV_", seq(1:ncol(seqtab)))

# ASSIGNING ASV NAMES TO SEQUENCES
seqs <- getSequences(seqtab)
names(seqs) <- seqtab_cols # This propagates to the tip labels of the tree
# Assuming 'seqs' is a character vector
seqs_df <- as.data.frame(seqs) # seqs_df <- data.frame(ASV=rownames(seqs), Sequence=seqs)   # writing the seqs with ASV name to df
seqs_df
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA)


## TREE
print("Building the tree")
# making the tree
phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)

## negative edges length changed to 0!
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                     rearrangement = "stochastic", control = pml.control(trace = 0))

## SAVE DATA
print("Saving data")
# Sequence file
write.table(seqs_df, file=snakemake@output@seq, row.names=TRUE)
# Tree file
write.tree(fitGTR$tree, file=snakemake@output$tree)
# RDS
saveRDS(fitGTR, file=snakemake@output$rds)
