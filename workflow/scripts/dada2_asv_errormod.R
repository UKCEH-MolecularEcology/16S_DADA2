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
print("Error models")

# model
err_forward <- dada2::learnErrors(
    fls=snakemake@input$r1,
    randomize=TRUE,
    multithread=THREADS,
    verbose=TRUE
)
err_reverse <- dada2::learnErrors(
    fls=snakemake@input$r2,
    randomize=TRUE,
    multithread=THREADS,
    verbose=TRUE
)

# plot
custom_theme <-
    theme_bw() + # black/white theme
    theme( # custom theme
        text=element_text(size=9),
        strip.background=element_rect(fill="white"),
        strip.text=element_text(color="black")
    )
pdf(snakemake@output$pdf, width=14, height=8)
dada2::plotErrors(err_forward, nominalQ=TRUE) + custom_theme
dada2::plotErrors(err_reverse, nominalQ=TRUE) + custom_theme
dev.off()

# save
saveRDS(list(forward=err_forward, reverse=err_reverse), file=snakemake@output$rds)