#!/usr/bin/Rscript

## NOTE
# Based on the following tutorial: https://astrobiomike.github.io/amplicon/dada2_workflow_ex#processing-with-dada2-in-r

## LOG FILE
sink(file=file(snakemake@log[[1]], open="wt"), type=c("output", "message"))

## IMPORT
suppressMessages(library(testit)) # checks, assertions
suppressMessages(library(dada2)) # DADA2 analysis
suppressMessages(library(ggplot2)) # plots

## UTILS
source(snakemake@params$utils)

## ARGS
testit::assert("Argument truncLen must be > 0", all(snakemake@config$dada2$fastq$truncLen >= 0))
testit::assert("Argument minLen must be > 0",   all(snakemake@config$dada2$fastq$minLen   >  0))
for(arg_name in c('truncLen')){
    testit::assert(sprintf("Argument %s must have 2 values", arg_name), length(snakemake@config$dada2$fastq[[arg_name]]) == 2)
}
for(arg_name in c("r1", "r2")){
    testit::assert(sprintf("Could not find file %s", arg_name), file.exists(snakemake@input[[arg_name]]))
}

## INFO
print(sessionInfo())
print(snakemake@input)
print(snakemake@output)
print(snakemake@params)
print(snakemake@config$dada2$fastq)

## FILTER & TRIM
processed <- dada2::filterAndTrim(
    fwd=snakemake@input$r1,
    filt=snakemake@output$r1,
    rev=snakemake@input$r2,
    filt.rev=snakemake@output$r2,
    compress=TRUE,
    truncLen=snakemake@config$dada2$fastq$truncLen,
    minLen=snakemake@config$dada2$fastq$minLen,
    rm.phix=TRUE,
    matchIDs=TRUE,
    multithread=ifelse(snakemake@threads==1, FALSE, snakemake@threads),
    qualityType="Auto",
    verbose=TRUE
)

## PLOTS
custom_theme <-
    theme_bw() + # black/white theme
    theme( # custom theme
        text=element_text(size=9),
        strip.background=element_rect(fill="white"),
        strip.text=element_text(color="black")
    )
pdf(snakemake@output$pdf, width=12, height=6)
dada2::plotQualityProfile(c(snakemake@input$r1,  snakemake@input$r2 )) + custom_theme
dada2::plotQualityProfile(c(snakemake@output$r1, snakemake@output$r2)) + custom_theme
dev.off()