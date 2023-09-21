##################################################
# Files

FASTQ_SUFFIX  <- "_R[1,2]\\.fastq\\.gz$"

#' Get sample ID from FASTQ file path/name
#' @param fname: File path/name
#' @param suffix: File suffix to be removed to get WGS-ID
#' @return WGS-ID extracted from file basename
sid_from_fastq <- function(fname, suffix=FASTQ_SUFFIX){
    fname <- basename(fname)
    fname <- sub(suffix, "", fname)
    return(fname)
}

#' Read sample metadata
#' @param fname: File path
#' @return data.frame
read_meta <- function(fname){
    print(sprintf("Reading in: %s", fname))
    tab <- read.csv(file=fname, sep='\t', header=TRUE, row.names=1, stringsAsFactors=FALSE, check.names=FALSE)
    return(tab)
}

#' Read in ASV count table
#' Requirements: samples x features; header, rownames, tab-sep.
#' @param fname: File path
#' @return data.frame
read_asv_counts <- function(fname){
    print(sprintf("Reading in: %s", fname))
    tab <- read.csv(file=fname, sep='\t', header=TRUE, row.names=1, stringsAsFactors=FALSE, check.names=FALSE)
    return(tab)
}

#' Read in ASV taxonomy table
#' Requirements: features x taxonomy ranks; header, rownames, tab-sep.
#' @param fname: File path
#' @return data.frame
read_asv_tax <- function(fname){
    print(sprintf("Reading in: %s", fname))
    tab <- read.csv(file=fname, sep='\t', header=TRUE, row.names=1, stringsAsFactors=FALSE, check.names=FALSE)
    return(tab)
}

#' Read in ASV stat.s table
#' @param fname: File path
#' @return data.frame
read_asv_stats <- function(fname){
    print(sprintf("Reading in: %s", fname))
    tab <- read.csv(file=fname, sep='\t', header=TRUE, row.names=1, stringsAsFactors=FALSE, check.names=FALSE)
    return(tab)
}

#' Read in MultiQC FastQC data
#' @param fname: File path
#' @param rm_r2: rm entries for *_R2 keep only *_R1, sample suffix (_R1|R2) will be removed
#' @return data.frame
read_multiqc_fastq <- function(fname, rm_r2=TRUE){
    print(sprintf("Reading in: %s", fname))
    tab <- read.csv(file=fname, sep='\t', header=TRUE, row.names=1, stringsAsFactors=FALSE, check.names=FALSE)
    if(rm_r2){
        tab <- tab[grepl("_R2$", rownames(tab)),]
        rownames(tab) <- sub("_R1$", "", rownames(tab))
    }
    return(tab)
}

FASTQ_ASV_STATS <- c(
    raw="Reads\n(raw)",
    trim1="Reads\n(trimmed I)",
    trim2="Reads\n(trimmed II)",
    Reads2ASVmerged="Reads\n(DADA2)",
    Reads2ASVnochim="Reads\n(DADA2, bimera-free)",
    ASV="ASVs",
    ASVnochim="ASVs\n(bimera-free)"
)

##################################################
# DADA2
#' Count number of reads
#' @param x: DADA2 ASV sample object from dada2::dada or dada2::mergePairs
#' @return number of reads in given sample object
dada2_count_reads <- function(x){
    sum(dada2::getUniques(x))
}

##################################################
# Plots
META_SHAPE <- list(
    Ctrl=c(Control=15, Sample=16),
    ConcLabel=c("< 1"=15, ">= 1"=16)
)
META_SIZE <- list(
    ConcLabel=c("< 1"=1, ">= 1"=3)
)