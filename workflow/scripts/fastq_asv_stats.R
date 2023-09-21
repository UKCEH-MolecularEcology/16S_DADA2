#!/usr/bin/Rscript

# Collect and plot FastQC/ASV stat.s

set.seed(42)

## LOG FILE
sink(file=file(snakemake@log[[1]], open="wt"), type=c("output","message"))

## IMPORT
suppressMessages(library(testit)) # checks, assertions
suppressMessages(library(dplyr)) # table manipulation
suppressMessages(library(reshape2)) # reshape data
suppressMessages(library(ggplot2)) # plots
suppressMessages(library(ggsci)) # colors

## UTILS
source(snakemake@params$utils)

## INFO
print(sessionInfo())
print(snakemake@input)
print(snakemake@output)
print(snakemake@params)

## DATA
tab <- read_asv_stats(snakemake@input$asv)
tab$raw   <- read_multiqc_fastq(snakemake@input$raw)[rownames(tab),"Total Sequences"]
tab$trim1 <- read_multiqc_fastq(snakemake@input$trim1)[rownames(tab),"Total Sequences"]
tab$trim2 <- read_multiqc_fastq(snakemake@input$trim2)[rownames(tab),"Total Sequences"]

## PLOTS
plots <- list()

custom_theme <-
    theme_bw() + # black/white theme
    theme(legend.position = "none") + # remove legend
    theme( # custom theme
        text=element_text(size=10),
        axis.text=element_text(color='black'),
        axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=6),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks=element_blank(),
        strip.background=element_rect(fill="#333333"),
        strip.text=element_text(color="white")
    )

# barplot
plot_df <- reshape2::melt(cbind(sample=rownames(tab), tab), id.vars=c("sample"))
plot_df$variable <- as.character(plot_df$variable)
plot_df$group <- sapply(plot_df$variable, function(x){ ifelse(grepl("^ASV", x), "ASVs", "Reads") })
plot_df$threshold <- sapply(plot_df$variable, function(x){ ifelse(grepl("^ASV", x), NA, 15000) })
plot_df$variable <- factor(
    x=FASTQ_ASV_STATS[plot_df$variable],
    levels=FASTQ_ASV_STATS,
    ordered=TRUE
)

plots$barplot <-
    ggplot(data=plot_df, aes(x=sample, y=value, fill=group)) +
    geom_col() +
    geom_hline(aes(yintercept=threshold), linetype="dotted") +
    scale_fill_manual(values=setNames(pal_jama("default")(7)[c(1,5)], c("ASVs", "Reads"))) +
    facet_grid(rows=vars(variable), scales="free_y") +
    labs(
        x="", y="Number of reads/ASVs"
    ) +
    custom_theme

## PDF
pdf(snakemake@output$pdf, width=16, height=9)
for(pp in plots){
    print(pp)
}
dev.off()

## TSV
write.table(cbind(Sample=rownames(tab), tab), file=snakemake@output$tsv, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)