#!/usr/bin/Rscript

# Based on https://vaulot.github.io/tutorials/Phyloseq_tutorial.html
#
# Basic analysis w/ phyloseq

set.seed(42)

## LOG FILE
sink(file=file(snakemake@log[[1]], open="wt"), type=c("output","message"))

## IMPORT
require(testit) # assertions
require(scales) # scaling
require(ggplot2) # plots
require(ggrepel) # add labels to plots
require(ggsci) # colors
require(viridis) # colors
require(phyloseq) # analysis

## UTILS
source(snakemake@params$utils)

## INFO
print(sessionInfo())
print(snakemake@input)
print(snakemake@output)
print(snakemake@params)

## PARAMS
ORD_METHOD <- snakemake@params$ord_method
ORD_DIST   <- snakemake@params$ord_dist

## DATA
asv_counts <- read_asv_counts(snakemake@input$counts)
asv_tax    <- read_asv_tax(snakemake@input$tax)
asv_tree <- read_tree(snakemake@input$tree) 
# sample metadata
meta       <- read_meta(snakemake@input$meta)
# meta$ConcLabel <- sapply(meta[,"ng/microl"], function(x){ ifelse(x=="< 1",x,">= 1") })
meta$Ctrl      <- sapply(rownames(meta), function(x){ ifelse(grepl("^VIP", x), "VIP", "VPAC") })
print(head(meta, 30))

# phyloseq object
asv_ph <- phyloseq::phyloseq(
    phyloseq::otu_table(as.matrix(asv_counts), taxa_are_rows=FALSE),
    phyloseq::tax_table(as.matrix(asv_tax)),
    phyloseq::sample_data(cbind(sampleID=rownames(meta), meta)),
    phyloseq::phy_tree(asv_tree)
)

# Removing specific sample (MB023-C12) due to very low read counts
asv_ph = subset_samples(asv_ph, sampleID != "MB023-C12")

# count normalization by median sequencing depth
counts_median <- median(phyloseq::sample_sums(asv_ph))
asv_ph_norm   <- phyloseq::transform_sample_counts(
    physeq=asv_ph,
    fun=function(x, t=counts_median){ round(t * (x / sum(x))) }
)

# ordination
asv_ph_norm_ord <- phyloseq::ordinate(physeq=asv_ph_norm, method=ORD_METHOD, distance=ORD_DIST)

## PLOTS
plots <- list()

custom_theme <-
    theme_bw() + # black/white theme
    theme( # custom theme
        text=element_text(size=10),
        axis.ticks=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        strip.background=element_rect(fill="#333333"),
        strip.text=element_text(color="white")
    )

# alpha div
plots$alphadiv <-
    phyloseq::plot_richness(
        physeq=asv_ph_norm,
        measures=c("Chao1", "Shannon", "Simpson", "InvSimpson"),
        nrow=4
    ) +
    labs(
        x=""
    ) +
    custom_theme + theme(
        axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=6)
    )


# barplots: 
# NOTE: skipped - too many samples and ASVs/taxa
carbom_group <- merge_samples(asv_ph_norm, "group")
plots$group_barplot <- plot_bar(carbom_group, fill = "Phylum") +
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack") +
  custom_theme + theme(
    axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=6)
  )

carbom_genotype <- merge_samples(asv_ph_norm, "genotype")
plots$genotype_barplot <- plot_bar(carbom_genotype, fill = "Phylum") +
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack") +
  custom_theme + theme(
    axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=6)
  )


# heatmap
plots$heatmap <-
    phyloseq::plot_heatmap(
        physeq=asv_ph_norm,
        method=ORD_METHOD, distance=ORD_DIST,
        trans=scales::log_trans(2) # not sure if relevant
    ) +
    scale_fill_viridis(
        begin=0, end =0.75,
        discrete=FALSE,
        trans=scales::log_trans(2),
        na.value="white",
        option="magma"
    ) +
    labs(
        title=sprintf("All ASVs (%s, %s)", ORD_DIST, ORD_METHOD),
        x="", y="ASVs",
        fill="Scaled\nabundance"
    ) +
    custom_theme + theme(
        axis.text.y=element_blank(),
        axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=6)
    )

# ordination
plots$ord_asvs <-
    phyloseq::plot_ordination(
        asv_ph_norm,
        asv_ph_norm_ord,
        type="taxa",
        color="Phylum",
        # shape= "Division", 
        title=sprintf("ASVs (all ASVs, %s, %s)", ORD_DIST, ORD_METHOD)
    ) +
    geom_point(size=2) +
    custom_theme + theme(
        axis.text=element_blank()
    )

plots$ord_asvs_split <-
    plots$ord_asvs +
    facet_wrap(~Phylum) +
    theme(
        legend.position="none"
    )

plots$ord_samples_site <-
    phyloseq::plot_ordination(
        asv_ph_norm,
        asv_ph_norm_ord,
        type="samples",
        color="group"
    ) +
    geom_point(size=2) +
    geom_text_repel(
        aes(label=sampleID),
        max.overlaps=10,
        min.segment.length=0,
        size=3,
        seed=42,
        show.legend=FALSE
    ) +
    scale_color_d3() +
    labs(
        title=sprintf("Samples (all ASVs, %s, %s)", ORD_DIST, ORD_METHOD),
        color=""
    ) +
    custom_theme + theme(
        axis.text=element_blank()
    )

plots$ord_samples_site_split <-
    plots$ord_samples_site +
    facet_wrap(~group) +
    theme(
        legend.position="none"
    )

plots$ord_samples_source <-
    phyloseq::plot_ordination(
        asv_ph_norm,
        asv_ph_norm_ord,
        type="samples",
        color="genotype"
    ) +
    geom_point(size=2) +
    geom_text_repel(
        aes(label=sampleID),
        max.overlaps=10,
        min.segment.length=0,
        size=3,
        seed=42,
        show.legend=FALSE
    ) +
    scale_color_d3() +
    labs(
        title=sprintf("Samples (all ASVs, %s, %s)", ORD_DIST, ORD_METHOD),
        color=""
    ) +
    custom_theme + theme(
        axis.text=element_blank()
    )

plots$ord_samples_source_split <-
    plots$ord_samples_source +
    facet_wrap(~genotype) +
    theme(
        legend.position="none"
    )

plot_df <- phyloseq::plot_ordination(
    asv_ph_norm,
    asv_ph_norm_ord,
    type="samples",
    justDF=TRUE
)
plots$ord_samples_source2 <-
    ggplot(data=plot_df, aes(x=DCA1, y=DCA2, color=group, shape=genotype)) +
    geom_point(size=2) +
#    scale_shape_manual(values=META_SHAPE$genotype) +
    scale_color_d3() +
    geom_text_repel(
        data=plot_df[plot_df$group == "VPAC2",],
        aes(label=sampleID),
        max.overlaps=Inf,
        min.segment.length=0,
        force=15,
        size=3,
        seed=42,
        show.legend=FALSE
    ) +
    labs(
        title=sprintf("Samples (all ASVs, %s, %s)", ORD_DIST, ORD_METHOD),
        color=""
    ) +
    custom_theme + theme(
        axis.text=element_blank()
    )

plots$ord_samples_source2_split <- plots$ord_samples_source2 + facet_wrap(~group)
# plots$ord_samples_source2_split <- plots$ord_samples_source2 

# Simple Network
carbom_abund <- filter_taxa(asv_ph_norm, function(x) sum(x > counts_median*0.20) > 0, TRUE)
plots$network <- plot_net(carbom_abund, distance = "(A+B-2*J)/(A+B)", type = "taxa",
         laymeth="circle", maxdist = 0.8, color="Phylum", point_label="Genus")

## PDF
pdf(snakemake@output$pdf, width=16, height=9)
for(pp in plots){
    print(pp)
}
dev.off()

## SAVING DATA
saveRDS(asv_ph_norm, file=snakemake@output$rds)
