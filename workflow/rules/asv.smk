##################################################
# ASVs

rule dada2_download_taxa:
    output:
        os.path.join(DBS_DIR, "dada2", os.path.basename(config["dada2"]["asv"]["tax"]))
    params:
        url=config["dada2"]["asv"]["tax"]
    log:
        os.path.join(DBS_DIR, "dada2", os.path.basename(config["dada2"]["asv"]["tax"]) + ".log")
    message:
        "DADA2: download taxa data"
    shell:
        "(date && wget -O {output} {params.url} && date) &> {log}"

rule dada2_asv_errormod:
    input:
        r1=expand(os.path.join(RESULTS_DIR, "fastq/cutadapt_adapters/{sid}_R1.fastq.gz"), sid=SAMPLES),
        r2=expand(os.path.join(RESULTS_DIR, "fastq/cutadapt_adapters/{sid}_R2.fastq.gz"), sid=SAMPLES)
    output:
        pdf=os.path.join(RESULTS_DIR, "dada2/errormod.pdf"),
        rds=os.path.join(RESULTS_DIR, "dada2/errormod.rds")
    log:
        os.path.join(RESULTS_DIR, "dada2/errormod.log")
    threads:
        config["dada2"]["asv"]["threads"]
    params:
        utils=os.path.join(SRC_DIR, "utils.R")
    conda:
        os.path.join(ENV_DIR, "r.yaml")
    message:
        "DADA2: error model"
    script:
        os.path.join(SRC_DIR, "dada2_asv_errormod.R")

rule dada2_asv:
    input:
        r1=os.path.join(RESULTS_DIR, "fastq/cutadapt_adapters/{sid}_R1.fastq.gz"),
        r2=os.path.join(RESULTS_DIR, "fastq/cutadapt_adapters/{sid}_R2.fastq.gz"),
        errormod=os.path.join(RESULTS_DIR, "dada2/errormod.rds")
    output:
        rds=temp(os.path.join(RESULTS_DIR, "dada2/{sid}.rds"))
    log:
        os.path.join(RESULTS_DIR, "dada2/{sid}.log")
    threads:
        config["dada2"]["asv"]["threads_single"]
    params:
        utils=os.path.join(SRC_DIR, "utils.R")
    conda:
        os.path.join(ENV_DIR, "r.yaml")
    message:
        "DADA2: ASVs: {input}"
    script:
        os.path.join(SRC_DIR, "dada2_asv.R")

rule dada2_merge:
    input:
        asv=expand(os.path.join(RESULTS_DIR, "dada2/{sid}.rds"), sid=SAMPLES),
        tax=os.path.join(DBS_DIR, "dada2", os.path.basename(config["dada2"]["asv"]["tax"]))
    output:
        counts=os.path.join(RESULTS_DIR, "dada2/ASV.counts.tsv"),
        fasta=os.path.join(RESULTS_DIR, "dada2/ASV.fasta"),
        tax=os.path.join(RESULTS_DIR, "dada2/ASV.tax.tsv"),
        stats=os.path.join(RESULTS_DIR, "dada2/ASV.stats.tsv"),
        rds=os.path.join(RESULTS_DIR, "dada2/ASV.rds")
    log:
        os.path.join(RESULTS_DIR, "dada2/ASV.log")
    threads:
        config["dada2"]["asv"]["threads"]
    params:
        utils=os.path.join(SRC_DIR, "utils.R")
    conda:
        os.path.join(ENV_DIR, "r.yaml")
    message:
        "DADA2: merge"
    script:
        os.path.join(SRC_DIR, "dada2_asv_merge.R")    

rule asv_tree:
    input:
        rds=os.path.join(RESULTS_DIR, "dada2/ASV.rds")
    output:
        rds=os.path.join(RESULTS_DIR, "dada2/tree.rds"),
        tree=os.path.join(RESULTS_DIR, "dada2/ASV.tree")
    log:
        os.path.join(RESULTS_DIR, "dada2/tree.log")
    threads:
        config["dada2"]["tree"]["threads"]
    params:
        utils=os.path.join(SRC_DIR, "utils.R")
    conda:
        os.path.join(ENV_DIR, "tree.yaml")
    message:
        "Building an ASV tree"
    script:
        os.path.join(SRC_DIR, "asv_tree.R")



##################################################
# Stats

rule fastq_asv_stats:
    input:
        raw=os.path.join(RESULTS_DIR, "multiqc/fastqc/raw/multiqc_data/multiqc_fastqc.txt"),
        trim1=os.path.join(RESULTS_DIR, "multiqc/fastqc/cutadapt_primers/multiqc_data/multiqc_fastqc.txt"),
        trim2=os.path.join(RESULTS_DIR, "multiqc/fastqc/cutadapt_adapters/multiqc_data/multiqc_fastqc.txt"),
        asv=os.path.join(RESULTS_DIR, "dada2/ASV.stats.tsv")
    output:
        tsv=os.path.join(RESULTS_DIR, "dada2/stats.tsv"),
        pdf=os.path.join(RESULTS_DIR, "dada2/stats.pdf")
    log:
        os.path.join(RESULTS_DIR, "dada2/stats.log")
    params:
        utils=os.path.join(SRC_DIR, "utils.R")
    conda:
        os.path.join(ENV_DIR, "r.yaml")
    message:
        "FastQ/ASV stats"
    script:
        os.path.join(SRC_DIR, "fastq_asv_stats.R")


##################################################
# Krona

# Krona text input: https://github.com/marbl/Krona/wiki/Importing-text-and-XML-data
rule krona_text:
    input:
        counts=os.path.join(RESULTS_DIR, "dada2/ASV.counts.tsv"),
        tax=os.path.join(RESULTS_DIR, "dada2/ASV.tax.tsv")
    output:
        temp(expand(os.path.join(RESULTS_DIR, "dada2/ASV.krona.{sid}.txt"), sid=SAMPLES))
    threads:
        1
    message:
        "KronaTools: TEXT input from {input}"
    run:
        import re
        from pandas import read_csv
        # read in data
        counts = read_csv(input.counts, sep='\t', header=0, index_col=0) # sample x asv
        tax    = read_csv(input.tax,    sep='\t', header=0, index_col=0) # asv x taxonomy
        # create output files
        for sid in counts.index:
            sid_re = re.compile(".*\.krona\.%s\.txt$" % sid) # pattern should match the output file pattern
            with open(list(filter(sid_re.match, output))[0], "w") as ofile: # matching output file
                for asv in counts.columns:
                    if counts.loc[sid, asv] > 0:
                        assert asv in tax.index
                        c = counts.loc[sid, asv] # count
                        t = "\t".join(tax.loc[asv].fillna(value="NA")) # taxonomy
                        t = re.sub("^NA(\tNA)+$", "Unknown", t) # replace completely unknown taxonomy
                        t = re.sub("(\tNA)+$", "", t) # remove trailing NAs
                        ofile.write("%d\t%s\n" % (c, t))

rule krona_textimport:
    input:
        expand(os.path.join(RESULTS_DIR, "dada2/ASV.krona.{sid}.txt"), sid=SAMPLES)
    output:
        os.path.join(RESULTS_DIR, "dada2/ASV.krona.html")
    threads:
        1
    conda:
        os.path.join(ENV_DIR, "krona.yaml")
    message:
        "KronaTools: plot from TEXT input: {input}"
    shell:
        "ktImportText {input} -o {output} && sed -i 's/ASV\.krona\.//g' {output}"

##################################################
# Phyloseq

rule phyloseq:
    input:
        counts=os.path.join(RESULTS_DIR, "dada2/ASV.counts.tsv"),
        tax=os.path.join(RESULTS_DIR, "dada2/ASV.tax.tsv"),
        tree=os.path.join(RESULTS_DIR, "dada2/ASV.tree"),
        meta=config["samples_meta"]
    output:
        pdf=os.path.join(RESULTS_DIR, "dada2/ASV.phyloseq.pdf"),
        rds=os.path.join(RESULTS_DIR, "dada2/phyloseq.rds")
    log:
        os.path.join(RESULTS_DIR, "dada2/ASV.phyloseq.log")
    params:
        utils=os.path.join(SRC_DIR, "utils.R"),
        ord_method="DCA",
        ord_dist="bray",
    conda:
        os.path.join(ENV_DIR, "r.yaml")
    message:
        "Phyloseq"
    script:
        os.path.join(SRC_DIR, "phyloseq.R") 
