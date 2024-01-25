import os

samples = [
    line.rstrip().replace(".hg38.cram", "")
    for line in open("/home/wdecoster/p200/1000G/all_cram_names_ftp.txt")
    if line
]

outdir = "/home/wdecoster/p200/1000G/"
ref = "/home/wdecoster/database/1KG_ONT_VIENNA_hg38.fa"


rule all:
    input:
        strdust=expand(
            os.path.join(outdir, "pathSTR_STRdust/{sample}.vcf.gz"),
            sample=samples,
        ),


rule samtools_view:
    output:
        os.path.join(outdir, "crams_pathSTR/{sample}.cram"),
    params:
        sample="{sample}",
        targets="/home/wdecoster/p200/data/hg38.STRchive-disease-loci.TRGT.bed",
        ref=ref,
    log:
        "logs/samtools_view_pathSTR/{sample}.log",
    conda:
        "/home/wdecoster/p200/1000G/envs/samtools.yml"
    shell:
        """samtools view -T {params.ref} \
        --region-file {params.targets} \
        -o {output} \
        https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1KG_ONT_VIENNA/hg38/{params.sample}.hg38.cram 2> {log}
        """


rule index_crams:
    input:
        os.path.join(outdir, "crams_pathSTR/{sample}.cram"),
    output:
        os.path.join(outdir, "crams_pathSTR/{sample}.cram.crai"),
    log:
        "logs/index_crams_pathSTR/{sample}.log",
    conda:
        "/home/wdecoster/p200/1000G/envs/samtools.yml"
    shell:
        "samtools index {input}"


rule strdust_golga8a_unphased:
    input:
        cram=os.path.join(outdir, "crams_pathSTR/{sample}.cram"),
        crai=os.path.join(outdir, "crams_pathSTR/{sample}.cram.crai"),
    output:
        os.path.join(outdir, "pathSTR_STRdust/{sample}.vcf.gz"),
    log:
        "logs/pathSTR_STRdust/{sample}.log",
    params:
        ref=ref,
        targets="/home/wdecoster/p200/data/hg38.STRchive-disease-loci.TRGT.bed",
        binary="/home/wdecoster/repositories/STRdust/target/release/STRdust",
    conda:
        "/home/wdecoster/p200/1000G/envs/samtools.yml"
    shell:
        """RUST_LOG=debug {params.binary} \
        -R {params.targets} \
        --support 2 \
        --unphased \
        --find-outliers \
        --somatic \
        {params.ref} \
        {input.cram} 2> {log} \
        | bgzip > {output} 2>> {log}
        """
