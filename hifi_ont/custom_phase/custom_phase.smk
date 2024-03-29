import os, sys
import pandas as pd

SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)

if not os.path.exists("log"):
    os.makedirs("log")


configfile: "hifiasm_custom.yaml"


MANIFEST = config.get("MANIFEST", "manifest_custom.tab")
VERSION = config.get("HIFIASM_VERSION", "0.19.6")
ASM_THREADS = config.get("ASM_THREADS", 16)

manifest_df = pd.read_csv(
    MANIFEST, sep="\t", header=0, index_col="sample", dtype=str, comment="#"
).fillna("NA")


def find_hifi_fofn(wildcards):
    return manifest_df.at[wildcards.sample, "hifi_fofn"]


def find_ont_fofn(wildcards):
    return manifest_df.at[wildcards.sample, "ont_fofn"]


def find_yak(wildcards):
    return [manifest_df.at[wildcards.sample, x] for x in ["hap1_yak", "hap2_yak"]]


def find_ont_files(wildcards):
    with open(manifest_df.at[wildcards.sample, "ont_fofn"], "r") as infile:
        ont_files = [line.rstrip() for line in infile]
    return ont_files


def find_hifi_files(wildcards):
    with open(manifest_df.at[wildcards.sample, "hifi_fofn"], "r") as infile:
        hifi_files = [line.rstrip() for line in infile]
    return hifi_files


def nanopore_string(wildcards):
    with open(manifest_df.at[wildcards.sample, "ont_fofn"], "r") as infile:
        ont_files = [line.rstrip() for line in infile]
    return ",".join(ont_files)


def find_final(wildcards):
    trio = expand(
        "{sample}/assemblies/hifiasm_ont/phase/{version}/{sample}.hifiasm.dip.{hap}.p_ctg.gfa.fasta",
        hap=["hap1", "hap2"],
        sample=manifest_df.index,
        version=VERSION,
    )
    n50_stats_dip = expand(
        "{sample}/assemblies/hifiasm_ont/{sample}-{phase}-custom_n50-{version}.txt",
        sample=manifest_df.index,
        phase="dip",
        version=VERSION,
    )
    return trio


# return trio + n50_stats_dip


wildcard_constraints:
    version=VERSION,
    sample="|".join(manifest_df.index),
    phase="dip|bp",


localrules:
    all,
    make_fai,


rule all:
    input:
        find_final,


rule gather_psuedo:
    input:
        expand(
            "{sample}/assemblies/hifiasm_ont/{version}/{sample}.hifiasm.bp.p_utg.gfa",
            sample=manifest_df.index,
            version=VERSION,
        ),


rule hifiasm_prim_gfa:
    input:
        hifi_fofn=find_hifi_fofn,
        ont_fofn=find_ont_fofn,
        hifi_in=find_hifi_files,
        ont_in=find_ont_files,
    output:
        gfa="{sample}/assemblies/hifiasm_ont/{version}/{sample}.hifiasm.bp.p_utg.gfa",
    threads: ASM_THREADS
    params:
        nanopore=nanopore_string,
    resources:
        mem=int(336 / ASM_THREADS),
        hrs=900,
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        f"hifiasm/{VERSION}",
    shell:
        """
        hifiasm -o {wildcards.sample}/assemblies/hifiasm_ont/{wildcards.version}/{wildcards.sample}.hifiasm -t {threads} --ul {params.nanopore} $( cat {input.hifi_fofn} )
        """


rule hifiasm_trio:
    input:
        yak=find_yak,
        asm=rules.hifiasm_prim_gfa.output.gfa,
    output:
        asm_pat="{sample}/assemblies/hifiasm_ont/{version}/{sample}.hifiasm.dip.hap1.p_ctg.gfa",
        asm_mat="{sample}/assemblies/hifiasm_ont/{version}/{sample}.hifiasm.dip.hap2.p_ctg.gfa",
    threads: 8
    resources:
        mem=int(336 / 8),
        hrs=128,
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        f"hifiasm/{VERSION}",
    shell:
        """
        hifiasm -o {wildcards.sample}/assemblies/hifiasm_ont/{wildcards.version}/{wildcards.sample}.hifiasm -t {threads} -1 {input.yak[0]} -2 {input.yak[1]} --ul /dev/null /dev/null
        """


rule make_fasta:
    input:
        asm_pat="{sample}/assemblies/hifiasm_ont/{version}/{sample}.hifiasm.{phase}.hap1.p_ctg.gfa",
        asm_mat="{sample}/assemblies/hifiasm_ont/{version}/{sample}.hifiasm.{phase}.hap2.p_ctg.gfa",
    output:
        fa_pat="{sample}/assemblies/hifiasm_ont/{version}/{sample}.hifiasm.{phase}.hap1.p_ctg.gfa.fasta",
        fa_mat="{sample}/assemblies/hifiasm_ont/{version}/{sample}.hifiasm.{phase}.hap2.p_ctg.gfa.fasta",
    threads: 1
    resources:
        mem=4,
        hrs=1,
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
    shell:
        """
        {SNAKEMAKE_DIR}/make_fasta.sh {input.asm_pat} > {output.fa_pat}
        {SNAKEMAKE_DIR}/make_fasta.sh {input.asm_mat} > {output.fa_mat}
        if [[ $( echo {wildcards.phase} ) == 'dip' ]]; then
            ln -s $( basename {output.fa_pat} ) $( echo {output.fa_pat} | sed 's/hap1/pat/' )
            ln -s $( basename {output.fa_mat} ) $( echo {output.fa_mat} | sed 's/hap2/mat/' )
        fi
        """


rule make_fai:
    input:
        fa_pat=rules.make_fasta.output.fa_pat,
        fa_mat=rules.make_fasta.output.fa_mat,
    output:
        fai_pat="{sample}/assemblies/hifiasm_ont/{version}/{sample}.hifiasm.{phase}.hap1.p_ctg.gfa.fasta.fai",
        fai_mat="{sample}/assemblies/hifiasm_ont/{version}/{sample}.hifiasm.{phase}.hap2.p_ctg.gfa.fasta.fai",
    threads: 1
    resources:
        mem=4,
        hrs=1,
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        "samtools/1.12",
    shell:
        """
        samtools faidx {input.fa_pat}
        samtools faidx {input.fa_mat}
        """


rule get_n50:
    input:
        fai_pat=rules.make_fai.output.fai_pat,
        fai_mat=rules.make_fai.output.fai_mat,
    output:
        n50_stats="{sample}/assemblies/hifiasm_ont/{sample}-{phase}_n50-{version}.txt",
    threads: 1
    resources:
        mem=1,
        hrs=1,
    shell:
        """
        script=/net/eichler/vol26/7200/software/pipelines/compteam_tools/n50
        # n50 stats for paternal haplotype
        echo $(basename {input.fai_pat}) > {output.n50_stats}
        $script {input.fai_pat} >> {output.n50_stats}

        # n50 stats for maternal haplotype
        echo $(basename {input.fai_mat}) >> {output.n50_stats}
        $script {input.fai_mat} >> {output.n50_stats}
        """


rule organize_trio:
    input:
        fa_pat="{sample}/assemblies/hifiasm_ont/{version}/{sample}.hifiasm.dip.hap1.p_ctg.gfa.fasta",
        fa_mat="{sample}/assemblies/hifiasm_ont/{version}/{sample}.hifiasm.dip.hap2.p_ctg.gfa.fasta",
        fai_pat="{sample}/assemblies/hifiasm_ont/{version}/{sample}.hifiasm.dip.hap1.p_ctg.gfa.fasta.fai",
        fai_mat="{sample}/assemblies/hifiasm_ont/{version}/{sample}.hifiasm.dip.hap2.p_ctg.gfa.fasta.fai",
    output:
        fa_pat="{sample}/assemblies/hifiasm_ont/phase/{version}/{sample}.hifiasm.dip.hap1.p_ctg.gfa.fasta",
        fa_mat="{sample}/assemblies/hifiasm_ont/phase/{version}/{sample}.hifiasm.dip.hap2.p_ctg.gfa.fasta",
        fai_pat="{sample}/assemblies/hifiasm_ont/phase/{version}/{sample}.hifiasm.dip.hap1.p_ctg.gfa.fasta.fai",
        fai_mat="{sample}/assemblies/hifiasm_ont/phase/{version}/{sample}.hifiasm.dip.hap2.p_ctg.gfa.fasta.fai",
    threads: 1
    resources:
        mem=1,
        hrs=1,
    shell:
        """
        cp -l {input.fa_pat} {output.fa_pat}
        cp -l {input.fa_mat} {output.fa_mat}
        cp -l {input.fai_pat} {output.fai_pat}
        cp -l {input.fai_mat} {output.fai_mat}
        ln -sf $( basename {output.fa_pat} ) $( echo {output.fa_pat} | sed 's/hap1/pat/' )
        ln -sf $( basename {output.fa_mat} ) $( echo {output.fa_mat} | sed 's/hap2/mat/' )
        """
