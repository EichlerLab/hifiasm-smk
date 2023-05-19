import os, sys
import pandas as pd

SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)

configfile: "hifiasm.yaml"

MANIFEST = config['MANIFEST']

manifest_df = pd.read_csv(MANIFEST, sep='\t', header=0, index_col='sample')

shell.prefix("source {SNAKEMAKE_DIR}/env.cfg; ")

localrules: cat_hapotypes, combined_fa_to_fai

if not os.path.exists("log"):
	os.makedirs("log")



def find_sample_fofn(wildcards):
    return manifest_df.at[wildcards.sample,'fofn']

rule all:
	input:
		expand("{sample}/{sample}.hifiasm.bp.combined.fa.fai", sample=manifest_df.index)

rule hifiasm:
    input:
        sample_fofn_hifi = find_sample_fofn
    output:
        a_contig = "{sample}/{sample}.hifiasm.bp.hap1.p_ctg.gfa",
        p_contig = "{sample}/{sample}.hifiasm.bp.hap2.p_ctg.gfa"
    threads: 16
    resources:
        mem = 10,
        hrs = 24
    shell:
        '''
        hifiasm -o {wildcards.sample}/{wildcards.sample}.hifiasm -t{threads} $(cat {input.sample_fofn_hifi})
        '''


rule gfa_to_fasta:
	input: "{sample}/{sample}.hifiasm.bp.hap{haplotype}.p_ctg.gfa"
	output: "{sample}/{sample}.hifiasm.bp.hap{haplotype}.p_ctg.gfa.fa"
	resources: mem = 10
	run:
		szCommand = SNAKEMAKE_DIR + "/gfa_to_fasta.sh " + str( input ) + " >" + str( output )
		print( "about to execute: " + szCommand )
		shell( szCommand )

rule cat_hapotypes:
	input: "{sample}/{sample}.hifiasm.bp.hap1.p_ctg.gfa.fa", "{sample}/{sample}.hifiasm.bp.hap2.p_ctg.gfa.fa"
	output: "{sample}/{sample}.hifiasm.bp.combined.fa"
	run:
		szCommand = "cat " + str( input[0] ) + " " + str( input[1] ) + " >" + str( output )
		print( "about to execute: " + szCommand )
		shell( szCommand )

rule combined_fa_to_fai:
	input: "{sample}/{sample}.hifiasm.bp.combined.fa"
	output: "{sample}/{sample}.hifiasm.bp.combined.fa.fai"
	run:
		szCommand = "module load samtools/1.9 && samtools faidx " + str( input )
		print( "about to execute: " + szCommand )
		shell( szCommand )
