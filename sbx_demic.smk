# -*- mode: Snakemake -*-

import os
import sys
from pathlib import Path
from sunbeamlib import samtools


TARGET_DEMIC = [
    str(MAPPING_FP/'demic'/'DEMIC_OUT'/'all_PTR.txt')
]
BINNED_DIR = str(ASSEMBLY_FP/'coassembly'/'max_bin')
CONTIGS_FASTA = BINNED_DIR + '/all_final_contigs.fa'

COASSEMBLY_DEMIC_FP = ASSEMBLY_FP / "coassembly_demic"

try:
    BENCHMARK_FP
except NameError:
    BENCHMARK_FP = output_subdir(Cfg, "benchmarks")
try:
    LOG_FP
except NameError:
    LOG_FP = output_subdir(Cfg, "logs")


def get_demic_path() -> str:
    for fp in sys.path:
        if fp.split("/")[-1] == "sbx_demic":
            return fp
    raise Error("Filepath for demic not found, are you sure it's installed under extensions/sbx_demic?")

rule all_demic:
    input:
        TARGET_DEMIC


def zip3l(l1, l2, l3):
    return list(zip(l1, l2, l3))


def coassembly_groups(fp, sample_list):
    if fp == "":
        K = ["all"] * (len(sample_list) * 2)
        V = list(sorted(sample_list)) * 2
        R = [1] * len(sample_list) + [2] * len(sample_list)
        return [K, V, R]
    groups = ruamel.yaml.safe_load(open(str(fp)).read())
    sorted_keys = sorted(groups.keys())
    K = []  # group
    V = []  # sample
    for k in sorted_keys:
        K += [k] * len(groups[k])
        V += groups[k]
    R = [1] * len(V) + [2] * len(V)
    return [K + K, V + V, R]


rule all_coassemble_demic:
    input:
        a=expand(
            str(COASSEMBLY_DEMIC_FP / "{group}_final_contigs.fa"),
            group=list(
                set(
                    coassembly_groups(
                        Cfg["coassembly_demic"]["group_file"], Samples.keys()
                    )[0]
                )
            ),
        ),
        b=expand(
            str(COASSEMBLY_DEMIC_FP / "agglomerate" / "{sample}_{group}_{rp}.fastq"),
            zip3l,
            group=coassembly_groups(
                Cfg["coassembly_demic"]["group_file"], Samples.keys()
            )[0],
            sample=coassembly_groups(
                Cfg["coassembly_demic"]["group_file"], Samples.keys()
            )[1],
            rp=coassembly_groups(
                Cfg["coassembly_demic"]["group_file"], Samples.keys()
            )[2],
        ),


rule prep_samples_for_concatenation_paired_demic:
    input:
        r1=str(QC_FP / "decontam" / "{sample}_1.fastq.gz"),
        r2=str(QC_FP / "decontam" / "{sample}_2.fastq.gz"),
    output:
        r1=temp(str(COASSEMBLY_DEMIC_FP / "agglomerate" / "{sample}_{group}_1.fastq")),
        r2=temp(str(COASSEMBLY_DEMIC_FP / "agglomerate" / "{sample}_{group}_2.fastq")),
    benchmark:
        BENCHMARK_FP / "prep_samples_for_concatenation_paired_{sample}_{group}.tsv"
    log:
        LOG_FP / "prep_samples_for_concatenation_paired_demic_{sample}_{group}.log",
    threads: Cfg["coassembly_demic"]["threads"]
    conda:
        "envs/sbx_demic_coassembly_env.yml"
    shell:
        """
        pigz -d -p {threads} -c {input.r1} > {output.r1}
        pigz -d -p {threads} -c {input.r2} > {output.r2}
        """


rule all_prep_paired:
    input:
        expand(str(ASSEMBLY_FP/'coassembly'/'agglomerate'/'{sample}_{group}_{rp}.fastq'), zip3l, group=coassembly_groups(Cfg['sbx_coassembly']['group_file'],Samples.keys())[0], sample=coassembly_groups(Cfg['sbx_coassembly']['group_file'],Samples.keys())[1], rp=coassembly_groups(Cfg['sbx_coassembly']['group_file'],Samples.keys())[2])
    output:
        touch(ASSEMBLY_FP/'coassembly'/'agglomerate'/'prepped.done')
    shell:
        'echo "samples ready to combine for co-assembly"'


rule combine_groups_paired_demic:
    input:
        rules.all_coassemble_demic.input.b,
    output:
        r1=str(COASSEMBLY_DEMIC_FP / "fastq" / "{group}_1.fastq.gz"),
        r2=str(COASSEMBLY_DEMIC_FP / "fastq" / "{group}_2.fastq.gz"),
    params:
        w1=str(str(COASSEMBLY_DEMIC_FP / "agglomerate") + str("/*{group}_1.fastq")),
        w2=str(str(COASSEMBLY_DEMIC_FP / "agglomerate") + str("/*{group}_2.fastq")),
    threads: Cfg["coassembly_demic"]["threads"]
    conda:
        "envs/sbx_demic_coassembly_env.yml"
    resources:
        runtime=120,
    shell:
        """
        cat {params.w1} | pigz -p {threads} > {output.r1}
        cat {params.w2} | pigz -p {threads} > {output.r2}
        """


rule coassemble_paired_demic:
    input:
        r1=str(COASSEMBLY_DEMIC_FP / "fastq" / "{group}_1.fastq.gz"),
        r2=str(COASSEMBLY_DEMIC_FP / "fastq" / "{group}_2.fastq.gz"),
    output:
        str(COASSEMBLY_DEMIC_FP / "{group}_final_contigs.fa"),
    benchmark:
        BENCHMARK_FP / "coassemble_paired_{group}.tsv"
    log:
        LOG_FP / "coassemble_paired_demic_{group}.log",
    params:
        assembly_dir=str(COASSEMBLY_DEMIC_FP / "{group}"),
    threads: Cfg["coassembly_demic"]["threads"]
    conda:
        "envs/sbx_demic_coassembly_env.yml"
    resources:
        mem_mb=20000,
        runtime=720,
    shell:
        """
        megahit -1 {input.r1} -2 {input.r2} -t {threads} -o {params.assembly_dir} 2>&1 | tee {log}
        mv {params.assembly_dir}/final.contigs.fa {output}
        """


rule maxbin:
    input:
        a = expand(str(ASSEMBLY_FP/'coassembly'/'{group}_final_contigs.fa'), group = list(set(coassembly_groups(Cfg['sbx_coassembly']['group_file'],Samples.keys())[0]))),
        b = rules.all_prep_paired.input
    output:
        str(Cfg['all']['output_fp']) + CONTIGS_FASTA
    benchmark:
        BENCHMARK_FP / "maxbin.tsv"
    log:
        LOG_FP / "maxbin.log",
    params:
        basename = str(Cfg['all']['output_fp']),
        binned_dir = str(Cfg['all']['output_fp']) + BINNED_DIR,
        contigs_fasta = str(Cfg['all']['output_fp']) + CONTIGS_FASTA,
        maxbin_dir=str(Path(get_demic_path()) / "MaxBin_2.2.7_scripts"),
        script=str(Path(get_demic_path()) / "MaxBin_2.2.7_scripts" / "run_MaxBin.pl"),
    conda:
        "envs/demic_bio_env.yml"
    shell:
        """
        find {params.basename}/qc/decontam -iname '*.fastq.gz' > {params.basename}/decontam_list
        mkdir -p {params.binned_dir}
        cp {params.basename}/assembly/coassembly/all_final_contigs.fa {output}
        
        if command -v MaxBin &> /dev/null
        then
            cd {params.maxbin_dir}
            {params.script} -thread 10 -contig {input.a} \
            -out {params.binned_dir} -reads_list {params.basename}/decontam_list \
            -verbose 2>&1 | tee {log}
        elif command -v run_MaxBin.pl &> /dev/null
        then
            run_MaxBin.pl -thread 10 -contig {input.a} \
            -out {params.binned_dir} -reads_list {params.basename}/decontam_list \
            -verbose 2>&1 | tee {log}
        else
            echo "Could not find MaxBin or run_MaxBin.pl in $PATH" > {log}
        fi
        """


rule bowtie2_build:
    input:
        str(Cfg['all']['output_fp']) + CONTIGS_FASTA
    params:
        basename = str(Cfg['all']['output_fp']) + CONTIGS_FASTA
    threads:
        Cfg['sbx_demic']['threads']
    output:
        touch(str(Cfg['all']['output_fp']) + CONTIGS_FASTA + '.1.bt2')
    conda:
        "envs/demic_bio_env.yml"
    shell:
        "bowtie2-build --threads {threads} {input} {params.basename}"

# Run bowtie2 with index
rule bowtie2:
    output:
        str(MAPPING_FP/'demic'/'raw'/'{sample}.sam')
    input:
        rules.bowtie2_build.output,
        reads = expand(
            str(QC_FP/'decontam'/'{sample}_{rp}.fastq.gz'),
            sample = Samples.keys(),
            rp = Pairs)
    threads:
        Cfg['sbx_demic']['threads']
    params:
        db_basename = str(Cfg['all']['output_fp']) + CONTIGS_FASTA
    conda:
        "envs/demic_bio_env.yml"
    shell:
        """
        bowtie2 -q -x {params.db_basename} \
        -1 {input.reads[0]} -2 {input.reads[1]} -p {threads} \
        -S {output}
        """

rule samtools_sort:
    input:
        str(MAPPING_FP/'demic'/'raw'/'{sample}.sam')
    output:
        temp_files = temp(str(MAPPING_FP/'demic'/'sorted'/'{sample}.bam')),
        sorted_files = str(MAPPING_FP/'demic'/'sorted'/'{sample}.sam')
    threads:
        Cfg['sbx_demic']['threads']
    conda:
        "envs/demic_bio_env.yml"
    log:
        str(MAPPING_FP/'demic'/'logs'/'samtools_{sample}.error')
    shell:
        """
        echo "converting to bam, sorting, and converting back to sam"
        samtools view -@ {threads} -bS {input} | samtools sort -@ {threads} - -o {output.temp_files} 2> {log}
        samtools view -@ {threads} -h {output.temp_files} > {output.sorted_files} 2>> {log}
        """

# TODO
# how to get the directory of this output:
#        str(MAPPING_FP/'demic'/'sorted'/'{sample}.sam')
# and how to get the directory of:
#       CONTIGS_FASTA
# because those are the inputs of the next rule
#
# Maybe this will work:
#
# os.path.dirname

rule run_demic:
    input:
        expand(str(MAPPING_FP/'demic'/'sorted'/'{sample}.sam'),
        sample = Samples.keys())
    output:
        str(MAPPING_FP/'demic'/'DEMIC_OUT'/'all_PTR.txt')
    params:
        r_installer = get_demic_path() + "/envs/install.R",
        demic = get_demic_path() + "/vendor_demic_v1.0.2/DEMIC.pl",
        sam_dir = str(MAPPING_FP/'demic'/'sorted'),
        fasta_dir = str(Cfg['all']['output_fp']) + BINNED_DIR,
        keep_all = Cfg['sbx_demic']['keepall'],
        extras = Cfg['sbx_demic']['extras'],
    threads:
        Cfg['sbx_demic']['threads']
    conda:
        "envs/demic_env.yml"
    log:
        str(MAPPING_FP/'demic'/'logs'/'demic.error')
    # Rscript {params.r_installer} && \
    shell:
        """
        
        {params.demic} --output_all {params.keep_all} {params.extras} \
        --thread_num {threads} \
        -S {params.sam_dir} -F {params.fasta_dir} \
        -O $(dirname {output}) 2> {log}
        """

rule aggregate_demic:
    input:
        expand(str(MAPPING_FP/'demic'/'DEMIC_OUT'/'{group}'/'all_PTR.txt'), group = list(set(coassembly_groups(Cfg['sbx_coassembly']['group_file'],Samples.keys())[0])))
    output:
        MAPPING_FP / "demic" / "DEMIC_OUT" / "all_PTR.txt"
    shell:
        "cat {input} > output"