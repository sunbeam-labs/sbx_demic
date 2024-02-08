try:
    BENCHMARK_FP
except NameError:
    BENCHMARK_FP = output_subdir(Cfg, "benchmarks")
try:
    LOG_FP
except NameError:
    LOG_FP = output_subdir(Cfg, "logs")


DEMIC_FP = MAPPING_FP / "demic"
DEMIC_REF_FP = DEMIC_FP / "ref"


def get_demic_path() -> Path:
    for fp in sys.path:
        if fp.split("/")[-1] == "sbx_demic":
            return Path(fp)
    raise Error(
        "Filepath for demic not found, are you sure it's installed under extensions/sbx_demic?"
    )


localrules:
    all_demic_ref,


rule all_demic_ref:
    input:
        all=DEMIC_REF_FP / "all_PTR.txt",
        contig=DEMIC_REF_FP / "contig_PTR.txt",
        sample=DEMIC_REF_FP / "sample_PTR.txt",
        hist=expand(
            DEMIC_REF_FP / "coverage" / "{sample}_001.txt", sample=Samples.keys()
        ),
        depth=expand(
            DEMIC_REF_FP / "coverage" / "{sample}_001.depth", sample=Samples.keys()
        ),
        tsv=expand(
            DEMIC_REF_FP / "coverage" / "{sample}_001.tsv", sample=Samples.keys()
        ),


rule spades_assemble_with_ref:
    input:
        r1=expand(QC_FP / "decontam" / "{sample}_1.fastq.gz", sample=Samples.keys()),
        r2=expand(QC_FP / "decontam" / "{sample}_2.fastq.gz", sample=Samples.keys()),
        contigs=Cfg["sbx_demic"]["ref_fp"],
    output:
        r1=temp(ASSEMBLY_FP / "demic_ref" / "contigs" / "spades.001.1.fastq"),
        r2=temp(ASSEMBLY_FP / "demic_ref" / "contigs" / "spades.001.2.fastq"),
        r1s=temp(expand(QC_FP / "decontam" / "{sample}_1.fastq", sample=Samples.keys())),
        r2s=temp(expand(QC_FP / "decontam" / "{sample}_2.fastq", sample=Samples.keys())),
        contigs=ASSEMBLY_FP / "demic_ref" / "contigs" / "spades.001.fasta",
    params:
        output_dir=str(ASSEMBLY_FP / "demic_ref"),
    threads: 8
    conda:
        "envs/demic_ref_env.yml"
    shell:
        """
        gzip -dk {input.r1}
        gzip -dk {input.r2}

        cat {output.r1s} > {output.r1}
        cat {output.r2s} > {output.r2}

        spades.py -1 {output.r1} -2 {output.r2} -o {params.output_dir} -t {threads} --trusted-contigs {input.contigs}
        mv {params.output_dir}/contigs.fasta {output.contigs}
        """


rule bowtie2_build_demic_ref:
    input:
        ASSEMBLY_FP / "demic_ref" / "contigs" / "spades.001.fasta",
    output:
        ASSEMBLY_FP / "demic_ref" / "contigs" / "spades.001.fasta.1.bt2",
    threads: Cfg["sbx_demic"]["demic_threads"]
    conda:
        "envs/demic_bio_env.yml"
    shell:
        """
        bowtie2-build --threads {threads} {input} {input}
        """


rule bowtie2_demic_ref:
    input:
        contigs=ASSEMBLY_FP / "demic_ref" / "contigs" / "spades.001.fasta",
        indexes=ASSEMBLY_FP / "demic_ref" / "contigs" / "spades.001.fasta.1.bt2",
        reads=expand(
            QC_FP / "decontam" / "{{sample}}_{rp}.fastq.gz",
            rp=Pairs,
        ),
    output:
        DEMIC_REF_FP / "raw" / "{sample}_001.sam",
    threads: Cfg["sbx_demic"]["demic_threads"]
    params:
        reads_dir=str(QC_FP / "decontam"),
    conda:
        "envs/demic_bio_env.yml"
    shell:
        """
        bowtie2 -q -x {input.contigs} -1 {input.reads[0]} -2 {input.reads[1]} -p {threads} -S {output}
        """


rule samtools_sort_demic_ref:
    input:
        DEMIC_REF_FP / "raw" / "{sample}_001.sam",
    output:
        temp_files=temp(DEMIC_REF_FP / "sorted" / "{sample}_001.bam"),
        sorted_files=DEMIC_REF_FP / "sorted" / "{sample}_001.sam",
    threads: Cfg["sbx_demic"]["demic_threads"]
    conda:
        "envs/demic_bio_env.yml"
    shell:
        """
        samtools view -@ {threads} -bS {input} | samtools sort -@ {threads} - -o {output.temp_files}
        samtools view -@ {threads} -h {output.temp_files} > {output.sorted_files}
        """


rule samtools_coverage_ref:
    input:
        DEMIC_REF_FP / "sorted" / "{sample}_001.sam",
    output:
        hist=DEMIC_REF_FP / "coverage" / "{sample}_001.txt",
        depth=DEMIC_REF_FP / "coverage" / "{sample}_001.depth",
        tsv=DEMIC_REF_FP / "coverage" / "{sample}_001.tsv",
    conda:
        "envs/demic_bio_env.yml"
    shell:
        """
        samtools coverage {input} -m -o {output.hist}
        samtools coverage {input} -D -o {output.depth}
        samtools coverage {input} -o {output.tsv}
        """


rule run_pycov3_ref:
    input:
        sams=expand(DEMIC_REF_FP / "sorted" / "{sample}_001.sam", sample=Samples.keys()),
        contigs=ASSEMBLY_FP / "demic_ref" / "contigs" / "spades.001.fasta",
    output:
        DEMIC_REF_FP / "pycov3" / "spades.001.cov3",
    params:
        sam_dir=str(DEMIC_REF_FP / "sorted"),
        fasta_dir=str(ASSEMBLY_FP / "demic_ref" / "contigs"),
        output_dir=str(DEMIC_REF_FP / "pycov3"),
        extras=Cfg["sbx_demic"]["extras"],
    threads: Cfg["sbx_demic"]["demic_threads"]
    resources:
        mem_mb=20000,
        runtime=720,
    conda:
        "envs/demic_bio_env.yml"
    log:
        LOG_FP / "run_pycov3_ref.log",
    benchmark:
        BENCHMARK_FP / "run_pycov3_ref.txt"
    shell:
        """
        pycov3 -S {params.sam_dir} -F {params.fasta_dir} -O {params.output_dir} -X {params.extras} 2>&1 | tee {log}
        """


rule run_demic_ref:
    input:
        input=DEMIC_REF_FP / "pycov3" / "spades.001.cov3",
        #installed=DEMIC_FP / ".installed",
    output:
        all=DEMIC_REF_FP / "DEMIC_OUT" / "spades.001.all.ptr",
        contig=DEMIC_REF_FP / "DEMIC_OUT" / "spades.001.contig.ptr",
        sample=DEMIC_REF_FP / "DEMIC_OUT" / "spades.001.sample.ptr",
    params:
        output_dir=str(DEMIC_REF_FP / "DEMIC_OUT"),
    threads: Cfg["sbx_demic"]["demic_threads"]
    resources:
        mem_mb=20000,
        runtime=720,
    conda:
        "envs/demic_env.yml"
    log:
        LOG_FP / "run_demic.log",
    script:
        "scripts/run_demic_ref.R"


rule aggregate_demic_ref:
    input:
        all=expand(
            DEMIC_REF_FP / "DEMIC_OUT" / "spades.001.all.ptr", sample=Samples.keys()
        ),
        contig=expand(
            DEMIC_REF_FP / "DEMIC_OUT" / "spades.001.contig.ptr", sample=Samples.keys()
        ),
        sample=expand(
            DEMIC_REF_FP / "DEMIC_OUT" / "spades.001.sample.ptr", sample=Samples.keys()
        ),
    output:
        all=DEMIC_REF_FP / "all_PTR.txt",
        contig=DEMIC_REF_FP / "contig_PTR.txt",
        sample=DEMIC_REF_FP / "sample_PTR.txt",
    params:
        dir=DEMIC_REF_FP / "DEMIC_OUT",
    shell:
        """
        echo "sample\testPTR\tcoefficient\tpValue\tcor\tcorrectY" | tee {output.all} {output.contig} {output.sample} > /dev/null

        tail -n +1 {params.dir}/*.all.ptr >> {output.all}

        tail -n +1 {params.dir}/*.contig.ptr >> {output.contig}

        tail -n +1 {params.dir}/*.sample.ptr >> {output.sample}
        """
