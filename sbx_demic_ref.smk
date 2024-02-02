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


localrules:
    all_demic_ref,


rule all_demic_ref:
    input:
        DEMIC_REF_FP / "all_PTR.txt"


rule install_opera_ms:
    output:
        DEMIC_REF_FP / ".installed"
    params:
        demic_fp=str(get_demic_path()),
    conda:
        "envs/demic_ref_env.yml"
    shell:
        """
        cd {params.demic_fp}
        git clone https://github.com/CSB5/OPERA-MS.git
        cd OPERA-MS
        make
        perl OPERA-MS.pl check-dependency
        """

rule run_opera_ms:
    input:
        DEMIC_REF_FP / ".installed",
        reads=QC_FP / "decontam" / "{sample}_{rp}.fastq.gz",
        contigs=Cfg["sbx_demic"]["ref_fp"],








rule run_pycov3_ref:
    input:
        DEMIC_FP / "sorted",
        COASSEMBLY_DEMIC_FP / "max_bin" / "max_bin",
    output:
        directory(DEMIC_FP / "pycov3"),
    params:
        sam_dir=str(DEMIC_FP / "sorted"),
        fasta_dir=str(COASSEMBLY_DEMIC_FP / "max_bin"),
        extras=Cfg["sbx_demic"]["extras"],
    threads: Cfg["sbx_demic"]["demic_threads"]
    resources:
        mem_mb=20000,
        runtime=720,
    conda:
        "envs/demic_bio_env.yml"
    log:
        LOG_FP / "run_pycov3.log",
    shell:
        """
        pycov3 -S {params.sam_dir} -F {params.fasta_dir} -O {output} -X {params.extras} 2>&1 | tee {log}
        """


rule run_demic_ref:
    input:
        input=DEMIC_FP / "pycov3",
        installed=DEMIC_FP / ".installed",
    output:
        out=directory(DEMIC_FP / "DEMIC_OUT"),
    threads: Cfg["sbx_demic"]["demic_threads"]
    resources:
        mem_mb=20000,
        runtime=720,
    conda:
        "envs/demic_env.yml"
    log:
        LOG_FP / "run_demic.log",
    script:
        "scripts/run_demic.R"


rule aggregate_demic_ref:
    input:
        DEMIC_FP / "DEMIC_OUT",
    output:
        DEMIC_REF_FP / "all_PTR.txt",
    shell:
        """
        echo "sample\testPTR\tcoefficient\tpValue\tcor\tcorrectY" > {output}
        tail -n +1 {input}/*.ptr >> {output}
        """