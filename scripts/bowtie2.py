import os
import subprocess as sp
from pathlib import Path

# Original bash command:
# bowtie2 -q -x {input.fasta} -1 {input.reads[0]} -2 {input.reads[1]} -p {threads} -S {output}

fasta_dir = Path(snakemake.params.base_dir)
fastq_dir = Path(snakemake.params.reads_dir)
output_dir = Path(snakemake.output[0])

os.makedirs(output_dir, exist_ok=True)

for fasta in [fn for fn in os.listdir(fasta_dir) if fn.endswith(".fasta")]:
    for fastq in [fn for fn in os.listdir(fastq_dir) if fn.endswith("_1.fastq.gz")]:
        sample = fastq.replace("_1.fastq.gz", "")
        bin = fasta.split(".")[1]

        args = [
            "bowtie2",
            "-q",
            "-x",
            str(fasta_dir / fasta),
            "-1",
            str(fastq_dir / fastq),
            "-2",
            str(fastq_dir / f"{sample}_2.fastq.gz"),
            "-p",
            str(snakemake.threads),
            "-S",
            str(output_dir / f"{sample}_{bin}.sam"),
        ]
        sp.run(args)
