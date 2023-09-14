import os
import subprocess as sp
from pathlib import Path

# Original bash command:
# bowtie2 -q -x {input.fasta} -1 {input.reads[0]} -2 {input.reads[1]} -p {threads} -S {output}

fasta_dir = Path(snakemake.params.base_dir)
fastq_dir = Path(snakemake.params.reads_dir)
output_dir = Path(snakemake.output[0])

os.makedirs(output_dir, exist_ok=True)

fastas = [fn for fn in os.listdir(fasta_dir) if fn.endswith(".fasta")]
fastqs = [fn for fn in os.listdir(fastq_dir) if fn.endswith("_1.fastq.gz")]
fastas_len = len(fastas)
fastqs_len = len(fastqs)

for i, fasta in enumerate(fastas):
    for j, fastq in enumerate(fastqs):
        print(
            f"Progress: {str(round(100 * (i / fastas_len + (j / fastqs_len) / fastas_len), 2))}%"
        )
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
