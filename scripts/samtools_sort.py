import os
import subprocess as sp
from pathlib import Path

# Original bash commands:
# samtools view -@ {threads} -bS {input} | samtools sort -@ {threads} - -o {output.temp_files} 2> {log}
# samtools view -@ {threads} -h {output.temp_files} > {output.sorted_files} 2>> {log}

input_dir = Path(snakemake.input[0])
output_dir = Path(snakemake.output[0])

os.makedirs(output_dir, exist_ok=True)

sams = [fn for fn in os.listdir(input_dir) if fn.endswith(".sam")]
sams_len = len(sams)

for i, sam in enumerate(sams):
    print(f"Progress: {str(round(100 * (i / sams_len), 2))}%")

    fn = sam.replace(".sam", "")

    args1 = [
        "samtools",
        "view",
        "-@",
        str(snakemake.threads),
        "-bS",
        str(input_dir / f"{fn}.sam"),
    ]
    args2 = [
        "samtools",
        "sort",
        "-@",
        str(snakemake.threads),
        "-",
        "-o",
        str(output_dir / f"{fn}.bam"),
    ]
    args3 = [
        "samtools",
        "view",
        "-@",
        str(snakemake.threads),
        "-h",
        str(output_dir / f"{fn}.bam"),
    ]

    ps = sp.run(args1, check=True, capture_output=True)
    bam = sp.run(args2, input=ps.stdout, capture_output=True)
    with open(output_dir / f"{fn}.sam", "w") as f_out:
        sp.run(args3, stdout=f_out)

for bam in [fn for fn in os.listdir(output_dir) if fn.endswith(".bam")]:
    os.remove(output_dir / bam)
