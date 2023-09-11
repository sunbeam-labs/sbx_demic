import os
import subprocess as sp
from pathlib import Path

input_dir = Path(snakemake.input[0])
output_dir = Path(snakemake.output[0])

os.makedirs(output_dir, exist_ok=True)

for sam in [fn for fn in os.listdir(input_dir) if fn.endswith(".sam")]:
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