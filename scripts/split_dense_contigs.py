import os
from pathlib import Path
from sunbeamlib.parse import parse_fasta

MAX_CONTIG_LEN = 10000000

contig_clusters = [Path(snakemake.params.base_dir) / fn for fn in os.listdir(snakemake.params.base_dir) if fn.endswith(".fasta")]

os.makedirs(snakemake.output[0], exist_ok=True)

for cluster in contig_clusters: # Handle files in dir from MaxBin
    with open(cluster) as f_in, open(Path(snakemake.output[0]) / cluster.name, "w") as f_out:
        for header, seq in parse_fasta(f_in): # Handle records in file
            contig_name = header.split(" ")[0].removeprefix(">")
            contig_len = len(seq)

            if contig_len < MAX_CONTIG_LEN:
                f_out.write(f">{contig_name.strip()}\n{seq.strip()}\n")
            else:
                count = 1
                while contig_len > MAX_CONTIG_LEN:
                    f_out.write(f">{contig_name.strip()}_{int(count)}\n{seq[:MAX_CONTIG_LEN].strip()}\n")
                    count += 1
                    seq = seq[MAX_CONTIG_LEN:]
                    contig_len = len(seq)
                f_out.write(f">{contig_name.strip()}_{int(count)}\n{seq.strip()}\n")