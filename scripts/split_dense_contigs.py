import os
import shutil
import textwrap
from pathlib import Path
from sunbeamlib.parse import parse_fasta

contig_clusters = [Path(snakemake.params.base_dir) / fn for fn in os.listdir(snakemake.params.base_dir) if fn.endswith(".fasta")]

for cluster in contig_clusters:
    with open(cluster) as f:
        contig_descs = [x[0] for x in parse_fasta(f)]
    
    num_contigs = len(contig_descs)
    contig_lens = {x.split(" ")[0].removeprefix(">"): int(x.split("len=")[1].split(" ")[0]) for x in contig_descs}
    os.makedirs(snakemake.output[0], exist_ok=True)

    if num_contigs < 25:
        with open(cluster) as f_in, open(Path(snakemake.output[0]) / cluster.name, "w") as f_out:
            for header, seq in parse_fasta(f_in):
                contig_name = header.split(" ")[0].removeprefix(">")

                for i, subseq in enumerate(textwrap.wrap(seq, int(len(seq) / 25))):
                    header_new_name = header.replace(contig_name, f'{contig_name}_{i}')
                    header_new_len = header_new_name.replace(str(contig_lens[contig_name]), str(len(subseq)))
                    f_out.write(f">{header_new_len}\n{subseq}\n")
    else:
        shutil.copyfile(cluster, Path(snakemake.output[0]) / cluster.name)