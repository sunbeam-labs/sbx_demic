import numpy as np
import matplotlib.pyplot as plt

def generate_reads(genome, amplitude, shift, read_length):
    """Generate reads over the genome with a sin wave distribution of coverage."""
    genome_length = len(genome)
    subsample_factor = 50
    # Use half the period for the sine wave
    coverage = (amplitude * np.sin(np.pi * np.arange(genome_length) / (genome_length))) + shift
    coverage = coverage[0::subsample_factor]

    print(f'Generating {sum(coverage) / 1000000} million reads')

    for i in range(int(len(coverage) / subsample_factor)):
        num_reads_at_position = int(coverage[i])
        for _ in range(num_reads_at_position):
            start = max(0, i - read_length // 2)
            end = min(genome_length, start + read_length)
            read = genome[start:end]
            yield read

def complement_read(read):
    """Generate the complement of a DNA sequence."""
    comp = ""
    for c in read:
        if c == "A":
            comp += "T"
        elif c == "T":
            comp += "A"
        elif c == "G":
            comp += "C"
        elif c == "C":
            comp += "G"
    return comp

def write_fastq(output_prefix, reads, index, id):
    """Write reads to paired FASTQ files as they are generated."""
    with open(f'{output_prefix}_{index}_1.fastq', 'a') as f1, open(f'{output_prefix}_{index}_2.fastq', 'a') as f2:
        for i, read in enumerate(reads):
            f1.write(f'@Read{i + 1}_{id}/1\n{read}\n+\n{"I" * len(read)}\n')

            # Generate reversed and complemented sequence for the second FASTQ file
            comp_read = complement_read(read[::-1])
            f2.write(f'@Read{i + 1}_{id}/2\n{comp_read}\n+\n{"I" * len(comp_read)}\n')

def plot_coverage(coverage):
    """Plot the coverage distribution."""
    plt.plot(coverage)
    plt.title('Coverage Distribution')
    plt.xlabel('Genome Position')
    plt.ylabel('Coverage')
    plt.show()

def generate_and_write_reads(output_prefix, genome, amplitudes, shift, read_length, id):
    """Generate and write n paired read files with specified amplitudes."""
    for amplitude in amplitudes:
        reads_generator = generate_reads(genome, amplitude, shift, read_length)
        write_fastq(output_prefix, reads_generator, amplitude, id)

# Parameters
amplitudes = [2, 3, 4]  # Specify the amplitudes for each set of reads
shift = 1
read_length = 150
output_prefix = 'triple'

# Generate synthetic genome
with open("reference/Ecoli.fasta") as f:
    f.readline()  # Info line
    ecoli = "".join([l.strip() for l in f.readlines()])
with open("reference/Bfragilis.fasta") as f:
    f.readline()  # Info line
    bfragilis = "".join([l.strip() for l in f.readlines()])
with open("reference/akk-genome.fasta") as f:
    f.readline()  # Info line
    akk = "".join([l.strip() for l in f.readlines()])

# Generate and write n paired read files with specified amplitudes
generate_and_write_reads(output_prefix, ecoli, amplitudes, shift, read_length, "ecoli")
generate_and_write_reads(output_prefix, bfragilis, amplitudes, shift, read_length, "bfragilis")
generate_and_write_reads(output_prefix, akk, amplitudes, shift, read_length, "akk")
