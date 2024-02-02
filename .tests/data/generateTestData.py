from io import TextIOWrapper
import random
import matplotlib.pyplot as plt
from collections import Counter
import gzip
import os
import math


def triangle(position: float, PTR: float, min_reads: int) -> int:
    assert PTR > 1
    assert 0 <= position <= 1
    if position < 0.5:
        return int((PTR - 1) * min_reads * 2 * position) + min_reads
    else:
        return int((PTR - 1) * min_reads * 2 * (1 - position)) + min_reads


def get_read(genome: str, start_index: int, read_length: int) -> str:
    genome_length = len(genome)
    start_index = start_index % genome_length  # Make sure start_index is on genome
    if start_index + read_length < genome_length:
        return genome[start_index : start_index + read_length]
    else:
        spillover = start_index + read_length - genome_length
        return genome[start_index:] + genome[:spillover]
    
def adjust_genome(genome: str, name_str: str) -> str:
    # Return genome str adjusted to have the replication origin in the middle if its location is known
    ref_dict = {"Bacteroides_fragilis_638R": 105766, "Bifidobacterium_animalis_lactis_B420": 1886, "Streptomyces_scabiei_87_22": 5138728}
    for species, rep_index in ref_dict.items():
        if species in name_str:
            start_at_rep = genome[rep_index:] + genome[:rep_index]
            return start_at_rep[len(genome) // 2:] + start_at_rep[:len(genome) // 2]
        
    return genome


# Generates a triangular distribution of reads from input sequence data
# mimicing the distribution of reads from a growing bacterial population
# @param genome is the path to the txt file containing only sequences (our test genome)
#        To manually prepare this file, just make a copy of your fasta and remove the
#        lines starting with >
# @param reads is the path to the output fastq file of generated reads WITHOUT the file extension
def generateTestData(genome_fp: str, reads: str, num_reads: int, PTR: float) -> Counter:
    genome: str = ""
    with open(genome_fp) as genome_f:
        for r in genome_f.readlines():
            if r[0] != ">":
                genome += r.strip()
    genome = adjust_genome(genome, genome_fp)

    with open(reads + "_R1.fastq", "w") as readsF_R1, open(
        reads + "_R2.fastq", "w"
    ) as readsF_R2:
        # Generate n sequences from the replication origin (assumed to be the start of the file)
        # of random length (0 to len(genome)) and then take a random substring of length m
        # from that to add to the fastq file
        c: Counter = Counter()
        genome_length: int = len(genome)
        read_length: int = 150
        min_triangle: int = 100

        def D(x):
            # Sine probability density function
            return (1 / genome_length) * (1 - ((1 - PTR) / (1 + PTR)) * math.cos((2 * math.pi * x) / genome_length))

        for i in range(num_reads):
            # Randomly sample sine probability distribution
            x = random.uniform(0, genome_length)
            y = random.uniform(0, D(0))
            while y > D(x):
                x = random.uniform(0, genome_length)
                y = random.uniform(0, D(0))
            
            start_index = int(x) - read_length
            forward_read = get_read(genome, start_index, read_length)
            reverse_read = get_read(genome, start_index + read_length, read_length)

            if random.random() < 0.05:
                writeReads(readsF_R1, forward_read, i)
                writeReads(readsF_R2, complementRead(reverse_read[::-1]), i)
        
        return c



        for initial_index in [0, 10, 20]:
            i = initial_index
            read_index = 0
            while i < genome_length:
                forward_read = get_read(genome, i, read_length)
                reverse_read = get_read(genome, i + read_length, read_length)

                writeReads(readsF_R1, forward_read, read_index)
                writeReads(readsF_R2, complementRead(reverse_read[::-1]), read_index)
                read_index += 1

                # Counter
                #for j in range(i, i + (2 * read_length)):
                #    c[j] += 1
                
                if i < genome_length / 2:
                    mirrored_index = i
                else:
                    mirrored_index = genome_length - i

                #increment = int(read_length - ((2 * read_length * (1 - 1 / PTR)) / genome_length) * mirrored_index)
                # linear increase, top gets chopped off too aggressively
                #increment = read_length + int(read_length / PTR) - int((4 * read_length / (genome_length ^ 2)) * (1 - 1 / PTR) * (mirrored_index ^ 2) + read_length / PTR)
                # quadratic increase
                increment = int(read_length / (math.cos(((2 * math.pi) / genome_length) * i) * ((PTR - 1) / 2) + (PTR + 1) / 2))
                # sine increase

                assert increment > 0
                i += 2 * increment # double because we're doing both forward and reverse reads
        
        return c


        # num_reads = triangle_area + underneath_rect_area
        # num_reads = 0.5 * (PTR * min_triangle - min_triangle) * num_steps + min_triangle * num_steps
        # num_reads = (0.5 * PTR + 1 - 0.5) * min_triangle * num_steps
        num_steps: int = int(num_reads / ((0.5 * PTR + 0.5) * min_triangle))
        step_range: int = genome_length / (2 * num_steps)

        total: int = 0
        for i in range(num_steps):
            reads_for_step = triangle(i / num_steps, PTR, min_triangle)
            c[i] = reads_for_step
            for j in range(reads_for_step):
                start_index = int(
                    random.random() * 2 * step_range
                    + genome_length * (i / num_steps)
                    - step_range
                )
                forward_read = get_read(genome, start_index, read_length)
                reverse_read = get_read(
                    genome, start_index + read_length + gap, read_length
                )

                writeReads(readsF_R1, forward_read, total)
                writeReads(readsF_R2, complementRead(reverse_read[::-1]), total)

                total += 1

        return c


# Write the given read to the given file
# @param readsF is the file to write to
# @param read is the read to write
# @param it is the iterator to append to the read label
def writeReads(readsF: TextIOWrapper, read: str, it: int) -> None:
    readsF.write("@TEST" + str(it) + "\n")
    readsF.write(read + "\n")
    readsF.write("+\n")
    # Create a random quality string
    #rand_qual = ""
    #for i in range(len(read)):
    #    rand_qual += chr(random.randint(33, 126))
    #readsF.write(rand_qual + "\n")
    readsF.write("~" * len(read) + "\n")


def complementRead(read: str) -> str:
    comp = ""
    for c in read:
        if c == "A":
            comp += "T"
        if c == "T":
            comp += "A"
        if c == "G":
            comp += "C"
        if c == "C":
            comp += "G"
    return comp


# Generate n sets of test data
# @param genome is the path to the txt file containing only sequences (our test genome)
#        To manually prepare this file, just make a copy of your fasta and remove the
#        lines starting with >
# @param reads is the path to the output fastq file of generated reads WITHOUT the file extension
# @param n is the number of read file pairs to create
def generateN(genome: str, reads: str, n: int, reads_per_sample: int = 50000):
    for i in range(n):
        ptr = i / 4 + 1.5
        ptr_str = str(ptr).replace(".", "-")
        c = generateTestData(genome, reads + ptr_str, reads_per_sample, ptr)
        print(f"Finished {reads + ptr_str}")
        #print(f"{reads + ptr_str}: {sorted(c.items())}")
        #plt.bar(c.keys(), c.values())
        #plt.show()

        with open(reads + ptr_str + "_R1.fastq", "rb") as r1, open(
            reads + ptr_str + "_R2.fastq", "rb"
        ) as r2, gzip.open(reads + ptr_str + "_R1.fastq.gz", "wb") as w1, gzip.open(
            reads + ptr_str + "_R2.fastq.gz", "wb"
        ) as w2:
            w1.writelines(r1)
            w2.writelines(r2)
        os.remove(reads + ptr_str + "_R1.fastq")
        os.remove(reads + ptr_str + "_R2.fastq")


# multi-reads
# generateN("reference/akk-genome.fasta", "multi-reads/Akk", 3)
# generateN("reference/Bfragilis.fasta", "multi-reads/Bfrag", 3)
# generateN("reference/Ecoli.fasta", "multi-reads/Ecoli", 3)

# reads
# generateN("reference/Bfragilis.fasta", "reads/Bfrag", 3, 25000)

# new reads
# generateN("reference/Bfragilis.fasta", "new/Bfrag", 5, 100000)

# new multi-reads
#generateN("reference/Bfragilis.fasta", "new-multi/Bfrag", 3, 25000)
#generateN("reference/Ecoli.fasta", "new-multi/Ecoli", 3, 25000)
#generateN("reference/akk-genome.fasta", "new-multi/Akk", 3, 25000)

# tiny-multi-reads
# generateN("tiny-ref/akk-genome.fasta", "tiny-multi-reads/Akk", 3, 1500)
# generateN("tiny-ref/Bfragilis.fasta", "tiny-multi-reads/Bfrag", 3, 1500)
# generateN("tiny-ref/Ecoli.fasta", "tiny-multi-reads/Ecoli", 3, 1500)

# source-forge multi-reads
#generateN("source-forge-ref/Bacteroides_fragilis_638R_uid84217__NC_016776.fna", "source-forge-multi/Bfragilis", 3, 500000)
#generateN("source-forge-ref/Bifidobacterium_animalis_lactis_B420_uid163691__NC_017866.fna", "source-forge-multi/Bifidobacterium", 3, 500000)
#generateN("source-forge-ref/Streptomyces_scabiei_87_22_uid46531__NC_013929.fna", "source-forge-multi/Streptomyces", 3, 500000)
        
# source-forge reads
#generateN("source-forge-ref/Bacteroides_fragilis_638R_uid84217__NC_016776.fna", "source-forge-reads/Bfragilis", 4, 50000)
        
# sparse reads
generateN("source-forge-ref/Bacteroides_fragilis_638R_uid84217__NC_016776.fna", "sparse/Bfragilis", 6, 500000)