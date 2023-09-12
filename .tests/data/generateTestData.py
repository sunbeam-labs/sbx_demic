from io import TextIOWrapper
from random import randrange
from typing import NoReturn
import math
from collections import Counter
import gzip
import sys
import os


def f(x: int, l: int, a: int) -> int:
    return int(100 * (-a * math.cos(x / l * 2 * math.pi) + a + 1))


def triangle(position: float, scale: float) -> int:
    # if position < 0.5:
    if True:
        return int(scale * position * 5)
    else:
        return int(scale * (1 - position))


# Generates a triangular distribution of reads from input sequence data
# mimicing the distribution of reads from a growing bacterial population
# @param genome is the path to the txt file containing only sequences (our test genome)
#        To manually prepare this file, just make a copy of your fasta and remove the
#        lines starting with >
# @param reads is the path to the output fastq file of generated reads WITHOUT the file extension
def generateTestData(genome: str, reads: str, num: int, PTR: float) -> Counter:
    fullGenome: str = ""
    with open(genome) as genomeF:
        for r in genomeF.readlines():
            if r[0] != ">":
                fullGenome += r.strip()

    with open(reads + "_R1.fastq", "w") as readsF_R1, open(
        reads + "_R2.fastq", "w"
    ) as readsF_R2:
        # Generate n sequences from the replication origin (assumed to be the start of the file)
        # of random length (0 to len(genome)) and then take a random substring of length m
        # from that to add to the fastq file
        c: Counter = Counter()
        n: int = 100
        m: int = 125
        l: int = len(fullGenome)
        b: int = int(l / n)
        for it in range(0, n):
            readIndex: int = it * b
            # size: int = randrange(2*m+2, l+1)
            # readIndex: int = randrange(size-(2*m+1))
            # for it2 in range(0, f(readIndex, l, num)):
            for it2 in range(0, triangle(it / n, PTR * 100)):
                # gap: int = randrange(250)
                gap: int = 0
                startIndex: int = readIndex + randrange(b - 2 * m - gap)
                fRead: str = fullGenome[startIndex : startIndex + m]
                rRead: str = fullGenome[startIndex + m + gap : startIndex + 2 * m + gap]
                c[round(startIndex / b, 0)] += 1

                writeReads(readsF_R1, fRead, it * 300 + it2)
                writeReads(readsF_R2, complementRead(rRead[::-1]), it * 300 + it2)

        return c


# Write the given read to the given file
# @param readsF is the file to write to
# @param read is the read to write
# @param it is the iterator to append to the read label
def writeReads(readsF: TextIOWrapper, read: str, it: int) -> None:
    readsF.write("@TEST" + str(it) + "\n")
    readsF.write(read + "\n")
    readsF.write("+\n")
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
def generateN(genome: str, reads: str, n: int):
    for i in range(n):
        c = generateTestData(genome, reads + str(i), n, n + i / n)

        with open(reads + str(i) + "_R1.fastq", "rb") as r1, open(
            reads + str(i) + "_R2.fastq", "rb"
        ) as r2, gzip.open(reads + str(i) + "_R1.fastq.gz", "wb") as w1, gzip.open(
            reads + str(i) + "_R2.fastq.gz", "wb"
        ) as w2:
            w1.writelines(r1)
            w2.writelines(r2)
        os.remove(reads + str(i) + "_R1.fastq")
        os.remove(reads + str(i) + "_R2.fastq")


generateN("reference/akk-genome.fasta", "multi-reads/Akk", 3)
generateN("reference/Bfragilis.fasta", "multi-reads/Bfrag", 3)
generateN("reference/Ecoli.fasta", "multi-reads/Ecoli", 3)
